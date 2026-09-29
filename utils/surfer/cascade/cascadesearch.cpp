// cascadesearch.cpp
//
// Goal-directed chained searches for one target knot or link. See README.md.
//
// Usage:
//   cascadesearch --target-pd '<PD>' --target-name NAME --work DIR
//                 --knot-table CSV --link-table CSV --census-db PATH
//                 [--knot-symmetry CSV] [--master-witnesses CSV]
//                 [--hop-mode process|child] [--verifyslicegenus PATH]
//                 [--goal-genus G] [--goal disjoint|connected]
//                 [--hop-surfaces N] [--max-hop-surfaces N] [--threads T]
//                 [--max-expansions K] [--cpu-budget SECONDS]
//                 [--max-crossings C] [--strategy best|dfs|bfs]
//                 [--literature | --constructive]
//                 [--witness-store CSV --run-name NAME [--dedupe-against CSV]...]
//   cascadesearch --sign-only --work DIR --witness-store CSV
//                 --knot-table CSV --link-table CSV [--dedupe-against CSV]...
//
// With --witness-store, every surface a hop keeps is also recorded for the
// atlas (keptstore.h): its hop appends it to <hop dir>/kept.csv at once, and
// the run's end signs the ones whose witness is new (to the store and to
// every --dedupe-against file) and appends them to the store, as the sweep
// would have recorded them. A hop's subject is the target's own name, a
// node's proved table name, or cascade:<run name>/<target>/n<node>, so no two
// runs can ever record different links under one name. --sign-only does the
// end-of-run step alone, for a run that was killed.
//
// Each expansion is one row searched on a node's diagram with the
// campaign's search shape: in this process (hoprunner.h; the default), or
// as a verifyslicegenus child (--hop-mode child, which needs
// --verifyslicegenus). Its surfaces become edges of the proof graph
// (hopedges.h), far sides become nodes (nodes.h), and bounds are relaxed to
// a fixed point (proofgraph.h) after every hop. The run stops when the
// target's goal has a proof, or the budget is spent.
//
// Writes, under DIR: hop_<k>_n<node>/ (each hop's row and log, and a child
// hop's witnesses), cascade.jsonl (one line per hop), and certificate.json
// when the goal is met.

#include <spawn.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <link/link.h>

#include "csvwriter.h"
#include "exactnaming/exactnamer.h"
#include "exactnaming/exacttables.h"
#include "farsidenaming.h"
#include "hopedges.h"
#include "hoprunner.h"
#include "identifycomplement.h"
#include "keptstore.h"
#include "leaves.h"
#include "nodes.h"
#include "proofgraph.h"
#include "witnesskey.h"
#include "witnessstore.h"

extern char **environ;

using exactnaming::GaussDiagram;
using namespace cascade;
namespace fs = std::filesystem;

namespace {

struct Config {
  std::string targetPD, targetName, work, verify, knotTable, linkTable,
      knotSymmetry, censusDb;
  int goalGenus = 0;
  bool goalDisjoint = false; // goal on the singleton partition (disjoint discs)
  long hopSurfaces = 100000, maxHopSurfaces = 1000000;
  int threads = 10;
  int maxExpansions = 20;
  double cpuBudget = 3600.0 * 4;
  size_t maxCrossings = 24;
  std::string strategy = "best";
  bool literature = true;
  /// The atlas's master cobordisms.csv (read only): a node that IS a table
  /// entry the atlas searched gets that row's witnesses as free edges.
  std::string masterWitnesses;
  /// "process": each hop searched in this process (hoprunner.h); "child": a
  /// verifyslicegenus row per hop, its witnesses read back from their pair
  /// signatures. Both search the same thickening with the same shape.
  std::string hopMode = "process";
  /// Each hop's search shape: the campaign's (hosts.conf) unless the
  /// --hop-* options change it. Either hop mode uses it.
  HopShape hopShape;
  bool verbose = false;
  /// The atlas-format witness store every kept surface is recorded in
  /// (keptstore.h); empty for none. Read-only stores to deduplicate against
  /// (the master, say), and the run's name for cascade: subjects.
  std::string witnessStore, runName;
  std::vector<std::string> dedupeAgainst;
  bool signOnly = false;
  /// Write lower_report.jsonl at the end (writeLowerReport()); the sources
  /// file (atlas data/lower_bound_sources.csv: name,...,special) marks which
  /// literature bounds are not Lipschitz.
  bool lowerReport = false;
  std::string lowerSources;
};

std::string jsonEscape(const std::string &s) {
  std::string o;
  for (char c : s) {
    if (c == '"' || c == '\\') o += '\\';
    if (c == '\n') { o += "\\n"; continue; }
    o += c;
  }
  return o;
}

// A row PD as knotbuilder reads it (labels from 0 or 1) as a Regina link.
regina::Link linkFromRowPD(const std::string &pd) {
  std::vector<int> nums;
  std::string cur;
  for (char c : pd) {
    if (std::isdigit(static_cast<unsigned char>(c))) cur += c;
    else if (!cur.empty()) { nums.push_back(std::stoi(cur)); cur.clear(); }
  }
  if (!cur.empty()) nums.push_back(std::stoi(cur));
  if (nums.empty() || nums.size() % 4)
    throw std::invalid_argument("target PD: not a list of 4-tuples");
  const int shift = *std::min_element(nums.begin(), nums.end()) == 0 ? 1 : 0;
  std::vector<std::array<int, 4>> xs;
  for (size_t i = 0; i < nums.size(); i += 4)
    xs.push_back({nums[i] + shift, nums[i + 1] + shift, nums[i + 2] + shift,
                  nums[i + 3] + shift});
  return regina::Link::fromPD(xs.begin(), xs.end());
}

GaussDiagram of(const regina::Link &l) {
  std::vector<size_t> origin(l.countComponents());
  std::iota(origin.begin(), origin.end(), 0);
  return GaussDiagram::of(l, origin);
}

struct Witness {
  std::string other, pairsig;
  int genus = 0, otherComponents = 0, layers = 2;
};

// Splits one cobordisms.csv line (quoted fields may hold commas).
std::vector<std::string> csvFields(const std::string &line) {
  std::vector<std::string> f;
  std::string cur;
  bool q = false;
  for (char c : line) {
    if (c == '"') q = !q;
    else if (c == ',' && !q) { f.push_back(cur); cur.clear(); }
    else cur += c;
  }
  f.push_back(cur);
  return f;
}

/**
 * The atlas's master witness file, read only: one scan records where each
 * subject's lines start; a subject's witnesses are read on demand. Columns
 * as in cobordisms.csv: kind,subject,subject_components,other,
 * other_candidates,other_components,genus,tubed,pairsig,source_row,
 * thicken_layers,max_faces[,resolved_vertices].
 */
class MasterIndex {
public:
  /// A table name's base: orientation tag and a knot's mirror prefix
  /// dropped. Used only to FIND candidate witnesses; every one found is
  /// redrawn and identified exactly before it means anything.
  static std::string base(std::string name) {
    if (auto b = name.find('{'); b != std::string::npos) name.resize(b);
    if (name.size() > 1 && name[0] == 'm' && std::isdigit(static_cast<unsigned char>(name[1])))
      name.erase(0, 1);
    return name;
  }

  explicit MasterIndex(const std::string &path) : path_(path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read " + path);
    std::string line;
    std::getline(in, line);
    std::streamoff at = in.tellg();
    while (std::getline(in, line)) {
      const auto a = line.find(','), b = line.find(',', a + 1);
      if (a != std::string::npos && b != std::string::npos) {
        offsets_[line.substr(a + 1, b - a - 1)].push_back(at);
        // field 3 (other) follows subject_components
        const auto c = line.find(',', b + 1), d = c == std::string::npos
                                                      ? std::string::npos
                                                      : line.find(',', c + 1);
        if (d != std::string::npos)
          byOther_[base(line.substr(c + 1, d - c - 1))].push_back(at);
      }
      at = in.tellg();
    }
  }
  bool has(const std::string &subject) const { return offsets_.count(subject) > 0; }
  /// Witnesses of the row `subject`.
  std::vector<std::pair<std::string, Witness>> rows(const std::string &subject) const {
    auto it = offsets_.find(subject);
    return it == offsets_.end() ? std::vector<std::pair<std::string, Witness>>{}
                                : read(it->second, 1u << 30);
  }
  /// Witnesses of OTHER rows whose recorded far side has this base name (a
  /// hint only), with their subjects; at most `cap`.
  std::vector<std::pair<std::string, Witness>> byFarSide(const std::string &b, size_t cap) const {
    auto it = byOther_.find(b);
    return it == byOther_.end() ? std::vector<std::pair<std::string, Witness>>{}
                                : read(it->second, cap);
  }
  size_t subjects() const { return offsets_.size(); }

private:
  std::vector<std::pair<std::string, Witness>> read(const std::vector<std::streamoff> &offs,
                                                    size_t cap) const {
    std::vector<std::pair<std::string, Witness>> out;
    std::ifstream in(path_);
    std::string line;
    for (std::streamoff off : offs) {
      if (out.size() >= cap) break;
      in.seekg(off);
      if (!std::getline(in, line)) continue;
      auto f = csvFields(line);
      if (f.size() < 11) continue;
      out.push_back({f[1], {f[3], f[8], std::stoi(f[6]), f[5].empty() ? 0 : std::stoi(f[5]),
                            f[10].empty() ? 2 : std::stoi(f[10])}});
    }
    return out;
  }
  std::string path_;
  std::unordered_map<std::string, std::vector<std::streamoff>> offsets_;
  std::unordered_map<std::string, std::vector<std::streamoff>> byOther_;
};

std::vector<Witness> readWitnesses(const std::string &path) {
  std::vector<Witness> out;
  std::ifstream in(path);
  if (!in) return out;
  std::string line;
  std::getline(in, line);
  while (std::getline(in, line)) {
    auto f = csvFields(line);
    if (f.size() < 11) continue; // torn last line
    out.push_back({f[3], f[8], std::stoi(f[6]), f[5].empty() ? 0 : std::stoi(f[5]),
                   f[10].empty() ? 2 : std::stoi(f[10])});
  }
  return out;
}

// name -> the table's PD string, as the atlas's rows were searched from it.
std::map<std::string, std::string> tablePDs(const std::vector<std::string> &files) {
  std::map<std::string, std::string> out;
  for (const std::string &file : files) {
    std::ifstream in(file);
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line)) {
      auto f = csvFields(line);
      if (f.size() >= 2) out[f[0]] = f[1];
    }
  }
  return out;
}

struct ChildRun {
  int status = -1;
  double wall = 0, cpu = 0;
};

// Runs argv with stdout/stderr to files; returns wall and child CPU time.
ChildRun runChild(const std::vector<std::string> &argv, const std::string &out,
                  const std::string &err) {
  std::vector<char *> a;
  for (const auto &s : argv) a.push_back(const_cast<char *>(s.c_str()));
  a.push_back(nullptr);
  posix_spawn_file_actions_t fa;
  posix_spawn_file_actions_init(&fa);
  posix_spawn_file_actions_addopen(&fa, 1, out.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
  posix_spawn_file_actions_addopen(&fa, 2, err.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
  rusage before{}, after{};
  getrusage(RUSAGE_CHILDREN, &before);
  const auto t0 = std::chrono::steady_clock::now();
  pid_t pid;
  ChildRun r;
  if (posix_spawn(&pid, a[0], &fa, nullptr, a.data(), environ) != 0) {
    posix_spawn_file_actions_destroy(&fa);
    return r;
  }
  int st = 0;
  waitpid(pid, &st, 0);
  posix_spawn_file_actions_destroy(&fa);
  r.wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  getrusage(RUSAGE_CHILDREN, &after);
  auto secs = [](const timeval &t) { return t.tv_sec + t.tv_usec / 1e6; };
  r.cpu = secs(after.ru_utime) - secs(before.ru_utime) + secs(after.ru_stime) -
          secs(before.ru_stime);
  r.status = WIFEXITED(st) ? WEXITSTATUS(st) : -1;
  return r;
}

class Cascade {
public:
  explicit Cascade(Config c)
      : cfg_(std::move(c)), reg_(g_),
        tables_(exactnaming::ExactTables::load(cfg_.knotTable, cfg_.linkTable,
                                               cfg_.knotSymmetry)),
        namer_(tables_, cheapLimits()) {}

  int run();

private:
  static exactnaming::NamerLimits cheapLimits() {
    exactnaming::NamerLimits l;
    l.simplifyTries = 4;
    l.exhaustiveHeight = 0;
    l.searchHeight = -1;
    l.maxDeepCrossings = 0;
    l.tableSideHeight = -1;
    return l;
  }

  Partition goalPartition(NodeId n) const {
    const int k = g_.node(n).components;
    return cfg_.goalDisjoint ? Partition::singletons(k) : Partition::coarsest(k);
  }
  static bool goalMetIn(const ProofGraph &g, NodeId t, const Partition &p, int genus) {
    auto b = g.best(t, p);
    return b && b->genus <= genus;
  }
  bool goalMet() const {
    return goalMetIn(g_, target_, goalPartition(target_), cfg_.goalGenus);
  }

  void onNewNode(NodeId n, int depth);
  /// Names every node in `ns` (in parallel), then records each at `depth`
  /// with its table name, literature leaf and lower bound (in order).
  void onNewNodes(const std::vector<NodeId> &ns, int depth);
  void applyName(NodeId n, const exactnaming::PieceName &pn);
  std::vector<NodeId> nodesSince(size_t first) const {
    std::vector<NodeId> ns;
    for (size_t m = first; m < g_.nodeCount(); ++m) ns.push_back(static_cast<NodeId>(m));
    return ns;
  }
  bool useful(NodeId n) const;
  std::optional<NodeId> choose(bool freeOnly = false);
  void expand(NodeId n, long surfaces);
  /// The name node n's hop records its witnesses under: the target's own
  /// name, a proved table name, or cascade:<run>/<target>/n<n>.
  std::string subjectName(NodeId n) const;
  /// With --witness-store: signs and stores every kept surface of the run,
  /// and writes <work>/nodes.csv for the cascade: subjects. Idempotent.
  void storeWitnesses();
  void printOutcome(const std::string &outcome) const;
  void printProfile() const;
  /// With --lower-report: <work>/lower_report.jsonl, what each tabulated
  /// node's lower bound carries to the target (README.md, "Lower bounds").
  void writeLowerReport() const;
  void writeCertificate() const;
  void log(const std::string &line) {
    std::ofstream(cfg_.work + "/cascade.jsonl", std::ios::app) << line << "\n";
  }

  Config cfg_;
  ProofGraph g_;
  NodeRegistry reg_;
  exactnaming::ExactTables tables_;
  exactnaming::ExactNamer namer_;
  NodeId target_ = -1;
  std::string targetCanonical_;
  std::map<NodeId, int> depth_;
  std::map<NodeId, std::string> tableName_;
  std::map<NodeId, std::vector<long>> expansions_;
  std::set<NodeId> refused_;
  // What a checker needs to replay each witness edge (certificate.json).
  struct EdgeInfo {
    std::string hopDir, rowPD, key;
    HopEdge he;
    int layers = 2;
    std::string pairsig; ///< inline for master witnesses (no hop directory)
    /// The row diagram's component i is the node's rowNodeMap[i]: identity
    /// for a hop on the node's own diagram, the registry's map for a master
    /// row (the table's diagram).
    std::vector<int> rowNodeMap;
    /// An in-process hop's surface, as triangles of the row's thickening,
    /// and that thickening's digest (WitnessRedrawer::buildChecksum()). A
    /// certificate carries both; the checker rebuilds the row, refuses a
    /// different digest, and rebuilds the surface from its faces. No pair
    /// signature is ever computed for it (that is only for a witness bound
    /// for the atlas; pairSigsOf()).
    std::vector<int> faces;
    std::string build;
  };
  static void writeSurface(std::ostream &c, const EdgeInfo &info);
  std::optional<farside::SignatureTable> signatures_;
  std::unique_ptr<HopSearcher> searcher_;
  std::unique_ptr<MasterIndex> master_;
  std::map<std::string, std::string> tablePD_;
  std::set<NodeId> masterDone_;
  bool masterRowsFor(NodeId n, std::vector<std::string> *rows = nullptr) const;
  void loadMaster(NodeId n);
  std::map<EdgeId, EdgeInfo> edgeInfo_;
  std::map<std::string, EdgeInfo> directInfo_; // by witness key
  double cpuSpent_ = 0, wallSpent_ = 0;
  int hops_ = 0;
  int invariantFailures_ = 0;
  std::map<NodeId, std::string> hopSubject_; ///< each searched node's subject name
  bool stored_ = false;
  size_t storedAppended_ = 0;             ///< witnesses the store gained
  std::string stopReason_ = "nothing-useful"; ///< why the loop ended short of the goal
};

std::string Cascade::subjectName(NodeId n) const {
  if (n == target_ && cfg_.targetName != "target") return cfg_.targetName;
  if (auto it = tableName_.find(n); it != tableName_.end() && tables_.entry(it->second))
    return it->second;
  return "cascade:" + cfg_.runName + "/" + cfg_.targetName + "/n" + std::to_string(n);
}

void Cascade::storeWitnesses() {
  if (cfg_.witnessStore.empty() || stored_) return;
  stored_ = true;
  // The cascade: subjects, as the atlas's results/cascade/nodes.csv lists
  // them (cascade_record.py), so a later identity can be attached to each.
  {
    std::ofstream nodes(cfg_.work + "/nodes.csv");
    nodes << "name,components,crossings,pd,signs,gauss,label\n";
    for (const auto &[n, name] : hopSubject_) {
      if (name.rfind("cascade:", 0) != 0) continue;
      const GaussDiagram &d = reg_.info(n).diagram;
      std::ostringstream signs, gauss;
      signs << '[';
      for (size_t i = 0; i < d.signs.size(); ++i) signs << (i ? ", " : "") << d.signs[i];
      signs << ']';
      gauss << '[';
      for (size_t c = 0; c < d.comps.size(); ++c) {
        gauss << (c ? ", [" : "[");
        for (size_t j = 0; j < d.comps[c].size(); ++j) gauss << (j ? ", " : "") << d.comps[c][j];
        gauss << ']';
      }
      gauss << ']';
      nodes << csvField(name) << ',' << d.components() << ',' << d.crossings() << ','
            << csvField(rowPD(d)) << ',' << csvField(signs.str()) << ','
            << csvField(gauss.str()) << ",node " << n << '\n';
    }
  }
  cobordismgraph::NameTable names;
  witnessstore::loadNameTable(cfg_.knotTable, names);
  witnessstore::loadNameTable(cfg_.linkTable, names);
  const StoreResult s = storeKept(readKept(cfg_.work), cfg_.witnessStore, cfg_.dedupeAgainst,
                                  names, static_cast<unsigned>(cfg_.threads));
  storedAppended_ = s.appended;
  std::cout << "[+] witness store: " << s.kept << " kept, " << s.fresh << " new, "
            << s.appended << " appended to " << cfg_.witnessStore << " (signed in "
            << std::fixed << std::setprecision(0) << s.signSeconds << " s)\n";
}

void Cascade::onNewNode(NodeId n, int depth) { onNewNodes({n}, depth); }

void Cascade::onNewNodes(const std::vector<NodeId> &ns, int depth) {
  // Naming is most of what a hop does outside its search, and each node's
  // name is independent of the others', so the names are found on a pool
  // (ExactNamer is safe to share: its caches are locked, and its SnapPea
  // calls serialised) and applied here in node order, as one at a time would.
  std::vector<std::optional<exactnaming::PieceName>> names(ns.size());
  std::atomic<size_t> next{0};
  auto work = [&] {
    for (size_t i; (i = next.fetch_add(1)) < ns.size();) {
      const NodeId n = ns[i];
      if (!reg_.known(n) || n == reg_.unknot()) continue;
      try {
        names[i] = namer_.identify(reg_.info(n).diagram);
      } catch (const std::exception &) {
      }
    }
  };
  std::vector<std::thread> pool;
  const size_t threads = std::min<size_t>(std::max(cfg_.threads, 1), ns.size());
  for (size_t t = 1; t < threads; ++t) pool.emplace_back(work);
  work();
  for (auto &t : pool) t.join();
  for (size_t i = 0; i < ns.size(); ++i) {
    depth_.emplace(ns[i], depth);
    if (names[i]) applyName(ns[i], *names[i]);
  }
}

void Cascade::applyName(NodeId n, const exactnaming::PieceName &pn) {
  if (pn.by == exactnaming::PieceName::By::untabulated || !pn.pinned() || pn.names.size() != 1)
    return;
  const std::string &name = pn.names.front();
  tableName_[n] = name;
  const exactnaming::TableEntry *e = tables_.entry(name);
  if (!e) return;
  auto g4 = parseTableG4(e->g4);
  if (!g4) return;
  const auto [lo, hi] = *g4;
  g_.setGenusLowerBound(n, lo, "literature " + name + " " + e->g4);
  // Never let the target's own literature value prove the target, even
  // through a duplicate node of it (README.md, "Leaf facts").
  if (!mayUseLiteratureUpperBound(tables_.canonical(name), targetCanonical_,
                                  cfg_.literature))
    return;
  g_.addLeaf(n, Partition::coarsest(g_.node(n).components), hi,
             "literature " + name + " " + e->g4);
}

bool Cascade::masterRowsFor(NodeId n, std::vector<std::string> *rows) const {
  if (!master_) return false;
  auto it = tableName_.find(n);
  if (it == tableName_.end()) return false;
  // Every table entry of this link's class (one oriented link up to mirror
  // and global reversal) is the same node; each is a row of its own.
  const std::string canon = tables_.canonical(it->second);
  const exactnaming::TableEntry *e = tables_.entry(it->second);
  if (!e) return false;
  bool any = !master_->byFarSide(MasterIndex::base(it->second), 1).empty();
  for (const exactnaming::TableEntry *v : tables_.variants(e->base))
    if (tables_.canonical(v->name) == canon && master_->has(v->name)) {
      any = true;
      if (rows) rows->push_back(v->name);
    }
  return any;
}

void Cascade::loadMaster(NodeId n) {
  masterDone_.insert(n);
  std::vector<std::string> own;
  if (!masterRowsFor(n, &own)) return;
  const size_t nodesBefore = g_.nodeCount();
  int assembled = 0, failed = 0, refusedRows = 0;
  // Witnesses by subject row: this node's own rows (every witness), and
  // other rows whose recorded far side names this node's base (a hint: the
  // far side is redrawn and identified exactly like any other).
  std::map<std::string, std::vector<Witness>> bySubject;
  std::set<std::string> ownRows(own.begin(), own.end());
  for (const std::string &name : own)
    for (auto &[subj, w] : master_->rows(name)) bySubject[subj].push_back(w);
  size_t reverse = 0;
  for (auto &[subj, w] : master_->byFarSide(MasterIndex::base(tableName_[n]), 300))
    if (!ownRows.count(subj)) {
      bySubject[subj].push_back(w);
      ++reverse;
    }
  for (const auto &[name, ws] : bySubject) {
    const exactnaming::TableEntry *e = tables_.entry(name);
    auto pdIt = tablePD_.find(name);
    if (!e || pdIt == tablePD_.end()) continue;
    // The row as the atlas searched it: the table's own diagram and PD, a
    // node by the registry's exact tests. An own row must be THIS node.
    GaussDiagram d = of(e->diagram);
    NodeMatch m = reg_.intern(simplifyKeepingComponents(d), "row " + name);
    if (ownRows.count(name) && m.node != n) {
      ++refusedRows;
      std::cout << "[!] master row " << name << " did not intern as node " << n
                << " (got " << m.node << "); not used\n";
      continue;
    }
    std::map<int, std::vector<const Witness *>> byLayers;
    for (const Witness &w : ws) byLayers[w.layers].push_back(&w);
    for (const auto &[layers, group] : byLayers) {
      HopRow row;
      row.node = m.node;
      row.diagram = d;
      row.nodeMap = m.componentMap;
      row.pd = pdIt->second;
      row.layers = layers;
      std::unique_ptr<HopAssembler> hop;
      try {
        hop = std::make_unique<HopAssembler>(g_, reg_, row);
      } catch (const std::exception &ex) {
        ++refusedRows;
        std::cout << "[!] master row " << name << " refused: " << ex.what() << "\n";
        continue;
      }
      for (const Witness *w : group) {
        const std::string key = witnesskey::witnessKey(w->pairsig);
        HopEdge he;
        try {
          he = hop->add({w->pairsig, w->genus, "master:" + key});
        } catch (const std::logic_error &ex) {
          ++invariantFailures_;
          std::cout << "[!!] master witness " << key << ": INVARIANT: " << ex.what() << "\n";
        } catch (const std::exception &ex) {
          he.why = ex.what();
        }
        if (!he.ok) { ++failed; continue; }
        ++assembled;
        EdgeInfo info{"master", row.pd, key, he, layers, w->pairsig, row.nodeMap};
        if (he.direct) directInfo_["master:" + key] = info;
        else edgeInfo_[he.edge] = info;
      }
    }
  }
  expansions_[n].push_back(0); // counts as this level's expansion
  onNewNodes(nodesSince(nodesBefore), depth_.at(n) + 1);
  g_.propagate();
  g_.propagateLower();
  auto best = g_.best(target_, goalPartition(target_));
  std::ostringstream o;
  o << "{\"master\":\"" << jsonEscape(tableName_[n]) << "\",\"node\":" << n
    << ",\"own_rows\":" << own.size() << ",\"reverse_witnesses\":" << reverse
    << ",\"subject_rows\":" << bySubject.size() << ",\"refused_rows\":" << refusedRows
    << ",\"assembled\":" << assembled << ",\"failed\":" << failed
    << ",\"nodes\":" << g_.nodeCount() << ",\"target_best\":"
    << (best ? std::to_string(best->genus) : "null") << "}";
  log(o.str());
  std::cout << "[+] master rows of node " << n << " (" << tableName_[n] << "): "
            << assembled << " witnesses assembled, " << failed << " failed; target best "
            << (best ? std::to_string(best->genus) : "none") << "\n";
}

bool Cascade::useful(NodeId n) const {
  // What-if: give n the best profile it could conceivably have (every
  // partition its linking numbers allow, at its proved lower bound) and see
  // whether the target's goal would follow over the edges found so far.
  ProofGraph what = g_;
  const Node &node = what.node(n);
  // Optimistic only down to what is proved impossible: each partition at
  // its propagated lower bound (the linking condition included).
  for (const Partition &p : allPartitions(node.components)) {
    const int lo = g_.lower(n, p);
    if (lo < ProofGraph::kNoSurface)
      what.addLeaf(n, p, lo, "what-if");
  }
  what.propagate();
  return goalMetIn(what, target_, goalPartition(target_), cfg_.goalGenus);
}

// freeOnly: only nodes whose master rows are not yet loaded (no search).
std::optional<NodeId> Cascade::choose(bool freeOnly) {
  struct Cand {
    NodeId n;
    std::tuple<int, double, int, size_t> key;
  };
  std::vector<Cand> cands;
  for (const auto &[n, dep] : depth_) {
    if (!reg_.known(n) || n == reg_.unknot() || refused_.count(n)) continue;
    const NodeInfo &ni = reg_.info(n);
    if (ni.diagram.crossings() > cfg_.maxCrossings) continue;
    if (freeOnly) {
      if (masterDone_.count(n) || !masterRowsFor(n)) continue;
    } else if (!expansions_[n].empty()) {
      continue; // one expansion per budget level
    }
    const bool use = n == target_ || useful(n);
    if (cfg_.verbose && !freeOnly) {
      auto t = tableName_.find(n);
      std::cout << "    candidate " << n << " (" << ni.diagram.crossings() << "x, "
                << ni.diagram.components() << "c, "
                << (t == tableName_.end() ? "untabulated" : t->second) << ", depth " << dep
                << "): " << (use ? "useful" : "not useful")
                << (masterRowsFor(n) && !masterDone_.count(n) ? ", master rows" : "") << "\n";
    }
    if (!use) continue;
    const int crossings = static_cast<int>(ni.diagram.crossings());
    const double vol = ni.hyperbolic ? ni.volume : 1e9;
    int order = 0;
    if (cfg_.strategy == "dfs") order = -dep;
    else if (cfg_.strategy == "bfs") order = dep;
    cands.push_back({n, {order, crossings + 0.0 * vol, dep, static_cast<size_t>(n)}});
    (void)vol;
  }
  if (cands.empty()) return std::nullopt;
  // best: fewest crossings, then lower volume, then shallower; dfs/bfs: by depth first.
  std::sort(cands.begin(), cands.end(), [&](const Cand &a, const Cand &b) {
    if (cfg_.strategy == "best") {
      const NodeInfo &ia = reg_.info(a.n), &ib = reg_.info(b.n);
      if (ia.diagram.crossings() != ib.diagram.crossings())
        return ia.diagram.crossings() < ib.diagram.crossings();
      const double va = ia.hyperbolic ? ia.volume : 1e9, vb = ib.hyperbolic ? ib.volume : 1e9;
      if (va != vb) return va < vb;
      return depth_.at(a.n) < depth_.at(b.n);
    }
    return a.key < b.key;
  });
  return cands.front().n;
}

void Cascade::expand(NodeId n, long surfaces) {
  const int k = hops_++;
  const std::string dir = cfg_.work + "/hop_" + std::to_string(k) + "_n" + std::to_string(n);
  fs::create_directories(dir);
  const GaussDiagram &d = reg_.info(n).diagram;
  HopRow row;
  row.node = n;
  row.diagram = d;
  row.nodeMap.resize(d.components());
  std::iota(row.nodeMap.begin(), row.nodeMap.end(), 0);
  row.pd = rowPD(d);
  row.layers = 2;
  using clock = std::chrono::steady_clock;
  auto seconds = [](clock::time_point a, clock::time_point b) {
    return std::chrono::duration<double>(b - a).count();
  };
  const auto tRow = clock::now();
  std::unique_ptr<HopAssembler> hop;
  try {
    hop = std::make_unique<HopAssembler>(g_, reg_, row);
  } catch (const std::exception &e) {
    refused_.insert(n);
    // The row as given, to diagnose the refusal: its PD, whether Regina can
    // recover every orientation from that PD, and the node's own diagram.
    std::ostringstream gauss;
    gauss << "{\"signs\":[";
    for (size_t i = 0; i < d.signs.size(); ++i) gauss << (i ? "," : "") << d.signs[i];
    gauss << "],\"gauss\":[";
    for (size_t c = 0; c < d.comps.size(); ++c) {
      gauss << (c ? ",[" : "[");
      for (size_t j = 0; j < d.comps[c].size(); ++j) gauss << (j ? "," : "") << d.comps[c][j];
      gauss << "]";
    }
    gauss << "]}";
    log("{\"hop\":" + std::to_string(k) + ",\"node\":" + std::to_string(n) +
        ",\"refused\":\"" + jsonEscape(e.what()) + "\",\"pd\":\"" + jsonEscape(row.pd) +
        "\",\"pd_ambiguous\":" + (d.link().pdAmbiguous() ? "true" : "false") +
        ",\"diagram\":" + gauss.str() + "}");
    std::cout << "[!] node " << n << " refused: " << e.what() << "\n";
    return;
  }
  // Where a hop's time goes, logged per hop: the row's build and
  // certification, the search's setup and the search itself, adding its
  // surfaces, naming new nodes, and relaxing the graph.
  double rowSeconds = seconds(tRow, clock::now()), setupSeconds = 0, searchSeconds = 0;
  std::string roundsJson = "[]";
  size_t drainTail = 0;
  double drainTailSeconds = 0;
  // The hop's subject: what its witnesses are recorded under, and what its
  // log lines are named by (subjectName()).
  const std::string rowName = subjectName(n);
  hopSubject_[n] = rowName;
  const size_t nodesBefore = g_.nodeCount();
  int assembled = 0, failed = 0;
  size_t witnesses = 0;
  // One witness (or kept surface) into the graph. Its edge's key is its
  // provenance: the witness key of a child hop's witness, or hop<k>#<i> for
  // an in-process one, which has no pair signature.
  const std::string build = searcher_ ? hop->redrawer().buildChecksum() : std::string();
  auto take = [&](const std::string &key, const std::string &label,
                  const std::function<HopEdge()> &add, std::vector<int> faces) {
    HopEdge e;
    try {
      e = add();
    } catch (const std::logic_error &ex) {
      // A broken invariant (e.g. simplify changed a linking number): never
      // silently. The witness is dropped, which is sound; the run goes on.
      ++invariantFailures_;
      e.why = std::string("INVARIANT: ") + ex.what();
      std::cout << "[!!] witness " << key << " (" << label << "): " << e.why << "\n";
    } catch (const std::exception &ex) {
      e.why = std::string("exception: ") + ex.what();
    }
    if (e.ok) {
      ++assembled;
      const bool inProcess = !faces.empty();
      EdgeInfo info{dir, row.pd, key, e, 2, "", row.nodeMap, std::move(faces),
                    inProcess ? build : std::string()};
      if (e.direct) directInfo_[key] = std::move(info);
      else edgeInfo_[e.edge] = std::move(info);
    } else {
      ++failed;
      std::cout << "[!] witness " << key << " (" << label << "): " << e.why << "\n";
    }
  };

  ChildRun r;
  std::chrono::steady_clock::time_point t0;
  if (searcher_) {
    HopRun run;
    try {
      run = searcher_->run(hop->redrawer(), rowName, surfaces, 7200);
    } catch (const std::exception &e) {
      refused_.insert(n);
      log("{\"hop\":" + std::to_string(k) + ",\"node\":" + std::to_string(n) +
          ",\"refused\":\"" + jsonEscape(e.what()) + "\"}");
      std::cout << "[!] node " << n << " refused: " << e.what() << "\n";
      return;
    }
    r.status = run.accountingFailure.empty() ? 0 : 1;
    r.wall = run.wall;
    r.cpu = run.cpu;
    setupSeconds = run.setup;
    searchSeconds = run.search;
    std::ostringstream rj;
    rj << std::fixed << std::setprecision(2) << '[';
    for (size_t i = 0; i < run.rounds.size(); ++i) rj << (i ? "," : "") << run.rounds[i];
    roundsJson = rj.str() + ']';
    drainTail = run.drainTail;
    drainTailSeconds = run.drainTailSeconds;
    std::ofstream(dir + "/log.txt")
        << "[+] " << rowName << " " << row.pd << "\n[+] " << rowName << ": "
        << run.kept.size() << " kept, outcome " << run.outcome << "\n[+] " << rowName
        << ": accounting: " << run.accounting << "\n";
    // Every hop's accounting in the driver log too, in verifyslicegenus's
    // shape after the hop number, so a campaign audits each hop as it
    // audits a row (tools/orchestrate/audit_rows.py).
    std::cout << "[+] hop " << k << " " << rowName << ": accounting: " << run.accounting
              << "\n";
    if (!run.accountingFailure.empty())
      std::cout << "[!!] hop " << k << ": surface accounting failed -- "
                << run.accountingFailure << " (completeness only: nothing unsound "
                << "is recorded)\n";
    if (!cfg_.witnessStore.empty()) {
      // Every kept surface, durably, before the graph takes its faces.
      std::vector<PendingWitness> pending;
      pending.reserve(run.kept.size());
      for (const KeptSurface &ks : run.kept) {
        PendingWitness p{ks.witness, row.pd, row.layers, ks.faces};
        p.witness.sourceRow = rowName;
        p.witness.thickenLayers = row.layers;
        p.witness.maxFaces = cfg_.hopShape.maxFaces;
        pending.push_back(std::move(p));
      }
      appendKept(dir, pending);
    }
    t0 = std::chrono::steady_clock::now();
    witnesses = run.kept.size();
    for (size_t i = 0; i < run.kept.size(); ++i) {
      KeptSurface &ks = run.kept[i];
      const std::string key = "hop" + std::to_string(k) + "#" + std::to_string(i);
      take(key, ks.farName, [&] { return hop->addRead(ks.link, ks.genus, key); },
           std::move(ks.faces));
    }
  } else {
    std::ofstream(dir + "/input.csv") << "Name,PD Notation,Genus-4D\n"
                                      << rowName << "," << row.pd << ",[0;99]\n";
    const HopShape &shape = cfg_.hopShape;
    std::vector<std::string> argv = {
        cfg_.verify, "--input", dir + "/input.csv", "--output", dir + "/out.csv",
        "--cobordisms", dir + "/cob.csv", "--census-db", cfg_.censusDb,
        "--no-census-updates", "--knot-table", cfg_.knotTable, "--link-table",
        cfg_.linkTable, "--max-crossings", "999", "--thicken-layers", "2",
        "--collar-layers", "2", "--max-faces", std::to_string(shape.maxFaces),
        "--iddfs-iterations", std::to_string(shape.iddfsIterations), "--iddfs-start",
        std::to_string(shape.iddfsStart), "--iddfs-step", std::to_string(shape.iddfsStep),
        "--root-budget-start", std::to_string(shape.rootBudgetStart),
        "--root-budget-growth", std::to_string(shape.rootBudgetGrowth), "--no-cone",
        "--harvest", "--boundary-condition", "proper", "--research-settled",
        "--per-knot-time-limit", "7200", "--threads", std::to_string(cfg_.threads),
        "--surface-target", std::to_string(surfaces), "--resolve-unlinked",
        "--exact-far-side-names", "--no-retriangulate-on-miss",
        "--pending-surface-cap", std::to_string(shape.pendingSurfaceCap),
        "--petal-cache-limit", std::to_string(shape.petalCacheLimit),
        "--recognition-cache-limit", std::to_string(shape.recognitionCacheLimit),
        "--boundary-signature-cache-limit",
        std::to_string(shape.boundarySignatureCacheLimit)};
    r = runChild(argv, dir + "/log.txt", dir + "/err.txt");
    t0 = std::chrono::steady_clock::now();
    std::vector<Witness> ws = readWitnesses(dir + "/cob.csv");
    witnesses = ws.size();
    for (const Witness &w : ws) {
      const std::string key = witnesskey::witnessKey(w.pairsig);
      take(key, w.other, [&] { return hop->add({w.pairsig, w.genus, key}); }, {});
    }
  }
  cpuSpent_ += r.cpu;
  wallSpent_ += r.wall;
  expansions_[n].push_back(surfaces);
  const auto tNodes = clock::now();
  const double addSeconds = seconds(t0, tNodes);
  onNewNodes(nodesSince(nodesBefore), depth_.at(n) + 1);
  const auto tProp = clock::now();
  g_.propagate();
  g_.propagateLower();
  const double nodeSeconds = seconds(tNodes, tProp), propagateSeconds = seconds(tProp, clock::now());
  const int targetLower = g_.lower(target_, goalPartition(target_));
  const double assemble = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - t0).count();
  auto best = g_.best(target_, goalPartition(target_));
  std::ostringstream o;
  o << "{\"hop\":" << k << ",\"node\":" << n << ",\"crossings\":" << d.crossings()
    << ",\"surfaces\":" << surfaces << ",\"status\":" << r.status
    << ",\"wall\":" << std::fixed << std::setprecision(1) << r.wall << ",\"cpu\":" << r.cpu
    << ",\"assemble\":" << assemble << ",\"row\":" << rowSeconds
    << ",\"setup\":" << setupSeconds << ",\"search\":" << searchSeconds
    << ",\"add\":" << addSeconds << ",\"name_nodes\":" << nodeSeconds
    << ",\"propagate\":" << propagateSeconds << ",\"rounds\":" << roundsJson
    << ",\"drain_tail\":" << drainTail << ",\"drain_tail_s\":" << drainTailSeconds
    << ",\"witnesses\":" << witnesses
    << ",\"assembled\":" << assembled << ",\"failed\":" << failed
    << ",\"nodes\":" << g_.nodeCount() << ",\"new_nodes\":" << (g_.nodeCount() - nodesBefore)
    << ",\"records\":" << g_.recordCount()
    << ",\"target_best\":" << (best ? std::to_string(best->genus) : "null")
    << ",\"target_lower\":" << targetLower
    << ",\"contradictions\":" << g_.contradictions().size()
    << ",\"invariant_failures\":" << invariantFailures_ << "}";
  log(o.str());
  std::cout << "[+] hop " << k << ": node " << n << " (" << d.crossings() << " crossings, "
            << d.components() << " components): " << witnesses << " witnesses, "
            << assembled << " assembled, " << (g_.nodeCount() - nodesBefore)
            << " new nodes; " << std::fixed << std::setprecision(0) << r.wall << " s wall, "
            << r.cpu << " s CPU; target best " << (best ? std::to_string(best->genus) : "none")
            << "\n";
}

// How a certificate finds a witness's surface: a master witness's pair
// signature inline; an in-process one's faces in its row's thickening, with
// that thickening's digest; a child hop's by its key in the hop directory's
// witness file (nothing here).
void Cascade::writeSurface(std::ostream &c, const EdgeInfo &info) {
  if (!info.pairsig.empty()) c << ",\"pairsig\":\"" << jsonEscape(info.pairsig) << "\"";
  if (!info.faces.empty()) {
    c << ",\"faces\":[";
    for (size_t i = 0; i < info.faces.size(); ++i) c << (i ? "," : "") << info.faces[i];
    c << "],\"build\":\"" << info.build << "\"";
  }
}

void Cascade::writeCertificate() const {
  auto best = g_.best(target_, goalPartition(target_));
  if (!best) return;
  std::ofstream c(cfg_.work + "/certificate.json");
  c << "{\"target\":\"" << jsonEscape(cfg_.targetName) << "\",\"target_pd\":\""
    << jsonEscape(cfg_.targetPD) << "\",\"goal_genus\":" << cfg_.goalGenus
    << ",\"goal\":\"" << (cfg_.goalDisjoint ? "disjoint" : "connected") << "\",\"genus\":"
    << best->genus << ",\"records\":[\n";
  bool first = true;
  std::set<NodeId> nodes;
  for (RecordId r : g_.proof(best->record)) {
    const Record &rec = g_.record(r);
    nodes.insert(rec.node);
    c << (first ? "" : ",\n") << "{\"id\":" << r << ",\"node\":" << rec.node
      << ",\"partition\":\"" << rec.partition.str() << "\",\"genus\":" << rec.genus
      << ",\"kind\":\"" << kindName(rec.kind) << "\",\"edge\":" << rec.edge
      << ",\"source\":\"" << jsonEscape(rec.source) << "\",\"children\":[";
    for (size_t i = 0; i < rec.children.size(); ++i)
      c << (i ? "," : "") << rec.children[i];
    c << "]";
    auto ints = [](const auto &v) {
      std::ostringstream o;
      o << '[';
      for (size_t i = 0; i < v.size(); ++i) o << (i ? "," : "") << v[i];
      return o.str() + ']';
    };
    if (rec.kind == RecordKind::witnessForward || rec.kind == RecordKind::witnessReverse) {
      const WitnessEdge &we = g_.witness(rec.edge);
      nodes.insert(we.in);
      nodes.insert(we.out);
      const auto it = edgeInfo_.find(rec.edge);
      c << ",\"witness\":\"" << jsonEscape(we.key) << "\"";
      c << ",\"in\":" << we.in
        << ",\"out\":" << we.out << ",\"shape\":{\"components\":" << we.shape.components
        << ",\"genus\":" << we.shape.genus << ",\"inComponent\":" << ints(we.shape.inComponent)
        << ",\"outComponent\":" << ints(we.shape.outComponent) << "},\"inMap\":"
        << ints(we.inMap) << ",\"outMap\":" << ints(we.outMap);
      if (it != edgeInfo_.end()) {
        const HopEdge &he = it->second.he;
        c << ",\"hop_dir\":\"" << jsonEscape(it->second.hopDir) << "\",\"row_pd\":\""
          << jsonEscape(it->second.rowPD) << "\",\"layers\":" << it->second.layers
          << ",\"row_node_map\":" << ints(it->second.rowNodeMap);
        writeSurface(c, it->second);
        c << ",\"split_edge\":" << he.splitEdge << ",\"farCurveEdges\":[";
        for (size_t j = 0; j < he.farCurveEdges.size(); ++j)
          c << (j ? "," : "") << ints(he.farCurveEdges[j]);
        c << "],\"pieces\":[";
        for (size_t k = 0; k < he.pieces.size(); ++k) {
          const NodeMatch &m = he.pieces[k];
          c << (k ? "," : "") << "{\"node\":" << m.node << ",\"method\":\"" << m.method
            << "\",\"componentMap\":" << ints(m.componentMap) << ",\"mirrored\":"
            << (m.mirrored ? "true" : "false") << ",\"reversed\":"
            << (m.reversed ? "true" : "false") << ",\"origins\":" << ints(he.pieceOrigins[k])
            << "}";
          nodes.insert(m.node);
        }
        c << "]";
      }
    }
    if (rec.kind == RecordKind::splitCombine || rec.kind == RecordKind::splitRestrict) {
      const SplitEdge &se = g_.split(rec.edge);
      nodes.insert(se.whole);
      nodes.insert(se.pieces.begin(), se.pieces.end());
      c << ",\"whole\":" << se.whole << ",\"pieces\":" << ints(se.pieces) << ",\"pieceMap\":[";
      for (size_t k = 0; k < se.pieceMap.size(); ++k)
        c << (k ? "," : "") << ints(se.pieceMap[k]);
      c << "]";
    }
    if (rec.kind == RecordKind::leaf && rec.source.rfind("direct witness ", 0) == 0) {
      const std::string key = rec.source.substr(15);
      if (auto it = directInfo_.find(key); it != directInfo_.end()) {
        c << ",\"witness\":\"" << jsonEscape(it->second.key)
          << "\",\"hop_dir\":\"" << jsonEscape(it->second.hopDir) << "\",\"row_pd\":\""
          << jsonEscape(it->second.rowPD) << "\",\"layers\":" << it->second.layers;
        writeSurface(c, it->second);
      }
    }
    c << "}";
    first = false;
  }
  c << "\n],\"nodes\":[\n";
  first = true;
  for (NodeId n : nodes) {
    c << (first ? "" : ",\n") << "{\"id\":" << n << ",\"label\":\""
      << jsonEscape(g_.node(n).label) << "\",\"components\":" << g_.node(n).components;
    if (reg_.known(n) && reg_.info(n).diagram.crossings() > 0) {
      // The node's own diagram, as signed Gauss data: component maps refer
      // to ITS component order, which a PD round trip need not keep.
      const GaussDiagram &d = reg_.info(n).diagram;
      c << ",\"pd\":\"" << jsonEscape(rowPD(d)) << "\",\"signs\":[";
      for (size_t k = 0; k < d.signs.size(); ++k) c << (k ? "," : "") << d.signs[k];
      c << "],\"gauss\":[";
      for (size_t i = 0; i < d.comps.size(); ++i) {
        c << (i ? ",[" : "[");
        for (size_t j = 0; j < d.comps[i].size(); ++j) c << (j ? "," : "") << d.comps[i][j];
        c << "]";
      }
      c << "]";
    }
    if (auto it = tableName_.find(n); it != tableName_.end())
      c << ",\"table\":\"" << jsonEscape(it->second) << "\"";
    c << "}";
    first = false;
  }
  c << "\n]}\n";
}

int Cascade::run() {
  fs::create_directories(cfg_.work);
  if (cfg_.hopMode == "process") {
    // As a child hop runs verifyslicegenus: a private census copy, never
    // written, and no Pachner searches on a census miss.
    if (!census::setCensusPath(cfg_.censusDb))
      std::cout << "[!] census not found at " << cfg_.censusDb << "\n";
    census::retriangulateOnMiss.store(false);
    identify::recognitionCacheLimit.store(cfg_.hopShape.recognitionCacheLimit);
    const auto t0 = std::chrono::steady_clock::now();
    signatures_ = farside::SignatureTable::fromTables(cfg_.knotTable, cfg_.linkTable);
    // The hops' far-side namers share the node namer's table caches.
    searcher_ = std::make_unique<HopSearcher>(*signatures_, &tables_, cfg_.hopShape,
                                              static_cast<unsigned>(cfg_.threads),
                                              namer_.caches());
    std::cout << "[+] hops in process: " << signatures_->knots() << " knot and "
              << signatures_->links() << " link diagram signatures ("
              << std::fixed << std::setprecision(1)
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()
              << " s)\n";
  } else if (cfg_.hopMode != "child") {
    throw std::invalid_argument("--hop-mode must be process or child");
  }
  {
    const HopShape &s = cfg_.hopShape;
    std::cout << "[+] hop shape: cap " << s.maxFaces << ", IDDFS " << s.iddfsIterations
              << " from " << s.iddfsStart << " step " << s.iddfsStep << ", root budget "
              << s.rootBudgetStart << " x" << s.rootBudgetGrowth << "\n";
  }
  printProfile();
  if (!cfg_.masterWitnesses.empty()) {
    const auto t0 = std::chrono::steady_clock::now();
    master_ = std::make_unique<MasterIndex>(cfg_.masterWitnesses);
    tablePD_ = tablePDs({cfg_.knotTable, cfg_.linkTable});
    std::cout << "[+] master witnesses: " << master_->subjects() << " subjects indexed in "
              << std::fixed << std::setprecision(0)
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()
              << " s (read only)\n";
  }
  GaussDiagram raw = of(linkFromRowPD(cfg_.targetPD));
  GaussDiagram simp = simplifyKeepingComponents(raw);
  if (exactnaming::splitPieces(simp).size() != 1)
    throw std::runtime_error("the target is a split diagram; give one piece");
  // The target's own table name, so its literature value never proves it.
  try {
    auto pn = namer_.identify(simp);
    if (pn.names.size() == 1) targetCanonical_ = tables_.canonical(pn.names.front());
  } catch (...) {
  }
  if (targetCanonical_.empty() && tables_.entry(cfg_.targetName))
    targetCanonical_ = tables_.canonical(cfg_.targetName);
  NodeMatch t = reg_.intern(simp, "target " + cfg_.targetName);
  target_ = t.node;
  onNewNode(target_, 0);
  std::cout << "[+] target " << cfg_.targetName << " = node " << target_ << " ("
            << simp.crossings() << " crossings, " << simp.components()
            << " components; table class '" << targetCanonical_ << "'); goal genus "
            << cfg_.goalGenus << (cfg_.goalDisjoint ? " (disjoint pieces)" : " (connected)")
            << "\n";
  const auto start = std::chrono::steady_clock::now();
  long budget = cfg_.hopSurfaces;
  while (!goalMet()) {
    if (!g_.contradictions().empty()) {
      for (const auto &c : g_.contradictions()) std::cout << "[!!] CONTRADICTION: " << c << "\n";
      // The surfaces are real whatever the contradiction's cause (a naming
      // or solver bug), so they are kept, as verifyslicegenus writes its
      // witnesses before its fatal-bug halt.
      storeWitnesses();
      printOutcome("contradiction");
      return 3;
    }
    // Free edges first (no search, so no budget): the target's own master
    // rows, if it is a table entry the atlas searched.
    if (master_ && !masterDone_.count(target_) && masterRowsFor(target_)) {
      loadMaster(target_);
      continue;
    }
    // Then any useful node the atlas already searched: its rows are free.
    if (master_)
      if (auto f = choose(/*freeOnly=*/true)) {
        loadMaster(*f);
        continue;
      }
    auto n = choose();
    if (hops_ >= cfg_.maxExpansions) {
      std::cout << "[-] expansion limit\n";
      stopReason_ = "expansion-limit";
      break;
    }
    if (cpuSpent_ >= cfg_.cpuBudget) {
      std::cout << "[-] CPU budget spent\n";
      stopReason_ = "cpu-budget";
      break;
    }
    if (!n) {
      // Every useful node searched at this budget: search them again deeper
      // (a surface-target round only covers a prefix of the roots).
      if (budget * 2 > cfg_.maxHopSurfaces) {
        std::cout << "[-] nothing useful left\n";
        stopReason_ = "nothing-useful";
        break;
      }
      budget *= 2;
      for (auto &[m, v] : expansions_) v.clear();
      std::cout << "[+] raising the hop budget to " << budget << " surfaces\n";
      continue;
    }
    expand(*n, budget);
  }
  const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  auto best = g_.best(target_, goalPartition(target_));
  std::cout << "[+] done: " << hops_ << " hops, " << g_.nodeCount() << " nodes, "
            << g_.recordCount() << " records; " << std::fixed << std::setprecision(0)
            << wall << " s wall, " << cpuSpent_ << " s search CPU. Target best: "
            << (best ? std::to_string(best->genus) : "none") << "\n";
  storeWitnesses();
  writeLowerReport();
  if (goalMet()) {
    writeCertificate();
    bool constructive = true;
    for (RecordId r : g_.proof(best->record))
      if (g_.record(r).kind == RecordKind::leaf &&
          g_.record(r).source.rfind("literature", 0) == 0)
        constructive = false;
    std::cout << "[+] GOAL MET: " << cfg_.targetName << " genus <= " << best->genus << " ("
              << (constructive ? "constructive" : "literature-assisted") << "); certificate "
              << cfg_.work << "/certificate.json\n";
    printOutcome("met");
    return 0;
  }
  printOutcome(stopReason_);
  return 1;
}

void Cascade::writeLowerReport() const {
  // For every tabulated node Y: the least charge of carrying a lower bound
  // from Y to the target, over every path the graph holds. Measured by
  // seeding Y alone at a large M in a copy and reading what reaches the
  // target (propagateLower() takes the maximum over sources, and M dwarfs
  // every real bound), so charge = M - lower(target). A literature bound
  // lo(Y) carries lo(Y) - charge; hi(Y) - charge is the most Y could ever
  // carry. README.md, "Lower bounds": only a charge-0 path to a source whose
  // bound is not Lipschitz (lower_bound_sources.csv, `special`) can beat the
  // target's own literature bound.
  if (!cfg_.lowerReport) return;
  std::map<std::string, bool> special;
  if (!cfg_.lowerSources.empty()) {
    std::ifstream in(cfg_.lowerSources);
    std::string line;
    std::getline(in, line);
    std::vector<std::string> head = parseCsvLine(line);
    const auto col = [&](const std::string &c) {
      return static_cast<size_t>(std::find(head.begin(), head.end(), c) - head.begin());
    };
    const size_t nameCol = col("name"), specialCol = col("special");
    while (std::getline(in, line)) {
      const std::vector<std::string> f = parseCsvLine(line);
      if (nameCol < f.size() && specialCol < f.size())
        special[f[nameCol]] = f[specialCol] == "1";
    }
  }
  const Partition goal = goalPartition(target_);
  const int targetLower = g_.lower(target_, goal);
  int litLo = -1;
  if (const exactnaming::TableEntry *e = tables_.entry(cfg_.targetName))
    if (auto g4 = parseTableG4(e->g4)) litLo = g4->first;
  std::ofstream out(cfg_.work + "/lower_report.jsonl");
  out << "{\"target\":\"" << jsonEscape(cfg_.targetName) << "\",\"target_lower\":" << targetLower
      << ",\"lit_lo\":" << litLo << ",\"nodes\":" << g_.nodeCount() << "}\n";
  // What node n's lower bound `seed` alone carries to the target: every other
  // lower bound forgotten (clearLowerBounds()), n seeded, relaxed. Only
  // consistent facts are ever seeded -- a value the node could really have,
  // at most its best proved genus -- or the split rules, which read proved
  // surfaces, would pump bounds without limit. -1 when the what-if itself
  // meets a contradiction (then nothing it says is used).
  auto carried = [&](NodeId n, int seed) {
    ProofGraph what = g_;
    what.clearLowerBounds();
    const size_t before = what.contradictions().size();
    what.setGenusLowerBound(n, seed, "lower-report what-if");
    what.propagateLower();
    if (what.contradictions().size() != before) return -1;
    return what.lower(target_, goal);
  };
  int bestCarry = 0, bestCould = 0;
  std::string bestName, bestCouldName;
  for (const auto &[n, name] : tableName_) {
    if (n == target_) continue;
    const exactnaming::TableEntry *e = tables_.entry(name);
    if (!e) continue;
    auto g4 = parseTableG4(e->g4);
    if (!g4) continue;
    // The most n could be: its literature upper end, or less if a surface
    // for it is already proved.
    int could = g4->second;
    if (auto b = g_.bestConnected(n)) could = std::min(could, b->genus);
    const int carries = g4->first > 0 ? carried(n, g4->first) : 0;
    const int couldCarry = could > 0 ? (could == g4->first ? carries : carried(n, could)) : 0;
    if (carries <= 0 && couldCarry <= 0) continue; // n reaches the target with nothing
    const auto sp = special.find(name);
    out << "{\"node\":" << n << ",\"name\":\"" << jsonEscape(name) << "\",\"lit_lo\":"
        << g4->first << ",\"lit_hi\":" << g4->second << ",\"could\":" << could
        << ",\"special\":" << (sp == special.end() ? "null" : sp->second ? "true" : "false")
        << ",\"carries\":" << carries << ",\"could_carry\":" << couldCarry << "}\n";
    if (carries > bestCarry) {
      bestCarry = carries;
      bestName = name;
    }
    if (couldCarry > bestCould) {
      bestCould = couldCarry;
      bestCouldName = name;
    }
  }
  std::cout << "[+] lower report: target lower " << targetLower << " (literature " << litLo
            << "); best carried " << (bestName.empty() ? std::string("none")
                                                       : std::to_string(bestCarry) + " from " + bestName)
            << "; most any tabulated node could carry "
            << (bestCouldName.empty() ? std::string("none")
                                      : std::to_string(bestCould) + " from " + bestCouldName)
            << "\n";
}

void Cascade::printOutcome(const std::string &outcome) const {
  // The line a campaign's runner parses, in verifyslicegenus's own shape
  // (dispatch.py RE_OUTCOME): witnesses newly recorded, and why the run ended.
  std::cout << "[+] " << cfg_.targetName << ": " << storedAppended_
            << " new witnesses, outcome " << outcome << "\n";
}

void Cascade::printProfile() const {
  // Everything that decides what a run covers, as key=value, so a campaign
  // records what actually ran rather than what its configuration asked for.
  const HopShape &s = cfg_.hopShape;
  std::cout << "[+] profile: goal=" << (cfg_.goalDisjoint ? "disjoint" : "connected")
            << " goal_genus=" << cfg_.goalGenus << " literature=" << (cfg_.literature ? 1 : 0)
            << " hop_mode=" << cfg_.hopMode << " hop_surfaces=" << cfg_.hopSurfaces
            << " max_hop_surfaces=" << cfg_.maxHopSurfaces
            << " max_expansions=" << cfg_.maxExpansions << " cpu_budget=" << cfg_.cpuBudget
            << " strategy=" << cfg_.strategy << " max_crossings=" << cfg_.maxCrossings
            << " threads=" << cfg_.threads << " max_faces=" << s.maxFaces
            << " iddfs_iterations=" << s.iddfsIterations << " iddfs_start=" << s.iddfsStart
            << " iddfs_step=" << s.iddfsStep << " root_budget_start=" << s.rootBudgetStart
            << " root_budget_growth=" << s.rootBudgetGrowth << " layers=2"
            << " resolve_unlinked=" << (s.resolveUnlinked ? 1 : 0)
            << " exact_far_side_names=1 pending_surface_cap=" << s.pendingSurfaceCap
            << " petal_cache_limit=" << s.petalCacheLimit
            << " boundary_signature_cache_limit=" << s.boundarySignatureCacheLimit
            << " recognition_cache_limit=" << s.recognitionCacheLimit
            << " master_witnesses=" << (cfg_.masterWitnesses.empty() ? "none" : cfg_.masterWitnesses)
            << " witness_store=" << (cfg_.witnessStore.empty() ? "none" : cfg_.witnessStore)
            << " run_name=" << (cfg_.runName.empty() ? "none" : cfg_.runName) << "\n";
}

} // namespace

int main(int argc, char **argv) {
  Config c;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) throw std::invalid_argument("missing value for " + a);
      return argv[++i];
    };
    if (a == "--target-pd") c.targetPD = next();
    else if (a == "--target-name") c.targetName = next();
    else if (a == "--work") c.work = next();
    else if (a == "--verifyslicegenus") c.verify = next();
    else if (a == "--knot-table") c.knotTable = next();
    else if (a == "--link-table") c.linkTable = next();
    else if (a == "--knot-symmetry") c.knotSymmetry = next();
    else if (a == "--census-db") c.censusDb = next();
    else if (a == "--goal-genus") c.goalGenus = std::stoi(next());
    else if (a == "--goal") c.goalDisjoint = next() == "disjoint";
    else if (a == "--hop-surfaces") c.hopSurfaces = std::stol(next());
    else if (a == "--max-hop-surfaces") c.maxHopSurfaces = std::stol(next());
    else if (a == "--threads") c.threads = std::stoi(next());
    else if (a == "--max-expansions") c.maxExpansions = std::stoi(next());
    else if (a == "--cpu-budget") c.cpuBudget = std::stod(next());
    else if (a == "--max-crossings") c.maxCrossings = std::stoul(next());
    else if (a == "--strategy") c.strategy = next();
    else if (a == "--literature") c.literature = true;
    else if (a == "--constructive") c.literature = false;
    else if (a == "--master-witnesses") c.masterWitnesses = next();
    else if (a == "--hop-mode") c.hopMode = next();
    else if (a == "--hop-max-faces") c.hopShape.maxFaces = std::stoll(next());
    else if (a == "--hop-iddfs-start") c.hopShape.iddfsStart = std::stoll(next());
    else if (a == "--hop-iddfs-iterations") c.hopShape.iddfsIterations = std::stoul(next());
    else if (a == "--hop-root-budget") c.hopShape.rootBudgetStart = std::stoll(next());
    else if (a == "--hop-iddfs-step") c.hopShape.iddfsStep = std::stoll(next());
    else if (a == "--hop-root-growth") c.hopShape.rootBudgetGrowth = std::stoll(next());
    else if (a == "--hop-pending-cap") c.hopShape.pendingSurfaceCap = std::stoull(next());
    else if (a == "--hop-petal-cache") c.hopShape.petalCacheLimit = std::stoull(next());
    else if (a == "--hop-boundary-cache")
      c.hopShape.boundarySignatureCacheLimit = std::stoull(next());
    else if (a == "--hop-recognition-cache")
      c.hopShape.recognitionCacheLimit = std::stoull(next());
    else if (a == "--verbose") c.verbose = true;
    else if (a == "--witness-store") c.witnessStore = next();
    else if (a == "--run-name") c.runName = next();
    else if (a == "--dedupe-against") c.dedupeAgainst.push_back(next());
    else if (a == "--sign-only") c.signOnly = true;
    else if (a == "--lower-report") c.lowerReport = true;
    else if (a == "--lower-sources") c.lowerSources = next();
    else {
      std::cerr << "unknown argument " << a << "\n";
      return 2;
    }
  }
  if (c.signOnly) {
    // The end-of-run store step alone, for a run that was killed after its
    // hops wrote kept.csv (the subjects are already in those lines).
    if (c.work.empty() || c.witnessStore.empty() || c.knotTable.empty() ||
        c.linkTable.empty()) {
      std::cerr << "cascadesearch --sign-only: --work, --witness-store, --knot-table and "
                   "--link-table are required\n";
      return 2;
    }
    try {
      cobordismgraph::NameTable names;
      witnessstore::loadNameTable(c.knotTable, names);
      witnessstore::loadNameTable(c.linkTable, names);
      const StoreResult s = storeKept(readKept(c.work), c.witnessStore, c.dedupeAgainst,
                                      names, static_cast<unsigned>(c.threads));
      std::cout << "[+] witness store: " << s.kept << " kept, " << s.fresh << " new, "
                << s.appended << " appended to " << c.witnessStore << "\n";
      return 0;
    } catch (const std::exception &e) {
      std::cerr << "cascadesearch: " << e.what() << "\n";
      return 2;
    }
  }
  if (!c.witnessStore.empty() && (c.runName.empty() || c.hopMode != "process")) {
    // A child hop writes its own witnesses (hop dir cob.csv), under its own
    // row name; only in-process hops are recorded through keptstore.h.
    std::cerr << "cascadesearch: --witness-store needs --run-name, and in-process hops\n";
    return 2;
  }
  if (c.targetPD.empty() || c.work.empty() || c.knotTable.empty() ||
      c.linkTable.empty() || c.censusDb.empty() ||
      (c.hopMode == "child" && c.verify.empty())) {
    std::cerr << "cascadesearch: --target-pd, --work, --knot-table, --link-table and "
                 "--census-db are required, and --verifyslicegenus with --hop-mode child\n";
    return 2;
  }
  if (c.targetName.empty()) c.targetName = "target";
  try {
    Cascade cascade(c);
    return cascade.run();
  } catch (const std::exception &e) {
    std::cerr << "cascadesearch: " << e.what() << "\n";
    return 2;
  }
}
