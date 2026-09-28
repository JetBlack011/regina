// cascadesearch.cpp
//
// Goal-directed chained searches for one target knot or link. See README.md.
//
// Usage:
//   cascadesearch --target-pd '<PD>' --target-name NAME --work DIR
//                 --verifyslicegenus PATH --knot-table CSV --link-table CSV
//                 [--knot-symmetry CSV] [--census-db PATH]
//                 [--goal-genus G] [--goal disjoint|connected]
//                 [--hop-surfaces N] [--max-hop-surfaces N] [--threads T]
//                 [--max-expansions K] [--cpu-budget SECONDS]
//                 [--max-crossings C] [--strategy best|dfs|bfs]
//                 [--literature | --constructive]
//
// Each expansion is one verifyslicegenus row (the campaign's search shape)
// run as a child process on a node's diagram; its witnesses become edges of
// the proof graph (hopedges.h), far sides become nodes (nodes.h), and bounds
// are relaxed to a fixed point (proofgraph.h) after every hop. The run stops
// when the target's goal has a proof, or the budget is spent.
//
// Writes, under DIR: hop_<k>_n<node>/ (each hop's row, witnesses and logs),
// cascade.jsonl (one line per hop), and certificate.json when the goal is met.

#include <spawn.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
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
#include <vector>

#include <link/link.h>

#include "exactnaming/exactnamer.h"
#include "exactnaming/exacttables.h"
#include "hopedges.h"
#include "leaves.h"
#include "nodes.h"
#include "proofgraph.h"
#include "witnesskey.h"

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
  int genus = 0, otherComponents = 0;
};

std::vector<Witness> readWitnesses(const std::string &path) {
  std::vector<Witness> out;
  std::ifstream in(path);
  if (!in) return out;
  std::string line;
  std::getline(in, line);
  while (std::getline(in, line)) {
    std::vector<std::string> f;
    std::string cur;
    bool q = false;
    for (char c : line) {
      if (c == '"') q = !q;
      else if (c == ',' && !q) { f.push_back(cur); cur.clear(); }
      else cur += c;
    }
    f.push_back(cur);
    if (f.size() < 9) continue; // torn last line
    out.push_back({f[3], f[8], std::stoi(f[6]), f[5].empty() ? 0 : std::stoi(f[5])});
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
  bool useful(NodeId n) const;
  std::optional<NodeId> choose();
  void expand(NodeId n, long surfaces);
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
  };
  std::map<EdgeId, EdgeInfo> edgeInfo_;
  std::map<std::string, EdgeInfo> directInfo_; // by witness key
  double cpuSpent_ = 0, wallSpent_ = 0;
  int hops_ = 0;
  int invariantFailures_ = 0;
};

void Cascade::onNewNode(NodeId n, int depth) {
  depth_.emplace(n, depth);
  if (!reg_.known(n) || n == reg_.unknot())
    return;
  const GaussDiagram &d = reg_.info(n).diagram;
  exactnaming::PieceName pn;
  try {
    pn = namer_.identify(d);
  } catch (const std::exception &e) {
    return;
  }
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

std::optional<NodeId> Cascade::choose() {
  struct Cand {
    NodeId n;
    std::tuple<int, double, int, size_t> key;
  };
  std::vector<Cand> cands;
  for (const auto &[n, dep] : depth_) {
    if (!reg_.known(n) || n == reg_.unknot() || refused_.count(n)) continue;
    const NodeInfo &ni = reg_.info(n);
    if (ni.diagram.crossings() > cfg_.maxCrossings) continue;
    if (!expansions_[n].empty()) continue; // one expansion per budget level
    if (n != target_ && !useful(n)) continue;
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
  std::unique_ptr<HopAssembler> hop;
  try {
    hop = std::make_unique<HopAssembler>(g_, reg_, row);
  } catch (const std::exception &e) {
    refused_.insert(n);
    log("{\"hop\":" + std::to_string(k) + ",\"node\":" + std::to_string(n) +
        ",\"refused\":\"" + jsonEscape(e.what()) + "\"}");
    std::cout << "[!] node " << n << " refused: " << e.what() << "\n";
    return;
  }
  const std::string rowName = "cascade_n" + std::to_string(n);
  std::ofstream(dir + "/input.csv") << "Name,PD Notation,Genus-4D\n"
                                    << rowName << "," << row.pd << ",[0;99]\n";
  std::vector<std::string> argv = {
      cfg_.verify, "--input", dir + "/input.csv", "--output", dir + "/out.csv",
      "--cobordisms", dir + "/cob.csv", "--census-db", cfg_.censusDb, "--no-census-updates",
      "--knot-table", cfg_.knotTable, "--link-table", cfg_.linkTable,
      "--max-crossings", "999", "--thicken-layers", "2", "--collar-layers", "2",
      "--max-faces", "5", "--iddfs-iterations", "2", "--iddfs-start", "4", "--iddfs-step", "1",
      "--root-budget-start", "50000", "--root-budget-growth", "2", "--no-cone", "--harvest",
      "--boundary-condition", "proper", "--research-settled", "--per-knot-time-limit", "7200",
      "--threads", std::to_string(cfg_.threads), "--surface-target", std::to_string(surfaces),
      "--resolve-unlinked", "--exact-far-side-names", "--no-retriangulate-on-miss"};
  ChildRun r = runChild(argv, dir + "/log.txt", dir + "/err.txt");
  cpuSpent_ += r.cpu;
  wallSpent_ += r.wall;
  expansions_[n].push_back(surfaces);

  const size_t nodesBefore = g_.nodeCount();
  const auto t0 = std::chrono::steady_clock::now();
  int assembled = 0, failed = 0;
  std::vector<Witness> ws = readWitnesses(dir + "/cob.csv");
  for (const Witness &w : ws) {
    const std::string key = witnesskey::witnessKey(w.pairsig);
    HopEdge e;
    try {
      e = hop->add({w.pairsig, w.genus, key});
    } catch (const std::logic_error &ex) {
      // A broken invariant (e.g. simplify changed a linking number): never
      // silently. The witness is dropped, which is sound; the run goes on.
      ++invariantFailures_;
      e.why = std::string("INVARIANT: ") + ex.what();
      std::cout << "[!!] witness " << key << " (" << w.other << "): " << e.why << "\n";
    } catch (const std::exception &ex) {
      e.why = std::string("exception: ") + ex.what();
    }
    if (e.ok) {
      ++assembled;
      if (e.direct) directInfo_[key] = {dir, row.pd, key, e};
      else edgeInfo_[e.edge] = {dir, row.pd, key, e};
    } else {
      ++failed;
      std::cout << "[!] witness " << witnesskey::witnessKey(w.pairsig) << " (" << w.other
                << "): " << e.why << "\n";
    }
  }
  for (size_t m = nodesBefore; m < g_.nodeCount(); ++m)
    onNewNode(static_cast<NodeId>(m), depth_.at(n) + 1);
  g_.propagate();
  g_.propagateLower();
  const int targetLower = g_.lower(target_, goalPartition(target_));
  const double assemble = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - t0).count();
  auto best = g_.best(target_, goalPartition(target_));
  std::ostringstream o;
  o << "{\"hop\":" << k << ",\"node\":" << n << ",\"crossings\":" << d.crossings()
    << ",\"surfaces\":" << surfaces << ",\"status\":" << r.status
    << ",\"wall\":" << std::fixed << std::setprecision(1) << r.wall << ",\"cpu\":" << r.cpu
    << ",\"assemble\":" << assemble << ",\"witnesses\":" << ws.size()
    << ",\"assembled\":" << assembled << ",\"failed\":" << failed
    << ",\"nodes\":" << g_.nodeCount() << ",\"new_nodes\":" << (g_.nodeCount() - nodesBefore)
    << ",\"records\":" << g_.recordCount()
    << ",\"target_best\":" << (best ? std::to_string(best->genus) : "null")
    << ",\"target_lower\":" << targetLower
    << ",\"contradictions\":" << g_.contradictions().size()
    << ",\"invariant_failures\":" << invariantFailures_ << "}";
  log(o.str());
  std::cout << "[+] hop " << k << ": node " << n << " (" << d.crossings() << " crossings, "
            << d.components() << " components): " << ws.size() << " witnesses, "
            << assembled << " assembled, " << (g_.nodeCount() - nodesBefore)
            << " new nodes; " << std::fixed << std::setprecision(0) << r.wall << " s wall, "
            << r.cpu << " s CPU; target best " << (best ? std::to_string(best->genus) : "none")
            << "\n";
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
      c << ",\"witness\":\"" << jsonEscape(we.key) << "\",\"in\":" << we.in
        << ",\"out\":" << we.out << ",\"shape\":{\"components\":" << we.shape.components
        << ",\"genus\":" << we.shape.genus << ",\"inComponent\":" << ints(we.shape.inComponent)
        << ",\"outComponent\":" << ints(we.shape.outComponent) << "},\"inMap\":"
        << ints(we.inMap) << ",\"outMap\":" << ints(we.outMap);
      if (auto it = edgeInfo_.find(rec.edge); it != edgeInfo_.end()) {
        const HopEdge &he = it->second.he;
        c << ",\"hop_dir\":\"" << jsonEscape(it->second.hopDir) << "\",\"row_pd\":\""
          << jsonEscape(it->second.rowPD) << "\",\"split_edge\":" << he.splitEdge
          << ",\"farCurveEdges\":[";
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
      if (auto it = directInfo_.find(key); it != directInfo_.end())
        c << ",\"witness\":\"" << jsonEscape(key) << "\",\"hop_dir\":\""
          << jsonEscape(it->second.hopDir) << "\",\"row_pd\":\""
          << jsonEscape(it->second.rowPD) << "\"";
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
    if (hops_ >= cfg_.maxExpansions) { std::cout << "[-] expansion limit\n"; break; }
    if (cpuSpent_ >= cfg_.cpuBudget) { std::cout << "[-] CPU budget spent\n"; break; }
    if (!g_.contradictions().empty()) {
      for (const auto &c : g_.contradictions()) std::cout << "[!!] CONTRADICTION: " << c << "\n";
      return 3;
    }
    auto n = choose();
    if (!n) {
      // Every useful node searched at this budget: search them again deeper
      // (a surface-target round only covers a prefix of the roots).
      if (budget * 2 > cfg_.maxHopSurfaces) { std::cout << "[-] nothing useful left\n"; break; }
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
    return 0;
  }
  return 1;
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
    else {
      std::cerr << "unknown argument " << a << "\n";
      return 2;
    }
  }
  if (c.targetPD.empty() || c.work.empty() || c.verify.empty() || c.knotTable.empty() ||
      c.linkTable.empty() || c.censusDb.empty()) {
    std::cerr << "cascadesearch: --target-pd, --work, --verifyslicegenus, --knot-table, "
                 "--link-table and --census-db are required\n";
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
