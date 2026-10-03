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
//                 (--resolve-unlinked | --no-resolve-unlinked)
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
#include <condition_variable>
#include <functional>
#include <mutex>
#include <limits>
#include <cmath>
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

#include "surfer/report/csvwriter.h"
#include "cobound/driver/timers.h"
#include "cobound/json.h"
#include "cobound/parallelfor.h"
#include "linknaming/linknamer.h"
#include "linknaming/tables.h"
#include "cobound/outgoing/outgoingnamer.h"
#include "cobound/bounds/searchcobordisms.h"
#include "cobound/outgoing/fromdatabase.h"
#include "cobound/cobordisms/pairsigner.h"
#include "cobound/search/search.h"
#include "linknaming/census/censusnaming.h"
#include "linknaming/complement/complementcache.h"
#include "linknaming/complement/unlinknaming.h"
#include "linknaming/names.h"
#include "cobound/cobordisms/pending.h"
#include "cobound/bounds/axioms.h"
#include "cobound/bounds/links.h"
#include "cobound/bounds/cobordismgraph.h"
#include "cobound/cobordisms/cobordismkey.h"
#include "cobound/cobordisms/database.h"
#include "cobound/solver/literature.h"

extern char **environ;

using exactnaming::GaussDiagram;
using namespace cascade;
namespace fs = std::filesystem;

namespace {

using timers::Clock;
using timers::secondsSince;

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
  std::string pairSigCache; ///< --pair-sig-cache: stored pair-signature contexts
  /// --read-back-cache: master witnesses' read-backs kept across runs
  /// (readbackcache.h); empty for none.
  std::string readBackCache;
  std::vector<std::string> dedupeAgainst;
  bool signOnly = false;
  /// Write lower_report.jsonl at the end (writeLowerReport()); the sources
  /// file (atlas data/lower_bound_sources.csv: name,...,special) marks which
  /// literature bounds are not Lipschitz.
  bool lowerReport = false;
  std::string lowerSources;
  /// --goal-lower G: also stop once lower(target, goal partition) >= G
  /// (README.md, "Lower-bound mode"); -1 for no lower goal. Needs
  /// --lower-sources: the special sources' largest literature lower bound
  /// caps what any node could ever carry, which is what makes the lower
  /// gate prune.
  int goalLower = -1;
  /// A node kept only for the lower goal is never expanded above this many
  /// crossings: the chain must come back to a table entry.
  size_t lowerMaxCrossings = 16;
  /// --master-loads lazy|eager: a node's stored rows loaded when it is
  /// about to be expanded (the default since 2026-09-30), or as soon as it
  /// is met and useful (every table node the graph reaches).
  bool lazyMasterLoads = true;
  /// Hub breadth (John, 2026-09-29): a node with at least hubDegree witness
  /// edges is expanded once at hubSurfaces or more, as verifyslicegenus's
  /// wide rows are, for many more first-level far sides. 0: off.
  size_t hubDegree = 0;
  long hubSurfaces = 0;
};

GaussDiagram of(const regina::Link &l) {
  std::vector<size_t> origin(l.countComponents());
  std::iota(origin.begin(), origin.end(), 0);
  return GaussDiagram::of(l, origin);
}

struct Witness {
  std::string other, pairsig;
  int genus = 0, otherComponents = 0, layers = 2;
  /// The PD its row was searched on, when the file's `.rows.csv` sidecar
  /// records it (a cascade hop's row is its node's simplified diagram, not
  /// the table's); empty for the table's PD.
  std::string rowPD;
};

// A database line as the cascade carries it (cobordisms/database's reader
// parses it): the outgoing name, the pair signature, genus, outgoing
// components, layers and the searched row's PD.
Witness carried(const witnessstore::StoredCobordism &s) {
  return {s.witness.other, s.witness.pairSig, s.witness.genus, s.witness.otherComponents,
          s.witness.thickenLayers, s.rowPD};
}

std::vector<Witness> readWitnesses(const std::string &path) {
  std::vector<Witness> out;
  for (cobordismgraph::Witness &w : witnessstore::readWitnesses(path))
    out.push_back({w.other, std::move(w.pairSig), w.genus, w.otherComponents, w.thickenLayers});
  return out;
}

// name -> the table's PD string, as the atlas's rows were searched from it.
std::map<std::string, std::string> tablePDs(const std::vector<std::string> &files) {
  std::map<std::string, std::string> out;
  for (const std::string &file : files)
    for (const exactnaming::TableRow &row : exactnaming::readTableRows(file))
      out[row.name] = row.pd;
  return out;
}

// The tables' names and literature bounds (the store step's candidate sets)
// and, with --knot-symmetry, the knots' symmetry types (the slice-composite
// anchors): one NameTable per process, as verifyslicegenus loads it.
cobordismgraph::NameTable loadTableNames(const std::string &knotTable,
                                         const std::string &linkTable,
                                         const std::string &knotSymmetry,
                                         size_t *symmetryTypes = nullptr) {
  cobordismgraph::NameTable names;
  witnessstore::loadNameTable(knotTable, names);
  witnessstore::loadNameTable(linkTable, names);
  if (!knotSymmetry.empty()) {
    const exactnaming::SymmetryTable types = exactnaming::readSymmetryTable(knotSymmetry);
    for (const auto &[knot, type] : types) names.setSymmetry(knot, type);
    if (symmetryTypes) *symmetryTypes = types.size();
  }
  return names;
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
        namer_(tables_, NodeAxioms::namerLimits()),
        axioms_(g_, reg_, tables_, namer_, names_.symmetries(),
                {.literature = cfg_.literature,
                 .classes = true,
                 .threads = static_cast<unsigned>(std::max(cfg_.threads, 1)),
                 .log = &std::cout}) {}

  int run();

private:
  Partition goalPartition(NodeId n) const {
    const int k = g_.node(n).components;
    return cfg_.goalDisjoint ? Partition::singletons(k) : Partition::coarsest(k);
  }
  static bool goalMetIn(const ProofGraph &g, NodeId t, const Partition &p, int genus) {
    auto b = g.best(t, p);
    return b && b->genus <= genus;
  }
  bool upperMet() const {
    return goalMetIn(g_, target_, goalPartition(target_), cfg_.goalGenus);
  }
  bool lowerMet() const {
    return cfg_.goalLower >= 0 &&
           g_.lower(target_, goalPartition(target_)) >= cfg_.goalLower;
  }
  bool goalMet() const { return upperMet() || lowerMet(); }

  void onNewNode(NodeId n, int depth) { onNewNodes({n}, depth); }
  /// Names every node in `ns` (in parallel), then records each at `depth`
  /// with its table name, literature leaf and lower bound (in order):
  /// NodeAxioms::name().
  void onNewNodes(const std::vector<NodeId> &ns, int depth) { axioms_.name(ns, depth); }
  cobordismgraph::NameTable names_; ///< table names and symmetry types
  std::vector<NodeId> nodesSince(size_t first) const {
    std::vector<NodeId> ns;
    for (size_t m = first; m < g_.nodeCount(); ++m) ns.push_back(static_cast<NodeId>(m));
    return ns;
  }
  bool useful(NodeId n) const;
  /// The lower gate (README.md, "Lower-bound mode"): whether n, given the
  /// best lower bounds it could ever have, would carry the lower goal to
  /// the target over the edges found so far. `slack` gets how much more
  /// than the goal it would carry (charge still affordable).
  bool usefulLower(NodeId n, int *slack = nullptr) const;
  /// What n could carry to the target at best, cached per graph version;
  /// computed for every node of `ns` on the run's threads.
  void lowerSlacks(const std::vector<NodeId> &ns) const;
  std::optional<int> lowerSlack(NodeId n) const;
  /// --lower-sources: which table names are special sources (their lower
  /// bound is not a Lipschitz invariant's), and the largest such bound.
  void loadLowerSources();
  /// Prints the proof of lower(n, q): its reason, then the reasons of what
  /// it read, indented, down to literature and linking leaves.
  void describeLower(std::ostream &o, NodeId n, const Partition &q, int indent) const;
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
  /// <work>/profiles.jsonl: every node's Pareto (partition, genus) entries
  /// and per-partition lower bounds (ProofGraph::profileFields()), with its
  /// name, depth, crossings and whether it was searched. Written at every
  /// exit that writes the run's other files, for the atlas page.
  void writeProfiles() const;
  void writeCertificate() const;
  /// <work>/lower_certificate.json when the lower goal is met: the proof
  /// of lower(target, goal) as a tree of facts (README.md, "Lower-bound
  /// mode"), for tools/cascade_check.py --lower.
  void writeLowerCertificate() const;
  /// <work>/node_bounds.jsonl at a run's end: every node's identity,
  /// diagram, proved profile entries and lower bounds (README.md,
  /// "Lower-bound mode"), for the atlas's table beyond the tables.
  void writeNodeBounds() const;
  void writeWitnessEdge(std::ostream &c, EdgeId e, std::set<NodeId> &nodes) const;
  void writeRecords(std::ostream &c, const std::vector<RecordId> &ids,
                    std::set<NodeId> &nodes) const;
  void writeNodes(std::ostream &c, const std::set<NodeId> &nodes) const;
  void log(const std::string &line) {
    std::ofstream(cfg_.work + "/cascade.jsonl", std::ios::app) << line << "\n";
  }

  Config cfg_;
  ProofGraph g_;
  NodeRegistry reg_;
  exactnaming::ExactTables tables_;
  exactnaming::ExactNamer namer_;
  /// Each link's name and outside facts (bounds/axioms.h), shared with a
  /// depth-0 search's own graph; the target, its class, each node's depth
  /// and table name are its.
  NodeAxioms axioms_;
  NodeId &target_ = axioms_.target;
  std::string &targetCanonical_ = axioms_.targetClass;
  std::map<NodeId, int> &depth_ = axioms_.depth;
  std::map<NodeId, std::string> &tableName_ = axioms_.tableName;
  std::map<NodeId, std::vector<long>> expansions_;
  std::set<NodeId> refused_;
  /// Why the run must halt (an impossible state, divergence 2); empty if not.
  std::string halt_;
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
  std::unique_ptr<witnessstore::DatabaseIndex> master_;
  std::map<std::string, std::string> tablePD_;
  std::set<NodeId> masterDone_;
  bool masterRowsFor(NodeId n, std::vector<std::string> *rows = nullptr) const;
  void loadMaster(NodeId n, bool countsAsExpansion = true);
  /// The class a table name stands for, as the namer names nodes: its base's
  /// variants that are one oriented link up to mirror and global reversal,
  /// by diagram OR by a meridian-carrying isometry (ExactNamer::
  /// canonicalName(), data/table_link_classes.csv). ExactTables::canonical()
  /// joins by diagram only, so it must never be compared with a node's name.
  std::string classOf(const std::string &name) const { return axioms_.classOf(name); }
  std::map<EdgeId, EdgeInfo> edgeInfo_;
  std::map<std::string, EdgeInfo> directInfo_; // by witness key
  double cpuSpent_ = 0, wallSpent_ = 0;
  int hops_ = 0;
  int invariantFailures_ = 0;
  std::map<NodeId, std::string> hopSubject_; ///< each searched node's subject name
  bool stored_ = false;
  size_t storedAppended_ = 0;             ///< witnesses the store gained
  std::set<NodeId> boosted_;              ///< hubs already expanded wide (--hub-degree)
  std::string stopReason_ = "nothing-useful"; ///< why the loop ended short of the goal
  std::map<std::string, bool> special_;   ///< --lower-sources: name -> special
  int lowerLMax_ = 0;                     ///< the largest special lower bound
  /// Lower what-ifs by node, valid for one (witnesses, records, lower
  /// version) of the graph: nullopt when the what-if was inconsistent.
  struct LowerCache {
    std::tuple<size_t, size_t, long> version;
    std::map<NodeId, std::optional<int>> carried;
  };
  mutable LowerCache lowerCache_;
  /// Each searched node's latest usable frontier (searchfrontier.h): its next
  /// hop carries on from there instead of searching the prefix again.
  std::map<NodeId, SearchFrontier> frontiers_;
  /// Nodes whose search ran to the end at the hop shape: nothing is left.
  std::set<NodeId> searchedOut_;
  /// The driver's own time, outside what a hop record's timers cover
  /// (README.md, "Where a run's time goes"): wall seconds, summed over the
  /// run. Each hop record carries what accrued since the previous one.
  struct DriverTimes {
    double choose = 0;     ///< choose(), whole
    double useful = 0;     ///< of which the upper gate's what-ifs
    double lowerSlack = 0; ///< of which the lower gate's what-ifs
    double master = 0;     ///< loadMaster(), whole
    double kept = 0;       ///< kept.csv (fsynced) and frontier.txt per hop
  };
  DriverTimes driver_, driverAtLastHop_;
  StoreResult storeResult_; ///< storeWitnesses()'s counts and times
  /// Writes the run's own record to cascade.jsonl: its whole wall and CPU,
  /// and where the time outside the hops went.
  void logRun(double wall, double cpu, double startup, double loop, double store,
              double lowerReport, double nodeBounds);
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
  const StoreResult s = storeKept(readKept(cfg_.work), cfg_.witnessStore, cfg_.dedupeAgainst,
                                  names_, static_cast<unsigned>(cfg_.threads),
                                  cfg_.pairSigCache);
  storedAppended_ = s.appended;
  storeResult_ = s;
  std::cout << "[+] witness store: " << s.kept << " kept, " << s.fresh << " new, "
            << s.appended << " appended to " << cfg_.witnessStore << " (signed in "
            << std::fixed << std::setprecision(0) << s.signSeconds << " s)\n";
}

bool Cascade::masterRowsFor(NodeId n, std::vector<std::string> *rows) const {
  if (!master_) return false;
  auto it = tableName_.find(n);
  if (it == tableName_.end()) return false;
  // Every table entry of this link's class (one oriented link up to mirror
  // and global reversal) is the same node; each is a row of its own.
  const std::string canon = classOf(it->second);
  const exactnaming::TableEntry *e = tables_.entry(it->second);
  if (!e) return false;
  bool any = !master_->byOutgoing(witnessstore::DatabaseIndex::base(it->second), 1).empty();
  for (const exactnaming::TableEntry *v : tables_.variants(e->base))
    if (classOf(v->name) == canon && master_->has(v->name)) {
      any = true;
      if (rows) rows->push_back(v->name);
    }
  return any;
}

void Cascade::loadMaster(NodeId n, bool countsAsExpansion) {
  const auto tLoad = Clock::now();
  masterDone_.insert(n);
  std::vector<std::string> own;
  if (!masterRowsFor(n, &own)) {
    driver_.master += secondsSince(tLoad);
    return;
  }
  double readSeconds = 0; // phase A: the assembly's waits for the readers
  const size_t nodesBefore = g_.nodeCount();
  int assembled = 0, failed = 0, refusedRows = 0;
  // A witness can be on an upper proof only if its genus is at most the
  // goal (glue() never lowers a genus), and on a lower proof only if the
  // charge it costs, at least its genus, is affordable: at most the largest
  // special source's bound minus the lower goal. Anything above both is
  // skipped before it is read back (the genus is in the CSV line).
  const int maxGenus = std::max(cfg_.goalGenus,
                                cfg_.goalLower >= 0 ? lowerLMax_ - cfg_.goalLower : -1);
  size_t skippedGenus = 0;
  auto keep = [&](const Witness &w) {
    if (w.genus <= maxGenus) return true;
    ++skippedGenus;
    return false;
  };
  // Witnesses by subject row: this node's own rows (every witness), and
  // other rows whose recorded far side names this node's base (a hint: the
  // far side is redrawn and identified exactly like any other).
  std::map<std::string, std::vector<Witness>> bySubject;
  std::set<std::string> ownRows(own.begin(), own.end());
  for (const std::string &name : own)
    for (const witnessstore::StoredCobordism &s : master_->rows(name))
      if (const Witness w = carried(s); keep(w)) bySubject[s.witness.subject].push_back(w);
  size_t reverse = 0;
  for (const witnessstore::StoredCobordism &s :
       master_->byOutgoing(witnessstore::DatabaseIndex::base(tableName_[n]), 300))
    if (const Witness w = carried(s); !ownRows.count(s.witness.subject) && keep(w)) {
      bySubject[s.witness.subject].push_back(w);
      ++reverse;
    }
  // Two phases. First every subject row's witnesses are read back, rows in
  // parallel: a read-back needs only the row's redrawer (no graph, no
  // registry), it was most of a run's single-threaded time, and reading one
  // enumerates isomorphisms onto the row's thickening -- kept across runs in
  // --read-back-cache. Then, serially and in the same order as before, each
  // row is interned, certified and its read-backs added: the graph is
  // exactly what the one-phase loop built.
  struct RowRead {
    std::string name, pd;
    int layers = 2;
    std::vector<const Witness *> ws;
    std::unique_ptr<farside::WitnessRedrawer> redraw;
    std::string buildError; // the redrawer could not be built
    std::vector<std::optional<farside::OutgoingLink>> links;
    std::vector<std::string> why;
    std::vector<char> invariant; // the read-back broke an invariant (logic_error)
    size_t cacheHits = 0;
  };
  std::vector<RowRead> reads;
  for (const auto &[name, ws] : bySubject) {
    const exactnaming::TableEntry *e = tables_.entry(name);
    auto pdIt = tablePD_.find(name);
    if (!e || pdIt == tablePD_.end()) continue;
    // One read per (row PD, layers): the table's PD unless the sidecar
    // recorded the row the witness was really searched on.
    std::map<std::pair<std::string, int>, std::vector<const Witness *>> byRow;
    for (const Witness &w : ws)
      byRow[{w.rowPD.empty() ? pdIt->second : w.rowPD, w.layers}].push_back(&w);
    for (const auto &[key, group] : byRow) {
      RowRead r;
      r.name = name;
      r.pd = key.first;
      r.layers = key.second;
      r.ws = group;
      reads.push_back(std::move(r));
    }
  }
  // Read by a pool of `threads` readers that lives for the whole load, in
  // row order, at most `window` rows ahead of the assembly below, which takes
  // each row as soon as it is read: so at most a window of rows' thickenings
  // is held at once (reading every row first held them all: 5.5 GB on
  // L10a174), and no reader waits at a batch's end for the batch's slowest
  // row (2026-10-02: a fresh pool per batch of 2 x threads rows, joined before
  // the batch was assembled, cost L10n112's cold load 6% of its wall).
  auto readRow = [&](size_t i) {
    RowRead &r = reads[i];
    try {
      r.redraw = std::make_unique<farside::WitnessRedrawer>(r.pd, r.layers);
    } catch (const std::exception &ex) {
      r.buildError = ex.what();
      return;
    }
    RowReadBacks cache(cfg_.readBackCache, r.pd, r.layers,
                       cfg_.readBackCache.empty() ? std::string()
                                                  : r.redraw->buildChecksum());
    r.links.resize(r.ws.size());
    r.why.resize(r.ws.size());
    r.invariant.assign(r.ws.size(), 0);
    for (size_t k = 0; k < r.ws.size(); ++k) {
      const std::string key = witnesskey::witnessKey(r.ws[k]->pairsig);
      if (const CachedReadBack *c = cache.get(key)) {
        r.links[k] = c->link;
        r.why[k] = c->why;
        continue;
      }
      std::string why;
      try {
        r.links[k] = r.redraw->outgoingLinkFast(r.ws[k]->pairsig, why);
        r.why[k] = why;
        cache.put(key, {r.links[k], why});
      } catch (const std::logic_error &ex) {
        r.invariant[k] = 1; // reported, never cached
        r.why[k] = ex.what();
      } catch (const std::exception &ex) {
        r.why[k] = ex.what(); // not cached either: it may be transient
      }
    }
    r.cacheHits = cache.hits();
    try {
      cache.flush();
    } catch (const std::exception &ex) {
      std::cerr << "[!] read-back cache not written for " << r.name << ": " << ex.what()
                << "\n";
    }
  };
  const size_t window = 2 * static_cast<size_t>(std::max(cfg_.threads, 1));
  std::mutex readMutex;
  std::condition_variable readCv;
  std::vector<char> readDone(reads.size(), 0);
  size_t assembling = 0;      // rows before this one are assembled (under readMutex)
  bool stopReading = false;   // the assembly left early (under readMutex)
  size_t nextRead = 0;        // the next row to read (under readMutex)
  auto reader = [&] {
    for (;;) {
      // A reader takes a row only once it is inside the window, so no
      // reader sits on a row it may not read yet while the rows before it
      // wait for a thread.
      size_t i;
      {
        std::unique_lock<std::mutex> lock(readMutex);
        readCv.wait(lock, [&] {
          return stopReading || nextRead >= reads.size() || nextRead < assembling + window;
        });
        if (stopReading || nextRead >= reads.size()) return;
        i = nextRead++;
      }
      readRow(i);
      {
        std::lock_guard<std::mutex> lock(readMutex);
        readDone[i] = 1;
      }
      readCv.notify_all();
    }
  };
  // Joined on every way out of the assembly, an exception included: the
  // readers stop at once, finishing only the rows they hold.
  struct ReaderPool {
    std::vector<std::thread> threads;
    std::function<void()> release;
    ~ReaderPool() {
      release();
      for (std::thread &t : threads) t.join();
    }
  } pool;
  pool.release = [&] {
    {
      std::lock_guard<std::mutex> lock(readMutex);
      stopReading = true;
    }
    readCv.notify_all();
  };
  for (size_t t = 0, k = std::min<size_t>(std::max(cfg_.threads, 1), reads.size()); t < k; ++t)
    pool.threads.emplace_back(reader);
  // Row i, read: what the assembly waited for it is read-back time.
  auto awaitRow = [&](size_t i) {
    const auto tWait = Clock::now();
    {
      std::unique_lock<std::mutex> lock(readMutex);
      assembling = i;
      readCv.notify_all();
      readCv.wait(lock, [&] { return readDone[i] != 0; });
    }
    readSeconds += secondsSince(tWait);
  };
  const auto tRows = Clock::now();
  size_t cacheHits = 0;
  // A row is interned once per PD it was searched on: the table's own
  // diagram (simplified, then the registry's exact tests), or the recorded
  // diagram of a cascade hop, which is a node's diagram of another run:
  // already reduced, so interned as it is. An own row must be THIS node.
  struct RowNode {
    NodeMatch match;
    GaussDiagram diagram; ///< the row's diagram, in the row PD's component order
  };
  std::map<std::string, RowNode> interned; // by row PD
  for (size_t i = 0; i < reads.size(); ++i) {
    awaitRow(i);
    RowRead &r = reads[i];
    cacheHits += r.cacheHits;
    auto seen = interned.find(r.pd);
    if (seen == interned.end()) {
      RowNode rn;
      if (r.pd == tablePD_.at(r.name)) {
        rn.diagram = of(tables_.entry(r.name)->diagram);
        rn.match = reg_.intern(simplifyKeepingComponents(rn.diagram), "row " + r.name);
      } else {
        rn.diagram = of(exactnaming::linkFromTablePD(r.pd));
        rn.match = reg_.intern(rn.diagram, "row " + r.name + " (recorded diagram)");
      }
      seen = interned.emplace(r.pd, std::move(rn)).first;
      if (ownRows.count(r.name) && seen->second.match.node != n) {
        ++refusedRows;
        std::cout << "[!] master row " << r.name << " did not intern as node " << n
                  << " (got " << seen->second.match.node << "); not used\n";
      }
    }
    const NodeMatch &m = seen->second.match;
    if (ownRows.count(r.name) && m.node != n) continue;
    HopRow row;
    row.node = m.node;
    row.diagram = seen->second.diagram;
    row.nodeMap = m.componentMap;
    row.pd = r.pd;
    row.layers = r.layers;
    std::unique_ptr<HopAssembler> hop;
    try {
      if (!r.buildError.empty()) throw std::runtime_error(r.buildError);
      hop = std::make_unique<HopAssembler>(g_, reg_, row, std::move(r.redraw));
    } catch (const std::exception &ex) {
      ++refusedRows;
      std::cout << "[!] master row " << r.name << " refused: " << ex.what() << "\n";
      continue;
    }
    for (size_t k = 0; k < r.ws.size(); ++k) {
      const Witness *w = r.ws[k];
      const std::string key = witnesskey::witnessKey(w->pairsig);
      HopEdge he;
      if (r.invariant[k]) {
        ++invariantFailures_;
        std::cout << "[!!] master witness " << key << ": INVARIANT: " << r.why[k] << "\n";
      } else if (!r.links[k]) {
        he.why = r.why[k];
      } else {
        try {
          he = hop->addRead(*r.links[k], w->genus, "master:" + key);
        } catch (const std::logic_error &ex) {
          ++invariantFailures_;
          std::cout << "[!!] master witness " << key << ": INVARIANT: " << ex.what() << "\n";
        } catch (const std::exception &ex) {
          he.why = ex.what();
        }
      }
      if (!he.ok) { ++failed; continue; }
      ++assembled;
      EdgeInfo info{"master", row.pd, key, he, r.layers, w->pairsig, row.nodeMap};
      if (he.direct) directInfo_["master:" + key] = info;
      else edgeInfo_[he.edge] = info;
    }
  }
  // The target's own rows stand in for its first hop at this budget level;
  // a node loaded lazily is expanded right after, so its load counts nothing.
  if (countsAsExpansion) expansions_[n].push_back(0);
  // Phase B, the serial assembly: the row loop less its read-backs.
  const double assembleSeconds = secondsSince(tRows) - readSeconds;
  const auto tName = Clock::now();
  onNewNodes(nodesSince(nodesBefore), depth_.at(n) + 1);
  const double nameSeconds = secondsSince(tName);
  const auto tProp = Clock::now();
  g_.propagate();
  g_.propagateLower();
  const double propagateSeconds = secondsSince(tProp);
  const double loadSeconds = secondsSince(tLoad);
  driver_.master += loadSeconds;
  auto best = g_.best(target_, goalPartition(target_));
  std::ostringstream o;
  o << std::fixed << std::setprecision(1);
  o << "{\"master\":\"" << json::escape(tableName_[n]) << "\",\"node\":" << n
    << ",\"wall\":" << loadSeconds << ",\"readback_s\":" << readSeconds
    << ",\"assemble_s\":" << assembleSeconds << ",\"name_s\":" << nameSeconds
    << ",\"propagate_s\":" << propagateSeconds
    << ",\"own_rows\":" << own.size() << ",\"reverse_witnesses\":" << reverse
    << ",\"subject_rows\":" << bySubject.size() << ",\"refused_rows\":" << refusedRows
    << ",\"skipped_genus\":" << skippedGenus
    << ",\"assembled\":" << assembled << ",\"failed\":" << failed
    << ",\"read_back_cache_hits\":" << cacheHits
    << ",\"nodes\":" << g_.nodeCount() << ",\"target_best\":"
    << (best ? std::to_string(best->genus) : "null") << "}";
  log(o.str());
  std::cout << "[+] master rows of node " << n << " (" << tableName_[n] << "): "
            << assembled << " witnesses assembled, " << failed << " failed, " << skippedGenus
            << " skipped by genus; target best " << (best ? std::to_string(best->genus) : "none")
            << "\n";
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

void Cascade::loadLowerSources() {
  special_.clear();
  lowerLMax_ = 0;
  if (cfg_.lowerSources.empty()) return;
  std::ifstream in(cfg_.lowerSources);
  if (!in) throw std::runtime_error("cannot read --lower-sources " + cfg_.lowerSources);
  std::string line;
  std::getline(in, line);
  std::vector<std::string> head = parseCsvLine(line);
  const auto col = [&](const std::string &c) {
    return static_cast<size_t>(std::find(head.begin(), head.end(), c) - head.begin());
  };
  const size_t nameCol = col("name"), specialCol = col("special"), loCol = col("lit_lo");
  while (std::getline(in, line)) {
    const std::vector<std::string> f = parseCsvLine(line);
    if (nameCol >= f.size() || specialCol >= f.size()) continue;
    const bool sp = f[specialCol] == "1";
    special_[f[nameCol]] = sp;
    if (sp && loCol < f.size())
      if (auto lo = exactnaming::parseTableG4(f[loCol])) lowerLMax_ = std::max(lowerLMax_, lo->first);
  }
}

void Cascade::lowerSlacks(const std::vector<NodeId> &ns) const {
  // One what-if per node (ProofGraph::lowerIf copies the graph and relaxes
  // it, so they are independent), on the run's threads, cached until the
  // graph changes.
  const std::tuple<size_t, size_t, long> version{g_.witnessCount(), g_.recordCount(),
                                                 g_.lowerVersion()};
  if (lowerCache_.version != version) {
    lowerCache_.version = version;
    lowerCache_.carried.clear();
  }
  std::vector<NodeId> todo;
  for (NodeId n : ns)
    if (!lowerCache_.carried.count(n)) todo.push_back(n);
  if (todo.empty()) return;
  const Partition goal = goalPartition(target_);
  std::vector<std::optional<int>> out(todo.size());
  parallelFor(todo.size(), static_cast<unsigned>(std::max(cfg_.threads, 1)), [&](size_t i) {
    const NodeId n = todo[i];
    const Node &node = g_.node(n);
    // The most n could ever have, per partition: any transported bound
    // is a literature seed minus charges, so at most the largest special
    // source's (lowerLMax_); a proved surface refining P caps P; the
    // literature upper bound is a connected surface, capping the
    // coarsest partition only. Below what is known already a seed is a
    // no-op, and above a proved surface it is refused (nullopt).
    int litHi = std::numeric_limits<int>::max();
    if (auto t = tableName_.find(n); t != tableName_.end())
      if (const exactnaming::TableEntry *e = tables_.entry(t->second))
        if (auto g4 = exactnaming::parseTableG4(e->g4)) litHi = g4->second;
    std::vector<ProofGraph::LowerSeed> seeds;
    for (const Partition &p : allPartitions(node.components)) {
      int cap = lowerLMax_;
      if (p.blocks() == 1) cap = std::min(cap, litHi);
      if (auto b = g_.best(n, p)) cap = std::min(cap, b->genus);
      seeds.push_back({n, p, cap});
    }
    out[i] = g_.lowerIf(seeds, target_, goal);
  });
  for (size_t i = 0; i < todo.size(); ++i) lowerCache_.carried[todo[i]] = out[i];
}

std::optional<int> Cascade::lowerSlack(NodeId n) const {
  lowerSlacks({n});
  const auto &c = lowerCache_.carried.at(n);
  if (!c) return std::nullopt;
  return *c - cfg_.goalLower;
}

bool Cascade::usefulLower(NodeId n, int *slack) const {
  if (cfg_.goalLower < 0 || n == target_) return false;
  if (g_.node(n).components > ProofGraph::kMaxLowerComponents) return false;
  if (reg_.info(n).diagram.crossings() > cfg_.lowerMaxCrossings) return false;
  auto s = lowerSlack(n);
  if (!s || *s < 0) return false;
  if (slack) *slack = *s;
  return true;
}

void Cascade::describeLower(std::ostream &o, NodeId n, const Partition &q,
                            int indent) const {
  using Kind = ProofGraph::LowerReason::Kind;
  const auto fact = g_.lowerWhy(n, q);
  const std::string pad(static_cast<size_t>(2 * indent + 4), ' ');
  auto name = [&](NodeId m) {
    auto t = tableName_.find(m);
    return "node " + std::to_string(m) +
           (t == tableName_.end() ? std::string(" (untabulated)") : " (" + t->second + ")");
  };
  o << pad << name(n) << " " << q.str() << " >= "
    << (fact.value >= ProofGraph::kNoSurface ? std::string("(no such surface)")
                                             : std::to_string(fact.value));
  if (!(fact.storedFor == q)) o << ", stored for " << fact.storedFor.str();
  if (indent > 40) {
    o << " ...\n";
    return;
  }
  switch (fact.reason.kind) {
  case Kind::none:
    o << ": nothing known\n";
    return;
  case Kind::literature:
    o << ": " << g_.node(n).lowerBoundSource << "\n";
    return;
  case Kind::linking:
    o << ": the linking numbers forbid this partition\n";
    return;
  case Kind::seed:
    o << ": what-if seed\n";
    return;
  case Kind::sumPiece: {
    const SumEdge &s = g_.sum(fact.reason.edge);
    const NodeId piece = s.pieces[static_cast<size_t>(fact.reason.piece)];
    o << ": a sum along components: summand " << name(piece) << " >= " << fact.reason.from
      << ", minus the other summands' connected genera plus components minus one ("
      << fact.reason.addition << "; records";
    for (RecordId r : fact.reason.records) o << " " << r;
    o << ")\n";
    describeLower(o, piece, Partition::coarsest(g_.node(piece).components), indent + 1);
    return;
  }
  case Kind::witness: {
    const WitnessEdge &e = g_.witness(fact.reason.edge);
    const NodeId other = fact.reason.toIsIn ? e.out : e.in;
    const Partition op = Partition::fromLabels(fact.reason.fromPartition);
    o << ": across witness " << e.key << " (genus " << e.shape.genus << ", "
      << e.shape.components << " pieces) from " << name(other) << " " << op.str() << " >= "
      << fact.reason.from << ", the cap adding " << fact.reason.addition << "\n";
    describeLower(o, other, op, indent + 1);
    return;
  }
  case Kind::splitWhole: {
    o << ": the split link's pieces, summed\n";
    for (EdgeId sid : g_.node(n).splitEdges) {
      const SplitEdge &s = g_.split(sid);
      if (s.whole != n) continue;
      for (size_t k = 0; k < s.pieces.size() && k < fact.reason.pieces.size(); ++k)
        describeLower(o, s.pieces[k], Partition::fromLabels(fact.reason.pieces[k]),
                      indent + 1);
      return;
    }
    return;
  }
  case Kind::splitPiece: {
    o << ": a piece of a split link whose whole is bounded, minus the proved genera ("
      << fact.reason.addition << ") of the other pieces (records";
    for (RecordId r : fact.reason.records) o << " " << r;
    o << ")\n";
    for (EdgeId sid : g_.node(n).splitEdges) {
      const SplitEdge &s = g_.split(sid);
      if (s.whole == n) continue;
      if (std::find(s.pieces.begin(), s.pieces.end(), n) == s.pieces.end()) continue;
      describeLower(o, s.whole, Partition::fromLabels(fact.reason.fromPartition), indent + 1);
      return;
    }
    return;
  }
  }
}

// freeOnly: only nodes whose master rows are not yet loaded (no search).
std::optional<NodeId> Cascade::choose(bool freeOnly) {
  const auto tChoose = Clock::now();
  struct ChooseTimer {
    DriverTimes &d;
    Clock::time_point t;
    ~ChooseTimer() { d.choose += secondsSince(t); }
  } chooseTimer{driver_, tChoose};
  struct Cand {
    NodeId n;
    bool lowerOnly = false; ///< kept by the lower gate alone
    int slack = 0;          ///< lower gate: charge still affordable
    std::tuple<int, int, int, size_t> key;
  };
  std::vector<Cand> cands;
  std::vector<NodeId> eligible;
  for (const auto &[n, dep] : depth_) {
    if (!reg_.known(n) || n == reg_.unknot() || refused_.count(n)) continue;
    const NodeInfo &ni = reg_.info(n);
    if (ni.diagram.crossings() > cfg_.maxCrossings) continue;
    if (!freeOnly && searchedOut_.count(n)) continue; // nothing left to search
    if (freeOnly) {
      if (masterDone_.count(n) || !masterRowsFor(n)) continue;
    } else if (!expansions_[n].empty()) {
      continue; // one expansion per budget level
    }
    eligible.push_back(n);
  }
  // The upper gate first; the lower what-ifs for every node it rejects run
  // as one parallel batch (each is a graph copy relaxed to a fixed point).
  std::map<NodeId, bool> upper;
  std::vector<NodeId> needLower;
  const auto tUseful = Clock::now();
  for (NodeId n : eligible) {
    upper[n] = n == target_ || useful(n);
    if (!upper[n] && cfg_.goalLower >= 0 && n != target_ &&
        g_.node(n).components <= ProofGraph::kMaxLowerComponents &&
        reg_.info(n).diagram.crossings() <= cfg_.lowerMaxCrossings)
      needLower.push_back(n);
  }
  driver_.useful += secondsSince(tUseful);
  const auto tLower = Clock::now();
  if (!needLower.empty()) lowerSlacks(needLower);
  driver_.lowerSlack += secondsSince(tLower);
  for (NodeId n : eligible) {
    const NodeInfo &ni = reg_.info(n);
    const int dep = depth_.at(n);
    int slack = 0;
    const bool lowerOnly = !upper[n] && usefulLower(n, &slack);
    const bool use = upper[n] || lowerOnly;
    if (cfg_.verbose && !freeOnly) {
      auto t = tableName_.find(n);
      std::cout << "    candidate " << n << " (" << ni.diagram.crossings() << "x, "
                << ni.diagram.components() << "c, "
                << (t == tableName_.end() ? "untabulated" : t->second) << ", depth " << dep
                << "): "
                << (upper[n] ? "useful" : lowerOnly ? "useful for the lower goal" : "not useful");
      if (lowerOnly) std::cout << " (slack " << slack << ")";
      std::cout << (masterRowsFor(n) && !masterDone_.count(n) ? ", master rows" : "") << "\n";
    }
    if (!use) continue;
    const int crossings = static_cast<int>(ni.diagram.crossings());
    int order = 0;
    if (cfg_.strategy == "dfs") order = -dep;
    else if (cfg_.strategy == "bfs") order = dep;
    cands.push_back({n, lowerOnly, slack, {order, crossings, dep, static_cast<size_t>(n)}});
  }
  if (cands.empty()) return std::nullopt;
  // best: fewest crossings; among equals a node the upper gate keeps before
  // one only the lower gate keeps, and among those the most slack (the
  // charge it can still spend reaches more sources); then lower volume,
  // then shallower, then the older node; dfs/bfs: by depth first. Volumes
  // are compared to 1e-6: SnapPea computes them on randomly retriangulated
  // complements, so two nodes of one volume (a link and its mirror, say)
  // differed in the last bits from run to run, and the order they were
  // chosen in -- and so the whole run -- was not reproducible (2026-09-30).
  auto roundedVolume = [&](NodeId n) {
    const NodeInfo &i = reg_.info(n);
    return i.hyperbolic ? std::llround(i.volume * 1e6) : std::numeric_limits<long long>::max();
  };
  std::sort(cands.begin(), cands.end(), [&](const Cand &a, const Cand &b) {
    if (cfg_.strategy == "best") {
      const NodeInfo &ia = reg_.info(a.n), &ib = reg_.info(b.n);
      if (ia.diagram.crossings() != ib.diagram.crossings())
        return ia.diagram.crossings() < ib.diagram.crossings();
      if (a.lowerOnly != b.lowerOnly) return !a.lowerOnly;
      if (a.lowerOnly && a.slack != b.slack) return a.slack > b.slack;
      const long long va = roundedVolume(a.n), vb = roundedVolume(b.n);
      if (va != vb) return va < vb;
      if (depth_.at(a.n) != depth_.at(b.n)) return depth_.at(a.n) < depth_.at(b.n);
      return a.n < b.n;
    }
    return a.key < b.key;
  });
  return cands.front().n;
}

void Cascade::expand(NodeId n, long surfaces) {
  // A node searched before carries on from where that search stopped (the
  // hop's surface target is its breadth, so it adds only what is new); one
  // already searched this far (a hub's wide hop, say) has nothing new at
  // this budget.
  const SearchFrontier *resume = nullptr;
  if (searcher_)
    if (auto f = frontiers_.find(n); f != frontiers_.end()) resume = &f->second;
  if (resume && resume->satisfying >= surfaces) {
    expansions_[n].push_back(surfaces);
    std::cout << "[+] node " << n << " already searched to " << resume->satisfying
              << " surfaces; nothing new at " << surfaces << "\n";
    return;
  }
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
  row.layers = cfg_.hopShape.layers;
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
        ",\"refused\":\"" + json::escape(e.what()) + "\",\"pd\":\"" + json::escape(row.pd) +
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
  std::string namingJson; // in-process hops only: the drain's naming times
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
      SearchRequest request = searcher_->hopRequest(hop->redrawer(), rowName, surfaces, 7200);
      request.resume = resume;
      if (tables_.entry(rowName)) request.censusName = cobordismgraph::baseName(rowName);
      // Every kept surface, durably, as the search runs (divergence 7): the
      // hop's pending file, signed at the run's end (storeWitnesses()).
      request.rowPD = row.pd;
      request.layers = row.layers;
      if (!cfg_.witnessStore.empty()) request.pending = dir + "/kept.csv";
      run = searcher_->run(hop->redrawer().rowBuild(), request);
    } catch (const SeedInvariantFailure &e) {
      // Divergence 2: an impossible state halts the run, once what it found
      // is written (run()).
      halt_ = "node " + std::to_string(n) + ": " + e.what();
      std::cout << "[!!] HALT: " << halt_ << "\n";
      return;
    } catch (const std::exception &e) {
      refused_.insert(n);
      log("{\"hop\":" + std::to_string(k) + ",\"node\":" + std::to_string(n) +
          ",\"refused\":\"" + json::escape(e.what()) + "\"}");
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
        << ": accounting: " << run.accounting << "\n[+] " << rowName
        << ": diagram naming: " << run.naming << "\n";
    namingJson = [&] {
      std::ostringstream o;
      o << std::fixed << std::setprecision(1) << ",\"naming_diagram_s\":"
        << run.namingDiagramSeconds << ",\"naming_fallback_s\":" << run.namingFallbackSeconds
        << ",\"naming_exact_s\":" << run.namingExactSeconds
        << ",\"naming_slowest_s\":" << run.namingSlowestSeconds;
      return o.str();
    }();
    // Every hop's accounting in the driver log too, in verifyslicegenus's
    // shape after the hop number, so a campaign audits each hop as it
    // audits a row (tools/orchestrate/audit_rows.py).
    std::cout << "[+] hop " << k << " " << rowName << ": accounting: " << run.accounting
              << "\n[+] hop " << k << " " << rowName << ": diagram naming: " << run.naming
              << "\n";
    if (run.impossible > 0) {
      // Divergence 2: a state that cannot occur halts the run, once this
      // hop's finds are recorded (below) and stored (run()), even if they
      // meet the goal.
      halt_ = "hop " + std::to_string(k) + ": surface accounting failed -- " +
              std::to_string(run.impossible) + " surfaces hit a state that cannot occur";
      std::cout << "[!!] HALT: " << halt_ << "\n";
    } else if (!run.accountingFailure.empty()) {
      std::cout << "[!!] hop " << k << ": surface accounting failed -- "
                << run.accountingFailure << " (completeness only: nothing unsound "
                << "is recorded)\n";
    }
    // The node's breadth so far, and where its next hop carries on from.
    std::cout << "[+] hop " << k << " " << rowName << ": breadth: "
              << (run.frontier ? run.frontier->summary() : std::string("not recorded"))
              << "; resumed "
              << (!resume ? std::string("none")
                  : run.resumed ? std::string("yes")
                                : "no: " + run.resumeRefusal)
              << "\n";
    const auto tKept = Clock::now();
    if (run.frontier) {
      if (run.frontier->complete) searchedOut_.insert(n);
      try {
        run.frontier->save(dir + "/frontier.txt");
      } catch (const std::exception &e) {
        std::cout << "[!] hop " << k << ": frontier not written: " << e.what() << "\n";
      }
      frontiers_[n] = std::move(*run.frontier);
    } else {
      frontiers_.erase(n); // not vouched for: the next search starts afresh
    }
    driver_.kept += secondsSince(tKept);
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
        cfg_.linkTable, "--max-crossings", "999", "--thicken-layers",
        std::to_string(shape.layers), "--collar-layers", std::to_string(shape.layers),
        "--max-faces", std::to_string(shape.maxFaces),
        "--iddfs-iterations", std::to_string(shape.iddfsIterations), "--iddfs-start",
        std::to_string(shape.iddfsStart), "--iddfs-step", std::to_string(shape.iddfsStep),
        "--root-budget-start", std::to_string(shape.rootBudgetStart),
        "--root-budget-growth", std::to_string(shape.rootBudgetGrowth), "--no-cone",
        "--harvest", "--boundary-condition", "proper", "--research-settled",
        "--per-knot-time-limit", "7200", "--threads", std::to_string(cfg_.threads),
        "--surface-target", std::to_string(surfaces),
        *shape.resolveUnlinked ? "--resolve-unlinked" : "--no-resolve-unlinked",
        "--exact-far-side-names", "--no-retriangulate-on-miss",
        "--pending-surface-cap", std::to_string(shape.pendingSurfaceCap),
        "--petal-cache-limit", std::to_string(shape.petalCacheLimit),
        "--recognition-cache-limit", std::to_string(shape.recognitionCacheLimit),
        "--boundary-signature-cache-limit",
        std::to_string(shape.boundarySignatureCacheLimit)};
    if (!cfg_.pairSigCache.empty()) {
      argv.push_back("--pair-sig-cache");
      argv.push_back(cfg_.pairSigCache);
    }
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
    << namingJson << ",\"witnesses\":" << witnesses
    << ",\"assembled\":" << assembled << ",\"failed\":" << failed
    << ",\"nodes\":" << g_.nodeCount() << ",\"new_nodes\":" << (g_.nodeCount() - nodesBefore)
    << ",\"records\":" << g_.recordCount()
    << ",\"target_best\":" << (best ? std::to_string(best->genus) : "null")
    << ",\"target_lower\":" << targetLower
    << ",\"contradictions\":" << g_.contradictions().size()
    << ",\"invariant_failures\":" << invariantFailures_
    // The driver's time since the previous hop record: choosing this node
    // (and the loads before it), and this hop's kept.csv and frontier.
    << ",\"choose_s\":" << driver_.choose - driverAtLastHop_.choose
    << ",\"useful_s\":" << driver_.useful - driverAtLastHop_.useful
    << ",\"lower_slack_s\":" << driver_.lowerSlack - driverAtLastHop_.lowerSlack
    << ",\"master_s\":" << driver_.master - driverAtLastHop_.master
    << ",\"kept_s\":" << driver_.kept - driverAtLastHop_.kept << "}";
  driverAtLastHop_ = driver_;
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
  if (!info.pairsig.empty()) c << ",\"pairsig\":\"" << json::escape(info.pairsig) << "\"";
  if (!info.faces.empty()) {
    c << ",\"faces\":[";
    for (size_t i = 0; i < info.faces.size(); ++i) c << (i ? "," : "") << info.faces[i];
    c << "],\"build\":\"" << info.build << "\"";
  }
}

namespace {
} // namespace

void Cascade::writeWitnessEdge(std::ostream &c, EdgeId eid, std::set<NodeId> &nodes) const {
  // A witness edge as a checker replays it: its key, ends, shape and maps,
  // and (for an edge with a hop or master row) the row, the surface (faces
  // and build digest, or pair signature) and each far-side piece's match.
  const WitnessEdge &we = g_.witness(eid);
  nodes.insert(we.in);
  nodes.insert(we.out);
  const auto it = edgeInfo_.find(eid);
  c << ",\"witness\":\"" << json::escape(we.key) << "\"";
  c << ",\"in\":" << we.in << ",\"out\":" << we.out
    << ",\"shape\":{\"components\":" << we.shape.components << ",\"genus\":" << we.shape.genus
    << ",\"inComponent\":" << json::array(we.shape.inComponent)
    << ",\"outComponent\":" << json::array(we.shape.outComponent) << "},\"inMap\":" << json::array(we.inMap)
    << ",\"outMap\":" << json::array(we.outMap);
  if (it == edgeInfo_.end()) return;
  const HopEdge &he = it->second.he;
  c << ",\"hop_dir\":\"" << json::escape(it->second.hopDir) << "\",\"row_pd\":\""
    << json::escape(it->second.rowPD) << "\",\"layers\":" << it->second.layers
    << ",\"row_node_map\":" << json::array(it->second.rowNodeMap);
  writeSurface(c, it->second);
  c << ",\"split_edge\":" << he.splitEdge << ",\"farCurveEdges\":[";
  for (size_t j = 0; j < he.farCurveEdges.size(); ++j)
    c << (j ? "," : "") << json::array(he.farCurveEdges[j]);
  c << "],\"pieces\":[";
  for (size_t k = 0; k < he.pieces.size(); ++k) {
    const NodeMatch &m = he.pieces[k];
    c << (k ? "," : "") << "{\"node\":" << m.node << ",\"method\":\"" << m.method
      << "\",\"componentMap\":" << json::array(m.componentMap) << ",\"mirrored\":"
      << (m.mirrored ? "true" : "false") << ",\"reversed\":" << (m.reversed ? "true" : "false")
      << ",\"origins\":" << json::array(he.pieceOrigins[k]) << "}";
    nodes.insert(m.node);
  }
  c << "]";
}

void Cascade::writeRecords(std::ostream &c, const std::vector<RecordId> &ids,
                           std::set<NodeId> &nodes) const {
  bool first = true;
  for (RecordId r : ids) {
    const Record &rec = g_.record(r);
    nodes.insert(rec.node);
    c << (first ? "" : ",\n") << "{\"id\":" << r << ",\"node\":" << rec.node
      << ",\"partition\":\"" << rec.partition.str() << "\",\"genus\":" << rec.genus
      << ",\"kind\":\"" << kindName(rec.kind) << "\",\"edge\":" << rec.edge
      << ",\"source\":\"" << json::escape(rec.source) << "\",\"children\":[";
    for (size_t i = 0; i < rec.children.size(); ++i)
      c << (i ? "," : "") << rec.children[i];
    c << "]";
    if (rec.kind == RecordKind::witnessForward || rec.kind == RecordKind::witnessReverse)
      writeWitnessEdge(c, rec.edge, nodes);
    if (rec.kind == RecordKind::splitCombine || rec.kind == RecordKind::splitRestrict) {
      const SplitEdge &se = g_.split(rec.edge);
      nodes.insert(se.whole);
      nodes.insert(se.pieces.begin(), se.pieces.end());
      c << ",\"whole\":" << se.whole << ",\"pieces\":" << json::array(se.pieces) << ",\"pieceMap\":[";
      for (size_t k = 0; k < se.pieceMap.size(); ++k)
        c << (k ? "," : "") << json::array(se.pieceMap[k]);
      c << "]";
    }
    if (rec.kind == RecordKind::sumCombine) {
      const SumEdge &se = g_.sum(rec.edge);
      nodes.insert(se.whole);
      nodes.insert(se.pieces.begin(), se.pieces.end());
      c << ",\"whole\":" << se.whole << ",\"pieces\":" << json::array(se.pieces) << ",\"pieceMap\":[";
      for (size_t k = 0; k < se.pieceMap.size(); ++k)
        c << (k ? "," : "") << json::array(se.pieceMap[k]);
      c << "]";
    }
    if (rec.kind == RecordKind::leaf && rec.source.rfind("direct witness ", 0) == 0) {
      const std::string key = rec.source.substr(15);
      if (auto it = directInfo_.find(key); it != directInfo_.end()) {
        c << ",\"witness\":\"" << json::escape(it->second.key)
          << "\",\"hop_dir\":\"" << json::escape(it->second.hopDir) << "\",\"row_pd\":\""
          << json::escape(it->second.rowPD) << "\",\"layers\":" << it->second.layers;
        writeSurface(c, it->second);
      }
    }
    c << "}";
    first = false;
  }
}

void Cascade::writeLowerCertificate() const {
  // The proof of lower(target, goal) as a tree of facts, children before
  // parents (README.md, "Lower-bound mode"): each fact is a node, a
  // partition, the value the bound holds there, and its reason. A witness
  // fact carries the edge exactly as an upper record does, so the checker
  // replays the surface the same way, then recomputes the cap's addition
  // and the other end's partition itself. Split facts carry the upper
  // records they subtract, with those records' own proofs.
  using Kind = ProofGraph::LowerReason::Kind;
  struct Fact {
    NodeId node;
    Partition q;
    ProofGraph::LowerFact f;
    long from = -1;          ///< witness / split-piece: the fact read
    std::vector<long> pieces; ///< split-whole: the pieces' facts
  };
  std::vector<Fact> facts;
  std::map<std::pair<NodeId, std::vector<int>>, long> ids;
  std::set<std::pair<NodeId, std::vector<int>>> onStack;
  std::set<RecordId> records;
  std::function<long(NodeId, const Partition &)> visit = [&](NodeId n, const Partition &q) -> long {
    const ProofGraph::LowerFact f = g_.lowerWhy(n, q);
    const auto key = std::make_pair(n, f.storedFor.labels());
    if (auto it = ids.find(key); it != ids.end()) return it->second;
    if (!onStack.insert(key).second)
      throw std::logic_error("a lower bound's reasons cycle: node " + std::to_string(n));
    Fact fact{n, f.storedFor, f};
    switch (f.reason.kind) {
    case Kind::witness: {
      const WitnessEdge &e = g_.witness(f.reason.edge);
      fact.from = visit(f.reason.toIsIn ? e.out : e.in,
                        Partition::fromLabels(f.reason.fromPartition));
      break;
    }
    case Kind::splitWhole:
      for (EdgeId sid : g_.node(n).splitEdges) {
        const SplitEdge &s = g_.split(sid);
        if (s.whole != n) continue;
        for (size_t k = 0; k < s.pieces.size() && k < f.reason.pieces.size(); ++k)
          fact.pieces.push_back(visit(s.pieces[k], Partition::fromLabels(f.reason.pieces[k])));
        break;
      }
      break;
    case Kind::splitPiece:
      for (EdgeId sid : g_.node(n).splitEdges) {
        const SplitEdge &s = g_.split(sid);
        if (s.whole == n || std::find(s.pieces.begin(), s.pieces.end(), n) == s.pieces.end())
          continue;
        fact.from = visit(s.whole, Partition::fromLabels(f.reason.fromPartition));
        break;
      }
      for (RecordId r : f.reason.records)
        for (RecordId p : g_.proof(r)) records.insert(p);
      break;
    case Kind::sumPiece: {
      const SumEdge &s = g_.sum(f.reason.edge);
      const NodeId piece = s.pieces[static_cast<size_t>(f.reason.piece)];
      fact.from = visit(piece, Partition::coarsest(g_.node(piece).components));
      for (RecordId r : f.reason.records)
        for (RecordId p : g_.proof(r)) records.insert(p);
      break;
    }
    default:
      break;
    }
    onStack.erase(key);
    facts.push_back(std::move(fact));
    const long id = static_cast<long>(facts.size()) - 1;
    ids[key] = id;
    return id;
  };
  const Partition goal = goalPartition(target_);
  const long top = visit(target_, goal);
  std::ofstream c(cfg_.work + "/lower_certificate.json");
  c << "{\"target\":\"" << json::escape(cfg_.targetName) << "\",\"target_pd\":\""
    << json::escape(cfg_.targetPD) << "\",\"goal_lower\":" << cfg_.goalLower << ",\"goal\":\""
    << (cfg_.goalDisjoint ? "disjoint" : "connected") << "\",\"lower\":"
    << g_.lower(target_, goal) << ",\"top\":" << top << ",\"facts\":[\n";
  std::set<NodeId> nodes;
  auto value = [](int v) {
    return v >= ProofGraph::kNoSurface ? std::string("\"inf\"") : std::to_string(v);
  };
  for (size_t i = 0; i < facts.size(); ++i) {
    const Fact &fact = facts[i];
    nodes.insert(fact.node);
    c << (i ? ",\n" : "") << "{\"id\":" << i << ",\"node\":" << fact.node << ",\"partition\":\""
      << fact.q.str() << "\",\"value\":" << value(fact.f.value);
    switch (fact.f.reason.kind) {
    case Kind::literature:
      c << ",\"kind\":\"literature\",\"source\":\"" << json::escape(g_.node(fact.node).lowerBoundSource)
        << "\"";
      break;
    case Kind::linking:
      c << ",\"kind\":\"linking\"";
      break;
    case Kind::witness: {
      const WitnessEdge &e = g_.witness(fact.f.reason.edge);
      c << ",\"kind\":\"witness\",\"to_is_in\":" << (fact.f.reason.toIsIn ? "true" : "false")
        << ",\"from\":" << fact.from << ",\"from_partition\":\""
        << Partition::fromLabels(fact.f.reason.fromPartition).str() << "\",\"from_value\":"
        << value(fact.f.reason.from) << ",\"addition\":" << fact.f.reason.addition;
      writeWitnessEdge(c, e.id, nodes);
      break;
    }
    case Kind::splitWhole: {
      c << ",\"kind\":\"split-whole\",\"pieces\":[";
      for (size_t k = 0; k < fact.pieces.size(); ++k) c << (k ? "," : "") << fact.pieces[k];
      c << "]";
      for (EdgeId sid : g_.node(fact.node).splitEdges) {
        const SplitEdge &se = g_.split(sid);
        if (se.whole != fact.node) continue;
        c << ",\"whole\":" << se.whole << ",\"piece_nodes\":" << json::array(se.pieces)
          << ",\"pieceMap\":[";
        for (size_t k = 0; k < se.pieceMap.size(); ++k) c << (k ? "," : "") << json::array(se.pieceMap[k]);
        c << "]";
        break;
      }
      break;
    }
    case Kind::splitPiece: {
      c << ",\"kind\":\"split-piece\",\"from\":" << fact.from << ",\"from_partition\":\""
        << Partition::fromLabels(fact.f.reason.fromPartition).str() << "\",\"from_value\":"
        << value(fact.f.reason.from) << ",\"subtracted\":" << fact.f.reason.addition
        << ",\"records\":[";
      for (size_t k = 0; k < fact.f.reason.records.size(); ++k)
        c << (k ? "," : "") << fact.f.reason.records[k];
      c << "]";
      break;
    }
    case Kind::sumPiece: {
      const SumEdge &se = g_.sum(fact.f.reason.edge);
      c << ",\"kind\":\"sum-piece\",\"from\":" << fact.from << ",\"piece\":" << fact.f.reason.piece
        << ",\"from_value\":" << value(fact.f.reason.from) << ",\"subtracted\":"
        << fact.f.reason.addition << ",\"records\":[";
      for (size_t k = 0; k < fact.f.reason.records.size(); ++k)
        c << (k ? "," : "") << fact.f.reason.records[k];
      c << "],\"whole\":" << se.whole << ",\"piece_nodes\":" << json::array(se.pieces) << ",\"pieceMap\":[";
      for (size_t k = 0; k < se.pieceMap.size(); ++k) c << (k ? "," : "") << json::array(se.pieceMap[k]);
      c << "]";
      for (NodeId p : se.pieces) nodes.insert(p);
      break;
    }
    default:
      c << ",\"kind\":\"none\"";
      break;
    }
    c << "}";
  }
  c << "\n],\"records\":[\n";
  writeRecords(c, std::vector<RecordId>(records.begin(), records.end()), nodes);
  c << "\n],\"nodes\":[\n";
  writeNodes(c, nodes);
  c << "\n]}\n";
}

void Cascade::writeNodeBounds() const {
  using Kind = ProofGraph::LowerReason::Kind;
  std::ofstream o(cfg_.work + "/node_bounds.jsonl");
  auto kindName = [](Kind k) {
    switch (k) {
    case Kind::literature: return "literature";
    case Kind::linking: return "linking";
    case Kind::witness: return "witness";
    case Kind::splitWhole: return "split-whole";
    case Kind::splitPiece: return "split-piece";
    case Kind::seed: return "seed";
    case Kind::sumPiece: return "sum-piece";
    default: return "none";
    }
  };
  for (size_t i = 0; i < g_.nodeCount(); ++i) {
    const NodeId n = static_cast<NodeId>(i);
    const Node &node = g_.node(n);
    o << "{\"node\":" << n << ",\"label\":\"" << json::escape(node.label) << "\",\"components\":"
      << node.components << ",\"target\":" << (n == target_ ? "true" : "false");
    if (auto t = tableName_.find(n); t != tableName_.end())
      o << ",\"table\":\"" << json::escape(t->second) << "\"";
    if (reg_.known(n)) {
      const NodeInfo &ni = reg_.info(n);
      const GaussDiagram &d = ni.diagram;
      o << ",\"crossings\":" << d.crossings() << ",\"hyperbolic\":"
        << (ni.hyperbolic ? "true" : "false");
      if (ni.hyperbolic) o << ",\"volume\":" << std::setprecision(12) << ni.volume;
      o << ",\"pd\":\"" << json::escape(rowPD(d)) << "\",\"signs\":" << json::array(d.signs) << ",\"gauss\":[";
      for (size_t c = 0; c < d.comps.size(); ++c) o << (c ? "," : "") << json::array(d.comps[c]);
      o << "],\"linking\":[";
      for (size_t a = 0; a < ni.linking.size(); ++a) o << (a ? "," : "") << json::array(ni.linking[a]);
      o << "]";
    }
    o << ",\"upper\":[";
    bool first = true;
    for (const ProfileEntry &e : node.profile.entries()) {
      bool constructive = true;
      for (RecordId r : g_.proof(e.record))
        if (g_.record(r).kind == RecordKind::leaf &&
            g_.record(r).source.rfind("literature", 0) == 0)
          constructive = false;
      o << (first ? "" : ",") << "{\"partition\":\"" << e.partition.str() << "\",\"genus\":"
        << e.genus << ",\"record\":" << e.record << ",\"constructive\":"
        << (constructive ? "true" : "false") << "}";
      first = false;
    }
    o << "],\"lower\":[";
    first = true;
    if (node.components <= ProofGraph::kMaxLowerComponents)
      for (const Partition &p : allPartitions(node.components)) {
        const auto f = g_.lowerWhy(n, p);
        if (f.value <= 0 || !(f.storedFor == p)) continue; // stored bounds only, once each
        o << (first ? "" : ",") << "{\"partition\":\"" << p.str() << "\",\"value\":"
          << (f.value >= ProofGraph::kNoSurface ? std::string("\"inf\"") : std::to_string(f.value))
          << ",\"kind\":\"" << kindName(f.reason.kind) << "\"";
        if (f.reason.kind == Kind::literature)
          o << ",\"source\":\"" << json::escape(node.lowerBoundSource) << "\"";
        o << "}";
        first = false;
      }
    o << "]}\n";
  }
}

void Cascade::writeCertificate() const {
  auto best = g_.best(target_, goalPartition(target_));
  if (!best) return;
  std::ofstream c(cfg_.work + "/certificate.json");
  c << "{\"target\":\"" << json::escape(cfg_.targetName) << "\",\"target_pd\":\""
    << json::escape(cfg_.targetPD) << "\",\"goal_genus\":" << cfg_.goalGenus
    << ",\"goal\":\"" << (cfg_.goalDisjoint ? "disjoint" : "connected") << "\",\"genus\":"
    << best->genus << ",\"records\":[\n";
  std::set<NodeId> nodes;
  writeRecords(c, g_.proof(best->record), nodes);
  c << "\n],\"nodes\":[\n";
  writeNodes(c, nodes);
  c << "\n]}\n";
}

void Cascade::writeNodes(std::ostream &c, const std::set<NodeId> &nodes) const {
  bool first = true;
  for (NodeId n : nodes) {
    c << (first ? "" : ",\n") << "{\"id\":" << n << ",\"label\":\""
      << json::escape(g_.node(n).label) << "\",\"components\":" << g_.node(n).components;
    if (reg_.known(n) && reg_.info(n).diagram.crossings() > 0) {
      // The node's own diagram, as signed Gauss data: component maps refer
      // to ITS component order, which a PD round trip need not keep.
      const GaussDiagram &d = reg_.info(n).diagram;
      c << ",\"pd\":\"" << json::escape(rowPD(d)) << "\",\"signs\":[";
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
      c << ",\"table\":\"" << json::escape(it->second) << "\"";
    c << "}";
    first = false;
  }
}

void Cascade::logRun(double wall, double cpu, double startup, double loop, double store,
                     double lowerReport, double nodeBounds) {
  std::ostringstream o;
  o << std::fixed << std::setprecision(1) << "{\"run\":\"" << json::escape(cfg_.targetName)
    << "\",\"threads\":" << cfg_.threads << ",\"wall\":" << wall << ",\"cpu\":" << cpu
    << ",\"cores\":" << (wall > 0 ? cpu / wall : 0.0) << ",\"startup_s\":" << startup
    << ",\"loop_s\":" << loop << ",\"hop_wall_s\":" << wallSpent_
    << ",\"hop_cpu_s\":" << cpuSpent_ << ",\"choose_s\":" << driver_.choose
    << ",\"useful_s\":" << driver_.useful << ",\"lower_slack_s\":" << driver_.lowerSlack
    << ",\"master_s\":" << driver_.master << ",\"kept_s\":" << driver_.kept
    << ",\"store_s\":" << store << ",\"store_dedupe_s\":" << storeResult_.dedupeSeconds
    << ",\"store_sign_s\":" << storeResult_.signSeconds
    << ",\"lower_report_s\":" << lowerReport << ",\"node_bounds_s\":" << nodeBounds
    << ",\"hops\":" << hops_ << ",\"nodes\":" << g_.nodeCount() << "}";
  log(o.str());
}

int Cascade::run() {
  const auto tRun = Clock::now();
  const double cpuRun = timers::processCpuSeconds();
  fs::create_directories(cfg_.work);
  if (cfg_.hopMode == "process") {
    // As a child hop runs verifyslicegenus: a private census copy, never
    // written, and no Pachner searches on a census miss.
    if (!census::setCensusPath(cfg_.censusDb))
      std::cout << "[!] census not found at " << cfg_.censusDb << "\n";
    census::retriangulateOnMiss.store(false);
    census::censusUpdates.store(false);
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
  loadLowerSources();
  // Table names and symmetry types, for the slice-composite anchors
  // (NodeAxioms) and the store step: as verifyslicegenus loads them.
  size_t symmetryTypes = 0;
  names_ = loadTableNames(cfg_.knotTable, cfg_.linkTable, cfg_.knotSymmetry, &symmetryTypes);
  if (!cfg_.knotSymmetry.empty())
    std::cout << "[+] knot symmetry: " << symmetryTypes << " types (slice composites beyond "
                 "3_1#m3_1 and 4_1#4_1 need them)\n";
  if (cfg_.goalLower >= 0)
    std::cout << "[+] lower goal " << cfg_.goalLower << ": " << special_.size()
              << " table names in --lower-sources, "
              << std::count_if(special_.begin(), special_.end(),
                               [](const auto &kv) { return kv.second; })
              << " special, largest special lower bound " << lowerLMax_
              << " (caps what any node could carry)\n";
  if (!cfg_.masterWitnesses.empty()) {
    const auto t0 = std::chrono::steady_clock::now();
    master_ = std::make_unique<witnessstore::DatabaseIndex>(cfg_.masterWitnesses);
    tablePD_ = tablePDs({cfg_.knotTable, cfg_.linkTable});
    std::cout << "[+] master witnesses: " << master_->subjects() << " subjects indexed in "
              << std::fixed << std::setprecision(0)
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()
              << " s (read only)\n";
  }
  GaussDiagram raw = of(exactnaming::linkFromTablePD(cfg_.targetPD));
  GaussDiagram simp = simplifyKeepingComponents(raw);
  if (exactnaming::splitPieces(simp).size() != 1)
    throw std::runtime_error("the target is a split diagram; give one piece");
  // The target's own table name, so its literature value never proves it.
  bool named = false;
  std::string composite;
  try {
    auto pn = namer_.identify(simp);
    if (pn.names.size() == 1 && pn.by != exactnaming::PieceName::By::untabulated) {
      targetCanonical_ = classOf(pn.names.front());
      named = true;
    } else if (simp.components() == 1) {
      // A composite target: its whole-diagram name (never an anchor for
      // itself: NodeAxioms skips the target), reported and recorded.
      exactnaming::FarSideName fs = namer_.name(simp.link());
      if (fs.exact && fs.pinned && fs.pieces.size() >= 2 &&
          fs.name.find('#') != std::string::npos)
        composite = fs.name;
    }
  } catch (...) {
  }
  if (tables_.entry(cfg_.targetName)) {
    // The PD must be the entry it is run as: its witnesses are recorded
    // under that name (--witness-store), so an alternative diagram, or a PD
    // copied wrongly, that the exact namer proves to be another link is
    // refused outright. A namer that cannot tell (cheap limits) is not a
    // refusal; the table's own PD needs no proof.
    const std::string claimed = classOf(cfg_.targetName);
    if (named && targetCanonical_ != claimed)
      throw std::runtime_error("the target PD is " + targetCanonical_ + ", not " +
                               cfg_.targetName + " (" + claimed + ")");
    targetCanonical_ = claimed;
  }
  NodeMatch t = reg_.intern(simp, "target " + cfg_.targetName);
  target_ = t.node;
  onNewNode(target_, 0);
  if (!composite.empty()) {
    tableName_[target_] = composite;
    std::cout << "[+] target is the composite " << composite
              << (exactnaming::isElementarySlice(composite, names_.symmetries())
                      ? " (a slice composite: its summands cancel in concordance)"
                      : "")
              << "\n";
  }
  std::cout << "[+] target " << cfg_.targetName << " = node " << target_ << " ("
            << simp.crossings() << " crossings, " << simp.components()
            << " components; table class '" << targetCanonical_ << "'); goal genus "
            << cfg_.goalGenus << (cfg_.goalDisjoint ? " (disjoint pieces)" : " (connected)")
            << "\n";
  const auto start = std::chrono::steady_clock::now();
  const double startupSeconds = secondsSince(tRun);
  long budget = cfg_.hopSurfaces;
  while (true) {
    // A halt first (divergence 2): even a goal met by the hop that found an
    // impossible state is not reported as met.
    if (!halt_.empty()) {
      // The surfaces found are real whatever broke, so they are kept.
      storeWitnesses();
      writeProfiles();
      printOutcome("halted");
      return 2;
    }
    if (!g_.contradictions().empty()) {
      for (const auto &c : g_.contradictions()) std::cout << "[!!] CONTRADICTION: " << c << "\n";
      // The surfaces are real whatever the contradiction's cause (a naming
      // or solver bug), so they are kept, as verifyslicegenus writes its
      // witnesses before its fatal-bug halt.
      storeWitnesses();
      writeProfiles();
      printOutcome("contradiction");
      return 3;
    }
    // Checked after the gates (divergence 6): a contradiction met together
    // with the goal still halts the run, never certifies through it.
    if (goalMet()) break;
    // Free edges first (no search, so no budget): the target's own master
    // rows, if it is a table entry the atlas searched.
    if (master_ && !masterDone_.count(target_) && masterRowsFor(target_)) {
      loadMaster(target_);
      continue;
    }
    // Any other node's stored rows are loaded only when that node is the
    // one about to be expanded (below): a proof can run through a node only
    // when the search picks it, and reading every table node's rows as it
    // was met cost more than the hops (2026-09-30, close1: 35-45 loads and
    // ~2,500 read-backs per row). --master-loads eager restores that.
    if (master_ && !cfg_.lazyMasterLoads)
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
    if (master_ && cfg_.lazyMasterLoads && !masterDone_.count(*n) && masterRowsFor(*n)) {
      // Its stored rows first: free edges, which may close the proof
      // without the hop, and are in any case what the hop would refind.
      loadMaster(*n, /*countsAsExpansion=*/false);
      if (goalMet() || !g_.contradictions().empty()) continue;
    }
    long surfaces = budget;
    if (cfg_.hubDegree > 0 && !boosted_.count(*n) &&
        g_.node(*n).witnessEdges.size() >= cfg_.hubDegree && cfg_.hubSurfaces > budget) {
      // A hub: many routes meet here, so one wide hop from it buys many
      // more first-level candidates than another narrow one elsewhere.
      surfaces = cfg_.hubSurfaces;
      boosted_.insert(*n);
      std::cout << "[+] hub: node " << *n << " has " << g_.node(*n).witnessEdges.size()
                << " witness edges; expanding it at " << surfaces << " surfaces\n";
    }
    expand(*n, surfaces);
  }
  const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  auto best = g_.best(target_, goalPartition(target_));
  std::cout << "[+] done: " << hops_ << " hops, " << g_.nodeCount() << " nodes, "
            << g_.recordCount() << " records; " << std::fixed << std::setprecision(0)
            << wall << " s wall, " << cpuSpent_ << " s search CPU. Target best: "
            << (best ? std::to_string(best->genus) : "none") << "\n";
  const auto tStore = Clock::now();
  storeWitnesses();
  const double storeSeconds = secondsSince(tStore);
  const auto tReport = Clock::now();
  writeLowerReport();
  const double reportSeconds = secondsSince(tReport);
  const auto tBounds = Clock::now();
  writeNodeBounds();
  writeProfiles();
  logRun(secondsSince(tRun), timers::processCpuSeconds() - cpuRun, startupSeconds, wall, storeSeconds,
         reportSeconds, secondsSince(tBounds));
  if (lowerMet() && !upperMet()) {
    // The bound's proof: the reasons from the target down to their leaves.
    const Partition goal = goalPartition(target_);
    std::cout << "[+] LOWER GOAL MET: " << cfg_.targetName << " genus >= "
              << g_.lower(target_, goal) << " (literature-assisted); the proof:\n";
    describeLower(std::cout, target_, goal, 0);
    writeLowerCertificate();
    std::cout << "[+] lower certificate " << cfg_.work << "/lower_certificate.json\n";
    printOutcome("met");
    return 0;
  }
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

void Cascade::writeProfiles() const {
  std::ofstream out(cfg_.work + "/profiles.jsonl");
  for (NodeId n = 0; n < static_cast<NodeId>(g_.nodeCount()); ++n) {
    // A crossingless node is named as the store names it (cascade_record.py):
    // it is never a hop's subject, so subjectName() has no better name.
    std::string name = subjectName(n);
    if (reg_.known(n) && reg_.info(n).diagram.signs.empty() && n != target_) {
      const int k = g_.node(n).components;
      name = identify::unlinkName(static_cast<size_t>(k));
    }
    out << "{\"node\":" << n << ",\"name\":\"" << json::escape(name) << '"';
    if (auto it = tableName_.find(n); it != tableName_.end())
      out << ",\"table\":\"" << json::escape(it->second) << '"';
    out << ",\"label\":\"" << json::escape(g_.node(n).label) << '"';
    if (auto it = depth_.find(n); it != depth_.end())
      out << ",\"depth\":" << it->second;
    // A split far side's whole is added to the graph, not the registry,
    // and has no diagram of its own.
    if (reg_.known(n))
      out << ",\"crossings\":" << reg_.info(n).diagram.signs.size();
    out << ",\"searched\":" << (hopSubject_.count(n) ? "true" : "false") << ','
        << g_.profileFields(n) << "}\n";
  }
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
  const std::map<std::string, bool> &special = special_;
  const Partition goal = goalPartition(target_);
  const int targetLower = g_.lower(target_, goal);
  int litLo = -1;
  if (const exactnaming::TableEntry *e = tables_.entry(cfg_.targetName))
    if (auto g4 = exactnaming::parseTableG4(e->g4)) litLo = g4->first;
  std::ofstream out(cfg_.work + "/lower_report.jsonl");
  out << "{\"target\":\"" << json::escape(cfg_.targetName) << "\",\"target_lower\":" << targetLower
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
  // Each what-if copies the graph and relaxes it alone, so they run on the
  // run's threads (serially they were a third of a run's driver time,
  // 2026-09-30); the lines are then written in the same order as before.
  struct Job {
    NodeId n;
    std::string name;
    int lo, hi, could;
    int carries = 0, couldCarry = 0;
  };
  std::vector<Job> jobs;
  for (const auto &[n, name] : tableName_) {
    if (n == target_) continue;
    const exactnaming::TableEntry *e = tables_.entry(name);
    if (!e) continue;
    auto g4 = exactnaming::parseTableG4(e->g4);
    if (!g4) continue;
    // The most n could be: its literature upper end, or less if a surface
    // for it is already proved.
    int could = g4->second;
    if (auto b = g_.bestConnected(n)) could = std::min(could, b->genus);
    jobs.push_back({n, name, g4->first, g4->second, could});
  }
  parallelFor(jobs.size(), static_cast<unsigned>(std::max(cfg_.threads, 1)), [&](size_t i) {
    Job &j = jobs[i];
    j.carries = j.lo > 0 ? carried(j.n, j.lo) : 0;
    j.couldCarry = j.could > 0 ? (j.could == j.lo ? j.carries : carried(j.n, j.could)) : 0;
  });
  int bestCarry = 0, bestCould = 0;
  std::string bestName, bestCouldName;
  for (const Job &j : jobs) {
    const std::string &name = j.name;
    const int carries = j.carries, couldCarry = j.couldCarry;
    if (carries <= 0 && couldCarry <= 0) continue; // n reaches the target with nothing
    const auto sp = special.find(name);
    out << "{\"node\":" << j.n << ",\"name\":\"" << json::escape(name) << "\",\"lit_lo\":"
        << j.lo << ",\"lit_hi\":" << j.hi << ",\"could\":" << j.could
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
            << " goal_genus=" << cfg_.goalGenus << " goal_lower=" << cfg_.goalLower
            << " lower_max_crossings=" << cfg_.lowerMaxCrossings
            << " master_loads=" << (cfg_.lazyMasterLoads ? "lazy" : "eager")
            << " literature=" << (cfg_.literature ? 1 : 0)
            << " hop_mode=" << cfg_.hopMode << " hop_surfaces=" << cfg_.hopSurfaces
            << " max_hop_surfaces=" << cfg_.maxHopSurfaces
            << " max_expansions=" << cfg_.maxExpansions << " cpu_budget=" << cfg_.cpuBudget
            << " strategy=" << cfg_.strategy << " max_crossings=" << cfg_.maxCrossings
            << " hub_degree=" << cfg_.hubDegree << " hub_surfaces=" << cfg_.hubSurfaces
            << " threads=" << cfg_.threads << " max_faces=" << s.maxFaces
            << " iddfs_iterations=" << s.iddfsIterations << " iddfs_start=" << s.iddfsStart
            << " iddfs_step=" << s.iddfsStep << " root_budget_start=" << s.rootBudgetStart
            << " root_budget_growth=" << s.rootBudgetGrowth << " layers=" << s.layers
            << " resolve_unlinked=" << (*s.resolveUnlinked ? 1 : 0)
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
    else if (a == "--resolve-unlinked") c.hopShape.resolveUnlinked = true;
    else if (a == "--no-resolve-unlinked") c.hopShape.resolveUnlinked = false;
    else if (a == "--verbose") c.verbose = true;
    else if (a == "--witness-store") c.witnessStore = next();
    else if (a == "--pair-sig-cache") c.pairSigCache = next();
    else if (a == "--read-back-cache") c.readBackCache = next();
    else if (a == "--run-name") c.runName = next();
    else if (a == "--dedupe-against") c.dedupeAgainst.push_back(next());
    else if (a == "--sign-only") c.signOnly = true;
    else if (a == "--hub-degree") c.hubDegree = std::stoul(next());
    else if (a == "--hub-surfaces") c.hubSurfaces = std::stol(next());
    else if (a == "--lower-report") c.lowerReport = true;
    else if (a == "--lower-sources") c.lowerSources = next();
    else if (a == "--goal-lower") c.goalLower = std::stoi(next());
    else if (a == "--master-loads") {
      const std::string v = next();
      if (v != "lazy" && v != "eager") throw std::invalid_argument("--master-loads lazy|eager");
      c.lazyMasterLoads = v == "lazy";
    }
    else if (a == "--lower-max-crossings") c.lowerMaxCrossings = std::stoul(next());
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
      const cobordismgraph::NameTable names = loadTableNames(c.knotTable, c.linkTable, "");
      const StoreResult s = storeKept(readKept(c.work), c.witnessStore, c.dedupeAgainst,
                                      names, static_cast<unsigned>(c.threads),
                                      c.pairSigCache);
      std::cout << "[+] witness store: " << s.kept << " kept, " << s.fresh << " new, "
                << s.appended << " appended to " << c.witnessStore << "\n";
      return 0;
    } catch (const std::exception &e) {
      std::cerr << "cascadesearch: " << e.what() << "\n";
      return 2;
    }
  }
  if (!c.hopShape.resolveUnlinked) {
    // Plan divergence 3: no default. It changes which surfaces a hop's
    // surface target counts, and every frontier's fingerprint.
    std::cerr << "cascadesearch: --resolve-unlinked or --no-resolve-unlinked is required "
                 "(resolve_unlinked has no default)\n";
    return 2;
  }
  if (!c.witnessStore.empty() && (c.runName.empty() || c.hopMode != "process")) {
    // A child hop writes its own witnesses (hop dir cob.csv), under its own
    // row name; only in-process hops are recorded through keptstore.h.
    std::cerr << "cascadesearch: --witness-store needs --run-name, and in-process hops\n";
    return 2;
  }
  if (c.goalLower >= 0 && c.lowerSources.empty()) {
    std::cerr << "cascadesearch: --goal-lower needs --lower-sources (the special sources' "
                 "largest bound is what makes the lower gate prune)\n";
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
