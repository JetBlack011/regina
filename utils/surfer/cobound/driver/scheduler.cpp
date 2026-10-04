//
//  scheduler.cpp
//
//  A goal run: chained searches for one target over its cobordism graph.
//  See ../README.md.
//
//  Each expansion is one search, in this process (search.h), on a node's
//  diagram with the run's search shape. Its surfaces become edges of the
//  cobordism graph (searchcobordisms.h), far sides become nodes (links.h),
//  and bounds are relaxed to a fixed point (cobordismgraph.h) after every
//  hop. The run stops when the target's goal has a proof, or the budget is
//  spent.

#include "cobound/driver/scheduler.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <link/link.h>

#include "surfer/report/csvwriter.h"
#include "cobound/driver/runrecords.h"
#include "cobound/driver/signals.h"
#include "cobound/driver/timers.h"
#include "cobound/json.h"
#include "cobound/parallelfor.h"
#include "linknaming/diagrams/simplification.h"
#include "linknaming/linknamer.h"
#include "linknaming/tables.h"
#include "cobound/outgoing/outgoingnamer.h"
#include "cobound/bounds/certificate.h"
#include "cobound/bounds/databasecobordisms.h"
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
#include "cobound/frozen.h"

namespace scheduler {

using bounds::CertificateGoal;
using bounds::CertificateWriter;
using bounds::DatabaseCobordisms;
using bounds::DatabaseLoad;
using bounds::CobordismSource;
using bounds::CobordismSources;
using bounds::CobordismAssembler;
using bounds::AddedCobordism;
using bounds::SearchedLink;
using bounds::GraphLink;
using bounds::LinkAxioms;
using bounds::LinkId;
using bounds::LinkInfo;
using bounds::LinkMatch;
using bounds::LinkRegistry;
using bounds::Partition;
using bounds::CobordismGraph;
using bounds::DerivationId;
using bounds::DerivationKind;
using bounds::allPartitions;
using bounds::diagramPD;
using cobordisms::SignResult;
using cobordisms::signPending;
using linknaming::simplifyKeepingComponents;
using search::SearchResult;
using search::Searcher;
using search::RunShape;
using search::KeptSurface;
using search::SearchRequest;
using search::SeedInvariantFailure;

using linknaming::GaussDiagram;
namespace fs = std::filesystem;

namespace {

using timers::Clock;
using timers::secondsSince;

class Scheduler {
public:
  explicit Scheduler(GoalOptions c)
      : cfg_(std::move(c)), reg_(g_),
        tables_(linknaming::Tables::load(cfg_.knotTable, cfg_.linkTable,
                                               cfg_.knotSymmetry)),
        namer_(tables_, LinkAxioms::namerLimits()),
        axioms_(g_, reg_, tables_, namer_, names_.symmetries(),
                {.literature = cfg_.literature,
                 .classes = true,
                 .threads = static_cast<unsigned>(std::max(cfg_.threads, 1)),
                 .log = &std::cout}) {}

  int run();

private:
  Partition goalPartition(LinkId n) const {
    const int k = g_.link(n).components;
    return cfg_.goalDisjoint ? Partition::singletons(k) : Partition::coarsest(k);
  }
  static bool goalMetIn(const CobordismGraph &g, LinkId t, const Partition &p, int genus) {
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

  void onNewLink(LinkId n, int depth) { onNewLinks({n}, depth); }
  /// Names every node in `ns` (in parallel), then records each at `depth`
  /// with its table name, literature leaf and lower bound (in order):
  /// NodeAxioms::name().
  void onNewLinks(const std::vector<LinkId> &ns, int depth) { axioms_.name(ns, depth); }
  solver::NameTable names_; ///< table names and symmetry types
  std::vector<LinkId> linksSince(size_t first) const {
    std::vector<LinkId> ns;
    for (size_t m = first; m < g_.linkCount(); ++m) ns.push_back(static_cast<LinkId>(m));
    return ns;
  }
  bool useful(LinkId n) const;
  /// The lower gate (README.md, "Lower-bound mode"): whether n, given the
  /// best lower bounds it could ever have, would carry the lower goal to
  /// the target over the edges found so far. `slack` gets how much more
  /// than the goal it would carry (charge still affordable).
  bool usefulLower(LinkId n, int *slack = nullptr) const;
  /// What n could carry to the target at best, cached per graph version;
  /// computed for every node of `ns` on the run's threads.
  void lowerSlacks(const std::vector<LinkId> &ns) const;
  std::optional<int> lowerSlack(LinkId n) const;
  /// --lower-sources: which table names are special sources (their lower
  /// bound is not a Lipschitz invariant's), and the largest such bound.
  void loadLowerSources();
  std::optional<LinkId> choose();
  void expand(LinkId n, long surfaces);
  /// The name node n's hop records its witnesses under: the target's own
  /// name, a proved table name, or cascade:<run>/<target>/n<n>.
  std::string subjectName(LinkId n) const;
  /// With --witness-store: signs and stores every kept surface of the run,
  /// and writes <work>/nodes.csv for the cascade: subjects. Idempotent.
  void signIntoDatabase();
  void printOutcome(const std::string &outcome) const;
  void printShape() const;
  /// The run's records' view of its graph.
  runrecords::GraphView view() const {
    return {g_, reg_, tableName_, depth_, target_, goalPartition(target_)};
  }
  /// <work>/profiles.jsonl, for the atlas page; written at every exit that
  /// writes the run's other files.
  void writePartitionGenera() const {
    runrecords::writePartitionGenera(
        cfg_.work, view(), [this](LinkId n) { return subjectName(n); },
        [this](LinkId n) { return searchSubject_.count(n) > 0; });
  }
  CertificateGoal certificateGoal() const {
    return {cfg_.targetName, cfg_.targetPD,  cfg_.goalGenus,         cfg_.goalLower,
            cfg_.goalDisjoint, target_,      goalPartition(target_)};
  }
  CertificateWriter certificates() const { return {g_, reg_, tableName_, sources_}; }
  void log(const std::string &line) { runrecords::append(cfg_.work, line); }

  GoalOptions cfg_;
  CobordismGraph g_;
  LinkRegistry reg_;
  linknaming::Tables tables_;
  linknaming::LinkNamer namer_;
  /// Each link's name and outside facts (bounds/axioms.h), shared with a
  /// depth-0 search's own graph; the target, its class, each node's depth
  /// and table name are its.
  LinkAxioms axioms_;
  LinkId &target_ = axioms_.target;
  std::string &targetCanonical_ = axioms_.targetClass;
  std::map<LinkId, int> &depth_ = axioms_.depth;
  std::map<LinkId, std::string> &tableName_ = axioms_.tableName;
  std::map<LinkId, std::vector<long>> expansions_;
  std::set<LinkId> refused_;
  /// Why the run must halt (an impossible state, divergence 2); empty if not.
  std::string halt_;
  /// What a checker needs to replay each witness edge (certificate.json).
  CobordismSources sources_;
  std::optional<linknaming::SignatureTable> signatures_;
  std::unique_ptr<Searcher> searcher_;
  /// The database's cobordisms as free edges (master_witnesses).
  std::unique_ptr<DatabaseCobordisms> database_;
  bool masterSubjectsFor(LinkId n) const {
    return database_ && database_->subjectsFor(n, axioms_, tables_);
  }
  bool masterDone(LinkId n) const { return database_ && database_->loaded(n); }
  void loadMaster(LinkId n, bool countsAsExpansion = true);
  /// The class a table name stands for, as the namer names nodes: its base's
  /// variants that are one oriented link up to mirror and global reversal,
  /// by diagram OR by a meridian-carrying isometry (ExactNamer::
  /// canonicalName(), data/table_link_classes.csv). ExactTables::canonical()
  /// joins by diagram only, so it must never be compared with a node's name.
  std::string classOf(const std::string &name) const { return axioms_.classOf(name); }
  double cpuSpent_ = 0, wallSpent_ = 0;
  int searches_ = 0;
  int invariantFailures_ = 0;
  std::map<LinkId, std::string> searchSubject_; ///< each searched node's subject name
  bool signed_ = false;
  size_t appended_ = 0;             ///< witnesses the store gained
  std::set<LinkId> boosted_;              ///< hubs already expanded wide (--hub-degree)
  std::string stopReason_ = "nothing-useful"; ///< why the loop ended short of the goal
  std::map<std::string, bool> special_;   ///< --lower-sources: name -> special
  int lowerLMax_ = 0;                     ///< the largest special lower bound
  /// Lower what-ifs by node, valid for one (witnesses, records, lower
  /// version) of the graph: nullopt when the what-if was inconsistent.
  struct LowerCache {
    std::tuple<size_t, size_t, long> version;
    std::map<LinkId, std::optional<int>> carried;
  };
  mutable LowerCache lowerCache_;
  /// Each searched node's latest usable frontier (searchfrontier.h): its next
  /// hop carries on from there instead of searching the prefix again.
  std::map<LinkId, SearchFrontier> frontiers_;
  /// Nodes whose search ran to the end at the hop shape: nothing is left.
  std::set<LinkId> searchedOut_;
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
  DriverTimes driver_, driverAtLastSearch_;
  SignResult signResult_; ///< storeWitnesses()'s counts and times
  /// Writes the run's own record to cascade.jsonl: its whole wall and CPU,
  /// and where the time outside the hops went.
  void logRun(double wall, double cpu, double startup, double loop, double databaseSeconds,
              double lowerReport, double linkBounds);
};

std::string Scheduler::subjectName(LinkId n) const {
  if (n == target_ && cfg_.targetName != "target") return cfg_.targetName;
  if (auto it = tableName_.find(n); it != tableName_.end() && tables_.entry(it->second))
    return it->second;
  return kFrozenCascadeSubjectPrefix + cfg_.runName + "/" + cfg_.targetName +
         kFrozenCascadeSubjectNodeMark + std::to_string(n);
}

void Scheduler::signIntoDatabase() {
  if (cfg_.cobordismsPath.empty() || signed_) return;
  signed_ = true;
  runrecords::writeLinksCsv(cfg_.work, searchSubject_, reg_);
  const SignResult s = signPending(cfg_.work, cfg_.cobordismsPath, cfg_.dedupeAgainst,
                                    names_, static_cast<unsigned>(cfg_.threads),
                                    cfg_.pairSigCache);
  appended_ = s.appended;
  signResult_ = s;
  std::cout << kFrozenWitnessStoreLine << s.kept << " kept, " << s.fresh << " new, "
            << s.appended << " appended to " << cfg_.cobordismsPath << " (signed in "
            << std::fixed << std::setprecision(0) << s.signSeconds << " s)\n";
}

void Scheduler::loadMaster(LinkId n, bool countsAsExpansion) {
  // A witness can be on an upper proof only if its genus is at most the
  // goal (glue() never lowers a genus), and on a lower proof only if the
  // charge it costs, at least its genus, is affordable: at most the largest
  // special source's bound minus the lower goal.
  DatabaseLoad load{g_,
                    reg_,
                    axioms_,
                    tables_,
                    sources_,
                    invariantFailures_,
                    std::max(cfg_.goalGenus, cfg_.goalLower >= 0 ? lowerLMax_ - cfg_.goalLower : -1),
                    static_cast<unsigned>(std::max(cfg_.threads, 1)),
                    cfg_.readBackCache,
                    target_,
                    goalPartition(target_),
                    [this](const std::string &line) { log(line); }};
  driver_.master += database_->load(n, load);
  // The target's own rows stand in for its first hop at this budget level;
  // a node loaded lazily is expanded right after, so its load counts nothing.
  if (countsAsExpansion && masterSubjectsFor(n)) expansions_[n].push_back(0);
}

bool Scheduler::useful(LinkId n) const {
  // What-if: give n the best profile it could conceivably have (every
  // partition its linking numbers allow, at its proved lower bound) and see
  // whether the target's goal would follow over the edges found so far.
  CobordismGraph what = g_;
  const GraphLink &link = what.link(n);
  // Optimistic only down to what is proved impossible: each partition at
  // its propagated lower bound (the linking condition included).
  for (const Partition &p : allPartitions(link.components)) {
    const int lo = g_.lower(n, p);
    if (lo < CobordismGraph::kNoSurface)
      what.addLeaf(n, p, lo, "what-if");
  }
  what.propagate();
  return goalMetIn(what, target_, goalPartition(target_), cfg_.goalGenus);
}

void Scheduler::loadLowerSources() {
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
      if (auto lo = linknaming::parseTableG4(f[loCol])) lowerLMax_ = std::max(lowerLMax_, lo->first);
  }
}

void Scheduler::lowerSlacks(const std::vector<LinkId> &ns) const {
  // One what-if per node (ProofGraph::lowerIf copies the graph and relaxes
  // it, so they are independent), on the run's threads, cached until the
  // graph changes.
  const std::tuple<size_t, size_t, long> version{g_.cobordismCount(), g_.derivationCount(),
                                                 g_.lowerVersion()};
  if (lowerCache_.version != version) {
    lowerCache_.version = version;
    lowerCache_.carried.clear();
  }
  std::vector<LinkId> todo;
  for (LinkId n : ns)
    if (!lowerCache_.carried.count(n)) todo.push_back(n);
  if (todo.empty()) return;
  const Partition goal = goalPartition(target_);
  std::vector<std::optional<int>> out(todo.size());
  parallelFor(todo.size(), static_cast<unsigned>(std::max(cfg_.threads, 1)), [&](size_t i) {
    const LinkId n = todo[i];
    const GraphLink &link = g_.link(n);
    // The most n could ever have, per partition: any transported bound
    // is a literature seed minus charges, so at most the largest special
    // source's (lowerLMax_); a proved surface refining P caps P; the
    // literature upper bound is a connected surface, capping the
    // coarsest partition only. Below what is known already a seed is a
    // no-op, and above a proved surface it is refused (nullopt).
    int litHi = std::numeric_limits<int>::max();
    if (auto t = tableName_.find(n); t != tableName_.end())
      if (const linknaming::TableEntry *e = tables_.entry(t->second))
        if (auto g4 = linknaming::parseTableG4(e->g4)) litHi = g4->second;
    std::vector<CobordismGraph::LowerSeed> seeds;
    for (const Partition &p : allPartitions(link.components)) {
      int cap = lowerLMax_;
      if (p.blocks() == 1) cap = std::min(cap, litHi);
      if (auto b = g_.best(n, p)) cap = std::min(cap, b->genus);
      seeds.push_back({n, p, cap});
    }
    out[i] = g_.lowerIf(seeds, target_, goal);
  });
  for (size_t i = 0; i < todo.size(); ++i) lowerCache_.carried[todo[i]] = out[i];
}

std::optional<int> Scheduler::lowerSlack(LinkId n) const {
  lowerSlacks({n});
  const auto &c = lowerCache_.carried.at(n);
  if (!c) return std::nullopt;
  return *c - cfg_.goalLower;
}

bool Scheduler::usefulLower(LinkId n, int *slack) const {
  if (cfg_.goalLower < 0 || n == target_) return false;
  if (g_.link(n).components > CobordismGraph::kMaxLowerComponents) return false;
  if (reg_.info(n).diagram.crossings() > cfg_.lowerMaxCrossings) return false;
  auto s = lowerSlack(n);
  if (!s || *s < 0) return false;
  if (slack) *slack = *s;
  return true;
}

std::optional<LinkId> Scheduler::choose() {
  const auto tChoose = Clock::now();
  struct ChooseTimer {
    DriverTimes &d;
    Clock::time_point t;
    ~ChooseTimer() { d.choose += secondsSince(t); }
  } chooseTimer{driver_, tChoose};
  struct Cand {
    LinkId n;
    bool lowerOnly = false; ///< kept by the lower gate alone
    int slack = 0;          ///< lower gate: charge still affordable
    std::tuple<int, int, int, size_t> key;
  };
  std::vector<Cand> cands;
  std::vector<LinkId> eligible;
  for (const auto &[n, dep] : depth_) {
    if (!reg_.known(n) || n == reg_.unknot() || refused_.count(n)) continue;
    const LinkInfo &ni = reg_.info(n);
    if (ni.diagram.crossings() > cfg_.maxCrossings) continue;
    if (searchedOut_.count(n)) continue; // nothing left to search
    if (!expansions_[n].empty()) continue; // one expansion per budget level
    eligible.push_back(n);
  }
  // The upper gate first; the lower what-ifs for every node it rejects run
  // as one parallel batch (each is a graph copy relaxed to a fixed point).
  std::map<LinkId, bool> upper;
  std::vector<LinkId> needLower;
  const auto tUseful = Clock::now();
  for (LinkId n : eligible) {
    upper[n] = n == target_ || useful(n);
    if (!upper[n] && cfg_.goalLower >= 0 && n != target_ &&
        g_.link(n).components <= CobordismGraph::kMaxLowerComponents &&
        reg_.info(n).diagram.crossings() <= cfg_.lowerMaxCrossings)
      needLower.push_back(n);
  }
  driver_.useful += secondsSince(tUseful);
  const auto tLower = Clock::now();
  if (!needLower.empty()) lowerSlacks(needLower);
  driver_.lowerSlack += secondsSince(tLower);
  for (LinkId n : eligible) {
    const LinkInfo &ni = reg_.info(n);
    const int dep = depth_.at(n);
    int slack = 0;
    const bool lowerOnly = !upper[n] && usefulLower(n, &slack);
    const bool use = upper[n] || lowerOnly;
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
  auto roundedVolume = [&](LinkId n) {
    const LinkInfo &i = reg_.info(n);
    return i.hyperbolic ? std::llround(i.volume * 1e6) : std::numeric_limits<long long>::max();
  };
  std::sort(cands.begin(), cands.end(), [&](const Cand &a, const Cand &b) {
    if (cfg_.strategy == "best") {
      const LinkInfo &ia = reg_.info(a.n), &ib = reg_.info(b.n);
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

void Scheduler::expand(LinkId n, long surfaces) {
  // A node searched before carries on from where that search stopped (the
  // hop's surface target is its breadth, so it adds only what is new); one
  // already searched this far (a hub's wide hop, say) has nothing new at
  // this budget.
  const SearchFrontier *resume = nullptr;
  if (auto f = frontiers_.find(n); f != frontiers_.end()) resume = &f->second;
  if (resume && resume->satisfying >= surfaces) {
    expansions_[n].push_back(surfaces);
    std::cout << "[+] node " << n << " already searched to " << resume->satisfying
              << " surfaces; nothing new at " << surfaces << "\n";
    return;
  }
  const int k = searches_++;
  const std::string dir = cfg_.work + "/" + kFrozenHopDirPrefix + std::to_string(k) +
                          kFrozenHopDirNodeMark + std::to_string(n);
  fs::create_directories(dir);
  const GaussDiagram &d = reg_.info(n).diagram;
  SearchedLink searched;
  searched.link = n;
  searched.diagram = d;
  searched.linkMap.resize(d.components());
  std::iota(searched.linkMap.begin(), searched.linkMap.end(), 0);
  searched.pd = diagramPD(d);
  searched.layers = *cfg_.runShape.layers;
  using clock = std::chrono::steady_clock;
  auto seconds = [](clock::time_point a, clock::time_point b) {
    return std::chrono::duration<double>(b - a).count();
  };
  const auto tBuild = clock::now();
  std::unique_ptr<CobordismAssembler> assembler;
  try {
    assembler = std::make_unique<CobordismAssembler>(g_, reg_, searched);
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
        ",\"refused\":\"" + json::escape(e.what()) + "\",\"pd\":\"" + json::escape(searched.pd) +
        "\",\"pd_ambiguous\":" + (d.link().pdAmbiguous() ? "true" : "false") +
        ",\"diagram\":" + gauss.str() + "}");
    std::cout << "[!] node " << n << " refused: " << e.what() << "\n";
    return;
  }
  // Where a hop's time goes, logged per hop: the row's build and
  // certification, the search's setup and the search itself, adding its
  // surfaces, naming new nodes, and relaxing the graph.
  double buildSeconds = seconds(tBuild, clock::now()), setupSeconds = 0, searchSeconds = 0;
  std::string roundsJson = "[]";
  size_t drainTail = 0;
  double drainTailSeconds = 0;
  std::string namingJson; // the drain's naming times
  // The hop's subject: what its witnesses are recorded under, and what its
  // log lines are named by (subjectName()).
  const std::string subject = subjectName(n);
  searchSubject_[n] = subject;
  const size_t linksBefore = g_.linkCount();
  int assembled = 0, failed = 0;
  size_t cobordisms = 0;
  // One kept surface into the graph. Its edge's key is its provenance:
  // hop<k>#<i> (it has no pair signature).
  const std::string build = assembler->redrawer().buildChecksum();
  auto take = [&](const std::string &key, const std::string &label,
                  const std::function<AddedCobordism()> &add, std::vector<int> faces) {
    AddedCobordism e;
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
      CobordismSource info{dir, searched.pd, key, e, 2, "", searched.linkMap, std::move(faces),
                    inProcess ? build : std::string()};
      if (e.direct) sources_.direct[key] = std::move(info);
      else sources_.byCobordism[e.cobordism] = std::move(info);
    } else {
      ++failed;
      std::cout << "[!] witness " << key << " (" << label << "): " << e.why << "\n";
    }
  };

  struct {
    int status = -1;
    double wall = 0, cpu = 0;
  } r;
  std::chrono::steady_clock::time_point t0;
  {
    SearchResult run;
    try {
      SearchRequest request =
          searcher_->requestFor(assembler->redrawer(), subject, surfaces, cfg_.searchSeconds);
      request.resume = resume;
      if (tables_.entry(subject)) request.censusName = linknaming::baseName(subject);
      // Every kept surface, durably, as the search runs (divergence 7): the
      // hop's pending file, signed at the run's end (storeWitnesses()).
      request.incomingPD = searched.pd;
      request.layers = searched.layers;
      if (!cfg_.cobordismsPath.empty()) request.pending = dir + "/kept.csv";
      request.runDirectory = cfg_.work;
      run = searcher_->run(assembler->redrawer().thickened(), request);
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
        << "[+] " << subject << " " << searched.pd << "\n[+] " << subject << ": "
        << run.kept.size() << " kept, outcome " << run.outcome << "\n[+] " << subject
        << ": accounting: " << run.accounting << "\n[+] " << subject
        << ": diagram naming: " << run.naming << "\n";
    namingJson = [&] {
      std::ostringstream o;
      o << std::fixed << std::setprecision(1) << ",\"naming_diagram_s\":"
        << run.namingDiagramSeconds << ",\"naming_fallback_s\":" << run.namingFallbackSeconds
        << ",\"naming_exact_s\":" << run.namingOrientedSeconds
        << ",\"naming_slowest_s\":" << run.namingSlowestSeconds;
      return o.str();
    }();
    // Every hop's accounting in the driver log too, in verifyslicegenus's
    // shape after the hop number, so a campaign audits each hop as it
    // audits a row (tools/orchestrate/audit_rows.py).
    std::cout << "[+] " << kFrozenHopLine << k << " " << subject << ": accounting: "
              << run.accounting << "\n[+] " << kFrozenHopLine << k << " " << subject
              << ": diagram naming: " << run.naming
              << "\n";
    if (run.impossible > 0) {
      // Divergence 2: a state that cannot occur halts the run, once this
      // hop's finds are recorded (below) and stored (run()), even if they
      // meet the goal.
      halt_ = kFrozenHopLine + std::to_string(k) + ": surface accounting failed -- " +
              std::to_string(run.impossible) + " surfaces hit a state that cannot occur";
      std::cout << "[!!] HALT: " << halt_ << "\n";
    } else if (!run.accountingFailure.empty()) {
      std::cout << "[!!] " << kFrozenHopLine << k << ": surface accounting failed -- "
                << run.accountingFailure << " (completeness only: nothing unsound "
                << "is recorded)\n";
    }
    // The node's breadth so far, and where its next hop carries on from.
    std::cout << "[+] " << kFrozenHopLine << k << " " << subject << ": breadth: "
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
        std::cout << "[!] " << kFrozenHopLine << k << ": frontier not written: " << e.what()
                  << "\n";
      }
      frontiers_[n] = std::move(*run.frontier);
    } else {
      frontiers_.erase(n); // not vouched for: the next search starts afresh
    }
    driver_.kept += secondsSince(tKept);
    t0 = std::chrono::steady_clock::now();
    cobordisms = run.kept.size();
    for (size_t i = 0; i < run.kept.size(); ++i) {
      KeptSurface &ks = run.kept[i];
      const std::string key = kFrozenHopKeyPrefix + std::to_string(k) + "#" + std::to_string(i);
      take(key, ks.outgoingName, [&] { return assembler->addRead(ks.link, ks.genus, key); },
           std::move(ks.faces));
    }
  }
  cpuSpent_ += r.cpu;
  wallSpent_ += r.wall;
  expansions_[n].push_back(surfaces);
  const auto tLinks = clock::now();
  const double addSeconds = seconds(t0, tLinks);
  onNewLinks(linksSince(linksBefore), depth_.at(n) + 1);
  const auto tProp = clock::now();
  g_.propagate();
  g_.propagateLower();
  const double linkSeconds = seconds(tLinks, tProp), propagateSeconds = seconds(tProp, clock::now());
  const int targetLower = g_.lower(target_, goalPartition(target_));
  const double assemble = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - t0).count();
  auto best = g_.best(target_, goalPartition(target_));
  std::ostringstream o;
  o << "{\"hop\":" << k << ",\"node\":" << n << ",\"crossings\":" << d.crossings()
    << ",\"surfaces\":" << surfaces << ",\"status\":" << r.status
    << ",\"wall\":" << std::fixed << std::setprecision(1) << r.wall << ",\"cpu\":" << r.cpu
    << ",\"assemble\":" << assemble << ",\"row\":" << buildSeconds
    << ",\"setup\":" << setupSeconds << ",\"search\":" << searchSeconds
    << ",\"add\":" << addSeconds << ",\"name_nodes\":" << linkSeconds
    << ",\"propagate\":" << propagateSeconds << ",\"rounds\":" << roundsJson
    << ",\"drain_tail\":" << drainTail << ",\"drain_tail_s\":" << drainTailSeconds
    << namingJson << ",\"witnesses\":" << cobordisms
    << ",\"assembled\":" << assembled << ",\"failed\":" << failed
    << ",\"nodes\":" << g_.linkCount() << ",\"new_nodes\":" << (g_.linkCount() - linksBefore)
    << ",\"records\":" << g_.derivationCount()
    << ",\"target_best\":" << (best ? std::to_string(best->genus) : "null")
    << ",\"target_lower\":" << targetLower
    << ",\"contradictions\":" << g_.contradictions().size()
    << ",\"invariant_failures\":" << invariantFailures_
    // The driver's time since the previous hop record: choosing this node
    // (and the loads before it), and this hop's kept.csv and frontier.
    << ",\"choose_s\":" << driver_.choose - driverAtLastSearch_.choose
    << ",\"useful_s\":" << driver_.useful - driverAtLastSearch_.useful
    << ",\"lower_slack_s\":" << driver_.lowerSlack - driverAtLastSearch_.lowerSlack
    << ",\"master_s\":" << driver_.master - driverAtLastSearch_.master
    << ",\"kept_s\":" << driver_.kept - driverAtLastSearch_.kept << "}";
  driverAtLastSearch_ = driver_;
  log(o.str());
  std::cout << "[+] " << kFrozenHopLine << k << ": node " << n << " (" << d.crossings()
            << " crossings, "
            << d.components() << " components): " << cobordisms << " witnesses, "
            << assembled << " assembled, " << (g_.linkCount() - linksBefore)
            << " new nodes; " << std::fixed << std::setprecision(0) << r.wall << " s wall, "
            << r.cpu << " s CPU; target best " << (best ? std::to_string(best->genus) : "none")
            << "\n";
}

void Scheduler::logRun(double wall, double cpu, double startup, double loop, double databaseSeconds,
                     double lowerReport, double linkBounds) {
  std::ostringstream o;
  o << std::fixed << std::setprecision(1) << "{\"run\":\"" << json::escape(cfg_.targetName)
    << "\",\"threads\":" << cfg_.threads << ",\"wall\":" << wall << ",\"cpu\":" << cpu
    << ",\"cores\":" << (wall > 0 ? cpu / wall : 0.0) << ",\"startup_s\":" << startup
    << ",\"loop_s\":" << loop << ",\"hop_wall_s\":" << wallSpent_
    << ",\"hop_cpu_s\":" << cpuSpent_ << ",\"choose_s\":" << driver_.choose
    << ",\"useful_s\":" << driver_.useful << ",\"lower_slack_s\":" << driver_.lowerSlack
    << ",\"master_s\":" << driver_.master << ",\"kept_s\":" << driver_.kept
    << ",\"store_s\":" << databaseSeconds << ",\"store_dedupe_s\":" << signResult_.dedupeSeconds
    << ",\"store_sign_s\":" << signResult_.signSeconds
    << ",\"lower_report_s\":" << lowerReport << ",\"node_bounds_s\":" << linkBounds
    << ",\"hops\":" << searches_ << ",\"nodes\":" << g_.linkCount() << "}";
  log(o.str());
}

int Scheduler::run() {
  const auto tRun = Clock::now();
  const double cpuRun = timers::processCpuSeconds();
  fs::create_directories(cfg_.work);
  // The process-wide settings (census, census writes, Pachner searches, the
  // complement cache) were set from the config before this run began
  // (setup::applyRunSettings()).
  if (!cfg_.censusLoaded)
    std::cout << "[!] census not found at " << cfg_.censusDb << "\n";
  {
    const auto t0 = std::chrono::steady_clock::now();
    // From the node namer's tables: one table load (phase 5).
    signatures_ = linknaming::SignatureTable::fromTables(tables_);
    // The searches' outgoing namers share the node namer's table caches.
    searcher_ = std::make_unique<Searcher>(*signatures_, &tables_, cfg_.runShape,
                                              static_cast<unsigned>(cfg_.threads),
                                              namer_.caches());
    std::cout << "[+] hops in process: " << signatures_->knots() << " knot and "
              << signatures_->links() << " link diagram signatures ("
              << std::fixed << std::setprecision(1)
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()
              << " s)\n";
  }
  {
    const RunShape &s = cfg_.runShape;
    std::cout << "[+] hop shape: cap " << *s.maxFaces << ", IDDFS " << *s.iddfsIterations
              << " from " << *s.iddfsStart << " step " << *s.iddfsStep << ", root budget "
              << *s.rootBudgetStart << " x" << *s.rootBudgetGrowth << "\n";
  }
  printShape();
  loadLowerSources();
  // Table names and symmetry types, for the slice-composite anchors
  // (NodeAxioms) and the store step: as verifyslicegenus loads them.
  size_t symmetryTypes = 0;
  names_ = solver::loadTableNames(cfg_.knotTable, cfg_.linkTable, cfg_.knotSymmetry,
                                        &symmetryTypes);
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
  if (!cfg_.masterCobordisms.empty()) {
    const auto t0 = std::chrono::steady_clock::now();
    database_ = std::make_unique<DatabaseCobordisms>(
        cfg_.masterCobordisms, std::vector<std::string>{cfg_.knotTable, cfg_.linkTable});
    std::cout << "[+] master witnesses: " << database_->subjects() << " subjects indexed in "
              << std::fixed << std::setprecision(0)
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()
              << " s (read only)\n";
  }
  GaussDiagram raw = GaussDiagram::of(linknaming::linkFromTablePD(cfg_.targetPD));
  GaussDiagram simp = simplifyKeepingComponents(raw);
  if (linknaming::splitPieces(simp).size() != 1)
    throw std::runtime_error("the target is a split diagram; give one piece");
  // The target's own table name, so its literature value never proves it.
  bool named = false;
  std::string composite;
  try {
    auto pn = namer_.namePiece(simp);
    if (pn.names.size() == 1 && pn.by != linknaming::PieceName::By::untabulated) {
      targetCanonical_ = classOf(pn.names.front());
      named = true;
    } else if (simp.components() == 1) {
      // A composite target: its whole-diagram name (never an anchor for
      // itself: NodeAxioms skips the target), reported and recorded.
      linknaming::LinkName fs = namer_.name(simp.link());
      if (fs.isName && fs.pinned && fs.pieces.size() >= 2 &&
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
  LinkMatch t = reg_.intern(simp, "target " + cfg_.targetName);
  target_ = t.link;
  onNewLink(target_, 0);
  if (!composite.empty()) {
    tableName_[target_] = composite;
    std::cout << "[+] target is the composite " << composite
              << (linknaming::isElementarySlice(composite, names_.symmetries())
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
  long budget = cfg_.surfaceTarget;
  while (true) {
    // A halt first (divergence 2): even a goal met by the hop that found an
    // impossible state is not reported as met.
    if (!halt_.empty()) {
      // The surfaces found are real whatever broke, so they are kept.
      signIntoDatabase();
      writePartitionGenera();
      printOutcome("halted");
      return 2;
    }
    if (!g_.contradictions().empty()) {
      for (const auto &c : g_.contradictions()) std::cout << "[!!] CONTRADICTION: " << c << "\n";
      // The surfaces are real whatever the contradiction's cause (a naming
      // or solver bug), so they are kept, as verifyslicegenus writes its
      // witnesses before its fatal-bug halt.
      signIntoDatabase();
      writePartitionGenera();
      printOutcome("contradiction");
      return 3;
    }
    // Checked after the gates (divergence 6): a contradiction met together
    // with the goal still halts the run, never certifies through it.
    if (goalMet()) break;
    // A first SIGINT or SIGTERM ended the search that was running, cleanly
    // (divergence 8): no further search.
    if (runsignals::interrupted()) {
      std::cout << "[!] " << runsignals::name() << ": the run searches no further\n";
      stopReason_ = "interrupted";
      break;
    }
    // Free edges first (no search, so no budget): the target's own master
    // rows, if it is a table entry the atlas searched.
    if (database_ && !masterDone(target_) && masterSubjectsFor(target_)) {
      loadMaster(target_);
      continue;
    }
    // Any other node's stored rows are loaded only when that node is the
    // one about to be expanded (below): a proof can run through a node only
    // when the search picks it, and reading every table node's rows as it
    // was met cost more than the hops (2026-09-30, close1: 35-45 loads and
    // ~2,500 read-backs per row).
    auto n = choose();
    if (searches_ >= cfg_.maxSearches) {
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
      if (budget * 2 > cfg_.maxSurfaceTarget) {
        std::cout << "[-] nothing useful left\n";
        stopReason_ = "nothing-useful";
        break;
      }
      budget *= 2;
      for (auto &[m, v] : expansions_) v.clear();
      std::cout << "[+] raising the hop budget to " << budget << " surfaces\n";
      continue;
    }
    if (database_ && !masterDone(*n) && masterSubjectsFor(*n)) {
      // Its stored rows first: free edges, which may close the proof
      // without the hop, and are in any case what the hop would refind.
      loadMaster(*n, /*countsAsExpansion=*/false);
      if (goalMet() || !g_.contradictions().empty()) continue;
    }
    long surfaces = budget;
    if (cfg_.hubDegree > 0 && !boosted_.count(*n) &&
        g_.link(*n).cobordisms.size() >= cfg_.hubDegree && cfg_.hubSurfaces > budget) {
      // A hub: many routes meet here, so one wide hop from it buys many
      // more first-level candidates than another narrow one elsewhere.
      surfaces = cfg_.hubSurfaces;
      boosted_.insert(*n);
      std::cout << "[+] hub: node " << *n << " has " << g_.link(*n).cobordisms.size()
                << " witness edges; expanding it at " << surfaces << " surfaces\n";
    }
    expand(*n, surfaces);
  }
  const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  auto best = g_.best(target_, goalPartition(target_));
  std::cout << "[+] done: " << searches_ << " hops, " << g_.linkCount() << " nodes, "
            << g_.derivationCount() << " records; " << std::fixed << std::setprecision(0)
            << wall << " s wall, " << cpuSpent_ << " s search CPU. Target best: "
            << (best ? std::to_string(best->genus) : "none") << "\n";
  const auto tDatabase = Clock::now();
  signIntoDatabase();
  const double databaseSeconds = secondsSince(tDatabase);
  const auto tReport = Clock::now();
  if (cfg_.lowerReport)
    runrecords::writeLowerReport(cfg_.work, view(), tables_, cfg_.targetName, special_,
                                 static_cast<unsigned>(std::max(cfg_.threads, 1)));
  const double reportSeconds = secondsSince(tReport);
  const auto tBounds = Clock::now();
  runrecords::writeLinkBounds(cfg_.work, view());
  writePartitionGenera();
  logRun(secondsSince(tRun), timers::processCpuSeconds() - cpuRun, startupSeconds, wall, databaseSeconds,
         reportSeconds, secondsSince(tBounds));
  if (lowerMet() && !upperMet()) {
    // The bound's proof: the reasons from the target down to their leaves.
    const Partition goal = goalPartition(target_);
    std::cout << "[+] LOWER GOAL MET: " << cfg_.targetName << " genus >= "
              << g_.lower(target_, goal) << " (literature-assisted); the proof:\n";
    certificates().describeLower(std::cout, target_, goal, 0);
    certificates().writeLower(cfg_.work + "/lower_certificate.json", certificateGoal());
    std::cout << "[+] lower certificate " << cfg_.work << "/lower_certificate.json\n";
    printOutcome("met");
    return 0;
  }
  if (goalMet()) {
    certificates().writeUpper(cfg_.work + "/certificate.json", certificateGoal());
    bool constructive = true;
    for (DerivationId r : g_.proof(best->derivation))
      if (g_.derivation(r).kind == DerivationKind::leaf &&
          g_.derivation(r).source.rfind("literature", 0) == 0)
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

void Scheduler::printOutcome(const std::string &outcome) const {
  // The line a campaign's runner parses, in verifyslicegenus's own shape
  // (dispatch.py RE_OUTCOME): witnesses newly recorded, and why the run ended.
  std::cout << "[+] " << cfg_.targetName << ": " << appended_
            << kFrozenNewWitnessesOutcome << outcome << "\n";
}

void Scheduler::printShape() const {
  // Everything that decides what a run covers, as key=value, so a campaign
  // records what actually ran rather than what its configuration asked for.
  // The keys are cascadesearch's option names (frozen: campaigns record
  // them), not the config's: hop_surfaces is surface_target,
  // max_hop_surfaces max_surface_target, max_expansions max_searches,
  // recognition_cache_limit complement_cache_limit, master_witnesses and
  // witness_store the databases read and signed into.
  const RunShape &s = cfg_.runShape;
  std::cout << "[+] profile: goal=" << (cfg_.goalDisjoint ? "disjoint" : "connected")
            << " goal_genus=" << cfg_.goalGenus << " goal_lower=" << cfg_.goalLower
            << " lower_max_crossings=" << cfg_.lowerMaxCrossings
            << " master_loads=lazy"
            << " literature=" << (cfg_.literature ? 1 : 0)
            << " hop_mode=process hop_surfaces=" << cfg_.surfaceTarget
            << " max_hop_surfaces=" << cfg_.maxSurfaceTarget
            << " max_expansions=" << cfg_.maxSearches << " cpu_budget=" << cfg_.cpuBudget
            << " strategy=" << cfg_.strategy << " max_crossings=" << cfg_.maxCrossings
            << " hub_degree=" << cfg_.hubDegree << " hub_surfaces=" << cfg_.hubSurfaces
            << " threads=" << cfg_.threads << " max_faces=" << *s.maxFaces
            << " iddfs_iterations=" << *s.iddfsIterations << " iddfs_start=" << *s.iddfsStart
            << " iddfs_step=" << *s.iddfsStep << " root_budget_start=" << *s.rootBudgetStart
            << " root_budget_growth=" << *s.rootBudgetGrowth << " layers=" << *s.layers
            << " resolve_unlinked=" << (*s.resolveUnlinked ? 1 : 0)
            << " exact_far_side_names=1 pending_surface_cap=" << *s.pendingSurfaceCap
            << " petal_cache_limit=" << *s.petalCacheLimit
            << " boundary_signature_cache_limit=" << *s.boundarySignatureCacheLimit
            << " recognition_cache_limit=" << cfg_.complementCacheLimit
            << " master_witnesses=" << (cfg_.masterCobordisms.empty() ? "none" : cfg_.masterCobordisms)
            << " witness_store=" << (cfg_.cobordismsPath.empty() ? "none" : cfg_.cobordismsPath)
            << " run_name=" << (cfg_.runName.empty() ? "none" : cfg_.runName) << "\n";
}

} // namespace

GoalOptions goalOptions(const config::Config &cfg) {
  GoalOptions o;
  o.targetPD = cfg.text("target_pd");
  o.targetName = cfg.text("target_name");
  o.work = cfg.text("work");
  o.knotTable = cfg.text("knot_table");
  o.linkTable = cfg.text("link_table");
  o.knotSymmetry = cfg.text("knot_symmetry");
  o.censusDb = cfg.text("census");
  o.goalGenus = static_cast<int>(cfg.integer("goal_genus"));
  o.goalDisjoint = cfg.text("goal_partition") == "disjoint";
  o.surfaceTarget = static_cast<long>(cfg.integer("surface_target"));
  o.maxSurfaceTarget = static_cast<long>(cfg.integer("max_surface_target"));
  o.threads = static_cast<int>(cfg.threads());
  o.maxSearches = static_cast<int>(cfg.integer("max_searches"));
  o.cpuBudget = cfg.real("cpu_budget");
  o.searchSeconds = cfg.real("search_seconds");
  o.maxCrossings = static_cast<size_t>(cfg.integer("max_crossings"));
  o.strategy = cfg.text("strategy");
  o.literature = cfg.flag("literature");
  o.masterCobordisms = cfg.text("master_witnesses");
  RunShape &h = o.runShape;
  h.layers = static_cast<int>(cfg.integer("layers"));
  h.maxFaces = cfg.integer("max_faces");
  h.iddfsIterations = static_cast<unsigned>(cfg.integer("iddfs_iterations"));
  h.iddfsStart = cfg.integer("iddfs_start");
  h.iddfsStep = cfg.integer("iddfs_step");
  h.rootBudgetStart = cfg.integer("root_budget_start");
  h.rootBudgetGrowth = cfg.integer("root_budget_growth");
  h.resolveUnlinked = cfg.flag("resolve_unlinked");
  h.pendingSurfaceCap = static_cast<size_t>(cfg.integer("pending_surface_cap"));
  h.petalCacheLimit = static_cast<size_t>(cfg.integer("petal_cache_limit"));
  h.boundarySignatureCacheLimit =
      static_cast<size_t>(cfg.integer("boundary_signature_cache_limit"));
  o.complementCacheLimit = static_cast<size_t>(cfg.integer("complement_cache_limit"));
  o.cobordismsPath = cfg.text("cobordisms");
  o.runName = cfg.text("run_name");
  o.pairSigCache = cfg.text("pair_sig_cache");
  o.readBackCache = cfg.text("read_back_cache");
  o.dedupeAgainst = cfg.paths("dedupe_against");
  o.lowerReport = cfg.flag("lower_report");
  o.lowerSources = cfg.text("lower_sources");
  o.goalLower = static_cast<int>(cfg.optionalInteger("goal_lower").value_or(-1));
  o.lowerMaxCrossings = static_cast<size_t>(cfg.integer("lower_max_crossings"));
  o.hubDegree = static_cast<size_t>(cfg.integer("hub_degree"));
  o.hubSurfaces = static_cast<long>(cfg.integer("hub_surfaces"));
  // What the retired cascadesearch refused, still refused.
  if (!cfg.flag("exact_far_side_names"))
    throw config::Error("exact_far_side_names cannot be 0 in a run with a goal (its "
                        "searches always name outgoing links exactly)");
  if (!o.cobordismsPath.empty() && o.runName.empty())
    throw config::Error("cobordisms needs run_name in a run with a goal");
  if (o.goalLower >= 0 && o.lowerSources.empty())
    throw config::Error("goal_lower needs lower_sources (the special sources' largest bound "
                        "is what makes the lower gate prune)");
  if (o.targetName.empty()) o.targetName = "target";
  return o;
}

int runToGoal(const GoalOptions &options) {
  Scheduler scheduler(options);
  return scheduler.run();
}

} // namespace scheduler
