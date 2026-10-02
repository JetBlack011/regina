//
//  search.h
//
//  One search from a link, as every driver runs it.
//

#ifndef SURFER_COBOUND_SEARCH_H
#define SURFER_COBOUND_SEARCH_H

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include "surfer/submanifold/submanifold.h"
#include "cobound/cobordisms/cobordism.h"
#include "cobound/cobordisms/pending.h"
#include "cobound/solver/solver.h"
#include "linknaming/tables.h"
#include "linknaming/complement/complementcache.h"
#include "cobound/outgoing/outgoinglink.h"
#include "cobound/outgoing/outgoingnamer.h"
#include "cobound/outgoing/fromdatabase.h"
#include "surfer/enumeration/searchstrategy.h"
#include "surfer/enumeration/surfacesearch.h"

/*! \file utils/surfer/cobound/search/search.h
 *  \brief THE search from a link (HopSearcher::run()), as every driver runs
 *  it: which BoundaryCondition it runs under (conditionFor()), the watchdog
 *  that ends it at its surface target or deadline (RowWatchdog), and what
 *  each caller still does its own way (SearchPolicy).
 *
 *  verifyslicegenus searches each table row through it (one witness per
 *  identity across the whole database, signed during the search), and the
 *  cascade runs each hop through it in process: the row's own thickening
 *  (the one its HopAssembler certified and reads from), the campaign's
 *  search shape, and the gates and accounting of preconditions.h, with the
 *  accepted surfaces read straight off the search instead of written out as
 *  pair signatures and read back.
 *
 *  What that saves, per hop: a process start and table load, a pair
 *  signature per kept surface (3.4 s each on a 10-crossing row, 2026-09-28:
 *  64 s of a 152 s hop), and decoding them again. A kept surface keeps its
 *  faces instead, from which its pair signature can be computed later
 *  (cobordisms/pairsigner.h), for the few witnesses a certificate needs.
 */

namespace rowsearch {

/** Which BoundaryCondition to search a row under; see conditionFor(). */
enum class BoundaryConditionMode { automatic, connected, proper };

/**
 * The BoundaryCondition for a row of `componentCount` components.
 *
 * `connected` (one curve per ambient boundary component) prunes far harder,
 * but it caps the curve count on the far side too, so a knot row searched
 * under it can never find a knot-to-link cobordism. It is honoured only for
 * a knot: a multi-component link can never meet it on its own search side.
 * `automatic` is `connected` for a knot and `proper` otherwise.
 */
BoundaryCondition conditionFor(BoundaryConditionMode mode, int componentCount);

/** When RowWatchdog ends a row's search. Unset limits never fire. */
struct WatchdogLimits {
    std::optional<long long> surfaceTarget;
    std::optional<double> rowSeconds;
    std::optional<double> sweepSeconds;
    std::chrono::steady_clock::time_point sweepStart{};
    std::optional<double> quiescenceSeconds;
    /** Milliseconds since the row last taught us something new; read only
     *  with `quiescenceSeconds`. */
    std::function<long long()> idleMillis;

    bool any() const {
        return surfaceTarget || rowSeconds || sweepSeconds || quiescenceSeconds;
    }
};

/**
 * Polls every 200 ms, while the search and its drain run, and calls
 * `endRow(why)` at most once: "surface-target", "timeout" or "quiescent".
 * stop() wakes it at once rather than waiting out its poll.
 * The surface target is checked before the clocks, so a row reaching it in
 * the same tick as a deadline records "surface-target": the two mean
 * different things to anyone later reading the negative.
 *
 * No thread is started when no limit is set.
 */
class RowWatchdog {
public:
    RowWatchdog(WatchdogLimits limits, std::function<void(const char *)> endRow);
    ~RowWatchdog();
    RowWatchdog(const RowWatchdog &) = delete;
    RowWatchdog &operator=(const RowWatchdog &) = delete;

    /** The live count of surfaces satisfying the boundary condition, from
     *  SearchCallbacks::onProgress (the only place it is handed out). */
    void publishSatisfying(long long count) {
        satisfying_.store(count, std::memory_order_relaxed);
    }
    /** Stops polling and joins; idempotent. Call once the search returns. */
    void stop();

private:
    WatchdogLimits limits_;
    std::function<void(const char *)> endRow_;
    std::atomic<long long> satisfying_{0};
    std::atomic<bool> done_{false};
    std::mutex wakeMutex_;
    std::condition_variable wake_;
    std::thread thread_;
};

} // namespace rowsearch

namespace cascade {

/// The campaign's search shape (atlas tools/orchestrate/hosts.conf
/// [campaign], root budget 840 from c5 on): the one place the cascade's
/// hops take it from, in either hop mode, and what its profile line prints.
struct HopShape {
  /// Thicken and collar layers of a hop's row (both): the in-process row,
  /// a child's --thicken-layers/--collar-layers, and the profile's layers=.
  /// A master row's layers are its witnesses' own.
  int layers = 2;
  long long maxFaces = 5;
  unsigned iddfsIterations = 2;
  long long iddfsStart = 4;
  long long iddfsStep = 1;
  long long rootBudgetStart = 840;
  long long rootBudgetGrowth = 2;
  /// Whether surfaces whose only self-intersections are unlinked count. No
  /// default (plan divergence 3): cascadesearch takes it from
  /// --resolve-unlinked or --no-resolve-unlinked, and a search refuses a
  /// shape without it.
  std::optional<bool> resolveUnlinked;
  /// hosts.conf's per-host limits, identical on every host. The pending
  /// cap in particular: at the binary's default (500,000) a search pauses
  /// to drain its queue, which a campaign row never does.
  size_t pendingSurfaceCap = 20'000'000;
  size_t petalCacheLimit = 12'000'000;
  size_t boundarySignatureCacheLimit = 1'000'000;
  /// Process-wide (identify::recognitionCacheLimit); set by the driver.
  size_t recognitionCacheLimit = 1'500'000;
};

/// How one search runs: the arguments of SurfaceSearch::search() that fix
/// its traversal and what it accepts, and the limits of its shared
/// structures. Every caller fills all of it (searchShape() for a hop).
struct SearchShape {
  BoundaryCondition condition = BoundaryCondition::proper;
  unsigned iddfsIterations = 0;
  long long iddfsStep = 0;
  std::optional<long long> iddfsStart;
  std::optional<unsigned> iddfsFinalThreads;
  std::optional<long long> maxFaces; ///< none: unbounded
  long long rootBudgetStart = 0;
  long long rootBudgetGrowth = 2;
  /// Accept surfaces whose only self-intersections are unlinked (paper
  /// §4.5). A required input of every search, with no default anywhere
  /// (plan divergence 3): it changes which surfaces count toward a surface
  /// target and is part of every frontier's fingerprint, so a caller that
  /// forgot it must not get one silently. HopSearcher::run() refuses a
  /// shape without it.
  std::optional<bool> resolveUnlinked;
  SurfaceSearchLimits limits;
};

/// A hop's SearchShape: `proper`, `shape`'s rounds, cap and budgets, its
/// limits, one name per multi-curve component and no pair signatures (a
/// hop keeps faces).
SearchShape searchShape(const HopShape &shape);

/**
 * Where one caller's searches still differ from the other's. Phase 4(a) put
 * verifyslicegenus's row loop and the cascade's in-process hops on one
 * search (HopSearcher::run()) without changing what either does; each
 * remaining difference is one field here, named after the plan's "Unified
 * divergences" entry that removes it in phase 4(b). The defaults are the
 * cascade's behaviour; verifyslicegenus sets its own. Fields read by the
 * drivers rather than the search are marked so: the drivers' own per-search
 * code (verifyslicegenus's row loop, the cascade's hop) reads them from
 * here too, so that every difference is listed in one place.
 *
 * Divergence 3 (defaults) is done: the search's inputs have one set of
 * defaults (layers 2/2, `proper`), and resolve_unlinked none. Divergence 8
 * (signals) has no field: it does not differ inside a search (both record
 * "interrupted" when the library's SIGINT scope stops one, and each driver
 * then carries on with its next search). Retriangulation on a census miss
 * is process-wide, set in each driver's main (verifyslicegenus: on unless
 * --no-retriangulate-on-miss; the cascade: off).
 */
struct SearchPolicy {
  /// Divergence 1, frontiers. HopRun::frontier is vouched for only when the
  /// accounting balanced and the drain completed; with this set, also only
  /// when something was examined (verifyslicegenus's rule; the cascade
  /// lacks that check).
  bool frontierNeedsExamined = false;

  /// Divergence 2, failures. `halt` (verifyslicegenus): a row whose diagram
  /// namer cannot be built is searched on the complement route, with a
  /// warning; a search whose accounting failed, or that examined nothing,
  /// records the outcome "unaccounted"; and the driver halts the process
  /// (exit 2) on a seed-invariant failure (SeedInvariantFailure), on an
  /// accounting failure and on a find below the literature (HopRun::fatal).
  /// `refuse` (the cascade): a namer that cannot be built, like the seed
  /// invariant, throws (the node is refused); the outcome is left as the
  /// search ended; and an accounting failure only marks the hop suspect.
  enum class Failures { halt, refuse };
  Failures failures = Failures::refuse;

  /// Divergences 6 and 10, who judges a find: the search itself, each new
  /// witness against the solver's bounds as of the last solve
  /// (cobordismgraph::upperBoundVia(); verifyslicegenus). It prints the
  /// CONSTRUCTIVE line and forces a checkpoint at the first constructive
  /// find, and a find below the row's literature lower bound ends the
  /// search with HopRun::fatal set. Needs SweepInputs.
  bool judgeInSearch = false;

  /// Divergence 7, signing. `duringSearch` (verifyslicegenus): one witness
  /// per identity across the whole database (SweepInputs::record), each new
  /// one signed off the drain threads (WitnessSigner) and appended at the
  /// 60 s checkpoints. `deferred` (the cascade): one kept surface per
  /// keptKey() within the search, its faces kept (HopRun::kept), signed
  /// later.
  enum class Signing { duringSearch, deferred };
  Signing signing = Signing::deferred;
};

/**
 * What a search that signs during the search, or judges its finds itself,
 * reads and writes beyond its row: verifyslicegenus's (SearchPolicy).
 * Pointers outlive the search.
 */
struct SweepInputs {
  /// Signing::duringSearch: where new witnesses are claimed and published.
  RecordedWitnesses *record = nullptr;
  /// The literature, for each witness's other_candidates (re-derived per
  /// witness, as every reader of the field does) and for judging.
  const cobordismgraph::NameTable *names = nullptr;
  /// judgeInSearch: the solver's bounds as of the last solve.
  const std::unordered_map<std::string, cobordismgraph::Bounds> *bounds = nullptr;
  /// The row's literature interval: judging, and the progress block.
  int literatureLo = 0, literatureHi = 0;
  /// The witnesses' thicken_layers (provenance).
  int thickenLayers = 1;
};

/// What a search writes besides what it returns: verifyslicegenus's
/// options. None by default.
struct SearchOutputs {
  /// The live progress block, and the drain's started/done lines, on
  /// stderr (rowsearch::progressBlock).
  bool progress = false;
  std::optional<std::string> surfaceStats; ///< --surface-stats, appended
  std::optional<std::string> surfaceLog;   ///< --surface-log, rewritten
  /// --rejection-sample-log: the first REJECTION_SAMPLES_PER_REASON
  /// surfaces each gate turns away. Outlives the search.
  std::ostream *rejectionSamples = nullptr;
  std::optional<std::string> selfIntersectionCensus; ///< appended
};

/// A surface a hop kept.
struct KeptSurface {
  farside::OutgoingLink link; ///< far side on T, oriented, and incoming side
  int genus = 0;              ///< tubed genus
  int resolvedVertices = 0;
  std::string farName;        ///< the namer's name; "" for no far side
  std::vector<int> faces;     ///< triangles of the row's thickening
  std::string key;            ///< its dedupe key (see HopSearcher::run())
  /// The witness verifyslicegenus would record for it, every column but the
  /// pair signature, other_candidates and the provenance (source row,
  /// layers, face cap): its subject is the hop's rowName.
  cobordismgraph::Witness witness;
};

/// One search's inputs besides its row.
struct SearchRequest {
  /// The incoming link's name: the subject of every witness, the name the
  /// incoming curves are primed with, and the name in every report line.
  std::string name;
  /// The row's redrawer: Signing::deferred orients each kept surface's
  /// outgoing link on it and keys it by its row components.
  const farside::WitnessRedrawer *row = nullptr;
  SearchShape shape;

  /// When to stop: at this many surfaces satisfying the condition, after
  /// `seconds` of this search, after `sweepSeconds` since `sweepStart`, or
  /// once no new witness has been recorded for `quiescenceSeconds`
  /// (SweepInputs::record). The drain then finishes, unless
  /// `skipDrainOnTimeout`.
  std::optional<long long> surfaceTarget;
  std::optional<double> seconds;
  std::optional<double> sweepSeconds;
  std::chrono::steady_clock::time_point sweepStart{};
  std::optional<double> quiescenceSeconds;
  bool skipDrainOnTimeout = false;

  /// An earlier search's frontier to carry on from (if it is this
  /// search's), and whether to record this one's.
  const SearchFrontier *resume = nullptr;
  bool recordFrontier = true;
  /// Where the pair-signature context is read from and kept
  /// (SurfaceSearch::setPairSigCacheDir()).
  std::optional<std::string> pairSigCacheDir;

  /// Name outgoing curves by their drawing when the searcher has signature
  /// tables (else every boundary by its complement).
  bool diagramNaming = true;
  /// The incoming knot's table name: after its search, its complement goes
  /// into the census under it, when census writes are on
  /// (census::censusUpdates). Unset, or a link: no insert.
  std::optional<std::string> censusName;
  /// Search a row without a collar seed (verifyslicegenus's --cone and
  /// --collar-layers 0, retired with coning); otherwise such a row is
  /// refused.
  bool unseeded = false;

  /// Signing::deferred: polled with each newly kept surface (from drain
  /// threads, one at a time); the search ends when it returns true.
  std::function<bool(const KeptSurface &)> stop;

  SweepInputs sweep;
  SearchOutputs outputs;
};

/// Thrown by HopSearcher::run() before searching when a searchable triangle
/// other than the seed has an edge on the incoming boundary, so found
/// surfaces could change the incoming link. what() is the cascade's
/// refusal; verifyslicegenus halts with its own message.
struct SeedInvariantFailure : std::runtime_error {
  explicit SeedInvariantFailure(size_t touching);
  size_t touching;
};

/// What a hop did.
struct HopRun {
  std::vector<KeptSurface> kept;
  long long accepted = 0;
  std::string accounting;        ///< the `accounting:` body (rowsearch.h)
  std::string accountingFailure; ///< empty iff every surface is accounted for
  std::string outcome;           ///< surface-target, exhausted, timeout, stopped
  double wall = 0;               ///< seconds
  double cpu = 0;                ///< process CPU seconds over the hop
  double setup = 0;              ///< wall before the search: namer, search, seed checks
  double search = 0;             ///< wall of the search itself, drain included
  std::vector<double> rounds;    ///< each IDDFS round's wall (SearchStats::Profile)
  size_t drainTail = 0;          ///< surfaces left for the drain when the search ended
  double drainTailSeconds = 0;   ///< and the wall it took to describe them
  /// The drain's far-side naming: the `diagram naming:` body
  /// (farside::NamingStats::summary()), its times by route, and the slowest
  /// single name, which is what can hold a drain's last thread alone.
  std::string naming;
  double namingDiagramSeconds = 0, namingFallbackSeconds = 0, namingExactSeconds = 0;
  double namingSlowestSeconds = 0;
  /// Where the search stopped (searchfrontier.h), cumulative over the
  /// frontier it resumed; unset when the hop cannot vouch for every surface
  /// in it (accounting failed, or its drain was cut short), so that a later
  /// hop never skips surfaces nobody examined.
  std::optional<SearchFrontier> frontier;
  bool resumed = false;          ///< carried on from the frontier it was given
  std::string resumeRefusal;     ///< why not, when given one it refused

  // What verifyslicegenus's row report reads besides (searchreport.h).
  SearchStats stats;               ///< the search's own, as it returned them
  long long described = 0;         ///< surfaces the drain described
  long long recorded = 0;          ///< new witnesses (or kept surfaces)
  long long otherOrientation = 0;  ///< rejected as another oriented variant
  bool drainSkipped = false;
  bool nothingExamined = false;    ///< RowAccounting::nothingExamined()
  /// The search's frontier as recorded, whether vouched for or not.
  std::optional<SearchFrontier> recordedFrontier;
  double frontierSeconds = 0;
  std::chrono::steady_clock::duration drainTailTime{};
  std::optional<PetalCache::Stats> petalsAtRoots; ///< as root filtering ended
  PetalCache::Stats petals;                        ///< at the search's end
  identify::BoundarySignatureCacheStats boundaryCache;
  identify::RecognitionCacheStats recognitionBefore; ///< as the search began
  identify::RecognitionCacheStats recognitionAfter;  ///< and as it ended
  std::pair<long long, long long> censusWritesBefore; ///< census::insertCounts()
  std::pair<long long, long long> censusWritesAfter;
  bool diagramNamed = false;       ///< a DiagramNamer named the outgoing curves
  long long nonPlanar = 0;         ///< its non-planar drawings
  long long pairSigsSigned = 0;    ///< Signing::duringSearch
  long long pairSigMillis = 0;
  /// The signer's context: when it was ready (seconds after the search
  /// began; 0 if never needed), whether the cache held it, and what
  /// waiting for the last signatures added after the drain.
  double pairSigContextSeconds = 0;
  bool pairSigContextLoaded = false;
  double pairSigFinishSeconds = 0;
  /// The `runs` of the frontier this search was offered to resume, if any.
  std::optional<unsigned> resumeOfferedRuns;
  bool linkingAudit = false;       ///< petal linking numbers were audited
  /// judgeInSearch: why the search found something impossible (a witness
  /// below the literature lower bound); empty when it did not.
  std::string fatal;
};

class HopSearcher {
public:
  /// `signatures` and `exact` name far sides as verifyslicegenus names them
  /// (farside::DiagramNamer), which is what surfaces are deduplicated by.
  /// Both must outlive the searcher. Every hop's exact names use one set of
  /// table caches (`exactCaches`, or the searcher's own when null), so what
  /// naming learns about the tables -- the HOMFLY index above all -- is
  /// built once, not once per hop. The cascade's: SearchPolicy's defaults.
  HopSearcher(const farside::SignatureTable &signatures,
              const exactnaming::ExactTables *exact, HopShape shape,
              unsigned threads,
              std::shared_ptr<exactnaming::TableCaches> exactCaches = nullptr);

  /// Any caller's: its own `policy`. Without `signatures`, every boundary is
  /// named by its complement. Exact names use this searcher's own caches.
  HopSearcher(const farside::SignatureTable *signatures,
              const exactnaming::ExactTables *exact, SearchPolicy policy,
              unsigned threads);

  /**
   * Searches `row`'s thickening, seeded with its collar, under `proper`,
   * until `surfaceTarget` surfaces satisfy it (or `seconds` pass; the drain
   * then finishes), and keeps one surface per key: verifyslicegenus's
   * witness identity (far-side name, its component count, genus, whether
   * tubed or resolved) together with which row components and how many
   * far-side curves each surface component carries. The second part is
   * what profiles read (partitions), and the witness identity alone
   * collapses it.
   *
   * `stop`, when given, is polled with each newly kept surface (from drain
   * threads, one at a time) and ends the search when it returns true.
   *
   * \throws std::runtime_error when the row's seed invariant fails.
   */
  /// `resume`, if set, is an earlier hop's frontier on this row: the search
  /// carries on from it (if it is this search's; see HopRun::resumeRefusal).
  /// `surfaceTarget` is the search's breadth, so a resumed hop adds only the
  /// surfaces beyond its frontier's (SearchCallbacks::surfaceTarget).
  /// `censusName`, if any, is SearchRequest::censusName.
  HopRun run(const farside::WitnessRedrawer &row, const std::string &rowName,
             long long surfaceTarget, double seconds,
             const std::function<bool(const KeptSurface &)> &stop = {},
             const SearchFrontier *resume = nullptr,
             std::optional<std::string> censusName = std::nullopt) const;

  /**
   * THE search: one search from `rb`'s incoming link, as `request` and this
   * searcher's SearchPolicy say. It builds the search over `rb`'s
   * thickening (seeded with its collar), names boundaries, checks the seed
   * invariant, runs the search with the watchdog, gates every surface the
   * drain describes (preconditions.h), records what passes (SearchPolicy::
   * Signing), and accounts for every surface. run(row, ...) above is the
   * cascade's hop through it.
   *
   * \throws SeedInvariantFailure when the seed invariant fails;
   * std::runtime_error for a row with no seed (unless `request.unseeded`),
   * and (SearchPolicy::Failures::refuse) when the diagram namer cannot be
   * built; a signer's failure (Signing::duringSearch).
   */
  HopRun run(const rowsearch::RowBuild &rb, const SearchRequest &request) const;

private:
  const farside::SignatureTable *signatures_;
  const exactnaming::ExactTables *exact_;
  HopShape shape_;
  SearchPolicy policy_;
  unsigned threads_;
  std::shared_ptr<exactnaming::TableCaches> exactCaches_;
};

} // namespace cascade

#endif // SURFER_COBOUND_SEARCH_H
