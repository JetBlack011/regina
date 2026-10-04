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
#include <unordered_set>
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
 *  that ends it at its surface target or deadline (RowWatchdog). Every
 *  caller's searches behave alike (phase 4(b) unified them).
 *
 *  verifyslicegenus searches each table row through it (its finds written to
 *  the row's pending file and signed at the run's end), and the
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

namespace search {

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
    std::optional<double> seconds;

    bool any() const { return surfaceTarget || seconds; }
};

/**
 * Polls every 200 ms, while the search and its drain run, and calls
 * `endRow(why)` at most once: "surface-target" or "timeout".
 * stop() wakes it at once rather than waiting out its poll.
 * The surface target is checked before the clocks, so a row reaching it in
 * the same tick as a deadline records "surface-target": the two mean
 * different things to anyone later reading the negative.
 *
 * No thread is started when no limit is set.
 */
class SearchWatchdog {
public:
    SearchWatchdog(WatchdogLimits limits, std::function<void(const char *)> endSearch);
    ~SearchWatchdog();
    SearchWatchdog(const SearchWatchdog &) = delete;
    SearchWatchdog &operator=(const SearchWatchdog &) = delete;

    /** The live count of surfaces satisfying the boundary condition, from
     *  SearchCallbacks::onProgress (the only place it is handed out). */
    void publishSatisfying(long long count) {
        satisfying_.store(count, std::memory_order_relaxed);
    }
    /** Stops polling and joins; idempotent. Call once the search returns. */
    void stop();

private:
    WatchdogLimits limits_;
    std::function<void(const char *)> endSearch_;
    std::atomic<long long> satisfying_{0};
    std::atomic<bool> done_{false};
    std::mutex wakeMutex_;
    std::condition_variable wake_;
    std::thread thread_;
};

} // namespace search

namespace search {

/// A goal run's search shape: what its profile line prints, and every
/// search of the run searches with (searchShape()). No field that decides
/// what a search visits or accepts has a default here: the config is the
/// one source of every value (plan, "search-shape defaults"; with a goal,
/// the campaign's shape: atlas tools/orchestrate/hosts.conf [campaign]), and
/// searchShape() refuses a shape with one unset.
struct RunShape {
  /// Thicken and collar layers of a hop's row (both), and the profile's
  /// layers=. A master row's layers are its witnesses' own.
  std::optional<int> layers;
  std::optional<long long> maxFaces;
  std::optional<unsigned> iddfsIterations;
  std::optional<long long> iddfsStart;
  std::optional<long long> iddfsStep;
  std::optional<long long> rootBudgetStart;
  std::optional<long long> rootBudgetGrowth;
  /// Whether surfaces whose only self-intersections are unlinked count (plan
  /// divergence 3).
  std::optional<bool> resolveUnlinked;
  /// The searches' resources (SurfaceSearchLimits): they decide no order,
  /// and unset is the library's default. With a goal the config sets
  /// hosts.conf's: at the library's pending cap (500,000) a search pauses to
  /// drain its queue, which a campaign row never does.
  std::optional<size_t> pendingSurfaceCap;
  std::optional<size_t> petalCacheLimit;
  std::optional<size_t> boundarySignatureCacheLimit;
};

/// How one search runs: the arguments of SurfaceSearch::search() that fix
/// its traversal and what it accepts, and the limits of its shared
/// structures. Every caller fills all of it (searchShape() for a hop).
struct SearchShape {
  BoundaryCondition condition = BoundaryCondition::proper;
  unsigned iddfsIterations = 0;
  long long iddfsStep = 0;
  std::optional<long long> iddfsStart;
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
SearchShape searchShape(const RunShape &shape);

/** The searched link's literature interval, as a search reports it: the
 *  CONSTRUCTIVE line and the progress block (verifyslicegenus's). */
struct LiteratureInterval {
  int literatureLo = 0, literatureHi = 0;
};

/// What a search writes besides what it returns: verifyslicegenus's
/// options. None by default.
struct SearchOutputs {
  /// The live progress block, and the drain's started/done lines, on
  /// stderr (search::progressBlock).
  bool progress = false;
  std::optional<std::string> surfaceStats; ///< --surface-stats, appended
  std::optional<std::string> surfaceLog;   ///< --surface-log, rewritten
  /// --rejection-sample-log: the first REJECTION_SAMPLES_PER_REASON
  /// surfaces each gate turns away. Outlives the search.
  std::ostream *rejectionSamples = nullptr;
  std::optional<std::string> selfIntersectionCensus; ///< appended
};

/// What the cobordism graph made of a find (SearchRequest::judge).
struct FindJudgement {
  /// A contradiction the find brought into the graph (its gates): the search
  /// ends, and HopRun::fatal says why. Empty, normally.
  std::string contradiction;
  /// The searched link's least genus, once a proof resting on no literature
  /// value reaches its literature lower bound (divergence 10): the first
  /// time, the search announces it (the CONSTRUCTIVE line, `-- ACHIEVED` in
  /// the progress block) and checkpoints at once.
  std::optional<int> constructive;
};

/// A surface a hop kept.
struct KeptSurface {
  outgoing::OutgoingLink link; ///< far side on T, oriented, and incoming side
  int genus = 0;              ///< tubed genus
  int resolvedVertices = 0;
  std::string outgoingName;        ///< the namer's name; "" for no far side
  std::vector<int> faces;     ///< triangles of the row's thickening
  std::string key;            ///< its dedupe key (see HopSearcher::run())
  /// The cobordism as the database will record it, every column but the
  /// pair signature and other_candidates (both made when it is signed): its
  /// subject is the search's name.
  cobordisms::Cobordism cobordism;
};

/// One search's inputs besides its row.
struct SearchRequest {
  /// The incoming link's name: the subject of every witness, the name the
  /// incoming curves are primed with, and the name in every report line.
  std::string name;
  /// The row's redrawer (required): each kept surface's outgoing link is
  /// oriented on it and keyed by its row components.
  const outgoing::OutgoingReader *reader = nullptr;
  SearchShape shape;

  /// When to stop: at this many surfaces satisfying the condition, or
  /// after `seconds` of this search (a wall-clock backstop). The drain then
  /// finishes: every surface the search found is examined.
  std::optional<long long> surfaceTarget;
  std::optional<double> seconds;

  /// An earlier search's frontier to carry on from (if it is this
  /// search's, and its finds are signed: see `runDirectory`), and whether
  /// to record this one's.
  const SearchFrontier *resume = nullptr;
  bool recordFrontier = true;
  /// The run's work directory: every pending file under it is signed by
  /// this run's sign step (or by `sign` after a kill), so a frontier whose
  /// pending file lies there may be resumed before it is signed. Any other
  /// frontier's pending file must be signed at least as far as the frontier
  /// recorded (plan divergence 1), or the resume is refused and says why.
  std::optional<std::string> runDirectory;
  /// Where the pair-signature context is read from and kept
  /// (SurfaceSearch::setPairSigCacheDir()).
  std::optional<std::string> pairSigCacheDir;

  /// The incoming knot's table name: after its search, its complement goes
  /// into the census under it, when census writes are on
  /// (census::censusUpdates). Unset, or a link: no insert.
  std::optional<std::string> censusName;

  /// The search's pending file (cobordisms/pending.h, PendingWriter): every
  /// surface it keeps is appended there with its faces -- fsynced once a
  /// minute while the search and its drain run, at the first constructive
  /// find, and the rest when the search ends -- for `sign` to sign into the
  /// database (plan divergence 7). Unset: nothing is written (a goal run
  /// without a witness store keeps its finds in memory only).
  std::optional<std::string> pending;
  /// The row's PD and layers, as a pending line and a cobordism's
  /// provenance record them (`sign` rebuilds the row from them).
  std::string incomingPD;
  int layers = 2;
  /// The identities of the database the run loaded (RecordedWitnesses): a
  /// find with one of them is a duplicate, one cobordism per identity
  /// across the database. Unset: no database loaded.
  const std::unordered_set<std::string> *knownIdentities = nullptr;

  /// Polled with each newly kept surface (from drain threads, one at a
  /// time); the search ends when it returns true.
  std::function<bool(const KeptSurface &)> stop;
  /// The cobordism graph's judgement of each find, as it is kept (plan
  /// divergences 6 and 10): called on a thread of its own, one find at a
  /// time, in the order they were kept, so the drain never waits for it;
  /// the search waits for the last one before it returns. A contradiction
  /// ends the search; the first constructive find is announced. Needs `row`
  /// (each find's outgoing link is read on it). Unset: the caller's graph
  /// judges the search's finds once it returns (a goal run's, which spans
  /// every search of the run).
  std::function<FindJudgement(const KeptSurface &)> judge;

  LiteratureInterval literature;
  SearchOutputs outputs;
};

/// Thrown by HopSearcher::run() before searching when a searchable triangle
/// other than the seed has an edge on the incoming boundary, so found
/// surfaces could change the incoming link: an impossible state, on which
/// every driver halts once what it found is written (plan divergence 2).
struct SeedInvariantFailure : std::runtime_error {
  explicit SeedInvariantFailure(size_t touching);
  size_t touching;
};

/// Thrown by HopSearcher::run() for a row it will not search: no collar
/// seed, or a far-side namer that cannot be built on it (plan divergence 2:
/// refused, never searched by a fallback route). The cascade refuses the
/// node; verifyslicegenus records the row as a build failure and goes on.
struct SearchRefused : std::runtime_error {
  using std::runtime_error::runtime_error;
};

/// What a hop did.
struct SearchResult {
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
  /// (linknaming::NamingStats::summary()), its times by route, and the slowest
  /// single name, which is what can hold a drain's last thread alone.
  std::string naming;
  double namingDiagramSeconds = 0, namingFallbackSeconds = 0, namingExactSeconds = 0;
  double namingSlowestSeconds = 0;
  /// Where the search stopped (searchfrontier.h), cumulative over the
  /// frontier it resumed, with its pending file and that file's fsynced
  /// length; unset when the search cannot vouch for every surface in it
  /// (its accounting failed, nothing was examined, or its drain was cut
  /// short: plan divergence 1), so that a later search never skips surfaces
  /// nobody examined.
  std::optional<SearchFrontier> frontier;
  bool resumed = false;          ///< carried on from the frontier it was given
  std::string resumeRefusal;     ///< why not, when given one it refused

  // What verifyslicegenus's row report reads besides (searchreport.h).
  SearchStats stats;               ///< the search's own, as it returned them
  long long described = 0;         ///< surfaces the drain described
  long long recorded = 0;          ///< kept surfaces
  /// Distinct cobordism identities among them: the outcome line's "N new
  /// witnesses" (the database gains these when they are signed).
  long long newCobordisms = 0;
  /// The pending file, and its fsynced length when the search ended; -1
  /// when none was written.
  std::string pendingPath;
  long long pendingBytes = -1;
  long long otherOrientation = 0;  ///< rejected as another oriented variant
  bool drainSkipped = false;
  bool nothingExamined = false;    ///< RowAccounting::nothingExamined()
  /// Surfaces in a state that cannot occur (RowAccounting::impossible()):
  /// unlike an imbalance, which ends only its search, a halt (divergence 2).
  long long impossible = 0;
  /// The search's frontier as recorded, whether vouched for or not.
  std::optional<SearchFrontier> recordedFrontier;
  double frontierSeconds = 0;
  std::chrono::steady_clock::duration drainTailTime{};
  std::optional<PetalCache::Stats> petalsAtRoots; ///< as root filtering ended
  PetalCache::Stats petals;                        ///< at the search's end
  namecache::BoundarySignatureCacheStats boundaryCache;
  complement::RecognitionCacheStats recognitionBefore; ///< as the search began
  complement::RecognitionCacheStats recognitionAfter;  ///< and as it ended
  std::pair<long long, long long> censusWritesBefore; ///< census::insertCounts()
  std::pair<long long, long long> censusWritesAfter;
  bool diagramNamed = false;       ///< a DiagramNamer named the outgoing curves
  long long nonPlanar = 0;         ///< its non-planar drawings
  /// Pair signatures made during the search: none since divergence 7 (they
  /// are made when the run signs its pending files).
  long long pairSigsSigned = 0;
  long long pairSigMillis = 0;
  /// The `runs` of the frontier this search was offered to resume, if any.
  std::optional<unsigned> resumeOfferedRuns;
  bool linkingAudit = false;       ///< petal linking numbers were audited
  /// Why the search found something impossible (the cobordism graph's
  /// gates: SearchRequest::judge); empty when it did not.
  std::string fatal;
};

class Searcher {
public:
  /// `signatures` and `exact` name far sides as verifyslicegenus names them
  /// (outgoing::DiagramNamer), which is what surfaces are deduplicated by.
  /// Both must outlive the searcher. Every hop's exact names use one set of
  /// table caches (`exactCaches`, or the searcher's own when null), so what
  /// naming learns about the tables -- the HOMFLY index above all -- is
  /// built once, not once per hop. The cascade's.
  Searcher(const linknaming::SignatureTable &signatures,
              const linknaming::ExactTables *exact, RunShape shape,
              unsigned threads,
              std::shared_ptr<linknaming::TableCaches> exactCaches = nullptr);

  /// Any caller's. Without `signatures`, every boundary is named by its
  /// complement. Exact names use this searcher's own caches.
  Searcher(const linknaming::SignatureTable *signatures,
              const linknaming::ExactTables *exact, unsigned threads);

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
  SearchResult run(const outgoing::OutgoingReader &reader, const std::string &subject,
             long long surfaceTarget, double seconds,
             const std::function<bool(const KeptSurface &)> &stop = {},
             const SearchFrontier *resume = nullptr,
             std::optional<std::string> censusName = std::nullopt) const;

  /// The request run(row, rowName, ...) makes: a hop's search of `row`
  /// (searchShape() of this searcher's HopShape, its frontier always
  /// recorded), for a caller that adds to it (its pending file, say).
  SearchRequest requestFor(const outgoing::OutgoingReader &reader, const std::string &subject,
                           long long surfaceTarget, double seconds) const;

  /**
   * THE search: one search from `rb`'s incoming link, as `request` says.
   * It builds the search over `rb`'s thickening (seeded with its collar),
   * names boundaries, checks the seed invariant, runs the search with the
   * watchdog, gates every surface the drain describes (preconditions.h),
   * keeps what passes (one per keptKey(), written to the pending file), and
   * accounts for every surface. run(row, ...) above is the cascade's hop
   * through it.
   *
   * \throws SeedInvariantFailure when the seed invariant fails;
   * SearchRefused for a row with no seed, or
   * whose diagram namer cannot be built; std::runtime_error when the
   * pending file cannot be written at the search's end.
   */
  SearchResult run(const search::IncomingThickening &rb, const SearchRequest &request) const;

private:
  const linknaming::SignatureTable *signatures_;
  const linknaming::ExactTables *exact_;
  RunShape shape_;
  unsigned threads_;
  std::shared_ptr<linknaming::TableCaches> exactCaches_;
};

} // namespace search

#endif // SURFER_COBOUND_SEARCH_H
