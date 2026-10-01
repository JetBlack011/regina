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
#include <string>
#include <thread>
#include <vector>

#include "surfer/submanifold/submanifold.h"
#include "cobound/cobordisms/cobordism.h"
#include "linknaming/tables.h"
#include "cobound/outgoing/outgoinglink.h"
#include "cobound/outgoing/outgoingnamer.h"
#include "cobound/outgoing/farsideredraw.h"
#include "surfer/enumeration/searchstrategy.h"

/*! \file utils/surfer/cobound/search/search.h
 *  \brief One search from a link: which BoundaryCondition it runs under
 *  (conditionFor()), and the watchdog that ends it at its surface target or
 *  deadline (RowWatchdog).
 *
 *  HopSearcher runs one hop of the cascade in process: the row's own
 *  thickening (the one its HopAssembler certified and reads from), the
 *  campaign's search shape, and the gates and accounting verifyslicegenus
 *  applies (preconditions.h), with the accepted surfaces read straight off
 *  the search instead of written out as pair signatures and read back.
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
/// [campaign], root budget 840 from c5 on). The layers are the row's.
struct HopShape {
  long long maxFaces = 5;
  unsigned iddfsIterations = 2;
  long long iddfsStart = 4;
  long long iddfsStep = 1;
  long long rootBudgetStart = 840;
  long long rootBudgetGrowth = 2;
  bool resolveUnlinked = true;
  /// hosts.conf's per-host limits, identical on every host. The pending
  /// cap in particular: at the binary's default (500,000) a search pauses
  /// to drain its queue, which a campaign row never does.
  size_t pendingSurfaceCap = 20'000'000;
  size_t petalCacheLimit = 12'000'000;
  size_t boundarySignatureCacheLimit = 1'000'000;
  /// Process-wide (identify::recognitionCacheLimit); set by the driver.
  size_t recognitionCacheLimit = 1'500'000;
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
};

class HopSearcher {
public:
  /// `signatures` and `exact` name far sides as verifyslicegenus names them
  /// (farside::DiagramNamer), which is what surfaces are deduplicated by.
  /// Both must outlive the searcher. Every hop's exact names use one set of
  /// table caches (`exactCaches`, or the searcher's own when null), so what
  /// naming learns about the tables -- the HOMFLY index above all -- is
  /// built once, not once per hop.
  HopSearcher(const farside::SignatureTable &signatures,
              const exactnaming::ExactTables *exact, HopShape shape,
              unsigned threads,
              std::shared_ptr<exactnaming::TableCaches> exactCaches = nullptr);

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
  HopRun run(const farside::WitnessRedrawer &row, const std::string &rowName,
             long long surfaceTarget, double seconds,
             const std::function<bool(const KeptSurface &)> &stop = {},
             const SearchFrontier *resume = nullptr) const;

private:
  const farside::SignatureTable &signatures_;
  const exactnaming::ExactTables *exact_;
  HopShape shape_;
  unsigned threads_;
  std::shared_ptr<exactnaming::TableCaches> exactCaches_;
};

} // namespace cascade

#endif // SURFER_COBOUND_SEARCH_H
