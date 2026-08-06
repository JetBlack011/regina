//
//  embeddingsearch.h
//
//  Created by John Teague on 07/15/2026.
//

#ifndef SURFACESEARCH_H

#define SURFACESEARCH_H

#include <map>
#include <mutex>
#include <thread>

#include "embeddingsearch.h"
#include "embeddedsubmanifold.h"
#include "identifycomplement.h"

/**
 * Everything needed to describe one found surface when no boundary-link
 * data is being tracked (BoundaryCondition::all/closed), handed to
 * SurfaceSearchCallbacks::onSurfaceFound.
 */
struct SurfaceFoundInfo {
  bool orientable; /**< Whether the surface is orientable. */
  int genus;
      /**< Meaningful ONLY when `connected` is true: chi = 2 - 2g - punctures
           (chi = 2 - g - punctures if non-orientable) is a statement about
           one connected surface, so this field carries whatever that
           formula evaluates to even when the surface is disconnected --
           e.g. two disjoint genus-0 1-punctured pieces (two discs) comes
           out as genus -1, not a real genus of anything. There's no single
           number that would correctly summarize a disconnected surface's
           genus for every purpose, so this deliberately isn't "fixed" to
           produce one; callers that need per-component detail for a
           disconnected find must decompose the surface themselves. Always
           check `connected` before reading this field. */
  int punctures;
      /**< The surface's total number of punctures (boundary components),
           across every connected piece if disconnected -- this count is
           well-defined regardless of `connected`. */
  bool connected; /**< Whether the surface is topologically connected -- see `genus`. */
  long long triangleCount; /**< The surface's number of triangles. */
  BoundaryCondition mostRestrictive;
      /**< The most restrictive BoundaryCondition this specific surface
           satisfies -- not necessarily the one the whole search was run
           under. Only isClosed()/isProper() are checked (both O(1)), not
           the O(ambient triangulation size)
           boundaryComponentsMapInjectively(), so this is never
           "connected". */
  std::function<std::string()> capturePairSig;
      /**< Computes (or returns the already-cached) reproducible
           isomorphism signature for this surface -- see
           EmbeddedSubmanifold::pairSig(). Only set (a non-empty
           std::function) when SurfaceSearchLimits::capturePairSig is
           true; call `static_cast<bool>(info.capturePairSig)` to check
           before calling.

           Deliberately lazy, not a precomputed string: pairSig() searches
           an automorphism group over the WHOLE ambient triangulation, so
           it's genuinely expensive (profiling verifyslicegenus showed
           this dominating >95% of its post-search boundary-processing
           time when every found surface's pairSig was captured
           unconditionally, even though almost none of them were ever
           read) -- callers should only invoke this for the rare surface
           they've actually decided is worth recording, never
           unconditionally for every found surface.

           \warning Only valid to call synchronously, within the callback
           this SurfaceFoundInfo/SurfaceBoundaryInfo was passed to: it
           closes over the embedding by reference, which the caller may
           mutate (removeFace()'s unwind) or reuse for the next entry
           immediately after the callback returns. Safe to call more than
           once within that window (pairSig() itself is cached). */
};

/**
 * SurfaceFoundInfo plus the boundary-link description, handed to
 * SurfaceSearchCallbacks::onSurfaceBoundaryProcessed once boundary links
 * are being tracked (BoundaryCondition::proper/connected).
 */
/**
 * One ambient boundary component's identification, as computed by
 * describeBoundary_(): the per-curve names of every one of the surface's
 * own boundary curves landing on that component, plus (when there's more
 * than one) the whole-link name of all of them drilled together.
 */
struct BoundaryComponentNames {
  size_t component; /**< Ambient boundary component index. */
  std::vector<std::string> curveNames;
      /**< The identify()'d name of each of this component's curves,
           individually, in describeBoundary_'s order. */
  std::optional<std::string> linkName;
      /**< identify()'d name of all of this component's curves drilled
           together (identify::identify() on the whole Link) -- the same
           value describeBoundary_ already folds into the parenthetical in
           boundaryDescription. Set if and only if curveNames.size() > 1;
           for a single curve, curveNames.front() already names it and no
           separate whole-link identification is meaningful. */
};

struct SurfaceBoundaryInfo : SurfaceFoundInfo {
  std::string boundaryDescription;
      /**< A formatted description of the surface's boundary; empty if
           the surface turned out closed. */
  std::vector<BoundaryComponentNames> boundaryComponents;
      /**< The same (ambient boundary component, curve names) grouping
           boundaryDescription is flattened from -- one entry per ambient
           boundary component the surface's boundary touches. Unlike
           boundaryDescription, safe to consume programmatically even when
           a component holds more than one curve (which boundaryDescription
           can't be reliably split back apart in that case, since curve
           names within a component and separate components are both
           comma-joined). */
  std::function<std::vector<std::pair<size_t, std::vector<OrientedCurve>>>()>
      captureOrientedBoundaryLinks;
      /**< Lazy, as capturePairSig above (and for the same reason: closes
           over the embedding by reference, only valid to call synchronously
           within the callback this SurfaceBoundaryInfo was passed to) --
           calls KnottedSurface::orientedBoundaryLinks(). Always populated
           (unlike capturePairSig, not gated behind a limits_ flag: no
           search currently runs often enough for this to need capping),
           but still a closure rather than a precomputed field, since most
           callers never need it (only the orientation-tracking gate in
           verifyslicegenus.cpp's own "search side, fully witnessed" branch
           does) and it costs a boundary-facet walk to compute. */
};

/**
 * SurfaceSearch's own callbacks: everything SearchCallbacks offers for
 * the main DFS phase, plus per-surface callbacks and three more for the
 * deferred boundary-link post-processing phase
 * (processRemainingSurfaceBoundaries) that runs after it.
 */
struct SurfaceSearchCallbacks : SearchCallbacks {
  std::function<void(const SurfaceFoundInfo &)> onSurfaceFound;
      /**< Fired once per found surface when boundary links are not
           being tracked (BoundaryCondition::all/closed), from the DFS
           hot path itself -- deliberately cheap (see
           SurfaceFoundInfo::mostRestrictive), since it can fire once per
           DFS node visited under BoundaryCondition::all. */
  std::function<void(const SurfaceBoundaryInfo &)> onSurfaceBoundaryProcessed;
      /**< Fired once per found surface when boundary links are being
           tracked (BoundaryCondition::proper/connected), from the
           deferred boundary-link processing phase. Mutually exclusive
           with onSurfaceFound: exactly one of the two fires per found
           surface. */
  std::function<void(size_t total, unsigned numThreads)>
      onBoundaryProcessingStarted;
      /**< Fired once, when boundary-link post-processing begins. */
  std::function<void(size_t processed, size_t total,
                     std::chrono::steady_clock::duration elapsed)>
      onBoundaryProcessingProgress;
      /**< Fired roughly once a second while boundary-link
           post-processing runs. */
  std::function<void(size_t total, std::chrono::steady_clock::duration elapsed)>
      onBoundaryProcessingComplete;
      /**< Fired once, when boundary-link post-processing finishes. `total`
           is always the full batch size (processBatchParallel_ always
           drains it completely once started), letting a caller report an
           overall average rate. */
  std::function<void(size_t queueSize, size_t cap)> onQueueDrainPause;
      /**< Fired when the live pending-surface queue (see
           SurfaceSearch::PendingSurfaceBatch) exceeds its cap and a DFS
           worker thread pauses the *entire* search -- every other worker
           blocks in its own next tryAdd() (see InterruptiblePredicate) --
           to drain the queue down to empty before any of them resume; see
           SurfaceSearch::ThreadHook::onFlush(). `queueSize`/`cap` are a
           snapshot taken right as the pause begins. */
  std::function<void(bool interrupted, size_t remaining)> onQueueDrainResume;
      /**< Fired once the drain triggered by onQueueDrainPause ends, either
           because the queue reached empty (`interrupted` false,
           `remaining` 0) or because Ctrl+C fired mid-drain (`interrupted`
           true, `remaining` the number of surfaces handed off to the
           normal post-search processing phase instead of being finished
           here). Search resumes immediately after this fires, unless
           `interrupted`. */
};

/**
 * A search over KnottedSurfaces that additionally recognizes each found
 * surface's boundary curves per ambient boundary component, tallying
 * which surface types have which per-component boundary identifications.
 *
 * This has its own search() (rather than delegating to
 * EmbeddingSearch<4,2>::search()) so that surface-type tallying and
 * boundary-link batching can be interleaved into the same DFS pass and
 * the same once-a-second progress report.
 */
/** Counters for how much recomputation the pending-surface queue's cap/backpressure is costing/avoiding. */
struct PendingSurfaceQueueStats {
  size_t currentSize = 0;
  size_t cap = 0;
  size_t peakSize = 0;
  /** How many onFlush() calls (across every DFS worker thread) found the queue over cap and helped drain it. */
  long long producerDrainEvents = 0;
  /** Total entries drained by a DFS worker thread helping out, as opposed to backgroundDrainLoop_'s own dedicated thread. */
  long long producerDrainedEntries = 0;
};

/**
 * Thresholds for every otherwise-unbounded structure SurfaceSearch shares
 * across its worker threads (see SurfaceSearch::configureLimits()) --
 * without these, a long-running search can grow memory without bound,
 * since candidates/finds vastly outnumber any small fixed set of distinct
 * outcomes. identify::recognitionCacheLimit (identifycomplement.h) is
 * deliberately not here: that cache is global/process-wide, not owned by any one
 * SurfaceSearch, so it's set directly (e.g. in surfer.cpp's main(), the
 * same way --no-simplify sets simplifyComplements).
 */
struct SurfaceSearchLimits {
  /** See SurfaceSearch::PendingSurfaceBatch::setCap(). */
  size_t pendingSurfaceCap = 500'000;
  /** See PetalCache::setClearThreshold(). */
  size_t petalCacheLimit = 2'000'000;
  /** See BoundarySignatureCache's clearThreshold constructor parameter. */
  size_t boundarySignatureCacheLimit = identify::BoundarySignatureCache::DEFAULT_CLEAR_THRESHOLD;
  /** See SurfaceSearch::LinkBoundaryTally::setCap(). */
  size_t boundaryTallyCap = 1'000'000;
  bool capturePairSig = false;
      /**< Whether onSurfaceFound/onSurfaceBoundaryProcessed populate
           SurfaceFoundInfo::pairSig (via the found embedding's own
           EmbeddedSubmanifold::pairSig()). Off by default so existing
           callers (surfer.cpp) pay no extra cost -- pairSig() runs a
           nontrivial automorphism search. */
};

class SurfaceSearch : public EmbeddingSearch<4, 2> {
public:
  /** See KnottedSurface::SurfaceTypeKey. */
  using SurfaceTypeKey = KnottedSurface::SurfaceTypeKey;

private:
  /** A thread-safe tally of how many surfaces of each SurfaceTypeKey were found. */
  class SurfaceTypeTally {
  private:
    mutable std::mutex mutex_;
    std::map<SurfaceTypeKey, long long> counts_;

  public:
    /** Merges `local`'s counts into this tally. */
    void merge(std::map<SurfaceTypeKey, long long> &local);

    /** Returns a human-readable summary of the current tally. */
    std::string summary() const;
  };

  /**
   * A thread-safe queue of surfaces (as face-index lists) awaiting
   * boundary-link processing.
   *
   * Bounded via cap_ (see setCap()): merge() itself never blocks or drops
   * anything (it's the data-structure primitive), but SurfaceSearch's
   * ThreadHook::onFlush() checks size() against cap() after merging and,
   * if over cap, has the calling DFS worker thread pop and process every
   * entry itself -- draining the queue down to empty, not just back under
   * cap -- before returning to the search. See onFlush() for why this is
   * safe (no double-processing, no new deadlock risk) and how it gives
   * elastic, multi-threaded draining exactly under pressure while
   * backgroundDrainLoop_'s single dedicated thread keeps running
   * regardless of load.
   */
  class PendingSurfaceBatch {
  private:
    mutable std::mutex mutex_;
    std::vector<std::vector<int>> pending_;
    size_t cap_ = 20'000;
    size_t peakSize_ = 0;
    long long producerDrainEvents_ = 0;
    long long producerDrainedEntries_ = 0;

  public:
    /** Merges `local`'s entries into this batch. */
    void merge(std::vector<std::vector<int>> &local);

    /** Removes and returns every currently pending entry. */
    std::vector<std::vector<int>> drain();

    /**
     * Removes and returns up to `maxCount` entries in one lock
     * acquisition (fewer if that's all that's available, none if
     * empty).
     *
     * Used by the background drain loop (see backgroundDrainLoop_) so
     * it never claims more than a small, bounded amount of work at
     * once, keeping it quick to hand the rest off once the search ends.
     */
    std::vector<std::vector<int>> popSome(size_t maxCount);

    /** The number of entries currently pending. */
    size_t size() const;

    /** The current cap; see setCap(). */
    size_t cap() const;

    /**
     * Sets the size this queue is considered "over cap" past -- see
     * SurfaceSearch::ThreadHook::onFlush() for what actually happens once
     * it is. Call before any worker thread starts searching.
     */
    void setCap(size_t cap);

    /** Records that a DFS worker thread (not backgroundDrainLoop_) drained `entriesDrained` entries itself while helping under backpressure. */
    void recordProducerDrain(size_t entriesDrained);

    /** A snapshot of this queue's current size/cap/peak/backpressure counters. */
    PendingSurfaceQueueStats stats() const;
  };

  /**
   * A thread-safe tally, per boundary descriptor (see describeBoundary_
   * -- one string capturing the identification of every curve in every
   * ambient boundary component a surface's boundary touches), of how
   * many surfaces of each type were found with that exact descriptor.
   */
  class LinkBoundaryTally {
  private:
    mutable std::mutex mutex_;
    std::unordered_map<std::string, std::map<SurfaceTypeKey, long long>>
        descriptorSurfaceTypes_;
    std::vector<std::string> order_;
        /**< Distinct descriptors, in the order record() first saw each one
             -- lets summary(maxRecent) show only the most recently
             encountered descriptors without an unbounded line count. */
    size_t descriptorCap_ = 100'000;
    long long refusedCount_ = 0;
        /**< Distinct descriptors seen but not recorded because
             descriptorSurfaceTypes_ was already at descriptorCap_. This is
             a reporting-only structure, so a destructive clear (as
             PetalCache/recognitionCache/BoundarySignatureCache use) would
             silently drop user-facing summary data -- refusing brand-new
             keys past the cap while still tallying into already-known ones
             is the safer trade-off here. */

  public:
    /** Records one surface of `type` found with boundary `descriptor`. */
    void record(const std::string &descriptor, const SurfaceTypeKey &type);

    /** Sets the distinct-descriptor cap past which record() stops admitting brand-new descriptors (still tallying into already-known ones). Call before any worker thread starts searching. */
    void setCap(size_t cap);

    /**
     * Returns a human-readable summary of the current tally.
     *
     * With `maxRecent` unset, every distinct descriptor is listed,
     * alphabetically -- suitable for a final, one-shot summary. With
     * `maxRecent` set, only the `maxRecent` most recently *first*
     * encountered descriptors are listed (in that discovery order, not
     * alphabetically) -- suitable for a live, repeatedly-redrawn progress
     * report, so its length stays bounded regardless of how many distinct
     * boundaries a long search eventually turns up.
     */
    std::string summary(std::optional<size_t> maxRecent = std::nullopt) const;
  };

  /**
   * Per-worker-thread hook passed to runSearch_(): tallies surface types
   * and (when boundary links are wanted) batches each find's face
   * indices for later link identification, merging into the owning
   * SurfaceSearch's shared accumulators whenever the harness flushes
   * this thread's counters.
   */
  class ThreadHook : public RunSearchThreadHook<4, 2> {
    SurfaceSearch &owner_;
    SurfaceTypeTally &tally_;
    bool wantLinks_;
    const SurfaceSearchCallbacks &callbacks_;
    std::map<SurfaceTypeKey, long long> localTypeCounts_;
    std::vector<std::vector<int>> localPending_;

    std::optional<KnottedSurface> helperEmbedding_;
        /**< Lazily constructed only if this thread ever has to help drain
             pendingSurfaces_ under backpressure (see onFlush()) -- kept
             alive and reused across every such event for this worker's
             lifetime, mirroring backgroundDrainLoop_'s own single reused
             KnottedSurface. Never the same object as this thread's own
             in-progress DFS embedding (owned by the `worker` lambda in
             runSearch_, entirely outside ThreadHook), so helping drain can
             never corrupt an in-flight DFS traversal. */

  public:
    ThreadHook(SurfaceSearch &owner, SurfaceTypeTally &tally, bool wantLinks,
              const SurfaceSearchCallbacks &callbacks)
        : owner_(owner), tally_(tally), wantLinks_(wantLinks),
          callbacks_(callbacks) {}

    void onFound(EmbeddedSubmanifold<4, 2> &embedding,
                const std::vector<int> &U, long long faceCount) override;

    void onFlush() override;
  };

  /**
   * Aux-thread hook passed to runSearch_(): continuously drains
   * pendingSurfaces_ into boundary links for as long as the search runs
   * (see backgroundDrainLoop_), then hands off whatever it hasn't gotten
   * to once the workers finish.
   */
  class AuxHooks : public RunSearchAuxHooks {
    SurfaceSearch &owner_;
    unsigned numThreads_;
    bool wantLinks_;
    const SurfaceSearchCallbacks &callbacks_;

  public:
    AuxHooks(SurfaceSearch &owner, unsigned numThreads, bool wantLinks,
             const SurfaceSearchCallbacks &callbacks)
        : owner_(owner), numThreads_(numThreads), wantLinks_(wantLinks),
          callbacks_(callbacks) {}

    std::thread spawn(std::atomic<bool> &workersFinished) override;

    void afterJoin() override;
  };

  PendingSurfaceBatch pendingSurfaces_; /**< Surfaces awaiting boundary-link processing. */
  LinkBoundaryTally linkTally_; /**< See linkTally(). */
  SurfaceTypeTally surfaceTypeTally_; /**< See surfaceTypeTally(). */

  /**
   * Deliberately separate from EmbeddingSearch::stopRequested_ (which
   * SIGINT also sets, via SigintScope): stopRequested_ stopping the final
   * drain too would mean a plain Ctrl+C silently drops whatever was still
   * queued but not yet identified, contradicting the "moving on to
   * whatever comes next with what was found so far" contract callers like
   * surfer.cpp's onInterrupted message promise (and requestStop()'s own
   * doc comment: it does NOT interrupt processRemainingSurfaceBoundaries()).
   * This flag exists for the narrower case a caller's own
   * onSurfaceBoundaryProcessed callback actually wants: "I've already
   * decided this search is satisfied, don't bother identifying the rest of
   * what's still queued" -- set only via skipRemainingBoundaryProcessing(),
   * never by SIGINT. See processBatchParallel_'s use of this.
   */
  std::atomic<bool> skipRemainingDrain_{false};

  /**
   * Shared across every KnottedSurface this search constructs (the main
   * DFS workers, the seed-validation/surfaceType probes, and the
   * boundary-link-processing workers) -- see PetalCache and
   * KnottedSurface's external-cache constructors. One instance regardless
   * of thread count, rather than one full cache per thread.
   */
  PetalCache petalCache_;

  /**
   * Set once via configureLimits(), before search() runs. Read by
   * ensureBoundarySigCaches_() when constructing each per-boundary-component
   * BoundarySignatureCache -- everything else configureLimits() controls is
   * applied directly to the relevant member at call time instead of being
   * read back out of this struct later.
   */
  SurfaceSearchLimits limits_;

  /**
   * mutable throughout: this whole group is a lazily-built cache (see
   * ensureBoundarySigCaches_()), populated on first use regardless of
   * whether that first use happens through a const accessor
   * (boundarySignatureCacheStats()) or describeBoundary_() itself.
   */
  mutable std::once_flag boundaryCachesOnce_;

  /**
   * Copies of the ambient triangulation's boundary components, built once,
   * lazily (see ensureBoundarySigCaches_()), and never modified again --
   * the exact same BoundaryComponent::build() call KnottedSurface's own
   * constructor uses for bdryComponents_, so boundarySigCaches_
   * canonicalizes marked edge sets against the identical edge indexing
   * boundaryLinks() itself produces them in.
   */
  mutable std::vector<regina::Triangulation<3>> boundaryComponentTris_;

  /**
   * One BoundarySignatureCache per entry of boundaryComponentTris_, shared
   * across every worker thread the same way petalCache_ is -- see
   * describeBoundary_(). unique_ptr since BoundarySignatureCache isn't
   * default-constructible (it tracks a fixed triangulation from
   * construction on) and must never be relocated once other threads may
   * be holding a reference to it.
   */
  mutable std::vector<std::unique_ptr<identify::BoundarySignatureCache>>
      boundarySigCaches_;

  /**
   * Builds boundaryComponentTris_/boundarySigCaches_ on first use,
   * regardless of which thread triggers it (std::once_flag) -- using the
   * same BoundaryComponent::build() call KnottedSurface's constructor
   * uses, so both stay in lockstep on edge indexing. const (like the
   * members it populates are mutable) since a first call may come through
   * either describeBoundary_() or the const stats accessors.
   */
  void ensureBoundarySigCaches_() const;

public:
  using EmbeddingSearch<4, 2>::EmbeddingSearch;

  /**
   * As the inherited seeded constructor, but additionally validates the
   * seed via a temporary KnottedSurface (not just a temporary
   * EmbeddedSubmanifold<4,2>), so a seed baked with an illegal
   * crossing/knot -- not just a facet-level collision -- is caught
   * eagerly too.
   *
   * `protectedBoundaryComponent` is forwarded straight through to the
   * base class constructor -- see EmbeddingSearch's own doc comment.
   */
  SurfaceSearch(
      const regina::Triangulation<4> &tri, const std::vector<int> &seedFaces,
      std::optional<size_t> protectedBoundaryComponent = std::nullopt);

  /**
   * Applies `limits` to every otherwise-unbounded structure this search
   * shares across its worker threads -- see SurfaceSearchLimits. Call
   * before search() (and before any const stats accessor that triggers
   * ensureBoundarySigCaches_(), e.g. boundarySignatureCacheStats()):
   * boundarySigCaches_ is built lazily, once, and each entry's clear
   * threshold is fixed at its own construction.
   */
  void configureLimits(const SurfaceSearchLimits &limits);

  /** Returns the current tally of found surfaces by boundary descriptor. */
  const LinkBoundaryTally &linkTally() const { return linkTally_; }

  /** Returns the current tally of found surfaces by SurfaceTypeKey. */
  const SurfaceTypeTally &surfaceTypeTally() const {
    return surfaceTypeTally_;
  }

  /**
   * Returns a snapshot of the shared PetalCache's current hit/miss
   * counters (see PetalCache::Stats) -- how much recomputation memoizing
   * addFace()'s local-flatness/transverse-self-intersection checks is
   * actually avoiding, aggregated across every thread since they all share
   * this one cache.
   */
  PetalCache::Stats petalCacheStats() const { return petalCache_.stats(); }

  /** As petalCacheStats(), but the shared cache's current distinct-petal count (since its last reset, if any). */
  size_t petalCacheSize() const { return petalCache_.size(); }

  /**
   * Returns a snapshot of pendingSurfaces_'s current size/cap/peak/
   * backpressure counters -- see PendingSurfaceQueueStats.
   */
  PendingSurfaceQueueStats pendingSurfaceQueueStats() const {
    return pendingSurfaces_.stats();
  }

  /**
   * Returns the aggregated hit/miss counters of every boundary component's
   * BoundarySignatureCache (see boundarySigCaches_) summed together -- how
   * much recomputation the pre-triangulation boundary-signature cache is
   * actually avoiding, across every search thread.
   */
  identify::BoundarySignatureCacheStats boundarySignatureCacheStats() const;

  /** As above, but the total number of distinct canonical boundary signatures seen across every boundary component. */
  size_t boundarySignatureCacheSize() const;

  /** As EmbeddingSearch::search(), reporting through `callbacks`. */
  SearchStats search(unsigned numThreads,
                     BoundaryCondition cond = BoundaryCondition::all,
                     const SurfaceSearchCallbacks &callbacks = {},
                     unsigned iddfsIterations = 0, long long iddfsStep = 0,
                     std::optional<long long> iddfsStart = std::nullopt,
                     std::optional<unsigned> finalThreads = std::nullopt,
                     bool orientableOnly = false);

  /**
   * Requests that processRemainingSurfaceBoundaries() abandon whatever is
   * still queued rather than identifying it, the next time (or currently,
   * if already running) it drains pendingSurfaces_. Unlike requestStop()
   * (which only stops the DFS itself, and is deliberately indistinguishable
   * from a real SIGINT -- see its own doc comment), this is never set by
   * SIGINT: call it only when a callback has independently decided this
   * search's purpose is already satisfied and the remaining queued
   * surfaces genuinely don't need identifying (e.g. verifyslicegenus,
   * once a resolving cobordism is found, has no use for identifying
   * whatever else happened to also queue up). Combine with requestStop()
   * to also stop the DFS from producing more work in the first place.
   */
  void skipRemainingBoundaryProcessing() {
    skipRemainingDrain_.store(true, std::memory_order_relaxed);
  }

  /**
   * Processes whatever boundary-link work is left after every DFS
   * worker thread (and the background drain loop) has finished, using
   * up to `numThreads` threads.
   *
   * Always announces itself via the three onBoundaryProcessing*
   * callbacks, since by construction this is the sole, complete pass
   * over anything not already handled by backgroundDrainLoop_.
   */
  void processRemainingSurfaceBoundaries(
      unsigned numThreads, const SurfaceSearchCallbacks &callbacks = {});

private:
  /**
   * Runs on the aux thread for as long as search() (with boundary links
   * wanted) is active: continuously pops and processes small batches
   * (via popSome(), never the whole queue), so pendingSurfaces_ never
   * grows unbounded and the once-a-second progress report keeps
   * reflecting linkTally_'s growth.
   *
   * Checks `workersFinished` between every individual entry, not just
   * between batches, so this stops within roughly one entry's
   * processing time once the DFS finishes or is interrupted -- handing
   * back any unprocessed remainder of its current batch to
   * pendingSurfaces_ rather than finishing it out.
   * processRemainingSurfaceBoundaries() then picks up everything left,
   * in one clean pass, once this has returned.
   */
  void backgroundDrainLoop_(const std::atomic<bool> &workersFinished,
                            const SurfaceSearchCallbacks &callbacks);

  /**
   * Shared parallel-processing core behind
   * processRemainingSurfaceBoundaries(): splits `batch` into CHUNK-sized
   * slices off a shared atomic cursor, processed by up to
   * min(numThreads, batch.size()) worker threads (each with its own
   * reused KnottedSurface, via processBatchRange_), reporting through
   * the onBoundaryProcessing* callbacks throughout.
   *
   * Each worker also checks skipRemainingDrain_ before claiming its next
   * chunk (see skipRemainingBoundaryProcessing()) -- deliberately not
   * stopRequested_/SIGINT, which stops the DFS but must NOT silently drop
   * whatever this drain has left to identify -- at CHUNK granularity, so
   * up to CHUNK already-claimed entries per thread still finish rather
   * than being abandoned mid-processing.
   */
  void processBatchParallel_(std::vector<std::vector<int>> batch,
                            unsigned numThreads,
                            const SurfaceSearchCallbacks &callbacks);

  /**
   * Runs the addFace/boundaryLinks/record/removeFace sequence for one
   * queued surface's face indices against `embedding` (reused by the
   * caller across many entries), recording one boundary descriptor (see
   * describeBoundary_) into linkTally_ if the surface has any boundary
   * at all, and firing `callbacks.onSurfaceBoundaryProcessed` with the
   * same descriptor.
   */
  void processEntry_(KnottedSurface &embedding,
                     const std::vector<int> &faceIndices,
                     const SurfaceSearchCallbacks &callbacks);

  /**
   * Replays batch[begin, end) via processEntry_, in the exact order the
   * DFS originally discovered each entry. `embedding` is reused across
   * the whole range so its boundary-component cache is only built once,
   * not once per queued surface.
   *
   * `processedCounter`, if given, is incremented by one after each
   * individual entry (not once for the whole range) -- per-entry cost
   * here is not uniform or necessarily cheap (identify()/pairSig() can
   * each cost real seconds), so onBoundaryProcessingProgress's once-a-
   * second reporting would otherwise sit frozen at a stale count for as
   * long as a single CHUNK-sized range takes a worker thread to finish,
   * rather than reflecting genuinely-live progress.
   */
  void processBatchRange_(KnottedSurface &embedding,
                          const std::vector<std::vector<int>> &batch,
                          size_t begin, size_t end,
                          const SurfaceSearchCallbacks &callbacks,
                          std::atomic<size_t> *processedCounter = nullptr);

  /**
   * Builds one human-readable descriptor of a surface's whole boundary
   * from `links` (already grouped and sorted by ambient boundary
   * component index): for each component, "<component+1>: " followed by
   * the comma-joined identify() of every curve of that component's
   * Link, followed by " (<link.identify()>)" only when that component
   * holds more than one curve -- with per-component entries joined by
   * ", ". E.g. "1: figure-eight, 2: unknot, unknot (Hopf link)".
   *
   * Each curve/link's identify() is resolved through this component's
   * BoundarySignatureCache (see boundarySigCaches_), so a boundary already
   * seen -- exactly, or up to a symmetry of its ambient boundary component
   * -- never reaches EdgeComplement::identify()'s
   * buildComplement()/simplify()/isoSig() at all. Not static (unlike
   * before) precisely because it now needs boundarySigCaches_.
   *
   * Returns both the flattened string and the same grouping structured (see
   * SurfaceBoundaryInfo::boundaryComponents) -- built together in one pass,
   * not two.
   */
  std::pair<std::string, std::vector<BoundaryComponentNames>>
  describeBoundary_(const std::vector<std::pair<size_t, Link>> &links);

  /**
   * Returns the most restrictive BoundaryCondition a surface satisfies,
   * given its already-computed `links` (see boundaryLinks()): closed (no
   * boundary at all), connected (every ambient boundary component
   * receives at most one of the surface's own boundary curves -- the
   * same condition boundaryComponentsMapInjectively() checks, but free
   * to compute here since `links` is already grouped by component), or
   * proper (otherwise).
   *
   * Shared by processEntry_ and search()'s onSeedFound, so the two never
   * drift apart on how this is computed.
   */
  static BoundaryCondition
  classifyByLinks_(const std::vector<std::pair<size_t, Link>> &links);

  /**
   * As classifyByLinks_(), but computed from `embedding` directly using
   * only isClosed()/isProper() (both O(1)) -- never "connected", since
   * that needs boundaryComponentsMapInjectively() (O(ambient
   * triangulation size)), deliberately not paid in the hot path this
   * feeds.
   *
   * Shared by ThreadHook::onFound and search()'s onSeedFound.
   */
  static BoundaryCondition
  classifyCheaply_(const EmbeddedSubmanifold<4, 2> &embedding);
};

#endif // SURFACESEARCH_H
