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

#include "embeddedsubmanifold.h"
#include "embeddingsearch.h"
#include "identifycomplement.h"
#include "pairsig.h"

/**
 * Everything needed to describe one found surface when no boundary-link
 * data is being tracked (BoundaryCondition::all/closed), handed to
 * SurfaceSearchCallbacks::onSurfaceFound.
 */
struct SurfaceFoundInfo {
    bool orientable; /**< Whether the surface is orientable. */
    int genus;
    /**< Only meaningful for orientable connected surfaces. Else use
       `tubedGenus` below. Always check `connected` before reading this field.
     */
    int tubedGenus;
    /**< The genus of the connected surface obtained by discarding this
         surface's closed components and tubing the rest together This is the
       number to use for a slice-genus bound: a disconnected find is still a
       legitimate witness, since tubing keeps the surface properly embedded with
       the same boundary. See KnottedSurface::tubedSurfaceType() for the
       derivation.

       \pre Only meaningful when `orientable` is true. Also 0 (with nothing to
       bound) in the degenerate case where *every* component is closed, i.e.
       `punctures == 0`; callers deducing a slice-genus bound are already gated
       on the surface's boundary matching what they were looking for, so they
       never reach this field then. */
    int closedComponents;
    /**< How many closed (boundary-free) components were discarded when
         computing `tubedGenus`. */
    int punctures;
    /**< The surface's total number of punctures (boundary components). */
    bool connected; /**< Whether the surface is topologically connected */
    long long triangleCount; /**< The surface's number of triangles. */
    BoundaryCondition mostRestrictive;
    /**< The most restrictive BoundaryCondition this specific surface
         satisfies */
    std::function<std::string()> capturePairSig;
    /**< Computes (or returns the already-cached) reproducible
         isomorphism signature for this surface

         \warning Only valid to call synchronously, within the callback
         this SurfaceFoundInfo/SurfaceBoundaryInfo was passed to: it
         closes over the embedding by reference, which the caller may
         mutate (removeFace()'s unwind) or reuse for the next entry
         immediately after the callback returns. Safe to call more than
         once within that window (pairSig() itself is cached). */
    int resolvedVertices = 0;
    /**< How many ambient vertices the surface meets itself at. Zero for an
         embedded surface; positive only for one accepted under
         KnottedSurface::SelfIntersectionOptions::resolveUnlinked, where
         every such vertex is an unlinked self-intersection (paper §4.5). */
};

/**
 * SurfaceFoundInfo plus the boundary-link description, handed to
 * SurfaceSearchCallbacks::onSurfaceBoundaryProcessed once boundary links
 * are being tracked (BoundaryCondition::proper/connected).
 */
struct BoundaryComponentNames {
    size_t component; /**< Ambient boundary component index. */
    std::vector<std::string> curveNames;
    /**< The identify()'d name of each of this component's curves,
         individually, in describeBoundary_'s order. */
    std::optional<std::string> linkName;
    /**< identify()'d name of all of this component's curves drilled
         together (identify::identify() on the whole Link) */
};

struct SurfaceBoundaryInfo : SurfaceFoundInfo {
    std::string boundaryDescription;
    /**< A formatted description of the surface's boundary; empty if
         the surface turned out closed. */
    std::vector<BoundaryComponentNames> boundaryComponents;
    /**< The same (ambient boundary component, curve names) grouping
         boundaryDescription is flattened from */
    std::function<std::vector<std::pair<size_t, std::vector<OrientedCurve>>>()>
        captureOrientedBoundaryLinks;
    /**< Lazy, as capturePairSig above. */
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
         being tracked (BoundaryCondition::all/closed). */
    std::function<void(const SurfaceBoundaryInfo &)> onSurfaceBoundaryProcessed;
    /**< Fired once per found surface when boundary links are being
         tracked (BoundaryCondition::proper/connected). Mutually exclusive with
       onSurfaceFound: exactly one of the two fires per found surface. */
    std::function<void(size_t total, unsigned numThreads)>
        onBoundaryProcessingStarted;
    /**< Fired once, when boundary-link post-processing begins. */
    std::function<void(size_t processed, size_t total,
                       std::chrono::steady_clock::duration elapsed)>
        onBoundaryProcessingProgress;
    /**< Fired roughly once a second while boundary-link
         post-processing runs. */
    std::function<void(size_t total,
                       std::chrono::steady_clock::duration elapsed)>
        onBoundaryProcessingComplete;
    /**< Fired once, when boundary-link post-processing finishes. */
    std::function<void(size_t queueSize, size_t cap)> onQueueDrainPause;
    /**< Fired when the live pending-surface queue (see
         SurfaceSearch::PendingSurfaceBatch) exceeds its cap and a DFS
         worker thread pauses the *entire* search. */
    std::function<void(bool interrupted, size_t remaining)> onQueueDrainResume;
    /**< Fired once the drain triggered by onQueueDrainPause ends, either
         because the queue reached empty or because Ctrl+C fired mid-drain. */
};

/** Counters for how much recomputation the pending-surface queue's
 * cap/backpressure is costing/avoiding. */
struct PendingSurfaceQueueStats {
    size_t currentSize = 0;
    size_t cap = 0;
    size_t peakSize = 0;
    /** How many onFlush() calls (across every DFS worker thread) found the
     * queue over cap and helped drain it. */
    long long producerDrainEvents = 0;
    /** Total entries drained by a DFS worker thread helping out */
    long long producerDrainedEntries = 0;
};

/**
 * Thresholds for every otherwise-unbounded structure SurfaceSearch shares
 * across its worker threads (see SurfaceSearch::configureLimits()).
 */
struct SurfaceSearchLimits {
    /** See SurfaceSearch::PendingSurfaceBatch::setCap(). */
    size_t pendingSurfaceCap = 500'000;
    /** See PetalCache::setClearThreshold(). */
    size_t petalCacheLimit = 2'000'000;
    /** See BoundarySignatureCache's clearThreshold constructor parameter. */
    size_t boundarySignatureCacheLimit =
        identify::BoundarySignatureCache::DEFAULT_CLEAR_THRESHOLD;
    /** See SurfaceSearch::LinkBoundaryTally::setCap(). */
    size_t boundaryTallyCap = 1'000'000;
    bool capturePairSig = false;
    /**< Whether onSurfaceFound/onSurfaceBoundaryProcessed populate
         SurfaceFoundInfo::pairSig (via the found embedding's own
         EmbeddedSubmanifold::pairSig()). */
};

class SurfaceSearch : public EmbeddingSearch<4, 2> {
  public:
    /** See KnottedSurface::SurfaceTypeKey. */
    using SurfaceTypeKey = KnottedSurface::SurfaceTypeKey;

  private:
    /** A thread-safe tally of how many surfaces of each SurfaceTypeKey were
     * found. */
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
         */
        std::vector<std::vector<int>> popSome(size_t maxCount);

        /** The number of entries currently pending. */
        size_t size() const;

        /** The current cap; see setCap(). */
        size_t cap() const;

        /**
         * Sets the size this queue is considered "over cap" past. Call before
         * any worker thread starts searching.
         */
        void setCap(size_t cap);

        /** Records that a DFS worker thread (not backgroundDrainLoop_) drained
         * `entriesDrained` entries itself while helping under backpressure. */
        void recordProducerDrain(size_t entriesDrained);

        /** A snapshot of this queue's current size/cap/peak/backpressure
         * counters. */
        PendingSurfaceQueueStats stats() const;
    };

    /**
     * A thread-safe tally, per boundary descriptor  of how many surfaces of
     * each type were found with that descriptor.
     */
    class LinkBoundaryTally {
      private:
        mutable std::mutex mutex_;
        std::unordered_map<std::string, std::map<SurfaceTypeKey, long long>>
            descriptorSurfaceTypes_;
        std::vector<std::string> order_;
        /**< Distinct descriptors, in the order record() first saw each one */
        size_t descriptorCap_ = 100'000;
        long long refusedCount_ = 0;
        /**< Distinct descriptors seen but not recorded because
             descriptorSurfaceTypes_ was already at descriptorCap_.  */

      public:
        /** Records one surface of `type` found with boundary `descriptor`. */
        void record(const std::string &descriptor, const SurfaceTypeKey &type);

        /** Sets the distinct-descriptor cap past which record() stops admitting
         * brand-new descriptors (still tallying into already-known ones). Call
         * before any worker thread starts searching. */
        void setCap(size_t cap);

        /**
         * Returns a human-readable summary of the current tally.
         */
        std::string
        summary(std::optional<size_t> maxRecent = std::nullopt) const;
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
             pendingSurfaces_ under backpressure (see onFlush()). */

      public:
        ThreadHook(SurfaceSearch &owner, SurfaceTypeTally &tally,
                   bool wantLinks, const SurfaceSearchCallbacks &callbacks)
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

    PendingSurfaceBatch
        pendingSurfaces_; /**< Surfaces awaiting boundary-link processing. */
    LinkBoundaryTally linkTally_;       /**< See linkTally(). */
    SurfaceTypeTally surfaceTypeTally_; /**< See surfaceTypeTally(). */

    std::atomic<bool> skipRemainingDrain_{false};

    /**
     * Shared across every KnottedSurface this search constructs.
     */
    PetalCache petalCache_;

    /**
     * Shared ambient pair-signature data, likewise handed to every
     * KnottedSurface this search constructs.
     *
     * The ambient here is the search's own cobordism -- fixed for the whole
     * search -- so everything pairSig() derives from it is the same for every
     * surface found. Computing it per surface made the drain's cost scale
     * with the number of witnesses rather than the amount of work (79% of a
     * whole run's CPU, measured). Built on first use, not here, so a search
     * that never signs anything never pays for it.
     */
    LazyPairSigContext<4, 2> pairSigCtx_{skeleton_.triangulation()};

    SurfaceSearchLimits limits_;

    /** Handed to every worker's KnottedSurface; see
     * configureSelfIntersections(). */
    KnottedSurface::SelfIntersectionOptions selfIntersections_;

    mutable std::once_flag boundaryCachesOnce_;

    mutable std::vector<regina::Triangulation<3>> boundaryComponentTris_;

    mutable std::vector<std::unique_ptr<identify::BoundarySignatureCache>>
        boundarySigCaches_;

    /**
     * Builds boundaryComponentTris_/boundarySigCaches_ on first use.
     */
    void ensureBoundarySigCaches_() const;

  public:
    using EmbeddingSearch<4, 2>::EmbeddingSearch;

    /**
     * As the inherited seeded constructor, but additionally validates the seed
     * via a temporary KnottedSurface, so a seed baked with an illegal
     * crossing/knot is caught.
     *
     * `protectedBoundaryComponent` is forwarded straight through to the
     * base class constructor; see EmbeddingSearch's own doc comment.
     */
    SurfaceSearch(
        const regina::Triangulation<4> &tri, const std::vector<int> &seedFaces,
        std::optional<size_t> protectedBoundaryComponent = std::nullopt);

    void configureLimits(const SurfaceSearchLimits &limits);

    /**
     * Sets how the search's workers treat self-intersections (see
     * KnottedSurface::SelfIntersectionOptions): whether resolvable ones are
     * accepted, and whether to record a SelfIntersectionCensus. Call before
     * search(). The defaults change nothing.
     */
    void configureSelfIntersections(
        const KnottedSurface::SelfIntersectionOptions &options) {
        selfIntersections_ = options;
    }

    /** Returns the current tally of found surfaces by boundary descriptor. */
    const LinkBoundaryTally &linkTally() const { return linkTally_; }

    /** Returns the current tally of found surfaces by SurfaceTypeKey. */
    const SurfaceTypeTally &surfaceTypeTally() const {
        return surfaceTypeTally_;
    }

    /**
     * Returns a snapshot of the shared PetalCache's current hit/miss
     * counters (see PetalCache::Stats).
     */
    PetalCache::Stats petalCacheStats() const { return petalCache_.stats(); }

    /** As petalCacheStats(), but the shared cache's current distinct-petal
     * count (since its last reset, if any). */
    size_t petalCacheSize() const { return petalCache_.size(); }

    /**
     * Returns a snapshot of pendingSurfaces_'s current size/cap/peak/
     * backpressure counters; see PendingSurfaceQueueStats.
     */
    PendingSurfaceQueueStats pendingSurfaceQueueStats() const {
        return pendingSurfaces_.stats();
    }

    /**
     * Returns the aggregated hit/miss counters of every boundary component's
     * BoundarySignatureCache (see boundarySigCaches_) summed together.
     */
    identify::BoundarySignatureCacheStats boundarySignatureCacheStats() const;

    /** As above, but the total number of distinct canonical boundary signatures
     * seen across every boundary component. */
    size_t boundarySignatureCacheSize() const;

    /** As EmbeddingSearch::search(), reporting through `callbacks`. */
    SearchStats search(unsigned numThreads,
                       BoundaryCondition cond = BoundaryCondition::all,
                       const SurfaceSearchCallbacks &callbacks = {},
                       unsigned iddfsIterations = 0, long long iddfsStep = 0,
                       std::optional<long long> iddfsStart = std::nullopt,
                       std::optional<unsigned> finalThreads = std::nullopt,
                       bool orientableOnly = false,
                       std::optional<long long> hardFaceCap = std::nullopt,
                       long long rootBudgetStart = 0,
                       long long rootBudgetGrowth = 2);

    void skipRemainingBoundaryProcessing() {
        skipRemainingDrain_.store(true, std::memory_order_relaxed);
    }

    /**
     * Processes whatever boundary-link work is left after every DFS
     * worker thread (and the background drain loop) has finished, using
     * up to `numThreads` threads.
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
     */
    void backgroundDrainLoop_(const std::atomic<bool> &workersFinished,
                              const SurfaceSearchCallbacks &callbacks);

    /**
     * Shared parallel-processing core behind
     * processRemainingSurfaceBoundaries().
     */
    void processBatchParallel_(std::vector<std::vector<int>> batch,
                               unsigned numThreads,
                               const SurfaceSearchCallbacks &callbacks);

    /**
     * Runs the addFace/boundaryLinks/record/removeFace sequence for one
     * queued surface's face indices against `embedding`, recording one boundary
     * descriptor into linkTally_ if the surface has any boundary at all, and
     * firing `callbacks.onSurfaceBoundaryProcessed` with the same descriptor.
     */
    void processEntry_(KnottedSurface &embedding,
                       const std::vector<int> &faceIndices,
                       const SurfaceSearchCallbacks &callbacks);

    /**
     * Replays batch[begin, end) via processEntry_, in the same order the
     * DFS originally discovered each entry. `embedding` is reused across
     * the whole range so its boundary-component cache is only built once,
     * not once per queued surface.
     */
    void processBatchRange_(KnottedSurface &embedding,
                            const std::vector<std::vector<int>> &batch,
                            size_t begin, size_t end,
                            const SurfaceSearchCallbacks &callbacks,
                            std::atomic<size_t> *processedCounter = nullptr);

    /**
     * Builds one human-readable descriptor of a surface's whole boundary
     * from `links`.
     *
     * Returns both the flattened string and the same grouping structured (see
     * SurfaceBoundaryInfo::boundaryComponents).
     */
    std::pair<std::string, std::vector<BoundaryComponentNames>>
    describeBoundary_(const std::vector<std::pair<size_t, Link>> &links);

    /**
     * Returns the most restrictive BoundaryCondition a surface satisfies,
     * given its already-computed `links`.
     *
     * Shared by processEntry_ and search()'s onSeedFound, so the two never
     * drift apart on how this is computed.
     */
    static BoundaryCondition
    classifyByLinks_(const std::vector<std::pair<size_t, Link>> &links);

    /**
     * As classifyByLinks_(), but computed from `embedding` directly using
     * only isClosed()/isProper() (both O(1)).
     *
     * Shared by ThreadHook::onFound and search()'s onSeedFound.
     */
    static BoundaryCondition
    classifyCheaply_(const EmbeddedSubmanifold<4, 2> &embedding);
};

#endif // SURFACESEARCH_H
