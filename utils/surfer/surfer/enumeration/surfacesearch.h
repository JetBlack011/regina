//
//  surfacesearch.h
//
//  Created by John Teague on 07/15/2026.
//

#ifndef SURFACESEARCH_H

#define SURFACESEARCH_H

#include <map>
#include <mutex>
#include <optional>
#include <thread>
#include <vector>

#include "surfer/submanifold/submanifold.h"
#include "surfer/enumeration/submanifoldsearch.h"
#include "surfer/enumeration/namecache.h"
#include "surfer/pairsig/pairsig.h"

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
       legitimate cobordism, since tubing keeps the surface properly embedded with
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
    /**< The BoundaryNamer's name of each of this component's curves,
         individually, in describeBoundary_'s order. */
    std::optional<std::string> linkName;
    /**< The BoundaryNamer's name of all of this component's curves
         together (BoundaryNamer::nameLink() on the whole Link) */
    std::vector<size_t> edgeIndices;
    /**< The component's boundary edges, sorted, as indices into that ambient
         boundary component's built triangulation -- the geometry the names
         above were computed from. Compared exactly; never approximated by
         a name. */
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
    std::function<std::map<const regina::Edge<3> *, size_t>()>
        captureBoundaryEdgeSurfaceComponent;
    /**< Lazy, as capturePairSig above: which connected component of the
         surface each boundary edge (the same Edge<3> objects as
         captureOrientedBoundaryLinks) belongs to. Needed because
         orientedBoundaryLinks() orients each component independently. */
    std::function<std::vector<int>()> captureFaces;
    /**< Lazy, as capturePairSig above: the surface's triangles, as indices
         into the ambient (EmbeddedSubmanifold::markedFaces()). With the
         ambient, all a pair signature is computed from (pairSig() in
         pairsig.h), so a caller can keep these and pay for signatures only
         of the surfaces it ends up needing. */
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
        namecache::BoundarySignatureCache::DEFAULT_CLEAR_THRESHOLD;
    /** See SurfaceSearch::LinkBoundaryTally::setCap(). */
    size_t boundaryTallyCap = 1'000'000;
    bool capturePairSig = false;
    /**< Whether onSurfaceFound/onSurfaceBoundaryProcessed populate
         SurfaceFoundInfo::pairSig (via the found embedding's own
         EmbeddedSubmanifold::pairSig()). */
    bool nameLinkCurves = true;
    /**< Whether a multi-curve boundary component's curves are also
         named one by one. Only their count is ever used downstream
         (splitBoundary(), classifyByLinks_()), so a caller that does not
         print them can switch this off; the curve names are then
         placeholders, and the count is unchanged. */
};

/**
 * Names the curves a search finds on its ambient boundary components. The
 * library has no naming route of its own: a search that tracks boundary
 * links (BoundaryCondition::proper or connected) must be given one
 * (SurfaceSearch::setBoundaryNamer()). UnlinkBoundaryNamer below names
 * without any census; cobound supplies the slice-genus search's own
 * (cobound/outgoing/outgoingnamer.h). Consulted only on a
 * BoundarySignatureCache miss (the cache still memoises the result per edge
 * set), from drain threads concurrently, so it must be thread-safe.
 *
 * A namer may say "Unknot" and "<n>-component unlink" only with a proof.
 * The search reads no other name (it counts curves), so any other name
 * answers to its reader. cobound's solver gates bounds by the curves'
 * count: one curve bears a bound whatever it is called, so a name of one
 * curve must denote one knot up to mirror (a name taken from its complement
 * does, by Gordon-Luecke), and 2 or more curves bear nothing unless they are
 * an unlink or proved by another input. See cobound/README.md, "Which
 * outgoing links bear a bound".
 */
class BoundaryNamer {
  public:
    virtual ~BoundaryNamer() = default;
    /** All the curves on boundary component `bc`, named together; a lone
     *  curve is a one-component link. */
    virtual std::string nameLink(size_t bc, const Link &curves) const = 0;
    /** One curve of a multi-curve boundary component `bc`, named on its own
     *  (wanted only with SurfaceSearchLimits::nameLinkCurves). */
    virtual std::string nameCurve(size_t bc, const Knot &curve) const = 0;
};

/**
 * Names boundary curves by their complements, through one pair of routes: a
 * lone curve (and each curve named on its own) by `knot`, several curves
 * together by `link`. UnlinkBoundaryNamer below is the census-free pair;
 * cobound's outgoing::ComplementNamer the census pair (census::nameComplement()).
 */
class ComplementBoundaryNamer : public BoundaryNamer {
  public:
    using KnotRoute = std::string (*)(const EdgeComplement &);
    using LinkRoute = std::string (*)(const Link &);

    ComplementBoundaryNamer(KnotRoute knot, LinkRoute link) : knot_(knot), link_(link) {}

    std::string nameLink(size_t bc, const Link &curves) const override;
    std::string nameCurve(size_t bc, const Knot &curve) const override;

  private:
    KnotRoute knot_;
    LinkRoute link_;
};

/**
 * Names boundary curves without any census, from their complements:
 * "Unknot", "<n>-component unlink", or else the complement's isoSig
 * (complement::unlinkNameOrIsoSig(), linknaming/complement/unlinknaming.h).
 */
class UnlinkBoundaryNamer : public ComplementBoundaryNamer {
  public:
    UnlinkBoundaryNamer();
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
     * indices for later link naming, merging into the owning
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
             pendingSurfaces_ under backpressure (see onFlush(), onPaused()). */

        /** Describes up to one batch of pendingSurfaces_ on this thread;
         *  returns how many it took (0 once the queue is empty). */
        size_t drainBatch_();

      public:
        ThreadHook(SurfaceSearch &owner, SurfaceTypeTally &tally,
                   bool wantLinks, const SurfaceSearchCallbacks &callbacks)
            : owner_(owner), tally_(tally), wantLinks_(wantLinks),
              callbacks_(callbacks) {}

        void onFound(EmbeddedSubmanifold<4, 2> &embedding,
                     const std::vector<int> &U, long long faceCount) override;

        void onFlush() override;

        /** While the queue's pause holds this worker, it drains too. */
        bool onPaused() override;
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
    /** See rebuildFailures(). */
    std::atomic<long long> rebuildFailures_{0};

    /**
     * The seed's own faces, when the seed is itself an accepted surface and
     * boundary links are wanted. backgroundDrainLoop_ describes it before
     * anything else, instead of the calling thread describing it before any
     * worker starts: describing it needs its pair signature, and the first
     * pair signature builds pairSigCtx_ (isoSigDetail of the whole
     * ambient, 21-26 s at 10 crossings), which used to hold every worker
     * back for that long. Written before the aux thread is spawned and
     * cleared by that thread, so it needs no lock.
     */
    std::optional<std::vector<int>> pendingSeed_;

    /**
     * The faces every drained surface shares: the seed's when seeded (every
     * surface a seeded search finds contains the seed), otherwise none.
     * Every drain embedding is built already holding them, and a queued
     * entry lists only the faces beyond them (see ThreadHook::onFound), so
     * describing a surface re-adds its few added faces rather than the whole
     * seed -- ~124 addFace() calls per surface before, ~4 now -- and the
     * queue holds a few ints per surface instead of the seed's ~120.
     *
     * Identical results: the embedding holds the same faces, added in the
     * same order (the seed's, then the entry's), and removing an entry's
     * faces in reverse returns it exactly to the seed-only state (rollback
     * union-find).
     */
    const std::vector<int> &residentFaces_() const;

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
     * with the number of cobordisms rather than the amount of work (79% of a
     * whole run's CPU, measured). Built on first use, not here, so a search
     * that never signs anything never pays for it.
     */
    LazyPairSigContext<4, 2> pairSigCtx_{skeleton_.triangulation()};

    SurfaceSearchLimits limits_;

    const BoundaryNamer *namer_ = nullptr; /**< See setBoundaryNamer(). */

    /** Handed to every worker's KnottedSurface; see
     * configureSelfIntersections(). */
    KnottedSurface::SelfIntersectionOptions selfIntersections_;

    mutable std::once_flag boundaryCachesOnce_;

    mutable std::vector<regina::Triangulation<3>> boundaryComponentTris_;

    mutable std::vector<std::unique_ptr<namecache::BoundarySignatureCache>>
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
    namecache::BoundarySignatureCacheStats boundarySignatureCacheStats() const;

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

    /** Whether skipRemainingBoundaryProcessing() was called, so some queued
     * surfaces may legitimately never reach onSurfaceBoundaryProcessed. */
    bool boundaryProcessingSkipped() const {
        return skipRemainingDrain_.load(std::memory_order_relaxed);
    }

    /** How many accepted surfaces failed to rebuild in the drain and were
     * therefore never described. Always 0 unless something is broken. */
    long long rebuildFailures() const {
        return rebuildFailures_.load(std::memory_order_relaxed);
    }

    /**
     * The pair-signature context of this search's ambient, built on the
     * first call (isoSigDetail of the whole ambient: tens of seconds at 10
     * crossings) and shared thereafter; concurrent callers wait for the one
     * build. What a caller signing surfaces off the drain threads needs:
     * context.sig(SurfaceBoundaryInfo::captureFaces()) is exactly the
     * surface's pairSig().
     */
    const PairSigContext<4, 2> &pairSigContext() const {
        return pairSigCtx_.get();
    }

    /**
     * Reads and writes the pair-signature context through `dir`
     * (PairSigContext::cached()), so a link searched again, or signed again
     * later, does not rebuild it. Call before search().
     */
    void setPairSigCacheDir(const std::string &dir) {
        pairSigCtx_.setCacheDir(dir);
    }

    /** Whether the context was read from the cache (see setPairSigCacheDir()). */
    bool pairSigContextLoaded() const { return pairSigCtx_.loaded(); }

    /**
     * Records `name` as boundary component `component`'s identity for the
     * edge set `edgeIndices` (sorted), so it is never named. For an edge
     * set known by construction: the incoming link, on the incoming boundary.
     */
    void primeBoundaryName(size_t component,
                           const std::vector<size_t> &edgeIndices,
                           const std::string &name);

    /**
     * Names every boundary curve the search describes through `namer`,
     * which must outlive the search. Required before a search that tracks
     * boundary links (BoundaryCondition::proper or connected): there is no
     * default namer, and search() refuses to start without one.
     */
    void setBoundaryNamer(const BoundaryNamer &namer) { namer_ = &namer; }

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
