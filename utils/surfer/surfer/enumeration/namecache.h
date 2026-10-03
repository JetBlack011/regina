//
//  namecache.h
//
//  The drain's names, per curve set up to automorphism.
//

#ifndef SURFER_NAMECACHE_H
#define SURFER_NAMECACHE_H

#include <functional>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

#include <triangulation/dim3.h>

/*! \file utils/surfer/surfer/enumeration/namecache.h
 *  \brief Memoizes the names of marked edge sets of one ambient boundary
 *  component, canonicalized against that component's automorphism group.
 */

namespace namecache {

/** Counters for how much recomputation BoundarySignatureCache is actually avoiding. */
struct BoundarySignatureCacheStats {
    long long checks = 0;
    long long hits = 0;
    long long cacheResets = 0; /**< How many times this cache has been fully cleared after exceeding its clear threshold. */
};

/**
 * Memoizes a BoundaryNamer's names for marked edge sets in one fixed ambient
 * boundary-component triangulation, canonicalized against that
 * triangulation's own automorphism group -- so a boundary curve already
 * seen (exactly, or up to a symmetry of the boundary component)
 * short-circuits before buildComplement()/simplify()/isoSig() ever runs,
 * rather than only being caught by the isoSig-keyed recognitionCache
 * *after* paying for that triangulation (see complementcache.h).
 *
 * recognitionCache is still needed on top of this: two curves can be the
 * same knot/link type without being combinatorially related by any
 * automorphism of the boundary component (e.g. two non-isomorphic edge
 * paths that happen to drill out to the same manifold), in which case only
 * recognitionCache -- keyed on the post-simplify isoSig -- catches the
 * duplicate. This cache is a cheaper pre-filter in front of that one, not a
 * replacement for it.
 *
 * The ambient boundary-component triangulation this tracks must stay fixed
 * for this cache's lifetime, and must be the exact triangulation
 * KnottedSurface::boundaryLinks() draws its edges from (i.e. built the same
 * way, via BoundaryComponent::build()) -- see SurfaceSearch, which owns one
 * instance per boundary component, shared across every worker thread's
 * KnottedSurface, the same way it shares one PetalCache.
 *
 * Thread-safe: the (lazily built, then read-only) automorphism group is
 * guarded by a std::once_flag, and the memoization table by its own mutex
 * -- so a cache hit never blocks behind another thread still computing the
 * group, and computing the group happens exactly once regardless of how
 * many threads race to trigger it.
 */
class BoundarySignatureCache {
  public:
    /** Default entry-count threshold; see the constructor's `clearThreshold` parameter. */
    static constexpr size_t DEFAULT_CLEAR_THRESHOLD = 200'000;

    /**
     * Tracks `boundary`, which must outlive this cache. Once this cache
     * holds `clearThreshold` distinct canonical signatures, it clears
     * itself entirely before admitting the next new one -- unlike
     * PetalCache, no epoch-tagging is needed here, since a cache_ entry is
     * keyed by its canonical signature (content), not a recyclable integer
     * id: a lookup racing a clear simply sees a clean miss, never a wrong
     * hit against an unrelated key.
     */
    explicit BoundarySignatureCache(
        const regina::Triangulation<3> &boundary,
        size_t clearThreshold = DEFAULT_CLEAR_THRESHOLD);

    /**
     * Returns the name of `edgeIndices` -- this boundary
     * component's marked edges, sorted, e.g. from
     * EdgeComplement::edgeIndices() -- computing it via `compute` on a
     * cache miss and memoizing the result. `compute` is only invoked on a
     * miss; typically a call to the search's BoundaryNamer.
     */
    std::string identifyCached(const std::vector<size_t> &edgeIndices,
                                const std::function<std::string()> &compute);

    /**
     * Records `name` as the answer for `edgeIndices` without computing
     * anything, for an edge set whose identity is known by construction (a
     * row's own link on its search side). Harmless if the cache is later
     * cleared: that entry is then just recomputed on demand.
     */
    void prime(const std::vector<size_t> &edgeIndices, const std::string &name);

    /**
     * Returns a snapshot of this cache's current hit/miss counters; see
     * BoundarySignatureCacheStats. Returned by value (not by reference)
     * since another thread may be concurrently updating the live counters.
     */
    BoundarySignatureCacheStats stats() const;

    /** This cache's current entry count (distinct canonical signatures seen). */
    size_t size() const;

  private:
    const regina::Triangulation<3> *boundary_;

    std::once_flag groupOnce_;
    std::vector<std::vector<size_t>> automorphismEdgeImages_;
        /**< Lazily built by ensureAutomorphismGroup_():
             automorphismEdgeImages_[g][e] is the image of edge `e` of
             *boundary_ under the g-th automorphism found. Read-only once
             built, so safe to read from any thread without further
             synchronization (std::once_flag already establishes
             happens-before for every caller, not just the one that ran
             the initializer). */

    size_t clearThreshold_;

    mutable std::mutex cacheMutex_; /**< Guards cache_ and stats_. */
    std::unordered_map<std::string, std::string> cache_;
    mutable BoundarySignatureCacheStats stats_;

    /** Builds automorphismEdgeImages_, exactly once, on first use. */
    void ensureAutomorphismGroup_();

    /**
     * Canonicalizes `edgeIndices` by minimizing its sorted image over
     * every automorphism in automorphismEdgeImages_ (triggering
     * ensureAutomorphismGroup_() first), returning the result as a
     * comma-joined string suitable as a cache_ key.
     */
    std::string canonicalKey_(const std::vector<size_t> &edgeIndices);
};

} // namespace namecache

#endif // SURFER_NAMECACHE_H
