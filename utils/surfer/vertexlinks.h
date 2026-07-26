//
//  vertexlinks.h
//
//  Created by John Teague on 07/26/2026.
//

#ifndef VERTEXLINKS_H

#define VERTEXLINKS_H

#include <cstdint>
#include <mutex>
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>

/*! \file utils/surfer/vertexlinks.h
 *  \brief Memoizes KnottedSurface::addFace()'s local-flatness and
 *  transverse-self-intersection checks.
 */

/**
 * Memoizes Knot::isUnknot() and Knot::linkingNumberWith() results for closed
 * petals, keyed by each petal's own combinatorial identity -- the exact set
 * of (ambient face, local vertex) "triangle-corner" pairs it's made of.
 *
 * That identity is a pure function of which triangles are currently present
 * at the ambient vertex (KnottedSurface::addFace()'s DSU-union logic only
 * ever unions triangle-corners that are both present and share a spoke, so
 * a petal's closure and content never depend on anything elsewhere in the
 * submanifold) -- so the same petal recurring across different DFS
 * branches/backtracks always has the same answer, and this cache turns
 * repeat occurrences into O(1) lookups instead of re-running Regina's
 * isSolidTorus()/HomologicalData machinery from scratch.
 *
 * Thread-safe (a single mutex guards every operation): intended to be
 * shared across every KnottedSurface instance in a search (one per worker
 * thread otherwise means one full copy of the cache per thread, which
 * multiplies with thread count -- see SurfaceSearch, which owns one
 * PetalCache and passes it by reference to each KnottedSurface it
 * constructs). KnottedSurface also supports owning a private, unshared
 * instance (its single-argument constructors) for standalone/test use where
 * sharing doesn't matter.
 *
 * The lock's critical sections are all tiny (a hash-map lookup/insert), so
 * a single mutex -- rather than sharding -- matches this codebase's own
 * convention for shared accumulators (see SurfaceSearch::SurfaceTypeTally
 * and friends in embeddingsearch.h) and keeps contention low relative to
 * what it replaces (Regina calls costing milliseconds or more per miss).
 */
class PetalCache {
public:
  /** A single (ambient face index, local vertex 0..2) triangle-corner. */
  using Corner = std::pair<int, int>;

  /**
   * Interns `corners` as a petal, returning its id -- the same corner set
   * always maps to the same id, regardless of insertion order (sorted
   * internally) or whether it was already known.
   */
  int internPetal(std::vector<Corner> corners);

  /** Returns the cached isUnknot() result for petal `id`, if known. */
  std::optional<bool> lookupUnknot(int id) const;

  /** Records `id`'s isUnknot() result. */
  void recordUnknot(int id, bool isUnknot);

  /**
   * Returns the cached (linkingNumberWith() != 0) result for the
   * unordered pair {a, b}, if known.
   */
  std::optional<bool> lookupLinksNonzero(int a, int b) const;

  /** Records the unordered pair {a, b}'s (linkingNumberWith() != 0) result. */
  void recordLinksNonzero(int a, int b, bool nonzero);

  /**
   * Records that KnottedSurface::addFace() rejected a face for closing a
   * knotted petal (non-local-flatness). Tracked here, alongside the
   * memoization counters, purely so every KnottedSurface sharing this
   * cache contributes to one aggregate total instead of an unshared
   * per-instance count that a multi-threaded search could never read back
   * in full.
   */
  void recordLocalFlatnessRejection();

  /** As recordLocalFlatnessRejection(), for a transverse self-intersection rejection. */
  void recordTransverseRejection();

  /** Counters for how much recomputation this cache is actually avoiding, and how often addFace()'s checks actually reject something. */
  struct Stats {
    long long unknotChecks = 0;
    long long unknotCacheHits = 0;
    long long linkingChecks = 0;
    long long linkingCacheHits = 0;
    long long localFlatnessRejections = 0;
    long long transverseRejections = 0;
  };

  /**
   * Returns a snapshot of this cache's current hit/miss counters; see
   * Stats. Returned by value (not by reference) since another thread may
   * be concurrently updating the live counters.
   */
  Stats stats() const;

private:
  struct CornersHash {
    size_t operator()(const std::vector<Corner> &corners) const;
  };

  /** Packs the unordered pair {a, b} into a single collision-free key. */
  static uint64_t linkKey_(int a, int b);

  mutable std::mutex mutex_; /**< Guards every member below. */

  std::unordered_map<std::vector<Corner>, int, CornersHash> byCorners_;
      /**< corners -> id; the interning table. */
  std::vector<std::optional<bool>> unknotById_;
      /**< isUnknot() result per petal id, indexed by id. */
  std::unordered_map<uint64_t, bool> linkingCache_;
      /**< (linkingNumberWith() != 0) result per unordered id pair. */

  mutable Stats stats_; /**< Updated even by the "read-only" lookup methods, hence mutable. */
};

#endif // VERTEXLINKS_H
