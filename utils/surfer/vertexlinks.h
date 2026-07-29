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
   * An interned petal's id, opaque outside this class: `(epoch << 32) |
   * localId`. The epoch lets internPetal() clear this cache outright once
   * it grows past setClearThreshold() (see that method) without risking a
   * stale id silently colliding with an unrelated, numerically-coincident
   * id minted after the clear -- every lookup/record call below checks the
   * embedded epoch first, and treats a mismatch as "unrelated to anything
   * currently cached" (a clean miss/no-op), never as a hit against the
   * wrong petal. Callers should treat this purely as an id, never assume
   * anything about its bit layout.
   */
  using PetalId = uint64_t;

  /**
   * Interns `corners` as a petal, returning its id -- the same corner set
   * always maps to the same id (within the current epoch; see PetalId),
   * regardless of insertion order (sorted internally) or whether it was
   * already known.
   *
   * If this cache has grown past its clear threshold (see
   * setClearThreshold()), this is where that's actually enforced: this is
   * the only method that ever mints a new id, so it's the only place that
   * needs to check.
   */
  PetalId internPetal(std::vector<Corner> corners);

  /** Returns the cached isUnknot() result for petal `id`, if known and not stale (see PetalId). */
  std::optional<bool> lookupUnknot(PetalId id) const;

  /** Records `id`'s isUnknot() result. A stale `id` (see PetalId) is silently ignored. */
  void recordUnknot(PetalId id, bool isUnknot);

  /**
   * Returns the cached (linkingNumberWith() != 0) result for the
   * unordered pair {a, b}, if known and neither is stale (see PetalId).
   */
  std::optional<bool> lookupLinksNonzero(PetalId a, PetalId b) const;

  /**
   * Records the unordered pair {a, b}'s (linkingNumberWith() != 0) result.
   * Silently ignored if either is stale (see PetalId).
   */
  void recordLinksNonzero(PetalId a, PetalId b, bool nonzero);

  /**
   * Sets the entry-count threshold (counted against unknotById_, i.e.
   * distinct petals interned since the last clear) past which
   * internPetal() clears this cache outright before minting the next new
   * id. Call before any worker thread starts using this cache; a call
   * racing with in-flight searches only risks applying one clear-check
   * late, which is harmless.
   */
  void setClearThreshold(size_t threshold);

  /** This cache's current distinct-petal count (since the last clear, if any). */
  size_t size() const;

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
    long long cacheResets = 0; /**< How many times this cache has been fully cleared after exceeding its clear threshold. */
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

  /** Packs the unordered pair of local ids {a, b} into a single collision-free key. */
  static uint64_t linkKey_(int a, int b);

  static PetalId makeId_(uint64_t epoch, int localId) {
    return (epoch << 32) | static_cast<uint32_t>(localId);
  }
  static uint64_t epochOf_(PetalId id) { return id >> 32; }
  static int localIdOf_(PetalId id) {
    return static_cast<int>(id & 0xFFFFFFFFu);
  }

  /** Default entry-count threshold; see setClearThreshold(). */
  static constexpr size_t DEFAULT_CLEAR_THRESHOLD = 2'000'000;

  mutable std::mutex mutex_; /**< Guards every member below. */

  std::unordered_map<std::vector<Corner>, int, CornersHash> byCorners_;
      /**< corners -> local id (current epoch only); the interning table. */
  std::vector<std::optional<bool>> unknotById_;
      /**< isUnknot() result per local id, indexed by local id (current epoch only). */
  std::unordered_map<uint64_t, bool> linkingCache_;
      /**< (linkingNumberWith() != 0) result per unordered local-id pair (current epoch only). */

  size_t clearThreshold_ = DEFAULT_CLEAR_THRESHOLD;
  uint64_t epoch_ = 0; /**< Bumped every time this cache is cleared; see PetalId. */

  mutable Stats stats_; /**< Updated even by the "read-only" lookup methods, hence mutable. */
};

#endif // VERTEXLINKS_H
