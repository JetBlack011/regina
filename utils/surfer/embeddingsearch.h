//
//  embeddingsearch.h
//
//  Created by John Teague on 07/15/2026.
//

#ifndef EMBEDDINGSEARCH_H

#define EMBEDDINGSEARCH_H

#include <atomic>
#include <chrono>
#include <functional>
#include <map>
#include <mutex>
#include <optional>
#include <thread>

#include <triangulation/dim2.h>
#include <unordered_map>

#include "embeddedsubmanifold.h"
#include "enumerate_cis.h"
#include "skeleton.h"

/*! \file utils/surfer/embeddingsearch.h
 *  \brief Parallel search over embedded subcomplexes of a triangulation.
 */

/** A graph's vertex count, together with its 1-indexed adjacency list. */
using AdjacencyList = std::pair<int, std::vector<std::vector<int>>>;

/** Formats a duration for progress/summary output (e.g. "1h23m45s"). */
std::string formatElapsed(std::chrono::steady_clock::duration d);

/** Returns the display name of `cond`. */
const char *boundaryConditionName(BoundaryCondition cond);

/**
 * The face-count cap for iterative-deepening pass `iter` (1-indexed): the
 * first pass is capped at `start`, and each subsequent pass grows the cap
 * by `step`. See EmbeddingSearch::search()'s iddfsStart/iddfsStep
 * parameters.
 */
inline long long iddfsCapForRound(unsigned iter, long long start,
                                  long long step) {
  return start + static_cast<long long>(iter - 1) * step;
}

/**
 * Converts a round's face-count cap into the DepthCappedPredicate maxDepth
 * that enforces it.
 *
 * `capFaces` counts faces *added by the search*, never a seed's own faces:
 * a seed (whatever its size) is committed via a single tryAdd() before any
 * root is explored (see ConnectedInducedSubgraphEnumerator::
 * seedFastForward_), so DepthCappedPredicate's nested-tryAdd() depth
 * charges it exactly one unit regardless of how many faces it contains --
 * every other commit corresponds to exactly one added face (see Graph's
 * doc comment). So a seeded round needs `capFaces + 1`: one depth unit for
 * that one-time seed commit, plus `capFaces` more for the search's own
 * single-face additions. Unseeded, every commit is already a single added
 * face, so the cap applies directly.
 */
inline long long iddfsMaxDepth(long long capFaces, bool isSeeded) {
  return isSeeded ? capFaces + 1 : capFaces;
}

/**
 * A snapshot of a search() call's progress/results, handed to
 * SearchCallbacks rather than printed directly.
 */
struct SearchStats {
  std::chrono::steady_clock::duration elapsed{}; /**< Time elapsed since the search started. */
  size_t rootsCompleted = 0; /**< The number of DFS roots fully explored so far. */
  size_t totalRoots = 0; /**< The total number of DFS roots this search will explore. */
  long long foundCount = 0; /**< Raw candidates visited, regardless of BoundaryCondition. */
  long long embeddedCount = 0; /**< Candidates satisfying isEmbedded(), regardless of BoundaryCondition. */
  long long satisfyingCount = 0; /**< Candidates satisfying isEmbedded() and the BoundaryCondition. */
  long long satisfyingFaceSum = 0; /**< The sum of face counts among satisfying finds. */
  long long largestSatisfying = 0; /**< The largest face count among satisfying finds. */

  unsigned iddfsRound = 1;
      /**< The iterative-deepening round currently running (1-indexed);
           see iddfsTotalRounds. */
  unsigned iddfsTotalRounds = 1;
      /**< The total number of iterative-deepening rounds this search will
           run: iddfsIterations capped passes, plus one final unbounded
           pass. 1 when iddfsIterations == 0 (the default, i.e. iterative
           deepening isn't in use). */
  bool iddfsCapped = false;
      /**< Whether the current round (see iddfsRound) is a capped pass, as
           opposed to the final unbounded one. */
  long long iddfsCap = 0;
      /**< The current round's face-count cap; only meaningful when
           iddfsCapped is true. */

  /** Returns the average face count among satisfying finds, or 0 if there are none. */
  double averageSatisfyingFaces() const {
    return satisfyingCount > 0
               ? static_cast<double>(satisfyingFaceSum) /
                     static_cast<double>(satisfyingCount)
               : 0.0;
  }
};

/**
 * Callbacks a caller of search() may supply to observe its progress, in
 * place of search() printing anything itself. Every field defaults to a
 * no-op, so supplying none of these costs nothing beyond the
 * once-a-second call overhead already inherent to the reporter.
 */
struct SearchCallbacks {
  std::function<void(const SearchStats &)> onProgress;
      /**< Fired roughly once a second while the search runs. */
  std::function<void()> onInterrupted;
      /**< Fired at most once, if Ctrl+C is caught mid-search. */
  std::function<void(const SearchStats &)> onSearchComplete;
      /**< Fired once, when the search has finished (the same
           SearchStats is also returned by search() itself). */
};

/**
 * Enumerates embedded subcomplexes (EmbeddedSubmanifold<dim, subdim>) of a
 * triangulation via a parallel depth-first search over its face-adjacency
 * graph (see Skeleton), reporting progress and results through
 * SearchCallbacks.
 *
 * \tparam dim the dimension of the ambient triangulation.
 * \tparam subdim the dimension of the faces being assembled into a
 * subcomplex.
 */
template <int dim, int subdim> class EmbeddingSearch {
protected:
  /**
   * Drives the DFS's embeddability check against an embedding object.
   *
   * \tparam EmbeddingT generic rather than fixed to
   * EmbeddedSubmanifold<dim,subdim>, so SurfaceSearch can supply
   * KnottedSurface instead (whose addFace()/removeFace() name-hide the
   * base class's to add transverse-self-intersection/local-flatness
   * checks) without EmbeddingSearch<dim,subdim> itself ever referencing
   * KnottedSurface. Since nothing here is virtual, resolving to the
   * override requires the predicate's *declared* member type to already
   * be the derived type -- template argument deduction means callers
   * still never need to name EmbeddingT explicitly.
   *
   * \note The members below are defined inline here, deliberately, rather
   * than out-of-line in embeddingsearch.cpp. Being a member template, its
   * specializations are *not* produced by that file's explicit
   * instantiations of EmbeddingSearch<dim,subdim> (explicitly instantiating
   * a class template does not instantiate its member templates), and the
   * trivial constructor is inlined at its call sites rather than emitted.
   * Any other translation unit constructing one -- as
   * tests/profile_driver.cpp does -- would then fail to link.
   */
  template <typename EmbeddingT>
  class EmbeddednessPredicate : public ConditionalPredicate {
    EmbeddingT &embedding_;
    const std::vector<std::vector<int>> &graphToSkel_;

  public:
    /** Drives `embedding` via the graph-to-skeleton-face mapping `graphToSkel` (see Graph). */
    EmbeddednessPredicate(EmbeddingT &embedding,
                          const std::vector<std::vector<int>> &graphToSkel)
        : embedding_(embedding), graphToSkel_(graphToSkel) {}

    bool tryAdd(int v) override {
      const auto &faces = graphToSkel_[v - 1];
      // Every non-seed node maps to exactly one face, where addFace() alone
      // is already correct (and this is the hot path -- called once per DFS
      // node visited) -- only a seeded graph's node 1 (the whole seed) can
      // have more than one, where addFaces()'s order search is needed. See
      // addFaces() for why: addFace() is order-sensitive, and the seed's
      // faces don't generally arrive in a directly-addable order.
      if (faces.size() == 1)
        return embedding_.addFace(faces[0]);
      return embedding_.addFaces(faces);
    }

    void undo(int v) override {
      const auto &faces = graphToSkel_[v - 1];
      for (auto it = faces.rbegin(); it != faces.rend(); ++it)
        embedding_.removeFace(*it);
    }
  };

  /** The face-adjacency graph searched over, together with its vertex-to-skeleton-face mapping. */
  struct Graph {
    AdjacencyList adjList; /**< The graph itself. */
    std::vector<std::vector<int>> graphToSkel;
        /**< graphToSkel[v-1]: the skeleton face indices graph vertex v
             corresponds to. Every vertex has a singleton list, except an
             unseeded graph's vertex 1, which has none -- in the seeded
             case, vertex 1 is the contracted seed and maps to every one
             of its skeleton faces at once. */
  };

  const Skeleton<dim, subdim> skeleton_; /**< The ambient skeleton. */
  const Graph graph_; /**< The face-adjacency graph searched over. */
  const bool isSeeded_ = false; /**< Whether this search was constructed seeded. */

private:
  std::vector<EmbeddedSubmanifold<dim, subdim>> valid_embeddings_; /**< Currently unused. */

public:
  /** Creates an unseeded search over `tri`'s skeleton of subdim-faces. */
  EmbeddingSearch(const regina::Triangulation<dim> &tri);

  /**
   * Seeded search, restricted to embeddings containing every face in
   * `seedFaces`.
   *
   * \pre `seedFaces` (skeleton face indices) is a connected, jointly
   * addable set of faces (e.g. the output of CollarBuilder, translated
   * to skeleton indices via Triangle::index()).
   *
   * \throws regina::InvalidArgument if the precondition is violated,
   * checked eagerly here rather than deferring discovery to search().
   */
  EmbeddingSearch(const regina::Triangulation<dim> &tri,
                  const std::vector<int> &seedFaces);

  /** Returns the number of subdim-faces eligible to appear in an embedding. */
  size_t numEmbeddableFaces() const { return graph_.graphToSkel.size(); }

  /**
   * Runs the search using `numThreads` worker threads, restricted to
   * embeddings satisfying `cond`, reporting through `callbacks`.
   *
   * \param iddfsIterations if greater than 0, runs this many capped
   * "fast sweep" passes over every root before the final unbounded pass,
   * with pass `i`'s cap set to `iddfsStart + (i - 1) * iddfsStep` faces
   * (see runSearch_). Defaults to 0 (a single unbounded pass, this
   * method's original behavior).
   * \param iddfsStep the face-count increment per capped pass after the
   * first; only consulted when `iddfsIterations > 0`.
   * \param iddfsStart the first capped pass's face-count cap; defaults to
   * `iddfsStep` when not given (i.e. pass `i`'s cap is `i * iddfsStep`,
   * this parameter's original behavior before iddfsStart existed). Only
   * consulted when `iddfsIterations > 0`.
   * \param finalThreads the number of threads to use for the final
   * unbounded pass; defaults to `numThreads` when not given.
   */
  SearchStats search(const unsigned numThreads,
                     BoundaryCondition cond = BoundaryCondition::all,
                     const SearchCallbacks &callbacks = {},
                     unsigned iddfsIterations = 0, long long iddfsStep = 0,
                     std::optional<long long> iddfsStart = std::nullopt,
                     std::optional<unsigned> finalThreads = std::nullopt);

protected:
  /**
   * Shared harness behind both EmbeddingSearch::search() and
   * SurfaceSearch::search(): the two differ only in a handful of
   * injection points, threaded through here so the bulk of either
   * function (threading, atomics, reporting) exists exactly once.
   *
   * \param numThreads the number of worker threads to search with.
   * \param cond restricts results to embeddings satisfying this
   * condition.
   * \param makeEmbedding factory called once per worker thread (plus
   * once more for the seeded-root prototype), returning the embedding
   * object the DFS predicate drives. Defaults to
   * EmbeddedSubmanifold<dim,subdim> (see EmbeddingSearch::search());
   * SurfaceSearch supplies KnottedSurface instead.
   * \param makeThreadHook factory called once per worker thread,
   * returning an object with `onFound(embedding, U, faceCount)` (called
   * for each candidate satisfying `cond`) and `onFlush()` (called
   * whenever that thread's counters flush to the shared accumulators,
   * including its final flush).
   * \param onSeedFound called at most once, if `isSeeded_` and the seed
   * alone already satisfies `cond`.
   * \param callbacks see SearchCallbacks.
   * \param auxHooks `spawn(workersFinished)` starts an optional thread
   * running alongside the workers (returning a default-constructed,
   * non-joinable std::thread if none is needed); `afterJoin()` runs once
   * that thread has been joined.
   * \param iddfsIterations see EmbeddingSearch::search(); 0 (the default)
   * means a single unbounded pass, this method's original behavior.
   * \param iddfsStep see EmbeddingSearch::search().
   * \param iddfsStart see EmbeddingSearch::search(); defaults to
   * `iddfsStep` when not given.
   * \param finalThreads see EmbeddingSearch::search(); defaults to
   * `numThreads` when not given.
   *
   * \note Each hot-path hook is a compile-time (lambda/functor) type,
   * not a virtual call, so the no-op hooks EmbeddingSearch::search()
   * passes cost nothing in the per-candidate hot path.
   *
   * \note SIGINT handling: for the duration of this call, Ctrl+C sets a
   * shared stop flag (see EmbeddednessPredicate) that prunes the DFS
   * everywhere, so worker threads wind down quickly and control falls
   * through to `auxHooks.afterJoin()` instead of continuing to search. A
   * second Ctrl+C -- at any point until this call returns, including
   * during afterJoin() -- terminates the process immediately.
   */
  template <typename EmbeddingFactory, typename ThreadHookFactory,
            typename OnSeedFound, typename AuxHooks>
  SearchStats runSearch_(unsigned numThreads, BoundaryCondition cond,
                         EmbeddingFactory makeEmbedding,
                         ThreadHookFactory makeThreadHook,
                         OnSeedFound onSeedFound,
                         const SearchCallbacks &callbacks, AuxHooks auxHooks,
                         unsigned iddfsIterations = 0, long long iddfsStep = 0,
                         std::optional<long long> iddfsStart = std::nullopt,
                         std::optional<unsigned> finalThreads = std::nullopt);

private:
  /** Builds the face-adjacency graph of `skeleton`'s subdim-faces (unseeded). */
  static Graph buildGraph_(const Skeleton<dim, subdim> &skeleton);

  /**
   * As buildGraph_(), then contracts `seedFaces` into a single vertex
   * (always ending up with the globally-smallest id, 1) via
   * ConnectedInducedSubgraphEnumerator::contractSeed().
   */
  static Graph buildSeededGraph_(const Skeleton<dim, subdim> &skeleton,
                                 const std::vector<int> &seedFaces);
};

extern template class EmbeddingSearch<3, 2>;
extern template class EmbeddingSearch<4, 2>;

/**
 * Everything needed to describe one found surface when no boundary-link
 * data is being tracked (BoundaryCondition::all/closed), handed to
 * SurfaceSearchCallbacks::onSurfaceFound.
 */
struct SurfaceFoundInfo {
  bool orientable; /**< Whether the surface is orientable. */
  int genus; /**< The surface's genus. */
  int punctures; /**< The surface's number of punctures (boundary components). */
  long long triangleCount; /**< The surface's number of triangles. */
  BoundaryCondition mostRestrictive;
      /**< The most restrictive BoundaryCondition this specific surface
           satisfies -- not necessarily the one the whole search was run
           under. Only isClosed()/isProper() are checked (both O(1)), not
           the O(ambient triangulation size)
           boundaryComponentsMapInjectively(), so this is never
           "connected". */
};

/**
 * SurfaceFoundInfo plus the boundary-link description, handed to
 * SurfaceSearchCallbacks::onSurfaceBoundaryProcessed once boundary links
 * are being tracked (BoundaryCondition::proper/connected).
 */
struct SurfaceBoundaryInfo : SurfaceFoundInfo {
  std::string boundaryDescription;
      /**< A formatted description of the surface's boundary; empty if
           the surface turned out closed. */
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
  std::function<void(std::chrono::steady_clock::duration elapsed)>
      onBoundaryProcessingComplete;
      /**< Fired once, when boundary-link post-processing finishes. */
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

  /** A thread-safe queue of surfaces (as face-index lists) awaiting boundary-link processing. */
  class PendingSurfaceBatch {
  private:
    mutable std::mutex mutex_;
    std::vector<std::vector<int>> pending_;

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

  public:
    /** Records one surface of `type` found with boundary `descriptor`. */
    void record(const std::string &descriptor, const SurfaceTypeKey &type);

    /** Returns a human-readable summary of the current tally. */
    std::string summary() const;
  };

  /**
   * Per-worker-thread hook passed to runSearch_(): tallies surface types
   * and (when boundary links are wanted) batches each find's face
   * indices for later link identification, merging into the owning
   * SurfaceSearch's shared accumulators whenever the harness flushes
   * this thread's counters.
   */
  class ThreadHook {
    SurfaceSearch &owner_;
    SurfaceTypeTally &tally_;
    bool wantLinks_;
    const SurfaceSearchCallbacks &callbacks_;
    std::map<SurfaceTypeKey, long long> localTypeCounts_;
    std::vector<std::vector<int>> localPending_;

  public:
    ThreadHook(SurfaceSearch &owner, SurfaceTypeTally &tally, bool wantLinks,
              const SurfaceSearchCallbacks &callbacks)
        : owner_(owner), tally_(tally), wantLinks_(wantLinks),
          callbacks_(callbacks) {}

    void onFound(EmbeddedSubmanifold<4, 2> &embedding,
                const std::vector<int> &U, long long faceCount);

    void onFlush();
  };

  /**
   * Aux-thread hook passed to runSearch_(): continuously drains
   * pendingSurfaces_ into boundary links for as long as the search runs
   * (see backgroundDrainLoop_), then hands off whatever it hasn't gotten
   * to once the workers finish.
   */
  class AuxHooks {
    SurfaceSearch &owner_;
    unsigned numThreads_;
    bool wantLinks_;
    const SurfaceSearchCallbacks &callbacks_;

  public:
    AuxHooks(SurfaceSearch &owner, unsigned numThreads, bool wantLinks,
             const SurfaceSearchCallbacks &callbacks)
        : owner_(owner), numThreads_(numThreads), wantLinks_(wantLinks),
          callbacks_(callbacks) {}

    std::thread spawn(std::atomic<bool> &workersFinished);

    void afterJoin();
  };

  PendingSurfaceBatch pendingSurfaces_; /**< Surfaces awaiting boundary-link processing. */
  LinkBoundaryTally linkTally_; /**< See linkTally(). */
  SurfaceTypeTally surfaceTypeTally_; /**< See surfaceTypeTally(). */

  /**
   * Shared across every KnottedSurface this search constructs (the main
   * DFS workers, the seed-validation/surfaceType probes, and the
   * boundary-link-processing workers) -- see PetalCache and
   * KnottedSurface's external-cache constructors. One instance regardless
   * of thread count, rather than one full cache per thread.
   */
  PetalCache petalCache_;

public:
  using EmbeddingSearch<4, 2>::EmbeddingSearch;

  /**
   * As the inherited seeded constructor, but additionally validates the
   * seed via a temporary KnottedSurface (not just a temporary
   * EmbeddedSubmanifold<4,2>), so a seed baked with an illegal
   * crossing/knot -- not just a facet-level collision -- is caught
   * eagerly too.
   */
  SurfaceSearch(const regina::Triangulation<4> &tri,
               const std::vector<int> &seedFaces);

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

  /** As EmbeddingSearch::search(), reporting through `callbacks`. */
  SearchStats search(unsigned numThreads,
                     BoundaryCondition cond = BoundaryCondition::all,
                     const SurfaceSearchCallbacks &callbacks = {},
                     unsigned iddfsIterations = 0, long long iddfsStep = 0,
                     std::optional<long long> iddfsStart = std::nullopt,
                     std::optional<unsigned> finalThreads = std::nullopt);

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
   */
  void processBatchRange_(KnottedSurface &embedding,
                          const std::vector<std::vector<int>> &batch,
                          size_t begin, size_t end,
                          const SurfaceSearchCallbacks &callbacks);

  /**
   * Builds one human-readable descriptor of a surface's whole boundary
   * from `links` (already grouped and sorted by ambient boundary
   * component index): for each component, "<component+1>: " followed by
   * the comma-joined identify() of every curve of that component's
   * Link, followed by " (<link.identify()>)" only when that component
   * holds more than one curve -- with per-component entries joined by
   * ", ". E.g. "1: figure-eight, 2: unknot, unknot (Hopf link)".
   */
  static std::string
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

#endif // EMBEDDINGSEARCH_H
