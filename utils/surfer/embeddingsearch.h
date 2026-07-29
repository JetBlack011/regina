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
#include <memory>
#include <optional>
#include <thread>
#include <vector>

#include <triangulation/dim2.h>

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
 * Per-worker-thread hook interface for EmbeddingSearch::runSearch_(): fired
 * on each candidate satisfying the search's BoundaryCondition (onFound) and
 * whenever a thread's counters flush to the shared accumulators, including
 * its final flush (onFlush).
 *
 * A genuine virtual interface, unlike the embedding type itself
 * (EmbeddednessPredicate's EmbeddingT, see below) -- onFound()/onFlush()
 * never touch anything beyond EmbeddedSubmanifold<dim,subdim>'s own base
 * methods (triangulation(), isClosed(), isProper()), regardless of whether
 * the search is actually driving a plain EmbeddedSubmanifold or a
 * KnottedSurface, so there's no correctness reason for this to be a
 * template parameter. Making it virtual instead is what lets runSearch_()
 * itself be templated on EmbeddingT alone (see its own doc comment) and
 * hence explicitly instantiated in embeddingsearch.cpp for the small, fixed
 * set of concrete embedding types actually used -- a real declaration in
 * the header, a real definition in the .cpp, no lambda-closure-type
 * visibility problem. The extra indirect call this costs only happens once
 * per *satisfying* candidate (see FLUSH_EVERY_BDRY), not once per raw DFS
 * candidate, so it's nowhere near the hot path EmbeddednessPredicate's
 * addFace()/removeFace() calls are on.
 */
template <int dim, int subdim>
class RunSearchThreadHook {
public:
  virtual void onFound(EmbeddedSubmanifold<dim, subdim> &embedding,
                       const std::vector<int> &U, long long faceCount) = 0;
  virtual void onFlush() = 0;
  virtual ~RunSearchThreadHook() = default;
};

/**
 * Aux-thread hook interface for EmbeddingSearch::runSearch_(): `spawn()`
 * starts an optional thread running alongside the workers (returning a
 * default-constructed, non-joinable std::thread if none is needed);
 * `afterJoin()` runs once that thread has been joined. Never
 * embedding-type-specific either (see RunSearchThreadHook for why that
 * matters for runSearch_()'s own separate-compilation story) -- SurfaceSearch's
 * own AuxHooks builds whatever KnottedSurface objects it needs internally,
 * independent of runSearch_()'s own EmbeddingT.
 */
class RunSearchAuxHooks {
public:
  virtual std::thread spawn(std::atomic<bool> &workersFinished) = 0;
  virtual void afterJoin() = 0;
  virtual ~RunSearchAuxHooks() = default;
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

  /**
   * Set for the duration of one runSearch_() call in response to SIGINT
   * (see runSearch_'s SIGINT-handling doc note), reset to false at that
   * call's start. A class member (not a runSearch_()-local variable) so
   * SurfaceSearch's ThreadHook::onFlush() -- which only has access to
   * *this via its owner_ reference, not to runSearch_'s own locals -- can
   * check it too, to cut a full-drain-to-empty pass (see pauseRequested_
   * below) short on interrupt rather than insisting on draining
   * everything first.
   */
  std::atomic<bool> stopRequested_{false};

  /**
   * Sits alongside stopRequested_ but with different semantics: set (see
   * SurfaceSearch::ThreadHook::onFlush()) while a full drain-to-empty pass
   * of a pending-surface queue is in progress, so every DFS worker's
   * tryAdd() *blocks* (not rejects -- see InterruptiblePredicate) until it
   * clears, pausing candidate production without losing any embeddings a
   * rejected tryAdd() would have permanently pruned. Always false, and
   * never set by anyone, for the base EmbeddingSearch<dim,subdim>::search()
   * itself -- only SurfaceSearch has a pending-surface queue to drain.
   */
  std::atomic<bool> pauseRequested_{false};

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
   * \tparam EmbeddingT the concrete embedding type the DFS predicate drives
   * (EmbeddedSubmanifold<dim,subdim> itself, see EmbeddingSearch::search(),
   * or KnottedSurface, see SurfaceSearch::search()) -- genuinely a template
   * parameter (not type-erased like the hook parameters below), since
   * EmbeddednessPredicate<EmbeddingT> calls addFace()/removeFace()
   * non-virtually; resolving to KnottedSurface's enhanced overrides
   * requires EmbeddingT to already be the derived type at compile time.
   * This is also runSearch_()'s *only* template parameter: declared here,
   * defined (and explicitly instantiated for every EmbeddingT actually
   * used) in embeddingsearch.cpp, exactly like any other member function
   * of an explicitly-instantiated class template -- unlike the hook
   * parameters below, EmbeddingT is a nameable, fixed set of concrete
   * types, so this doesn't have the lambda-closure-visibility problem
   * EmbeddednessPredicate's own doc comment describes.
   *
   * \param numThreads the number of worker threads to search with.
   * \param cond restricts results to embeddings satisfying this
   * condition.
   * \param makeEmbedding factory called once per worker thread (plus
   * once more for the seeded-root prototype), returning the embedding
   * object the DFS predicate drives.
   * \param makeThreadHook factory called once per worker thread,
   * returning a RunSearchThreadHook.
   * \param onSeedFound called at most once, if `isSeeded_` and the seed
   * alone already satisfies `cond`.
   * \param callbacks see SearchCallbacks.
   * \param auxHooks see RunSearchAuxHooks; must outlive this call.
   * \param iddfsIterations see EmbeddingSearch::search(); 0 (the default)
   * means a single unbounded pass, this method's original behavior.
   * \param iddfsStep see EmbeddingSearch::search().
   * \param iddfsStart see EmbeddingSearch::search(); defaults to
   * `iddfsStep` when not given.
   * \param finalThreads see EmbeddingSearch::search(); defaults to
   * `numThreads` when not given.
   *
   * \note SIGINT handling: for the duration of this call, Ctrl+C sets a
   * shared stop flag (see EmbeddednessPredicate) that prunes the DFS
   * everywhere, so worker threads wind down quickly and control falls
   * through to `auxHooks.afterJoin()` instead of continuing to search. A
   * second Ctrl+C -- at any point until this call returns, including
   * during afterJoin() -- terminates the process immediately.
   */
  template <typename EmbeddingT>
  SearchStats runSearch_(
      unsigned numThreads, BoundaryCondition cond,
      std::function<EmbeddingT()> makeEmbedding,
      std::function<std::unique_ptr<RunSearchThreadHook<dim, subdim>>()>
          makeThreadHook,
      std::function<void(const std::vector<int> &)> onSeedFound,
      const SearchCallbacks &callbacks, RunSearchAuxHooks &auxHooks,
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

// Explicit instantiations of runSearch_(), one per concrete EmbeddingT
// actually used (see runSearch_()'s own doc comment for why this is
// possible despite it being a member template): defined in
// embeddingsearch.cpp, so any other translation unit calling search() --
// e.g. surfacesearch.cpp, for SurfaceSearch<KnottedSurface> -- just needs
// these declarations, not the ~350-line definition itself.
extern template SearchStats
EmbeddingSearch<3, 2>::runSearch_<EmbeddedSubmanifold<3, 2>>(
    unsigned, BoundaryCondition, std::function<EmbeddedSubmanifold<3, 2>()>,
    std::function<std::unique_ptr<RunSearchThreadHook<3, 2>>()>,
    std::function<void(const std::vector<int> &)>, const SearchCallbacks &,
    RunSearchAuxHooks &, unsigned, long long, std::optional<long long>,
    std::optional<unsigned>);
extern template SearchStats
EmbeddingSearch<4, 2>::runSearch_<EmbeddedSubmanifold<4, 2>>(
    unsigned, BoundaryCondition, std::function<EmbeddedSubmanifold<4, 2>()>,
    std::function<std::unique_ptr<RunSearchThreadHook<4, 2>>()>,
    std::function<void(const std::vector<int> &)>, const SearchCallbacks &,
    RunSearchAuxHooks &, unsigned, long long, std::optional<long long>,
    std::optional<unsigned>);
extern template SearchStats
EmbeddingSearch<4, 2>::runSearch_<KnottedSurface>(
    unsigned, BoundaryCondition, std::function<KnottedSurface()>,
    std::function<std::unique_ptr<RunSearchThreadHook<4, 2>>()>,
    std::function<void(const std::vector<int> &)>, const SearchCallbacks &,
    RunSearchAuxHooks &, unsigned, long long, std::optional<long long>,
    std::optional<unsigned>);

#endif // EMBEDDINGSEARCH_H
