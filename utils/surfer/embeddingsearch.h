#ifndef EMBEDDINGSEARCH_H

#define EMBEDDINGSEARCH_H

#include <atomic>
#include <chrono>
#include <functional>
#include <map>
#include <mutex>
#include <thread>

#include <triangulation/dim2.h>
#include <unordered_map>

#include "embeddedsubmanifold.h"
#include "enumerate_cis.h"
#include "skeleton.h"

using AdjacencyList = std::pair<int, std::vector<std::vector<int>>>;

std::string formatElapsed(std::chrono::steady_clock::duration d);

const char *boundaryConditionName(BoundaryCondition cond);

// A snapshot of a search() call's progress/results, handed to SearchCallbacks
// rather than printed directly -- EmbeddingSearch has no opinion on how (or
// whether) this gets displayed.
struct SearchStats {
  std::chrono::steady_clock::duration elapsed{};
  size_t rootsCompleted = 0;
  size_t totalRoots = 0;
  long long foundCount = 0;         // raw candidates visited, any cond
  long long satisfyingCount = 0;    // satisfying the BoundaryCondition
  long long satisfyingFaceSum = 0;
  long long largestSatisfying = 0;  // max faces among satisfying finds

  double averageSatisfyingFaces() const {
    return satisfyingCount > 0
               ? static_cast<double>(satisfyingFaceSum) /
                     static_cast<double>(satisfyingCount)
               : 0.0;
  }
};

// Callbacks a caller of search() may supply to observe its progress instead
// of search() printing anything itself. Every field defaults to an empty
// std::function (a safe no-op), so passing none of these costs nothing
// beyond the once-a-second call overhead already inherent to the reporter.
struct SearchCallbacks {
  // Fired roughly once a second while the search runs.
  std::function<void(const SearchStats &)> onProgress;
  // Fired at most once, if Ctrl+C is caught mid-search.
  std::function<void()> onInterrupted;
  // Fired once, when the search has finished (the same SearchStats is also
  // returned by search() itself).
  std::function<void(const SearchStats &)> onSearchComplete;
};

template <int dim, int subdim> class EmbeddingSearch {
protected:
  // Generic over the embedding type so that SurfaceSearch can supply
  // KnottedSurface (whose addFace()/addFaces()/removeFace() override --
  // name-hide, non-virtually -- the base EmbeddedSubmanifold<4,2>'s to add
  // transverse-self-intersection/local-flatness checks) while
  // EmbeddingSearch<dim,subdim> itself never has to reference
  // KnottedSurface: since nothing here is virtual, resolving to the
  // override requires the predicate's *declared* member type to already
  // be the derived type, not just the object underneath it. C++17 class
  // template argument deduction means every existing construction call
  // site (`EmbeddednessPredicate predicate(embedding, ...)`) still needs
  // no explicit template argument -- EmbeddingT is deduced from whatever
  // makeEmbedding() (see runSearch_ below) returned.
  template <typename EmbeddingT>
  class EmbeddednessPredicate : public ConditionalPredicate {
    EmbeddingT &embedding_;
    const std::vector<std::vector<int>> &graphToSkel_;

  public:
    EmbeddednessPredicate(EmbeddingT &embedding,
                          const std::vector<std::vector<int>> &graphToSkel);

    bool tryAdd(int v) override;

    void undo(int v) override;
  };

  struct Graph {
    AdjacencyList adjList;
    // graphToSkel[v-1] is the list of skeleton face indices graph vertex v
    // corresponds to. Every vertex has a singleton list, except an
    // unseeded graph's vertex 1 has none: in the seeded case (see
    // buildSeededGraph_ below), vertex 1 is the contracted seed and maps to
    // every one of its skeleton faces at once.
    std::vector<std::vector<int>> graphToSkel;
  };

  const Skeleton<dim, subdim> skeleton_;
  const Graph graph_;
  const bool isSeeded_ = false;

private:
  std::vector<EmbeddedSubmanifold<dim, subdim>> valid_embeddings_;

public:
  EmbeddingSearch(const regina::Triangulation<dim> &tri);

  // Seeded search: seedFaces (skeleton face indices) must be a connected,
  // jointly-addable set of faces (e.g. the output of CollarBuilder,
  // translated to skeleton indices via Triangle::index()) -- checked
  // eagerly here, throwing regina::InvalidArgument on an invalid seed
  // rather than deferring discovery to search().
  EmbeddingSearch(const regina::Triangulation<dim> &tri,
                  const std::vector<int> &seedFaces);

  size_t numEmbeddableFaces() const { return graph_.graphToSkel.size(); }

  SearchStats search(const unsigned numThreads,
                     BoundaryCondition cond = BoundaryCondition::all,
                     const SearchCallbacks &callbacks = {});

protected:
  // Shared harness behind both EmbeddingSearch::search() and
  // SurfaceSearch::search(): the two differ only in a handful of injection
  // points, all threaded through here so the threading/atomics/reporting
  // code (the bulk of either function) exists exactly once. Each hot-path
  // hook parameter is a compile-time (lambda/functor) type, not a virtual
  // call, so the no-op hooks EmbeddingSearch::search() passes cost nothing
  // in the per-candidate hot path; callbacks (fired at most once a second)
  // are plain std::functions, since their overhead is irrelevant at that
  // rate.
  //
  //   makeEmbedding()   -- factory called once per worker thread (plus once
  //                        more for the seeded-root prototype below),
  //                        returning the embedding object the DFS predicate
  //                        drives. Defaults to EmbeddedSubmanifold<dim,subdim>
  //                        (see EmbeddingSearch::search()); SurfaceSearch
  //                        supplies KnottedSurface instead, alongside its
  //                        other hooks below -- this is the only place
  //                        EmbeddingSearch<dim,subdim> would otherwise need
  //                        to know KnottedSurface exists, so it's injected
  //                        the same way as the other three hooks rather
  //                        than referenced directly.
  //   makeThreadHook()  -- factory called once per worker thread, returning
  //                        an object with onFound(embedding, U, faceCount)
  //                        (called for each candidate satisfying cond) and
  //                        onFlush() (called whenever that thread's counters
  //                        flush to the globals, including its final
  //                        flush -- the place to merge thread-local state
  //                        into shared accumulators).
  //   onSeedFound(seedFaces) -- called at most once, if isSeeded_ and the
  //                        seed alone already satisfies cond.
  //   callbacks           -- see SearchCallbacks: onProgress fires roughly
  //                        once a second, onInterrupted fires at most once
  //                        on Ctrl+C, onSearchComplete fires once at the
  //                        end with the same SearchStats this returns.
  //   auxHooks            -- spawn(workersFinished) starts an optional
  //                        thread running alongside the workers (returning
  //                        a default-constructed, non-joinable
  //                        std::thread if none is needed); afterJoin() runs
  //                        once that thread has been joined.
  //
  // SIGINT handling: for the duration of this call, Ctrl+C sets a shared
  // stop flag (see EmbeddednessPredicate) that prunes the DFS everywhere,
  // so the worker threads wind down quickly and control falls through to
  // whatever auxHooks.afterJoin() does next (e.g. SurfaceSearch's boundary
  // batch processing) instead of continuing to search. A second Ctrl+C --
  // at any point until this call returns, including during afterJoin() --
  // terminates the process immediately.
  template <typename EmbeddingFactory, typename ThreadHookFactory,
            typename OnSeedFound, typename AuxHooks>
  SearchStats runSearch_(unsigned numThreads, BoundaryCondition cond,
                         EmbeddingFactory makeEmbedding,
                         ThreadHookFactory makeThreadHook,
                         OnSeedFound onSeedFound,
                         const SearchCallbacks &callbacks, AuxHooks auxHooks);

private:
  static Graph buildGraph_(const Skeleton<dim, subdim> &skeleton);

  // Builds the graph as buildGraph_ does, then contracts seedFaces into a
  // single vertex (always ending up with the globally-smallest id, 1) via
  // ConnectedInducedSubgraphEnumerator::contractSeed(), composing its id
  // mapping with buildGraph_'s own graphToSkel to produce a Graph whose
  // vertex 1 maps to the whole seed and every other vertex maps to a single
  // skeleton face, exactly as buildGraph_ would have.
  static Graph buildSeededGraph_(const Skeleton<dim, subdim> &skeleton,
                                 const std::vector<int> &seedFaces);
};

extern template class EmbeddingSearch<3, 2>;
extern template class EmbeddingSearch<4, 2>;

// Per-surface data handed to SurfaceSearchCallbacks::onSurfaceFound --
// everything needed to describe one found surface when no boundary-link
// data is being tracked (BoundaryCondition::all/closed).
struct SurfaceFoundInfo {
  bool orientable;
  int genus;
  int punctures;
  long long triangleCount;
  // The most restrictive BoundaryCondition this specific surface satisfies
  // -- not necessarily the one the whole search was run under. Only
  // isClosed() and isProper() (both O(1)) are checked here, not
  // boundaryComponentsMapInjectively() (O(size of the whole ambient
  // triangulation)) -- so this is never "connected", only all/closed/proper
  // -- see onSurfaceFound's doc comment below for why.
  BoundaryCondition mostRestrictive;
};

// Per-surface data handed to SurfaceSearchCallbacks::onSurfaceBoundaryProcessed
// -- SurfaceFoundInfo plus the boundary-link description, once boundary
// links are being tracked (BoundaryCondition::proper/connected).
struct SurfaceBoundaryInfo : SurfaceFoundInfo {
  // SurfaceSearch::describeBoundary_()'s formatted text -- "" when this
  // surface turned out closed (no boundary at all).
  std::string boundaryDescription;
};

// SurfaceSearch's own callbacks: everything SearchCallbacks offers for the
// main DFS phase, plus per-surface callbacks and three more for the
// boundary-link post-processing phase (processRemainingSurfaceBoundaries)
// that runs after it -- not covered by EmbeddingSearch's runSearch_, since
// they're specific to SurfaceSearch.
struct SurfaceSearchCallbacks : SearchCallbacks {
  // Fired once per found surface when boundary links are NOT being tracked
  // (BoundaryCondition::all/closed) -- i.e. from the DFS hot path itself.
  // Deliberately cheap (see SurfaceFoundInfo::mostRestrictive): this can
  // fire once per DFS node visited under BoundaryCondition::all, since
  // satisfies(all) doesn't prune anything.
  std::function<void(const SurfaceFoundInfo &)> onSurfaceFound;
  // Fired once per found surface when boundary links ARE being tracked
  // (BoundaryCondition::proper/connected) -- from the deferred boundary-link
  // processing phase, mutually exclusive with onSurfaceFound (exactly one of
  // the two fires per found surface, depending on which BoundaryCondition
  // the whole search was run under).
  std::function<void(const SurfaceBoundaryInfo &)> onSurfaceBoundaryProcessed;
  // Fired once, when boundary-link post-processing begins.
  std::function<void(size_t total, unsigned numThreads)>
      onBoundaryProcessingStarted;
  // Fired roughly once a second while boundary-link post-processing runs.
  std::function<void(size_t processed, size_t total,
                     std::chrono::steady_clock::duration elapsed)>
      onBoundaryProcessingProgress;
  // Fired once, when boundary-link post-processing finishes.
  std::function<void(std::chrono::steady_clock::duration elapsed)>
      onBoundaryProcessingComplete;
};

// A search over KnottedSurfaces (EmbeddedSubmanifold<4, 2>) that
// additionally recognizes each found surface's boundary curves, per ambient
// boundary component they lie in (individually, plus the whole link within
// that component when it holds more than one curve), tallying which surface
// types have which per-component boundary identifications. This has its
// own search() (rather than delegating to EmbeddingSearch<4,2>::search())
// so that surface-type tallying and boundary-link batching can be
// interleaved into the same DFS pass and the same once-a-second progress
// report that EmbeddingSearch<4,2>::search() would otherwise handle alone --
// all the surface-specific bookkeeping lives here, not in EmbeddingSearch.
class SurfaceSearch : public EmbeddingSearch<4, 2> {
public:
  using SurfaceTypeKey = KnottedSurface::SurfaceTypeKey;

private:
  class SurfaceTypeTally {
  private:
    mutable std::mutex mutex_;
    std::map<SurfaceTypeKey, long long> counts_;

  public:
    void merge(std::map<SurfaceTypeKey, long long> &local);

    std::string summary() const;
  };

  class PendingSurfaceBatch {
  private:
    mutable std::mutex mutex_;
    std::vector<std::vector<int>> pending_;

  public:
    void merge(std::vector<std::vector<int>> &local);

    std::vector<std::vector<int>> drain();

    // Pops up to maxCount entries in one lock acquisition (fewer if that's
    // all that's available, none if empty) -- used by the continuous
    // background worker (see backgroundDrainLoop_) so it never claims more
    // than a small, bounded amount of work out of the shared queue at
    // once, keeping it quick to hand the rest off the moment the search
    // ends.
    std::vector<std::vector<int>> popSome(size_t maxCount);
  };

  // Tallies, per boundary descriptor (see describeBoundary_ below -- one
  // string capturing the identification of every curve in every ambient
  // boundary component a surface's boundary touches), how many surfaces of
  // each type were found with that exact descriptor.
  class LinkBoundaryTally {
  private:
    mutable std::mutex mutex_;
    std::unordered_map<std::string, std::map<SurfaceTypeKey, long long>>
        descriptorSurfaceTypes_;

  public:
    void record(const std::string &descriptor, const SurfaceTypeKey &type);

    std::string summary() const;
  };

  // Per-worker-thread hook passed to runSearch_(): tallies surface types
  // and (when boundary links are wanted) batches each find's face indices
  // for later link identification, merging into the owning SurfaceSearch's
  // shared accumulators whenever the harness flushes this thread's
  // counters.
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

  // Aux-thread hook passed to runSearch_(): continuously drains
  // pendingSurfaces_ into boundary links, in small increments, for as long
  // as the search runs (see backgroundDrainLoop_), then hands off whatever
  // it hasn't gotten to once the workers finish.
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

  PendingSurfaceBatch pendingSurfaces_;
  LinkBoundaryTally linkTally_;
  SurfaceTypeTally surfaceTypeTally_;

public:
  using EmbeddingSearch<4, 2>::EmbeddingSearch;

  // Shadows (hides) the inherited seeded constructor of the same
  // signature: EmbeddingSearch<4,2>'s own seeded constructor only
  // validates the seed at the facet level, via a temporary
  // EmbeddedSubmanifold<4,2> -- this additionally validates it via a
  // temporary KnottedSurface, so a seed baked with an illegal crossing/
  // knot (not just a facet-level collision) is caught eagerly too.
  SurfaceSearch(const regina::Triangulation<4> &tri,
               const std::vector<int> &seedFaces);

  const LinkBoundaryTally &linkTally() const { return linkTally_; }

  const SurfaceTypeTally &surfaceTypeTally() const {
    return surfaceTypeTally_;
  }

  SearchStats search(unsigned numThreads,
                     BoundaryCondition cond = BoundaryCondition::all,
                     const SurfaceSearchCallbacks &callbacks = {});

  // Called once after every DFS worker thread (and the background drain
  // loop below) has finished, to process whatever's left across up to
  // numThreads worker threads. Always announces itself (via the three
  // onBoundaryProcessing* callbacks) since by construction this is the
  // sole, complete pass over anything not already handled by
  // backgroundDrainLoop_ -- see that method's doc comment for why there's
  // no longer any question of this racing a still-in-flight earlier batch.
  void processRemainingSurfaceBoundaries(
      unsigned numThreads, const SurfaceSearchCallbacks &callbacks = {});

private:
  // Runs on the aux thread for as long as search() (with boundary links
  // wanted) is active: continuously pops and processes small batches (via
  // popSome(), never the whole queue) so pendingSurfaces_ never grows
  // unbounded, silently since the search's own once-a-second report
  // already reflects linkTally_'s growth. Checks workersFinished between
  // every individual entry (not just between batches), so once the DFS
  // itself finishes or is interrupted, this stops within roughly one
  // entry's processing time -- handing back any unprocessed remainder of
  // its current small batch to pendingSurfaces_ rather than finishing it
  // out -- instead of a large already-claimed batch silently grinding on
  // for however long it takes (the earlier failure mode this replaces).
  // processRemainingSurfaceBoundaries() then picks up everything left, in
  // one clean pass, once this has returned.
  void backgroundDrainLoop_(const std::atomic<bool> &workersFinished,
                            const SurfaceSearchCallbacks &callbacks);

  // Shared parallel-processing core behind processRemainingSurfaceBoundaries():
  // splits `batch` into CHUNK-sized slices off a shared atomic cursor,
  // processed by up to min(numThreads, batch.size()) worker threads (each
  // with its own reused KnottedSurface, via processBatchRange_), reporting
  // through the onBoundaryProcessing* callbacks throughout.
  void processBatchParallel_(std::vector<std::vector<int>> batch,
                            unsigned numThreads,
                            const SurfaceSearchCallbacks &callbacks);

  // Runs the addFace/boundaryLinks/record/removeFace sequence for one
  // queued surface's face indices against `embedding` (reused by the
  // caller across many entries), recording one boundary descriptor (see
  // describeBoundary_) into linkTally_ if the surface has any boundary at
  // all (closed surfaces contribute nothing here), and firing
  // callbacks.onSurfaceBoundaryProcessed with the same descriptor.
  void processEntry_(KnottedSurface &embedding,
                     const std::vector<int> &faceIndices,
                     const SurfaceSearchCallbacks &callbacks);

  // Replays batch[begin, end) via processEntry_, in the exact order the
  // DFS originally discovered each entry -- deterministic, since it's the
  // same skeleton_. embedding is reused across the whole range so its
  // bdryComponents_ cache (see KnottedSurface's constructor) is only built
  // once, not once per queued surface.
  void processBatchRange_(KnottedSurface &embedding,
                          const std::vector<std::vector<int>> &batch,
                          size_t begin, size_t end,
                          const SurfaceSearchCallbacks &callbacks);

  // Builds one human-readable descriptor of a surface's whole boundary from
  // boundaryLinks() (already grouped and sorted by ambient boundary
  // component index): for each component, "<component+1>: " followed by
  // the comma-joined identify() of every curve of that component's Link,
  // followed by " (<link.identify()>)" only when that component holds more
  // than one curve -- and the per-component entries themselves joined by
  // ", ". E.g. "1: figure-eight, 2: unknot, unknot (Hopf link)".
  static std::string
  describeBoundary_(const std::vector<std::pair<size_t, Link>> &links);

  // The most restrictive BoundaryCondition a surface satisfies, given its
  // already-computed boundaryLinks() -- closed (no boundary at all, i.e.
  // links is empty), connected (every ambient boundary component receives
  // at most one of the surface's own boundary curves -- exactly what
  // boundaryComponentsMapInjectively() checks, but free here since
  // boundaryLinks() already groups curves by ambient component), or proper
  // (otherwise). Shared by processEntry_ and search()'s onSeedFound so the
  // two never drift apart on how this is computed.
  static BoundaryCondition
  classifyByLinks_(const std::vector<std::pair<size_t, Link>> &links);

  // The most restrictive BoundaryCondition a surface satisfies, using only
  // isClosed()/isProper() (both O(1)) -- never "connected", since that
  // needs boundaryComponentsMapInjectively() (O(size of the whole ambient
  // triangulation), deliberately not paid in the hot path this feeds).
  // Shared by ThreadHook::onFound and search()'s onSeedFound.
  static BoundaryCondition
  classifyCheaply_(const EmbeddedSubmanifold<4, 2> &embedding);
};

#endif // EMBEDDINGSEARCH_H
