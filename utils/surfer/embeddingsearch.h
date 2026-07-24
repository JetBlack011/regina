#ifndef EMBEDDINGSEARCH_H

#define EMBEDDINGSEARCH_H

#include <atomic>
#include <chrono>
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

template <int dim, int subdim> class EmbeddingSearch {
protected:
  class EmbeddednessPredicate : public ConditionalPredicate {
    EmbeddedSubmanifold<dim, subdim> &embedding_;
    const std::vector<std::vector<int>> &graphToSkel_;

  public:
    EmbeddednessPredicate(EmbeddedSubmanifold<dim, subdim> &embedding,
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

  void search(const unsigned numThreads,
              BoundaryCondition cond = BoundaryCondition::all);

protected:
  // Shared harness behind both EmbeddingSearch::search() and
  // SurfaceSearch::search(): the two differ only in a handful of injection
  // points, all threaded through here so the threading/atomics/reporting
  // code (the bulk of either function) exists exactly once. Each hook
  // parameter is a compile-time (lambda/functor) type, not a virtual call,
  // so the no-op hooks EmbeddingSearch::search() passes cost nothing in the
  // per-candidate hot path.
  //
  //   makeThreadHook()  -- factory called once per worker thread, returning
  //                        an object with onFound(embedding, U, faceCount)
  //                        (called for each candidate satisfying cond) and
  //                        onFlush() (called whenever that thread's counters
  //                        flush to the globals, including its final
  //                        flush -- the place to merge thread-local state
  //                        into shared accumulators).
  //   onSeedFound(seedFaces) -- called at most once, if isSeeded_ and the
  //                        seed alone already satisfies cond.
  //   extraReportLines()  -- extra text appended to both the once-a-second
  //                        progress report and the final summary.
  //   auxHooks            -- spawn(workersFinished) starts an optional
  //                        thread running alongside the workers (returning
  //                        a default-constructed, non-joinable
  //                        std::thread if none is needed); afterJoin() runs
  //                        once that thread has been joined.
  template <typename ThreadHookFactory, typename OnSeedFound,
            typename ExtraReportLines, typename AuxHooks>
  void runSearch_(unsigned numThreads, BoundaryCondition cond,
                  ThreadHookFactory makeThreadHook, OnSeedFound onSeedFound,
                  ExtraReportLines extraReportLines, AuxHooks auxHooks);

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
    std::map<SurfaceTypeKey, long long> localTypeCounts_;
    std::vector<std::vector<int>> localPending_;

  public:
    ThreadHook(SurfaceSearch &owner, SurfaceTypeTally &tally, bool wantLinks)
        : owner_(owner), tally_(tally), wantLinks_(wantLinks) {}

    void onFound(EmbeddedSubmanifold<4, 2> &embedding,
                const std::vector<int> &U, long long faceCount);

    void onFlush();
  };

  // Aux-thread hook passed to runSearch_(): periodically drains
  // pendingSurfaces_ into boundary links while the search runs, then
  // processes whatever is left once the workers finish.
  class AuxHooks {
    SurfaceSearch &owner_;
    unsigned numThreads_;
    bool wantLinks_;

  public:
    AuxHooks(SurfaceSearch &owner, unsigned numThreads, bool wantLinks)
        : owner_(owner), numThreads_(numThreads), wantLinks_(wantLinks) {}

    std::thread spawn(std::atomic<bool> &workersFinished);

    void afterJoin();
  };

  PendingSurfaceBatch pendingSurfaces_;
  LinkBoundaryTally linkTally_;

public:
  using EmbeddingSearch<4, 2>::EmbeddingSearch;

  const LinkBoundaryTally &linkTally() const { return linkTally_; }

  void search(unsigned numThreads,
              BoundaryCondition cond = BoundaryCondition::all);

  void processSurfaceBoundaries();

  void processRemainingSurfaceBoundaries(unsigned numThreads);

private:
  // Replays batch[begin, end) into embedding (via addFace(), in the exact
  // order the DFS originally discovered each entry -- deterministic, since
  // it's the same skeleton_), extracts its boundary Link(s), and records
  // one boundary descriptor (see describeBoundary_) per surface, tagged
  // with its surface type, into linkTally_. embedding is reused across the
  // whole range so its bdryComponents_ cache (see KnottedSurface's
  // constructor) is only built once, not once per queued surface.
  void processBatchRange_(KnottedSurface &embedding,
                          const std::vector<std::vector<int>> &batch,
                          size_t begin, size_t end);

  // Builds one human-readable descriptor of a surface's whole boundary from
  // boundaryLinks() (already grouped and sorted by ambient boundary
  // component index): for each component, "<component+1>: " followed by
  // the comma-joined identify() of every curve of that component's Link,
  // followed by " (<link.identify()>)" only when that component holds more
  // than one curve -- and the per-component entries themselves joined by
  // ", ". E.g. "1: figure-eight, 2: unknot, unknot (Hopf link)".
  static std::string
  describeBoundary_(const std::vector<std::pair<size_t, Link>> &links);
};

#endif // EMBEDDINGSEARCH_H
