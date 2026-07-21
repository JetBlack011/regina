#ifndef EMBEDDINGSEARCH_H

#define EMBEDDINGSEARCH_H

#include <chrono>
#include <map>
#include <mutex>

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
    const std::vector<int> &graphToSkel_;

  public:
    EmbeddednessPredicate(EmbeddedSubmanifold<dim, subdim> &embedding,
                          const std::vector<int> &graphToSkel);

    bool tryAdd(int v) override;

    void undo(int v) override;
  };

  struct Graph {
    AdjacencyList adjList;
    std::vector<int> graphToSkel;
  };

  const Skeleton<dim, subdim> skeleton_;
  const Graph graph_;

private:
  std::vector<EmbeddedSubmanifold<dim, subdim>> valid_embeddings_;

public:
  EmbeddingSearch(const regina::Triangulation<dim> &tri);

  size_t numEmbeddableFaces() const { return graph_.graphToSkel.size(); }

  void search(const unsigned numThreads,
              BoundaryCondition cond = BoundaryCondition::all);

private:
  static Graph buildGraph_(const Skeleton<dim, subdim> &skeleton);
};

extern template class EmbeddingSearch<3, 2>;
extern template class EmbeddingSearch<4, 2>;

// A search over KnottedSurfaces (EmbeddedSubmanifold<4, 2>) that
// additionally recognizes the link bounded by each found surface's
// boundary, tallying which surface types bound which links. This has its
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

  class LinkBoundaryTally {
  private:
    mutable std::mutex mutex_;
    std::unordered_map<std::string, std::map<SurfaceTypeKey, long long>>
        linkSurfaceTypes_;

  public:
    void record(const std::string &linkName, const SurfaceTypeKey &type);

    std::string summary() const;
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
  // what surface type each recognized link bounds into linkTally_.
  // embedding is reused across the whole range so its bdryComponents_
  // cache (see KnottedSurface's constructor) is only built once, not once
  // per queued surface.
  void processBatchRange_(KnottedSurface &embedding,
                          const std::vector<std::vector<int>> &batch,
                          size_t begin, size_t end);
};

#endif // EMBEDDINGSEARCH_H
