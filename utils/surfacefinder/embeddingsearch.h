#ifndef EMBEDDINGSEARCH_H

#define EMBEDDINGSEARCH_H

#include <chrono>
#include <map>
#include <mutex>
#include <tuple>

#include <triangulation/dim2.h>
#include <unordered_map>

#include "embeddedsubmanifold.h"
#include "enumerate_cis.h"
#include "skeleton.h"

using AdjacencyList = std::pair<int, std::vector<std::vector<int>>>;

using SurfaceTypeKey = std::tuple<bool, int, int>;

SurfaceTypeKey surfaceTypeKey(const regina::Triangulation<2> &surface);

std::string formatSurfaceType(const SurfaceTypeKey &key);

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

std::string formatElapsed(std::chrono::steady_clock::duration d);

const char *boundaryConditionName(BoundaryCondition cond);

template <int dim, int subdim>
class EmbeddednessPredicate : public ConditionalPredicate {
  EmbeddedSubmanifold<dim, subdim> &embedding;
  const std::vector<int> &graphToSkel;

public:
  EmbeddednessPredicate(EmbeddedSubmanifold<dim, subdim> &embedding,
                        const std::vector<int> &denseToOriginal);

  bool tryAdd(int v) override;

  void undo(int v) override;
};

template <int dim, int subdim> class EmbeddingSearch {
private:
  struct Graph {
    AdjacencyList adjList;
    std::vector<int> graphToSkel;
  };

  const Skeleton<dim, subdim> skeleton_;
  const Graph graph_;
  std::vector<EmbeddedSubmanifold<dim, subdim>> valid_embeddings_;

  PendingSurfaceBatch pendingSurfaces_;
  LinkBoundaryTally linkTally_;

public:
  EmbeddingSearch(const regina::Triangulation<dim> &tri);

  size_t numEmbeddableFaces() const { return graph_.graphToSkel.size(); }

  const LinkBoundaryTally &linkTally() const { return linkTally_; }

  void processSurfaceBoundaries();

  // Drains whatever is left in pendingSurfaces_ once (nothing else is
  // producing into it once the DFS workers have joined) and works through
  // it with numThreads worker threads -- the same threads search() used
  // for DFS, now idle. Meant to be called right after search()'s DFS
  // workers finish, as the (potentially very large) remaining backlog of
  // queued surfaces can dwarf the amount processed incrementally by
  // search()'s periodic single-threaded boundaryProcessor thread. Prints
  // linkTally_'s incremental summary once a second while it works, since
  // this can take a long time and would otherwise look like a hang.
  void processRemainingSurfaceBoundaries(unsigned numThreads);

  void search(const unsigned numThreads,
              BoundaryCondition cond = BoundaryCondition::all);

private:
  // Replays batch[begin, end) into embedding (via addFace(), in the exact
  // order the DFS originally discovered each entry -- deterministic,
  // since it's the same skeleton_/graph_), extracts its boundary Link(s),
  // and records what surface type each recognized link bounds into
  // linkTally_. embedding is reused across the whole range so its
  // bdryComponents_ cache (see EmbeddedSubmanifold's constructor) is only
  // built once, not once per queued surface.
  //
  // Only ever called from processSurfaceBoundaries()/
  // processRemainingSurfaceBoundaries() (both dim == 4 && subdim == 2
  // only, via static_assert), since it unconditionally calls
  // EmbeddedSubmanifold::boundaryLinks(), which itself is only defined
  // (explicitly instantiated) for <4, 2> -- see embeddedsubmanifold.cpp.
  void processBatchRange_(EmbeddedSubmanifold<dim, subdim> &embedding,
                          const std::vector<std::vector<int>> &batch,
                          size_t begin, size_t end);

  static Graph buildGraph_(const Skeleton<dim, subdim> &skeleton);
};

extern template class EmbeddednessPredicate<3, 2>;
extern template class EmbeddednessPredicate<4, 2>;

extern template class EmbeddingSearch<4, 2>;

#endif // EMBEDDINGSEARCH_H
