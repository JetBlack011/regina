#ifndef EMBEDDINGSEARCH_H

#define EMBEDDINGSEARCH_H

#include <algorithm>
#include <atomic>
#include <chrono>
#include <iomanip>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <thread>
#include <tuple>

#include <triangulation/dim2.h>
#include <unordered_map>

#include "embeddedsubmanifold.h"
#include "enumerate_cis.h"
#include "linkcomplement.h"
#include "skeleton.h"

using AdjacencyList = std::pair<int, std::vector<std::vector<int>>>;

using SurfaceTypeKey = std::tuple<bool, int, int>;

inline SurfaceTypeKey surfaceTypeKey(const regina::Triangulation<2> &surface) {
  bool isOrientable = surface.isOrientable();
  int punctures = surface.countBoundaryComponents();
  int genus = isOrientable ? (2 - surface.eulerChar() - punctures) / 2
                           : 2 - surface.eulerChar() - punctures;
  return {isOrientable, genus, punctures};
}

inline std::string formatSurfaceType(const SurfaceTypeKey &key) {
  auto [isOrientable, genus, punctures] = key;
  std::ostringstream ans;

  if (isOrientable) {
    if (genus == 0 && punctures == 1)
      ans << "Disc";
    else if (genus == 0 && punctures == 2)
      ans << "Annulus";
    else {
      if (genus == 0)
        ans << "Sphere";
      else if (genus == 1)
        ans << "Torus";
      else
        ans << "Orientable genus " << genus << " surface";

      if (punctures == 1)
        ans << ", 1 puncture";
      else if (punctures > 1)
        ans << ", " << punctures << " punctures";
    }
  } else {
    if (genus == 1 && punctures == 1)
      ans << "Möbius band";
    else {
      if (genus == 1)
        ans << "Projective plane";
      else if (genus == 2)
        ans << "Klein bottle";
      else
        ans << "Non-orientable genus " << genus << " surface";

      if (punctures == 1)
        ans << ", 1 puncture";
      else if (punctures > 1)
        ans << ", " << punctures << " punctures";
    }
  }

  return ans.str();
}

class SurfaceTypeTally {
private:
  mutable std::mutex mutex_;
  std::map<SurfaceTypeKey, long long> counts_;

public:
  void merge(std::map<SurfaceTypeKey, long long> &local) {
    if (local.empty())
      return;
    std::lock_guard<std::mutex> lock(mutex_);
    for (const auto &[key, n] : local)
      counts_[key] += n;
    local.clear();
  }

  std::string summary() const {
    std::lock_guard<std::mutex> lock(mutex_);
    std::ostringstream out;
    out << "Number of surfaces found:\n";
    if (counts_.empty()) {
      out << "  (none yet)\n";
    } else {
      for (const auto &[key, n] : counts_)
        out << "  " << formatSurfaceType(key) << " = " << n << "\n";
    }
    return out.str();
  }
};

class PendingSurfaceBatch {
private:
  mutable std::mutex mutex_;
  std::vector<std::vector<int>> pending_;

public:
  void merge(std::vector<std::vector<int>> &local) {
    if (local.empty())
      return;
    std::lock_guard<std::mutex> lock(mutex_);
    for (auto &faceIndices : local)
      pending_.push_back(std::move(faceIndices));
    local.clear();
  }

  std::vector<std::vector<int>> drain() {
    std::lock_guard<std::mutex> lock(mutex_);
    std::vector<std::vector<int>> out;
    std::swap(out, pending_);
    return out;
  }
};

class LinkBoundaryTally {
private:
  mutable std::mutex mutex_;
  std::unordered_map<std::string, std::map<SurfaceTypeKey, long long>>
      linkSurfaceTypes_;

public:
  void record(const std::string &linkName, const SurfaceTypeKey &type) {
    std::lock_guard<std::mutex> lock(mutex_);
    ++linkSurfaceTypes_[linkName][type];
  }

  std::string summary() const {
    std::lock_guard<std::mutex> lock(mutex_);
    std::ostringstream out;
    out << "Links found in surface boundaries:\n";
    if (linkSurfaceTypes_.empty()) {
      out << "  (none yet)\n";
    } else {
      std::vector<std::string> names;
      names.reserve(linkSurfaceTypes_.size());
      for (const auto &[name, _] : linkSurfaceTypes_)
        names.push_back(name);
      std::sort(names.begin(), names.end());

      for (const std::string &name : names) {
        out << "  " << name << " bounds: ";
        bool first = true;
        for (const auto &[type, count] : linkSurfaceTypes_.at(name)) {
          if (!first)
            out << ", ";
          out << formatSurfaceType(type) << " (" << count << ")";
          first = false;
        }
        out << "\n";
      }
    }
    return out.str();
  }
};

inline std::string formatElapsed(std::chrono::steady_clock::duration d) {
  using namespace std::chrono;
  long long totalSeconds = duration_cast<seconds>(d).count();
  long long hours = totalSeconds / 3600;
  long long minutes = (totalSeconds % 3600) / 60;
  long long seconds_ = totalSeconds % 60;

  std::ostringstream out;
  out << std::setfill('0') << std::setw(2) << hours << ":" << std::setw(2)
      << minutes << ":" << std::setw(2) << seconds_;
  return out.str();
}

inline const char *boundaryConditionName(BoundaryCondition cond) {
  switch (cond) {
  case BoundaryCondition::all:
    return "all";
  case BoundaryCondition::closed:
    return "closed";
  case BoundaryCondition::proper:
    return "proper";
  case BoundaryCondition::connected:
    return "connected";
  }
  return "unknown";
}

template <int dim, int subdim>
class EmbeddednessPredicate : public ConditionalPredicate {
  EmbeddedSubmanifold<dim, subdim> &embedding;
  const std::vector<int> &graphToSkel;

public:
  EmbeddednessPredicate(EmbeddedSubmanifold<dim, subdim> &embedding,
                        const std::vector<int> &denseToOriginal)
      : embedding(embedding), graphToSkel(denseToOriginal) {}

  bool tryAdd(int v) override { return embedding.addFace(graphToSkel[v - 1]); }

  void undo(int v) override { embedding.removeFace(graphToSkel[v - 1]); }
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
  EmbeddingSearch(const regina::Triangulation<dim> &tri)
      : skeleton_(tri), graph_(buildGraph_(skeleton_)) {}

  size_t numEmbeddableFaces() const { return graph_.graphToSkel.size(); }

  const LinkBoundaryTally &linkTally() const { return linkTally_; }

  void processSurfaceBoundaries() {
    static_assert(dim == 4 && subdim == 2,
                  "processSurfaceBoundaries() only makes sense for "
                  "surfaces (dim == 4, subdim == 2)");

    auto batch = pendingSurfaces_.drain();
    if (batch.empty())
      return;

    EmbeddedSubmanifold<dim, subdim> embedding(skeleton_);
    processBatchRange_(embedding, batch, 0, batch.size());
  }

  // Drains whatever is left in pendingSurfaces_ once (nothing else is
  // producing into it once the DFS workers have joined) and works through
  // it with numThreads worker threads -- the same threads search() used
  // for DFS, now idle. Meant to be called right after search()'s DFS
  // workers finish, as the (potentially very large) remaining backlog of
  // queued surfaces can dwarf the amount processed incrementally by
  // search()'s periodic single-threaded boundaryProcessor thread. Prints
  // linkTally_'s incremental summary once a second while it works, since
  // this can take a long time and would otherwise look like a hang.
  void processRemainingSurfaceBoundaries(unsigned numThreads) {
    static_assert(dim == 4 && subdim == 2,
                  "processRemainingSurfaceBoundaries() only makes sense "
                  "for surfaces (dim == 4, subdim == 2)");

    auto batch = pendingSurfaces_.drain();
    if (batch.empty())
      return;

    const auto phaseStart = std::chrono::steady_clock::now();
    const size_t total = batch.size();
    const unsigned workerCount =
        static_cast<unsigned>(std::min<size_t>(numThreads, total));

    std::cerr << "\n[*] Search complete; processing " << total
              << " remaining surface boundaries with " << workerCount
              << " threads...\n";

    constexpr size_t CHUNK = 64;
    std::atomic<size_t> nextIndex{0};
    std::atomic<size_t> processedCount{0};
    std::atomic<bool> done{false};

    std::thread reporter([&]() {
      using namespace std::chrono_literals;
      size_t prevLines = 0;
      while (!done.load(std::memory_order_relaxed)) {
        std::this_thread::sleep_for(1s);

        auto elapsed = std::chrono::steady_clock::now() - phaseStart;
        long long elapsedSeconds =
            std::chrono::duration_cast<std::chrono::seconds>(elapsed).count();
        size_t processed = processedCount.load();

        std::ostringstream report;
        report << "[+] elapsed: " << formatElapsed(elapsed) << "\n";
        report << "[+] boundaries processed: " << processed << "/" << total
               << " (" << std::fixed << std::setprecision(2)
               << (total > 0 ? 100.0 * static_cast<double>(processed) /
                                   static_cast<double>(total)
                             : 0.0)
               << "%)\n";
        if (processed > 0 && processed < total && elapsedSeconds > 0) {
          long long etaSeconds = elapsedSeconds *
                                 static_cast<long long>(total - processed) /
                                 static_cast<long long>(processed);
          report << "[+] ETA: "
                 << formatElapsed(std::chrono::seconds(etaSeconds)) << "\n";
        }
        report << linkTally_.summary();
        std::string text = report.str();

        if (prevLines > 0)
          std::cerr << "\x1b[" << prevLines << "F\x1b[0J";
        std::cerr << text;
        prevLines = std::count(text.begin(), text.end(), '\n');
      }
    });

    auto worker = [&]() {
      EmbeddedSubmanifold<dim, subdim> embedding(skeleton_);
      while (true) {
        size_t begin = nextIndex.fetch_add(CHUNK, std::memory_order_relaxed);
        if (begin >= total)
          break;
        size_t end = std::min(begin + CHUNK, total);
        processBatchRange_(embedding, batch, begin, end);
        processedCount.fetch_add(end - begin, std::memory_order_relaxed);
      }
    };

    std::vector<std::thread> threads;
    threads.reserve(workerCount);
    for (unsigned t = 0; t < workerCount; ++t)
      threads.emplace_back(worker);
    for (auto &th : threads)
      th.join();

    done.store(true, std::memory_order_relaxed);
    reporter.join();

    std::cerr << "\n[+] Finished processing remaining surface "
                 "boundaries in "
              << formatElapsed(std::chrono::steady_clock::now() - phaseStart)
              << "\n";
  }

  void search(const unsigned numThreads,
              BoundaryCondition cond = BoundaryCondition::all) {
    const auto searchStart = std::chrono::steady_clock::now();
    std::atomic<int> nextRoot{1}; // shared dynamic work queue over s = 1..n
    std::vector<long long> perThreadCount(numThreads, 0);

    std::atomic<long long> globalSubgraphCount{0};
    std::atomic<int> rootsCompleted{0};
    std::atomic<bool> workersFinished{false};
    SurfaceTypeTally surfaceTally;

    // batch size for progress report
    constexpr long long FLUSH_EVERY = 1000000;

    auto worker = [&](unsigned tid) {
      ConnectedInducedSubgraphEnumerator localEnumerator(graph_.adjList.first,
                                                         graph_.adjList.second);
      EmbeddedSubmanifold<dim, subdim> embedding(skeleton_);
      EmbeddednessPredicate<dim, subdim> predicate(embedding,
                                                   graph_.graphToSkel);
      // std::ofstream out("subgraphs_thread_" + std::to_string(tid) +
      //                   ".txt");
      long long localCount = 0;
      long long sinceFlush = 0;
      std::map<SurfaceTypeKey, long long> localTypeCounts;
      std::vector<std::vector<int>> localPending;

      while (true) {
        int s = nextRoot.fetch_add(1);
        if (s > graph_.adjList.first)
          break;
        localEnumerator.enumerateFromRootFiltered(
            s,
            [&](const std::vector<int> &U) {
              if (embedding.satisfies(cond)) {
                ++localCount;
                if constexpr (subdim == 2) {
                  ++localTypeCounts[surfaceTypeKey(embedding.triangulation())];
                }
                if constexpr (dim == 4 && subdim == 2) {
                  if (cond == BoundaryCondition::proper ||
                      cond == BoundaryCondition::connected) {
                    std::vector<int> faceIndices;
                    faceIndices.reserve(U.size());
                    for (int v : U)
                      faceIndices.push_back(graph_.graphToSkel[v - 1]);
                    localPending.push_back(std::move(faceIndices));
                  }
                }
                if (++sinceFlush >= FLUSH_EVERY) {
                  globalSubgraphCount.fetch_add(sinceFlush,
                                                std::memory_order_relaxed);
                  sinceFlush = 0;
                  if constexpr (subdim == 2) {
                    surfaceTally.merge(localTypeCounts);
                  }
                  if constexpr (dim == 4 && subdim == 2) {
                    pendingSurfaces_.merge(localPending);
                  }
                }
              }
              // for (size_t i = 0; i < U.size(); ++i)
              //     out << U[i] << (i + 1 < U.size() ? ' ' :
              //     '\n');
            },
            predicate);
        rootsCompleted.fetch_add(1, std::memory_order_relaxed);
      }
      // flush this thread's remainder so the global count ends up
      // exact
      if (sinceFlush > 0)
        globalSubgraphCount.fetch_add(sinceFlush, std::memory_order_relaxed);
      if constexpr (subdim == 2) {
        surfaceTally.merge(localTypeCounts);
      }
      if constexpr (dim == 4 && subdim == 2) {
        pendingSurfaces_.merge(localPending);
      }
      perThreadCount[tid] = localCount;
    };

    std::thread reporter([&]() {
      using namespace std::chrono_literals;
      size_t prevLines = 0;
      while (!workersFinished.load(std::memory_order_relaxed)) {
        std::this_thread::sleep_for(1s);

        std::ostringstream report;
        report << "[+] elapsed: "
               << formatElapsed(std::chrono::steady_clock::now() - searchStart)
               << "\n";
        report << "[+] roots completed: " << rootsCompleted.load() << "/"
               << graph_.adjList.first
               << "  | million embedded submanifolds found so far: "
               << (globalSubgraphCount.load() / FLUSH_EVERY) << "\n";
        if constexpr (subdim == 2) {
          report << surfaceTally.summary();
        }
        if constexpr (dim == 4 && subdim == 2) {
          if (cond == BoundaryCondition::proper ||
              cond == BoundaryCondition::connected)
            report << linkTally_.summary();
        }
        std::string text = report.str();

        if (prevLines > 0)
          std::cerr << "\x1b[" << prevLines << "F\x1b[0J";
        std::cerr << text;
        prevLines = std::count(text.begin(), text.end(), '\n');
      }
    });

    std::thread boundaryProcessor;
    if constexpr (dim == 4 && subdim == 2) {
      if (cond == BoundaryCondition::proper ||
          cond == BoundaryCondition::connected) {
        boundaryProcessor = std::thread([&]() {
          using namespace std::chrono_literals;
          auto lastProcessed = std::chrono::steady_clock::now();
          while (!workersFinished.load(std::memory_order_relaxed)) {
            std::this_thread::sleep_for(100ms);
            if (std::chrono::steady_clock::now() - lastProcessed < 10s)
              continue;
            processSurfaceBoundaries();
            lastProcessed = std::chrono::steady_clock::now();
          }
        });
      }
    }

    std::vector<std::thread> threads;
    threads.reserve(numThreads);
    for (unsigned t = 0; t < numThreads; ++t)
      threads.emplace_back(worker, t);
    for (auto &th : threads)
      th.join();

    workersFinished.store(true, std::memory_order_relaxed);
    reporter.join();
    if constexpr (dim == 4 && subdim == 2) {
      if (boundaryProcessor.joinable()) {
        boundaryProcessor.join();
        processRemainingSurfaceBoundaries(numThreads);
      }
    }

    long long total = 0;
    for (long long c : perThreadCount)
      total += c;
    std::cerr << "\n\nTotal elapsed: "
              << formatElapsed(std::chrono::steady_clock::now() - searchStart)
              << "\n";
    std::cerr << "Total embedded submanifolds (" << boundaryConditionName(cond)
              << "): " << total << " (across " << numThreads << " threads)\n";
    if constexpr (subdim == 2) {
      std::cerr << surfaceTally.summary() << "\n";
    }
    if constexpr (dim == 4 && subdim == 2) {
      if (cond == BoundaryCondition::proper ||
          cond == BoundaryCondition::connected)
        std::cerr << linkTally_.summary() << "\n";
    }

    assert(globalSubgraphCount.load() == total);
  }

private:
  // Replays batch[begin, end) into embedding (via addFace(), in the exact
  // order the DFS originally discovered each entry -- deterministic,
  // since it's the same skeleton_/graph_), extracts its boundary Link(s),
  // and records what surface type each recognized link bounds into
  // linkTally_. embedding is reused across the whole range so its
  // bdryComponents_ cache (see EmbeddedSubmanifold's constructor) is only
  // built once, not once per queued surface.
  void processBatchRange_(EmbeddedSubmanifold<dim, subdim> &embedding,
                          const std::vector<std::vector<int>> &batch,
                          size_t begin, size_t end) {
    for (size_t i = begin; i < end; ++i) {
      const auto &faceIndices = batch[i];
      for (int idx : faceIndices)
        embedding.addFace(idx);

      SurfaceTypeKey type = surfaceTypeKey(embedding.triangulation());
      for (auto &[component, link] : embedding.boundaryLinks())
        linkTally_.record(link.identify(), type);

      // Reverse order, mirroring how the DFS itself would back out --
      // resets embedding to empty for the next item in the batch.
      for (auto it = faceIndices.rbegin(); it != faceIndices.rend(); ++it)
        embedding.removeFace(*it);
    }
  }

  static Graph buildGraph_(const Skeleton<dim, subdim> &skeleton) {
    const auto &nodes = skeleton.getNodes();

    std::vector<int> graphToSkel;
    std::vector<int> skelToGraph(nodes.size(), -1);
    for (size_t i = 0; i < nodes.size(); ++i) {
      if (hasIrreparableSelfGluing<dim, subdim>(nodes[i].gluings))
        continue;
      skelToGraph[i] = static_cast<int>(graphToSkel.size());
      graphToSkel.push_back(static_cast<int>(i));
    }

    int n = static_cast<int>(graphToSkel.size());
    // 1-indexed, as the enumerator expects
    std::vector<std::vector<int>> adj(n + 1);
    for (int graphIdx = 0; graphIdx < n; ++graphIdx) {
      int i = graphToSkel[graphIdx];
      std::set<int> neighbors; // dedupes parallel gluings
      for (const auto &g : nodes[i].gluings) {
        assert(g.srcIndex == static_cast<size_t>(i));
        assert(g.dstIndex < nodes.size());
        int u = static_cast<int>(g.dstIndex);
        if (u == i)
          continue; // drop self-gluings
        int graphUdx = skelToGraph[u];
        if (graphUdx != -1) // neighbor also survived exclusion
          neighbors.insert(graphUdx);
      }
      for (int denseU : neighbors)
        adj[graphIdx + 1].push_back(denseU + 1);
    }
    return {AdjacencyList{n, std::move(adj)}, std::move(graphToSkel)};
  }
};

#endif // EMBEDDINGSEARCH_H
