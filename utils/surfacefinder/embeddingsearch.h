#ifndef EMBEDDINGSEARCH_H

#define EMBEDDINGSEARCH_H

#include <atomic>
#include <thread>

#include <triangulation/dim4.h>
#include <unordered_set>

#include "enumerate_cis.h"
#include "skeleton.h"

using AdjacencyList = std::pair<int, std::vector<std::vector<int>>>;

template <int dim, int subdim>
class EmbeddedSubmanifold {
  private:
    using Face = const typename Skeleton<dim, subdim>::Face;

    const Skeleton<dim, subdim> &skeleton_;
    regina::Triangulation<subdim> subtri_;
    std::vector<std::pair<regina::Simplex<subdim> *, std::vector<int>>> faces_;
    std::vector<bool> triVertsUsed_;

  public:
    EmbeddedSubmanifold(const Skeleton<dim, subdim> &skeleton)
        : skeleton_(skeleton), faces_(skeleton.numFaces()),
          triVertsUsed_(skeleton.numVertices()) {}

    bool addFace(int f) {
        const auto &node = skeleton_.getNodes()[f];

        // It's worth it to first check if adding this face would result in an
        // embedding before actually modifying our state

        // Facets of this simplex which are glued to this simplex
        std::vector<int> selfGluedFacets;
        // Facets of this simplex which are glued to a different simplex
        std::vector<int> adjGluedFacets;
        // Vertices appearing in facets of this simplex glued to different
        // simplices
        std::unordered_set<int> usedVerts;
        // Gluings to simplices in subtri_ (or the current simplex)
        std::vector<typename Skeleton<dim, subdim>::Gluing> inducedGluings;

        // First check facet conditions
        for (const auto &g : node.gluings) {
            if (g.srcIndex == g.dstIndex) {
                inducedGluings.push_back(g);
                selfGluedFacets.push_back(g.srcFacet);
                continue; // Check the self-gluings later
            }

            const auto *dst = faces_[g.dstIndex].first;

            if (dst == nullptr)
                continue;

            // If something else is glued to dstFacet, it's over
            size_t dstFacet = g.gluing[g.srcFacet];
            if (dst->adjacentSimplex(dstFacet) != nullptr)
                return false;

            const auto *srcFacet =
                node.face->template face<subdim - 1>(g.srcFacet);
            adjGluedFacets.push_back(g.srcFacet);
            for (int k = 0; k < subdim; ++k) {
                usedVerts.insert(srcFacet->vertex(k)->index());
            }

            inducedGluings.push_back(g);
        }

        for (int fi : selfGluedFacets) {
            for (int fj : adjGluedFacets) {
                // A facet cannot be glued to both this simplex and another
                // simplex
                if (fi == fj)
                    return false;
            }
        }

        // Next check vertex conditions. It suffices to check that "new"
        // vertices (i.e. vertices not already used by faces in the embedding)
        // admit regular neighborhoods.
        std::vector<int> newVerts;
        for (int i = 0; i <= subdim; ++i) {
            int v = node.face->vertex(i)->index();
            if (!usedVerts.contains(v)) {
                newVerts.push_back(v);
                if (triVertsUsed_[v])
                    return false;
            }
        }

        // Now we know it's safe to add this face and modify our state
        for (int v : newVerts) {
            faces_[f].second.push_back(v);
            triVertsUsed_[v] = true;
        }

        // TODO: I *think* it's more efficient to build the triangulation as
        // we go. Technically, we could modify this class so that this isn't
        // the case, but hopefully Regina's triangulation is efficient
        // enough (i.e. doesn't recomputeSkeleton() on every addFace() call)
        // for this to work.
        auto *src = subtri_.newSimplex();
        faces_[f].first = src;

        // Perform gluings
        for (const auto &g : inducedGluings) {
            if (src->adjacentSimplex(g.srcFacet) != nullptr) {
                // This should never happen unless we're dealing with a
                // self-gluing
                assert(g.srcIndex == g.dstIndex);
                continue;
            }
            src->join(g.srcFacet, faces_[g.dstIndex].first, g.gluing);
        }

        return true;
    }

    void removeFace(int f) {
        auto *face = faces_[f].first;
        // This should go through without issue (assuming the argument is
        // valid)
        if (face == nullptr)
            throw regina::InvalidArgument(
                "EmbeddedSubmanifold::removeFace(): Was asked to remove a face "
                "which is not in the embedded submanifold.");

        // Update triVertsUsed_
        for (int v : faces_[f].second)
            triVertsUsed_[v] = false;

        // Update faces_
        faces_[f].first = nullptr;
        faces_[f].second = {};

        // Update subtri_
        subtri_.removeSimplex(face);
    }
};

template <int dim, int subdim>
class EmbeddednessPredicate : public ConditionalPredicate {
    EmbeddedSubmanifold<dim, subdim> &embedding;

  public:
    explicit EmbeddednessPredicate(EmbeddedSubmanifold<dim, subdim> &embedding)
        : embedding(embedding) {}

    bool tryAdd(int v) override { return embedding.addFace(v - 1); }

    void undo(int v) override {
        // Reverse exactly what tryAdd(v) committed.
        embedding.removeFace(v - 1);
    }
};

template <int dim, int subdim>
class EmbeddingSearch {
  private:
    const Skeleton<dim, subdim> skeleton_;
    const AdjacencyList adjList_;
    std::vector<EmbeddedSubmanifold<dim, subdim>> valid_embeddings_;

  public:
    EmbeddingSearch(const regina::Triangulation<dim> &tri)
        : skeleton_(tri), adjList_(toAdjacencyList_(skeleton_)) {}

    void search(const unsigned numThreads) {
        std::atomic<int> nextRoot{1}; // shared dynamic work queue over s = 1..n
        std::vector<long long> perThreadCount(numThreads, 0);

        std::atomic<long long> globalSubgraphCount{0};
        std::atomic<int> rootsCompleted{0};
        std::atomic<bool> workersFinished{false};

        // batch size for progress report
        constexpr long long FLUSH_EVERY = 1000000;

        auto worker = [&](unsigned tid) {
            ConnectedInducedSubgraphEnumerator localEnumerator(adjList_.first,
                                                               adjList_.second);
            EmbeddedSubmanifold<dim, subdim> embedding(skeleton_);
            EmbeddednessPredicate<dim, subdim> predicate(embedding);
            // std::ofstream out("subgraphs_thread_" + std::to_string(tid) +
            //                   ".txt");
            long long localCount = 0;
            long long sinceFlush = 0;

            while (true) {
                int s = nextRoot.fetch_add(1);
                if (s > adjList_.first)
                    break;
                localEnumerator.enumerateFromRootFiltered(
                    s,
                    [&](const std::vector<int> &U) {
                        ++localCount;
                        if (++sinceFlush >= FLUSH_EVERY) {
                            globalSubgraphCount.fetch_add(
                                sinceFlush, std::memory_order_relaxed);
                            sinceFlush = 0;
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
                globalSubgraphCount.fetch_add(sinceFlush,
                                              std::memory_order_relaxed);
            perThreadCount[tid] = localCount;
        };

        // Prints total progress every second until all workers are done.
        std::thread reporter([&]() {
            using namespace std::chrono_literals;
            while (!workersFinished.load(std::memory_order_relaxed)) {
                std::this_thread::sleep_for(1s);
                std::cerr << "[+] roots completed: " << rootsCompleted.load()
                          << "/" << adjList_.first
                          << "  | million embedded submanifolds found so far: "
                          << (globalSubgraphCount.load() / FLUSH_EVERY) << "\n";
            }
        });

        std::vector<std::thread> threads;
        threads.reserve(numThreads);
        for (unsigned t = 0; t < numThreads; ++t)
            threads.emplace_back(worker, t);
        for (auto &th : threads)
            th.join();

        workersFinished.store(true, std::memory_order_relaxed);
        reporter.join();

        long long total = 0;
        for (long long c : perThreadCount)
            total += c;
        std::cerr << "Total embedded submanifolds: " << total
                  << " (across " << numThreads << " threads)\n";

        assert(globalSubgraphCount.load() == total);
    }

  private:
    AdjacencyList toAdjacencyList_(const Skeleton<dim, subdim> &skeleton) {
        const auto &nodes = skeleton.getNodes();
        int n = static_cast<int>(nodes.size());

        std::vector<std::vector<int>> adj(
            n + 1); // 1-indexed, as the enumerator expects
        for (int i = 0; i < n; ++i) {
            std::set<int> neighbors; // dedupes parallel gluings
            for (const auto &g : nodes[i].gluings) {
                // sanity: matches its own position
                assert(g.srcIndex == static_cast<size_t>(i));
                // sanity: in range
                assert(g.dstIndex < nodes.size());
                int u = static_cast<int>(g.dstIndex);
                if (u != i)
                    neighbors.insert(u); // drop self-gluings
            }
            for (int u : neighbors)
                adj[i + 1].push_back(u + 1);
        }
        return {n, std::move(adj)};
    }
};

#endif // EMBEDDINGSEARCH_H
