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
#include <triangulation/dim4.h>
#include <unordered_map>
#include <unordered_set>

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

// Thread-safe accumulator of discovered face-index sets (each one uniquely
// identifying a found embedding, exactly the `U` the DFS reports it under),
// staged for later, off-hot-path processing by
// EmbeddingSearch::processSurfaceBoundaries(). Mirrors SurfaceTypeTally's
// merge-under-lock/periodic-flush shape so the same FLUSH_EVERY cadence in
// search()'s worker loop can drive both.
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

// Thread-safe record of every distinct boundary link found so far (keyed by
// Link::identify()'s result -- a census name, "Unknot", or a bare isoSig
// fallback) together with how many times each surface type was found
// bounding it.
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
class EmbeddedSubmanifold {
  private:
    using Face = const typename Skeleton<dim, subdim>::Face;

    const Skeleton<dim, subdim> &skeleton_;
    regina::Triangulation<subdim> subtri_;
    std::vector<std::pair<regina::Simplex<subdim> *, std::vector<int>>> faces_;
    std::vector<bool> triVertsUsed_;

    // One triangulation per ambient boundary component, built once (via
    // BoundaryComponent<dim>::build()) and reused for the lifetime of this
    // instance. Only ever populated for dim == 4, subdim == 2 -- the only
    // case boundaryLinks() below makes sense for -- so it stays an empty,
    // harmless vector for every other instantiation of this class template.
    std::vector<regina::Triangulation<dim - 1>> bdryComponents_;

  public:
    EmbeddedSubmanifold(const Skeleton<dim, subdim> &skeleton)
        : skeleton_(skeleton), faces_(skeleton.numFaces()),
          triVertsUsed_(skeleton.numVertices()) {
        if constexpr (dim == 4 && subdim == 2) {
            const auto &tri = skeleton_.triangulation();
            bdryComponents_.reserve(tri.countBoundaryComponents());
            for (size_t c = 0; c < tri.countBoundaryComponents(); ++c)
                bdryComponents_.push_back(tri.boundaryComponent(c)->build());
        }
    }

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

    const regina::Triangulation<subdim> &triangulation() const {
        return subtri_;
    }

    // Unfortunately, Regina::Triangulation::isClosed() is only valid for
    // triangulations of dimensions 2, 3, and 4
    bool isClosed() const { return !subtri_.hasBoundaryFacets(); }

    bool isProper() const {
        for (size_t f = 0; f < faces_.size(); ++f) {
            const auto *simplex = faces_[f].first;
            if (simplex == nullptr)
                continue; // face f is not part of this embedding

            for (int i = 0; i <= subdim; ++i) {
                if (simplex->adjacentSimplex(i) != nullptr)
                    continue; // internal facet of subtri_

                const auto *ambientFacet =
                    skeleton_.getNodes()[f].face->template face<subdim - 1>(i);
                if (!ambientFacet->isBoundary())
                    return false;
            }
        }
        return true;
    }

    bool boundaryComponentsMapInjectively() const {
        std::unordered_map<size_t, const regina::BoundaryComponent<dim> *>
            subToAmbient;
        std::unordered_set<const regina::BoundaryComponent<dim> *> usedAmbient;

        for (size_t f = 0; f < faces_.size(); ++f) {
            const auto *simplex = faces_[f].first;
            if (simplex == nullptr)
                continue;

            for (int i = 0; i <= subdim; ++i) {
                if (simplex->adjacentSimplex(i) != nullptr)
                    continue;

                const auto *ambientFacet =
                    skeleton_.getNodes()[f].face->template face<subdim - 1>(i);
                const auto *ambientBC = ambientFacet->boundaryComponent();
                if (ambientBC == nullptr)
                    return false;

                const auto *subFacet = simplex->template face<subdim - 1>(i);
                const auto *subBC = subFacet->boundaryComponent();
                assert(subBC != nullptr);

                size_t subIdx = subBC->index();
                auto it = subToAmbient.find(subIdx);
                if (it == subToAmbient.end()) {
                    if (usedAmbient.contains(ambientBC))
                        return false;
                    subToAmbient.emplace(subIdx, ambientBC);
                    usedAmbient.insert(ambientBC);
                } else {
                    assert(it->second == ambientBC);
                    if (it->second != ambientBC)
                        return false;
                }
            }
        }
        return true;
    }

    // Returns one Link per ambient boundary component this embedding's
    // boundary actually touches, paired with that component's index. Only
    // meaningful when dim == 4 and subdim == 2 (a surface's boundary is a
    // 1-manifold living in a 3-manifold boundary component) -- guarded by a
    // static_assert rather than `if constexpr` since this is a template
    // member function, so it's simply never instantiated for any other
    // (dim, subdim) as long as callers only invoke it under their own
    // `if constexpr (dim == 4 && subdim == 2)`.
    //
    // Deliberately reuses the same un-glued-facet scan as isProper()/
    // boundaryComponentsMapInjectively() above (rather than
    // subtri_.boundaryComponents(), which KnottedSurface::boundary() used to
    // rely on) so this doesn't force a 2-manifold skeleton recompute of its
    // own -- surfaceTypeKey() already pays for that separately when it's
    // needed.
    std::vector<std::pair<size_t, Link>> boundaryLinks() const {
        static_assert(dim == 4 && subdim == 2,
                      "boundaryLinks() only makes sense for surfaces (dim == "
                      "4, subdim == 2)");

        std::vector<std::unordered_set<const regina::Edge<3> *>>
            edgesByComponent(bdryComponents_.size());

        for (size_t f = 0; f < faces_.size(); ++f) {
            const auto *simplex = faces_[f].first;
            if (simplex == nullptr)
                continue;

            for (int i = 0; i <= subdim; ++i) {
                if (simplex->adjacentSimplex(i) != nullptr)
                    continue; // internal facet of subtri_

                const auto *ambientFacet =
                    skeleton_.getNodes()[f].face->template face<subdim - 1>(i);
                const auto *ambientBC = ambientFacet->boundaryComponent();
                if (ambientBC == nullptr)
                    continue; // shouldn't happen when cond is proper/connected

                size_t c = ambientBC->index();

                // No O(1) "my index within my boundary component" accessor
                // exists in Regina, but BoundaryComponent<dim>::build()
                // preserves index-for-index correspondence with the
                // component's own edge list one dimension down -- same idiom
                // KnottedSurface::boundary() used.
                for (int k = 0; k < ambientBC->countEdges(); ++k) {
                    if (ambientBC->edge(k) == ambientFacet) {
                        edgesByComponent[c].insert(
                            bdryComponents_[c].edge(k));
                        break;
                    }
                }
            }
        }

        std::vector<std::pair<size_t, Link>> result;
        for (size_t c = 0; c < edgesByComponent.size(); ++c) {
            if (edgesByComponent[c].empty())
                continue;

            std::vector<const regina::Edge<3> *> edgeList(
                edgesByComponent[c].begin(), edgesByComponent[c].end());
            result.emplace_back(c, Link(bdryComponents_[c], edgeList));
        }
        return result;
    }

    bool satisfies(BoundaryCondition cond) const {
        switch (cond) {
        case BoundaryCondition::all:
            return true;
        case BoundaryCondition::closed:
            return isClosed();
        case BoundaryCondition::proper:
            return isProper();
        case BoundaryCondition::connected:
            return boundaryComponentsMapInjectively();
        }

        throw regina::InvalidArgument(
            "EmbeddedSubmanifold::satisfies(): Invalid BoundaryCondition");
    }
};

template <int dim, int subdim>
class EmbeddednessPredicate : public ConditionalPredicate {
    EmbeddedSubmanifold<dim, subdim> &embedding;

  public:
    explicit EmbeddednessPredicate(EmbeddedSubmanifold<dim, subdim> &embedding)
        : embedding(embedding) {}

    bool tryAdd(int v) override { return embedding.addFace(v - 1); }

    void undo(int v) override { embedding.removeFace(v - 1); }
};

template <int dim, int subdim>
class EmbeddingSearch {
  private:
    const Skeleton<dim, subdim> skeleton_;
    const AdjacencyList adjList_;
    std::vector<EmbeddedSubmanifold<dim, subdim>> valid_embeddings_;

    // Only ever populated/consumed when dim == 4, subdim == 2, and cond is
    // proper or connected -- see search()'s dedicated background thread and
    // processSurfaceBoundaries() below. Harmless empty members otherwise.
    PendingSurfaceBatch pendingSurfaces_;
    LinkBoundaryTally linkTally_;

  public:
    EmbeddingSearch(const regina::Triangulation<dim> &tri)
        : skeleton_(tri), adjList_(toAdjacencyList_(skeleton_)) {}

    const LinkBoundaryTally &linkTally() const { return linkTally_; }

    // Drains pendingSurfaces_ and, for each queued face-index set, replays
    // it (via addFace(), in the exact order the DFS originally discovered
    // it -- deterministic, since it's the same skeleton_/adjList_) into a
    // single reused EmbeddedSubmanifold, extracts its boundary Link(s), and
    // records what surface type each recognized link bounds into
    // linkTally_. This is the expensive step (pinch + simplify + census
    // lookup per surface) deliberately kept off the DFS hot path.
    void processSurfaceBoundaries() {
        static_assert(dim == 4 && subdim == 2,
                      "processSurfaceBoundaries() only makes sense for "
                      "surfaces (dim == 4, subdim == 2)");

        auto batch = pendingSurfaces_.drain();
        if (batch.empty())
            return;

        // Reused across every item in the batch so its bdryComponents_
        // cache (see EmbeddedSubmanifold's constructor) is only built once
        // per call, not once per queued surface.
        EmbeddedSubmanifold<dim, subdim> embedding(skeleton_);

        for (const auto &faceIndices : batch) {
            for (int idx : faceIndices)
                embedding.addFace(idx);

            SurfaceTypeKey type = surfaceTypeKey(embedding.triangulation());
            for (auto &[component, link] : embedding.boundaryLinks())
                linkTally_.record(link.identify(), type);

            // Reverse order, mirroring how the DFS itself would back out --
            // resets embedding to empty for the next item in the batch.
            for (auto it = faceIndices.rbegin(); it != faceIndices.rend();
                 ++it)
                embedding.removeFace(*it);
        }
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
            ConnectedInducedSubgraphEnumerator localEnumerator(adjList_.first,
                                                               adjList_.second);
            EmbeddedSubmanifold<dim, subdim> embedding(skeleton_);
            EmbeddednessPredicate<dim, subdim> predicate(embedding);
            // std::ofstream out("subgraphs_thread_" + std::to_string(tid) +
            //                   ".txt");
            long long localCount = 0;
            long long sinceFlush = 0;
            std::map<SurfaceTypeKey, long long> localTypeCounts;
            std::vector<std::vector<int>> localPending;

            while (true) {
                int s = nextRoot.fetch_add(1);
                if (s > adjList_.first)
                    break;
                localEnumerator.enumerateFromRootFiltered(
                    s,
                    [&](const std::vector<int> &U) {
                        if (embedding.satisfies(cond)) {
                            ++localCount;
                            if constexpr (subdim == 2) {
                                ++localTypeCounts[surfaceTypeKey(
                                    embedding.triangulation())];
                            }
                            if constexpr (dim == 4 && subdim == 2) {
                                if (cond == BoundaryCondition::proper ||
                                    cond == BoundaryCondition::connected) {
                                    std::vector<int> faceIndices;
                                    faceIndices.reserve(U.size());
                                    for (int v : U)
                                        faceIndices.push_back(v - 1);
                                    localPending.push_back(
                                        std::move(faceIndices));
                                }
                            }
                            if (++sinceFlush >= FLUSH_EVERY) {
                                globalSubgraphCount.fetch_add(
                                    sinceFlush, std::memory_order_relaxed);
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
                globalSubgraphCount.fetch_add(sinceFlush,
                                              std::memory_order_relaxed);
            if constexpr (subdim == 2) {
                surfaceTally.merge(localTypeCounts);
            }
            if constexpr (dim == 4 && subdim == 2) {
                pendingSurfaces_.merge(localPending);
            }
            perThreadCount[tid] = localCount;
        };

        // Prints total progress every second until all workers are done,
        // redrawing in place (rather than scrolling) so the tally doesn't
        // eat up the terminal.
        std::thread reporter([&]() {
            using namespace std::chrono_literals;
            size_t prevLines = 0;
            while (!workersFinished.load(std::memory_order_relaxed)) {
                std::this_thread::sleep_for(1s);

                std::ostringstream report;
                report << "[+] elapsed: "
                       << formatElapsed(std::chrono::steady_clock::now() -
                                        searchStart)
                       << "\n";
                report << "[+] roots completed: " << rootsCompleted.load()
                       << "/" << adjList_.first
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

                // Move the cursor up to the start of the previous report
                // and clear everything from there down, so a shorter report
                // (e.g. fewer distinct surface types) doesn't leave stale
                // lines behind.
                if (prevLines > 0)
                    std::cerr << "\x1b[" << prevLines << "F\x1b[0J";
                std::cerr << text;
                prevLines = std::count(text.begin(), text.end(), '\n');
            }
        });

        // Drains pendingSurfaces_ roughly every 10 seconds on its own
        // dedicated thread -- separate from `reporter` above, since
        // processSurfaceBoundaries() does genuinely heavy work (pinch +
        // simplify + census lookup per surface) that must not delay the
        // lightweight 1s progress display. Only started when there's
        // actually something for it to do.
        std::thread boundaryProcessor;
        if constexpr (dim == 4 && subdim == 2) {
            if (cond == BoundaryCondition::proper ||
                cond == BoundaryCondition::connected) {
                boundaryProcessor = std::thread([&]() {
                    using namespace std::chrono_literals;
                    // Polled in short (100ms) increments, rather than one
                    // long sleep_for(10s), so that when the DFS finishes
                    // (workersFinished flips true) this thread notices
                    // promptly instead of leaving search() blocked on
                    // join() for up to the remainder of a 10s sleep --
                    // otherwise even a trivially small search would always
                    // take >= 10s wall-clock once this thread is started.
                    auto lastProcessed = std::chrono::steady_clock::now();
                    while (!workersFinished.load(std::memory_order_relaxed)) {
                        std::this_thread::sleep_for(100ms);
                        if (std::chrono::steady_clock::now() - lastProcessed <
                            10s)
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
                // Flush whatever landed in the final partial window -- the
                // thread above only processes on its own ~10s cadence, so
                // without this, surfaces found in the last (<10s) stretch
                // would never get processed.
                processSurfaceBoundaries();
            }
        }

        long long total = 0;
        for (long long c : perThreadCount)
            total += c;
        std::cerr << "\n\nTotal elapsed: "
                  << formatElapsed(std::chrono::steady_clock::now() -
                                   searchStart)
                  << "\n";
        std::cerr << "Total embedded submanifolds ("
                  << boundaryConditionName(cond) << "): " << total
                  << " (across " << numThreads << " threads)\n";
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
