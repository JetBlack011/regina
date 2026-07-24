//
//  embeddingsearch.cpp
//

#include "embeddingsearch.h"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <thread>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

// satisfying-cond finds are rare relative to raw finds, so flush every one --
// batching these hid live progress behind a threshold that was rarely
// reached, making the rolling report look stuck at 0
#define FLUSH_EVERY_BDRY 1
// the found-count fires on every callback invocation (not just those
// satisfying cond), so it's flushed to the atomic less often to keep
// contention down
#define FLUSH_EVERY_FOUND 1000000

std::string formatElapsed(std::chrono::steady_clock::duration d) {
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

const char *boundaryConditionName(BoundaryCondition cond) {
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

// Labels a count that's been divided by `divisor` (a progress-report flush
// batch size) with the scale actually applied, e.g. divisor=100000 gives
// " (x10^5)". Returns "" for divisor <= 1, where the reported count is exact.
static std::string scaleLabel(long long divisor) {
    if (divisor <= 1)
        return "";
    int exponent =
        static_cast<int>(std::llround(std::log10(static_cast<double>(divisor))));
    std::ostringstream out;
    out << " (x10^" << exponent << ")";
    return out.str();
}

template <int dim, int subdim>
EmbeddingSearch<dim, subdim>::EmbeddednessPredicate::EmbeddednessPredicate(
    EmbeddedSubmanifold<dim, subdim> &embedding,
    const std::vector<std::vector<int>> &graphToSkel)
    : embedding_(embedding), graphToSkel_(graphToSkel) {}

template <int dim, int subdim>
bool EmbeddingSearch<dim, subdim>::EmbeddednessPredicate::tryAdd(int v) {
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

template <int dim, int subdim>
void EmbeddingSearch<dim, subdim>::EmbeddednessPredicate::undo(int v) {
    const auto &faces = graphToSkel_[v - 1];
    for (auto it = faces.rbegin(); it != faces.rend(); ++it)
        embedding_.removeFace(*it);
}

template <int dim, int subdim>
EmbeddingSearch<dim, subdim>::EmbeddingSearch(
    const regina::Triangulation<dim> &tri)
    : skeleton_(tri), graph_(buildGraph_(skeleton_)) {}

template <int dim, int subdim>
EmbeddingSearch<dim, subdim>::EmbeddingSearch(
    const regina::Triangulation<dim> &tri, const std::vector<int> &seedFaces)
    : skeleton_(tri), graph_(buildSeededGraph_(skeleton_, seedFaces)),
      isSeeded_(true) {
    // Eagerly validate that seedFaces is jointly addable in *some* order
    // (not just individually embeddable -- buildSeededGraph_ already
    // checked that) from an empty submanifold. Throws
    // regina::InvalidArgument on failure. Discarded immediately -- this is
    // purely a fail-fast check so a bad seed is caught here, at
    // construction, rather than later inside a worker thread during
    // search().
    EmbeddedSubmanifold<dim, subdim>(skeleton_, seedFaces);
}

template <int dim, int subdim>
void EmbeddingSearch<dim, subdim>::search(const unsigned numThreads,
                                          BoundaryCondition cond) {
    const auto searchStart = std::chrono::steady_clock::now();

    // Shared dynamic work queue over roots. Unseeded: every graph vertex,
    // unconditionally (matches the old s = 1..n sweep exactly). Seeded:
    // every sibling of the seed surviving the predicate, via one prototype
    // embedding/enumerator built once here on the calling thread -- which
    // also lets us handle the "seed alone, no additional faces" result
    // here, since no root's subtree can ever produce it (every root
    // requires descending into at least one more face).
    std::vector<int> roots;
    long long seedFoundCount = 0;
    long long seedSubgraphCount = 0;
    long long seedMaxFaces = 0;
    long long seedFaceSum = 0;
    if (isSeeded_) {
        EmbeddedSubmanifold<dim, subdim> protoEmbedding(skeleton_);
        EmbeddednessPredicate protoPredicate(protoEmbedding,
                                             graph_.graphToSkel);
        ConnectedInducedSubgraphEnumerator protoEnumerator(
            graph_.adjList.first, graph_.adjList.second, true, protoPredicate);
        roots = protoEnumerator.getRoots();

        seedFoundCount = 1;
        if (protoEmbedding.isEmbedded() && protoEmbedding.satisfies(cond)) {
            seedSubgraphCount = 1;
            seedMaxFaces = static_cast<long long>(
                protoEmbedding.triangulation().size());
            seedFaceSum = seedMaxFaces;
        }
    } else {
        roots.resize(graph_.adjList.first);
        std::iota(roots.begin(), roots.end(), 1);
    }
    const size_t totalRoots = roots.size();

    std::atomic<size_t> nextRootIdx{0};
    std::vector<long long> perThreadCount(numThreads, 0);

    std::atomic<long long> globalFoundCount{seedFoundCount};
    std::atomic<long long> globalSubgraphCount{seedSubgraphCount};
    std::atomic<long long> globalMaxFaces{seedMaxFaces};
    std::atomic<long long> globalFaceSum{seedFaceSum};
    std::atomic<size_t> rootsCompleted{0};
    std::atomic<bool> workersFinished{false};
    std::vector<long long> perThreadFoundCount(numThreads, 0);
    std::vector<long long> perThreadFaceSum(numThreads, 0);

    auto worker = [&](unsigned tid) {
        EmbeddedSubmanifold<dim, subdim> embedding(skeleton_);
        EmbeddednessPredicate predicate(embedding, graph_.graphToSkel);
        // ConnectedInducedSubgraphEnumerator holds a reference member, so
        // it isn't assignable -- construct it in place via optional rather
        // than choosing between two constructor calls via assignment.
        std::optional<ConnectedInducedSubgraphEnumerator> localEnumeratorOpt;
        if (isSeeded_)
            localEnumeratorOpt.emplace(graph_.adjList.first,
                                       graph_.adjList.second, true, predicate);
        else
            localEnumeratorOpt.emplace(graph_.adjList.first,
                                       graph_.adjList.second);
        auto &localEnumerator = *localEnumeratorOpt;
        long long localCount = 0;
        long long sinceFlush = 0;
        long long localFoundCount = 0;
        long long sinceFlushFound = 0;
        long long localFaceSum = 0;
        long long sinceFlushFaceSum = 0;

        while (true) {
            size_t idx = nextRootIdx.fetch_add(1);
            if (idx >= totalRoots)
                break;
            int s = roots[idx];
            localEnumerator.enumerateFromRootFiltered(
                s,
                [&](const std::vector<int> &) {
                    ++localFoundCount;
                    if (++sinceFlushFound >= FLUSH_EVERY_FOUND) {
                        globalFoundCount.fetch_add(sinceFlushFound,
                                                   std::memory_order_relaxed);
                        sinceFlushFound = 0;
                    }
                    if (embedding.isEmbedded() && embedding.satisfies(cond)) {
                        ++localCount;
                        long long faceCount = static_cast<long long>(
                            embedding.triangulation().size());
                        localFaceSum += faceCount;
                        sinceFlushFaceSum += faceCount;
                        if (++sinceFlush >= FLUSH_EVERY_BDRY) {
                            globalSubgraphCount.fetch_add(
                                sinceFlush, std::memory_order_relaxed);
                            globalFaceSum.fetch_add(
                                sinceFlushFaceSum, std::memory_order_relaxed);
                            sinceFlush = 0;
                            sinceFlushFaceSum = 0;
                        }
                        auto prevMax =
                            globalMaxFaces.load(std::memory_order_relaxed);
                        while (faceCount > prevMax &&
                               !globalMaxFaces.compare_exchange_weak(
                                   prevMax, faceCount,
                                   std::memory_order_relaxed))
                            ;
                    }
                },
                predicate);
            rootsCompleted.fetch_add(1, std::memory_order_relaxed);
        }
        // flush this thread's remainder so the global count ends up
        // exact
        if (sinceFlush > 0) {
            globalSubgraphCount.fetch_add(sinceFlush,
                                          std::memory_order_relaxed);
            globalFaceSum.fetch_add(sinceFlushFaceSum,
                                    std::memory_order_relaxed);
        }
        if (sinceFlushFound > 0)
            globalFoundCount.fetch_add(sinceFlushFound,
                                       std::memory_order_relaxed);
        perThreadCount[tid] = localCount;
        perThreadFoundCount[tid] = localFoundCount;
        perThreadFaceSum[tid] = localFaceSum;
    };

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
            report << "[+] roots completed: " << rootsCompleted.load() << "/"
                   << totalRoots
                   << "  | embedded submanifolds found so far"
                   << scaleLabel(FLUSH_EVERY_FOUND) << ": "
                   << (globalFoundCount.load() / FLUSH_EVERY_FOUND) << "\n";
            report << "[+] embedded submanifolds satisfying boundary "
                      "condition ("
                   << boundaryConditionName(cond) << ") found so far"
                   << scaleLabel(FLUSH_EVERY_BDRY) << ": "
                   << (globalSubgraphCount.load() / FLUSH_EVERY_BDRY) << "\n";
            report << "[+] largest embedded submanifold satisfying boundary "
                      "condition ("
                   << boundaryConditionName(cond)
                   << ") found so far: " << globalMaxFaces.load()
                   << " faces\n";
            {
                long long n = globalSubgraphCount.load();
                double avg = n > 0 ? static_cast<double>(globalFaceSum.load()) /
                                         static_cast<double>(n)
                                   : 0.0;
                report << "[+] average faces per embedded submanifold "
                          "satisfying boundary condition ("
                       << boundaryConditionName(cond)
                       << ") found so far: " << std::fixed
                       << std::setprecision(2) << avg << "\n";
            }
            std::string text = report.str();

            if (prevLines > 0)
                std::cerr << "\x1b[" << prevLines << "F\x1b[0J";
            std::cerr << text;
            prevLines = std::count(text.begin(), text.end(), '\n');
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

    long long total = seedSubgraphCount;
    for (long long c : perThreadCount)
        total += c;
    long long totalFound = seedFoundCount;
    for (long long c : perThreadFoundCount)
        totalFound += c;
    long long totalFaceSum = seedFaceSum;
    for (long long c : perThreadFaceSum)
        totalFaceSum += c;
    std::cerr << "\n\nTotal elapsed: "
              << formatElapsed(std::chrono::steady_clock::now() - searchStart)
              << "\n";
    std::cerr << "Total embedded submanifolds found: " << totalFound
              << " (across " << numThreads << " threads)\n";
    std::cerr << "Total embedded submanifolds (" << boundaryConditionName(cond)
              << "): " << total << " (across " << numThreads << " threads)\n";
    std::cerr << "Largest embedded submanifold (" << boundaryConditionName(cond)
              << "): " << globalMaxFaces.load() << " faces\n";
    std::cerr << "Average faces per embedded submanifold ("
              << boundaryConditionName(cond) << "): " << std::fixed
              << std::setprecision(2)
              << (total > 0 ? static_cast<double>(totalFaceSum) /
                                  static_cast<double>(total)
                            : 0.0)
              << "\n";

    assert(globalFoundCount.load() == totalFound);
    assert(globalSubgraphCount.load() == total);
    assert(globalFaceSum.load() == totalFaceSum);
}

template <int dim, int subdim>
typename EmbeddingSearch<dim, subdim>::Graph
EmbeddingSearch<dim, subdim>::buildGraph_(
    const Skeleton<dim, subdim> &skeleton) {
    const auto &nodes = skeleton.getNodes();

    std::vector<int> skelOf; // dense graph index -> skeleton index
    std::vector<int> skelToGraph(nodes.size(), -1);
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (EmbeddedSubmanifold<dim, subdim>::hasIrreparableSelfGluing(
                nodes[i].gluings)
            // hasUnexplainedSelfCollision() filters codimension >= 2
            // (vertex-level) self-collisions -- disabled along with Phase 2
            // in addFace() (embeddedsubmanifold.cpp), per the conjecture
            // that these are always resolvable cusp intersections. Only
            // the codimension-1 (facet-level) exclusion above remains.
            //
            // || EmbeddedSubmanifold<dim, subdim>::hasUnexplainedSelfCollision(
            //        nodes[i].face, nodes[i].gluings)
        )
            continue;
        skelToGraph[i] = static_cast<int>(skelOf.size());
        skelOf.push_back(static_cast<int>(i));
    }

    int n = static_cast<int>(skelOf.size());
    // 1-indexed, as the enumerator expects
    std::vector<std::vector<int>> adj(n + 1);
    for (int graphIdx = 0; graphIdx < n; ++graphIdx) {
        int i = skelOf[graphIdx];
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

    std::vector<std::vector<int>> graphToSkel;
    graphToSkel.reserve(n);
    for (int skelIdx : skelOf)
        graphToSkel.push_back({skelIdx});

    return {AdjacencyList{n, std::move(adj)}, std::move(graphToSkel)};
}

template <int dim, int subdim>
typename EmbeddingSearch<dim, subdim>::Graph
EmbeddingSearch<dim, subdim>::buildSeededGraph_(
    const Skeleton<dim, subdim> &skeleton, const std::vector<int> &seedFaces) {
    Graph base = buildGraph_(skeleton);

    // Invert base.graphToSkel: skeleton index -> graph id (1-indexed).
    std::vector<int> skelToGraph(skeleton.numFaces(), -1);
    for (int i = 0; i < base.adjList.first; ++i)
        skelToGraph[base.graphToSkel[i][0]] = i + 1;

    std::vector<int> seedGraphIds;
    seedGraphIds.reserve(seedFaces.size());
    for (int f : seedFaces) {
        int gid = skelToGraph[f];
        if (gid == -1) {
            const auto &node = skeleton.getNodes()[f];
            std::ostringstream reason;
            if (EmbeddedSubmanifold<dim, subdim>::hasIrreparableSelfGluing(
                    node.gluings))
                reason << "an irreparable self-gluing";
            // hasUnexplainedSelfCollision() usage disabled along with
            // Phase 2 in addFace() -- see buildGraph_() above.
            //
            // else if (EmbeddedSubmanifold<dim, subdim>::
            //              hasUnexplainedSelfCollision(node.face, node.gluings))
            //     reason << "an unexplained self-collision";
            else
                reason << "an unknown reason"; // shouldn't happen
            throw regina::InvalidArgument(
                "EmbeddingSearch: seed face " + std::to_string(f) +
                " is not embeddable (excluded by buildGraph_'s pre-filter: " +
                reason.str() + ")");
        }
        seedGraphIds.push_back(gid);
    }

    ConnectedInducedSubgraphEnumerator::SeededGraph sg =
        ConnectedInducedSubgraphEnumerator::contractSeed(
            base.adjList.first, base.adjList.second, seedGraphIds);

    Graph result;
    result.adjList = AdjacencyList{sg.n, std::move(sg.adj)};
    result.graphToSkel.resize(sg.n);
    result.graphToSkel[0] = seedFaces; // vertex 1 -> the whole seed
    for (int i = 1; i < sg.n; ++i) {
        int oldGraphId = sg.originalOf[i + 1];
        result.graphToSkel[i] = base.graphToSkel[oldGraphId - 1];
    }
    return result;
}

template class EmbeddingSearch<3, 2>;
template class EmbeddingSearch<4, 2>;

void SurfaceSearch::SurfaceTypeTally::merge(
    std::map<SurfaceTypeKey, long long> &local) {
    if (local.empty())
        return;
    std::lock_guard<std::mutex> lock(mutex_);
    for (const auto &[key, n] : local)
        counts_[key] += n;
    local.clear();
}

std::string SurfaceSearch::SurfaceTypeTally::summary() const {
    std::lock_guard<std::mutex> lock(mutex_);
    std::ostringstream out;
    out << "Number of surfaces found:\n";
    if (counts_.empty()) {
        out << "  (none yet)\n";
    } else {
        for (const auto &[key, n] : counts_)
            out << "  " << KnottedSurface::formatSurfaceType(key) << " = " << n
                << "\n";
    }
    return out.str();
}

void SurfaceSearch::PendingSurfaceBatch::merge(
    std::vector<std::vector<int>> &local) {
    if (local.empty())
        return;
    std::lock_guard<std::mutex> lock(mutex_);
    for (auto &faceIndices : local)
        pending_.push_back(std::move(faceIndices));
    local.clear();
}

std::vector<std::vector<int>> SurfaceSearch::PendingSurfaceBatch::drain() {
    std::lock_guard<std::mutex> lock(mutex_);
    std::vector<std::vector<int>> out;
    std::swap(out, pending_);
    return out;
}

void SurfaceSearch::LinkBoundaryTally::record(const std::string &linkName,
                                              const SurfaceTypeKey &type) {
    std::lock_guard<std::mutex> lock(mutex_);
    ++linkSurfaceTypes_[linkName][type];
}

std::string SurfaceSearch::LinkBoundaryTally::summary() const {
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
                out << KnottedSurface::formatSurfaceType(type) << " (" << count
                    << ")";
                first = false;
            }
            out << "\n";
        }
    }
    return out.str();
}

void SurfaceSearch::processSurfaceBoundaries() {
    auto batch = pendingSurfaces_.drain();
    if (batch.empty())
        return;

    KnottedSurface embedding(skeleton_);
    processBatchRange_(embedding, batch, 0, batch.size());
}

void SurfaceSearch::processRemainingSurfaceBoundaries(unsigned numThreads) {
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
                std::chrono::duration_cast<std::chrono::seconds>(elapsed)
                    .count();
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
                long long etaSeconds =
                    elapsedSeconds * static_cast<long long>(total - processed) /
                    static_cast<long long>(processed);
                report << "[+] ETA: "
                       << formatElapsed(std::chrono::seconds(etaSeconds))
                       << "\n";
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
        KnottedSurface embedding(skeleton_);
        while (true) {
            size_t begin =
                nextIndex.fetch_add(CHUNK, std::memory_order_relaxed);
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

void SurfaceSearch::processBatchRange_(
    KnottedSurface &embedding, const std::vector<std::vector<int>> &batch,
    size_t begin, size_t end) {
    for (size_t i = begin; i < end; ++i) {
        const auto &faceIndices = batch[i];
        for (int idx : faceIndices)
            embedding.addFace(idx);

        SurfaceTypeKey type = embedding.surfaceType();
        for (auto &[component, link] : embedding.boundaryLinks())
            linkTally_.record(link.identify(), type);

        // Reverse order, mirroring how the DFS itself would back out --
        // resets embedding to empty for the next item in the batch.
        for (auto it = faceIndices.rbegin(); it != faceIndices.rend(); ++it)
            embedding.removeFace(*it);
    }
}

void SurfaceSearch::search(unsigned numThreads, BoundaryCondition cond) {
    const auto searchStart = std::chrono::steady_clock::now();

    const bool wantLinks = cond == BoundaryCondition::proper ||
                           cond == BoundaryCondition::connected;

    // See EmbeddingSearch<dim,subdim>::search() for the rationale -- shared
    // root work-list built once, up front, also handling the "seed alone"
    // result here since no root's subtree can ever produce it.
    std::vector<int> roots;
    long long seedFoundCount = 0;
    long long seedSubgraphCount = 0;
    long long seedMaxFaces = 0;
    long long seedFaceSum = 0;
    std::map<SurfaceTypeKey, long long> seedTypeCounts;
    std::vector<std::pair<std::string, SurfaceTypeKey>> seedLinks;
    if (isSeeded_) {
        KnottedSurface protoEmbedding(skeleton_);
        EmbeddednessPredicate protoPredicate(protoEmbedding,
                                             graph_.graphToSkel);
        ConnectedInducedSubgraphEnumerator protoEnumerator(
            graph_.adjList.first, graph_.adjList.second, true, protoPredicate);
        roots = protoEnumerator.getRoots();

        seedFoundCount = 1;
        if (protoEmbedding.isEmbedded() && protoEmbedding.satisfies(cond)) {
            seedSubgraphCount = 1;
            seedMaxFaces = static_cast<long long>(
                protoEmbedding.triangulation().size());
            seedFaceSum = seedMaxFaces;
            SurfaceTypeKey type = protoEmbedding.surfaceType();
            ++seedTypeCounts[type];
            if (wantLinks)
                for (auto &[component, link] : protoEmbedding.boundaryLinks())
                    seedLinks.emplace_back(link.identify(), type);
        }
    } else {
        roots.resize(graph_.adjList.first);
        std::iota(roots.begin(), roots.end(), 1);
    }
    const size_t totalRoots = roots.size();

    std::atomic<size_t> nextRootIdx{0};
    std::vector<long long> perThreadCount(numThreads, 0);

    std::atomic<long long> globalFoundCount{seedFoundCount};
    std::atomic<long long> globalSubgraphCount{seedSubgraphCount};
    std::atomic<long long> globalMaxFaces{seedMaxFaces};
    std::atomic<long long> globalFaceSum{seedFaceSum};
    std::atomic<size_t> rootsCompleted{0};
    std::atomic<bool> workersFinished{false};
    SurfaceTypeTally surfaceTally;
    surfaceTally.merge(seedTypeCounts);
    for (auto &[linkName, type] : seedLinks)
        linkTally_.record(linkName, type);
    std::vector<long long> perThreadFoundCount(numThreads, 0);
    std::vector<long long> perThreadFaceSum(numThreads, 0);

    auto worker = [&](unsigned tid) {
        EmbeddedSubmanifold<4, 2> embedding(skeleton_);
        EmbeddednessPredicate predicate(embedding, graph_.graphToSkel);
        std::optional<ConnectedInducedSubgraphEnumerator> localEnumeratorOpt;
        if (isSeeded_)
            localEnumeratorOpt.emplace(graph_.adjList.first,
                                       graph_.adjList.second, true, predicate);
        else
            localEnumeratorOpt.emplace(graph_.adjList.first,
                                       graph_.adjList.second);
        auto &localEnumerator = *localEnumeratorOpt;
        long long localCount = 0;
        long long sinceFlush = 0;
        long long localFoundCount = 0;
        long long sinceFlushFound = 0;
        long long localFaceSum = 0;
        long long sinceFlushFaceSum = 0;
        std::map<SurfaceTypeKey, long long> localTypeCounts;
        std::vector<std::vector<int>> localPending;

        while (true) {
            size_t idx = nextRootIdx.fetch_add(1);
            if (idx >= totalRoots)
                break;
            int s = roots[idx];
            localEnumerator.enumerateFromRootFiltered(
                s,
                [&](const std::vector<int> &U) {
                    ++localFoundCount;
                    if (++sinceFlushFound >= FLUSH_EVERY_FOUND) {
                        globalFoundCount.fetch_add(sinceFlushFound,
                                                   std::memory_order_relaxed);
                        sinceFlushFound = 0;
                    }
                    if (embedding.isEmbedded() && embedding.satisfies(cond)) {
                        ++localCount;
                        ++localTypeCounts[KnottedSurface::surfaceTypeKey(
                            embedding.triangulation())];
                        if (wantLinks) {
                            std::vector<int> faceIndices;
                            for (int v : U)
                                for (int f : graph_.graphToSkel[v - 1])
                                    faceIndices.push_back(f);
                            localPending.push_back(std::move(faceIndices));
                        }
                        long long faceCount = static_cast<long long>(
                            embedding.triangulation().size());
                        localFaceSum += faceCount;
                        sinceFlushFaceSum += faceCount;
                        if (++sinceFlush >= FLUSH_EVERY_BDRY) {
                            globalSubgraphCount.fetch_add(
                                sinceFlush, std::memory_order_relaxed);
                            globalFaceSum.fetch_add(
                                sinceFlushFaceSum, std::memory_order_relaxed);
                            sinceFlush = 0;
                            sinceFlushFaceSum = 0;
                            surfaceTally.merge(localTypeCounts);
                            if (wantLinks)
                                pendingSurfaces_.merge(localPending);
                        }
                        auto prevMax =
                            globalMaxFaces.load(std::memory_order_relaxed);
                        while (faceCount > prevMax &&
                               !globalMaxFaces.compare_exchange_weak(
                                   prevMax, faceCount,
                                   std::memory_order_relaxed))
                            ;
                    }
                },
                predicate);
            rootsCompleted.fetch_add(1, std::memory_order_relaxed);
        }
        // flush this thread's remainder so the global count ends up
        // exact
        if (sinceFlush > 0) {
            globalSubgraphCount.fetch_add(sinceFlush,
                                          std::memory_order_relaxed);
            globalFaceSum.fetch_add(sinceFlushFaceSum,
                                    std::memory_order_relaxed);
        }
        if (sinceFlushFound > 0)
            globalFoundCount.fetch_add(sinceFlushFound,
                                       std::memory_order_relaxed);
        surfaceTally.merge(localTypeCounts);
        if (wantLinks)
            pendingSurfaces_.merge(localPending);
        perThreadCount[tid] = localCount;
        perThreadFoundCount[tid] = localFoundCount;
        perThreadFaceSum[tid] = localFaceSum;
    };

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
            report << "[+] roots completed: " << rootsCompleted.load() << "/"
                   << totalRoots
                   << "  | embedded submanifolds found so far"
                   << scaleLabel(FLUSH_EVERY_FOUND) << ": "
                   << (globalFoundCount.load() / FLUSH_EVERY_FOUND) << "\n";
            report << "[+] embedded submanifolds satisfying boundary "
                      "condition ("
                   << boundaryConditionName(cond) << ") found so far"
                   << scaleLabel(FLUSH_EVERY_BDRY) << ": "
                   << (globalSubgraphCount.load() / FLUSH_EVERY_BDRY) << "\n";
            report << "[+] largest embedded submanifold satisfying boundary "
                      "condition ("
                   << boundaryConditionName(cond)
                   << ") found so far: " << globalMaxFaces.load()
                   << " faces\n";
            {
                long long n = globalSubgraphCount.load();
                double avg = n > 0 ? static_cast<double>(globalFaceSum.load()) /
                                         static_cast<double>(n)
                                   : 0.0;
                report << "[+] average faces per embedded submanifold "
                          "satisfying boundary condition ("
                       << boundaryConditionName(cond)
                       << ") found so far: " << std::fixed
                       << std::setprecision(2) << avg << "\n";
            }
            report << surfaceTally.summary();
            if (wantLinks)
                report << linkTally_.summary();
            std::string text = report.str();

            if (prevLines > 0)
                std::cerr << "\x1b[" << prevLines << "F\x1b[0J";
            std::cerr << text;
            prevLines = std::count(text.begin(), text.end(), '\n');
        }
    });

    std::thread boundaryProcessor;
    if (wantLinks) {
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

    std::vector<std::thread> threads;
    threads.reserve(numThreads);
    for (unsigned t = 0; t < numThreads; ++t)
        threads.emplace_back(worker, t);
    for (auto &th : threads)
        th.join();

    workersFinished.store(true, std::memory_order_relaxed);
    reporter.join();
    if (boundaryProcessor.joinable()) {
        boundaryProcessor.join();
        processRemainingSurfaceBoundaries(numThreads);
    }

    long long total = seedSubgraphCount;
    for (long long c : perThreadCount)
        total += c;
    long long totalFound = seedFoundCount;
    for (long long c : perThreadFoundCount)
        totalFound += c;
    long long totalFaceSum = seedFaceSum;
    for (long long c : perThreadFaceSum)
        totalFaceSum += c;
    std::cerr << "\n\nTotal elapsed: "
              << formatElapsed(std::chrono::steady_clock::now() - searchStart)
              << "\n";
    std::cerr << "Total embedded submanifolds found: " << totalFound
              << " (across " << numThreads << " threads)\n";
    std::cerr << "Total embedded submanifolds (" << boundaryConditionName(cond)
              << "): " << total << " (across " << numThreads << " threads)\n";
    std::cerr << "Largest embedded submanifold (" << boundaryConditionName(cond)
              << "): " << globalMaxFaces.load() << " faces\n";
    std::cerr << "Average faces per embedded submanifold ("
              << boundaryConditionName(cond) << "): " << std::fixed
              << std::setprecision(2)
              << (total > 0 ? static_cast<double>(totalFaceSum) /
                                  static_cast<double>(total)
                            : 0.0)
              << "\n";
    std::cerr << surfaceTally.summary() << "\n";
    if (wantLinks)
        std::cerr << linkTally_.summary() << "\n";

    assert(globalFoundCount.load() == totalFound);
    assert(globalSubgraphCount.load() == total);
    assert(globalFaceSum.load() == totalFaceSum);
}
