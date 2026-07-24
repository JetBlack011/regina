//
//  embeddingsearch.cpp
//

#include "embeddingsearch.h"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <csignal>
#include <cmath>
#include <cstdlib>
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
#define FLUSH_EVERY_FOUND 10'000

namespace {
// Ctrl+C handling for runSearch_() (embeddingsearch.h's doc comment on that
// method has the user-facing contract). g_stopRequested/g_sigintCount are
// process-wide, not per-search -- signals are inherently process-wide too,
// and this tool only ever runs one search at a time, so that's not a
// practical limitation here.
//
// handleSigint only calls std::atomic<bool>::store (on an
// always-lock-free-in-practice atomic) and std::_Exit -- per
// [support.signal], these are among the few operations the standard
// actually guarantees are safe to perform from within a signal handler.
std::atomic<bool> *g_stopRequested = nullptr;
std::atomic<int> g_sigintCount{0};

void handleSigint(int) {
    if (g_sigintCount.fetch_add(1, std::memory_order_relaxed) > 0)
        std::_Exit(130); // 128 + SIGINT, the usual shell convention
    if (g_stopRequested)
        g_stopRequested->store(true, std::memory_order_relaxed);
}

// RAII: installs handleSigint for SIGINT for the lifetime of one
// runSearch_() call (constructed/destroyed on the calling thread only,
// before any worker thread exists/after they've all joined), restoring
// whatever handler was previously in place afterward.
class SigintScope {
    using Handler = void (*)(int);
    Handler previous_;

  public:
    explicit SigintScope(std::atomic<bool> &flag) {
        g_stopRequested = &flag;
        g_sigintCount.store(0, std::memory_order_relaxed);
        previous_ = std::signal(SIGINT, handleSigint);
    }

    ~SigintScope() {
        std::signal(SIGINT, previous_);
        g_stopRequested = nullptr;
    }
};
} // namespace

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

    EmbeddedSubmanifold<dim, subdim>(skeleton_, seedFaces);
}

template <int dim, int subdim>
template <typename ThreadHookFactory, typename OnSeedFound,
          typename ExtraReportLines, typename AuxHooks>
void EmbeddingSearch<dim, subdim>::runSearch_(
    unsigned numThreads, BoundaryCondition cond,
    ThreadHookFactory makeThreadHook, OnSeedFound onSeedFound,
    ExtraReportLines extraReportLines, AuxHooks auxHooks) {
    const auto searchStart = std::chrono::steady_clock::now();

    // See this method's doc comment (embeddingsearch.h) for the SIGINT
    // contract. Constructed before the seeded proto-embedding check below
    // and destroyed only once this whole function returns, so a Ctrl+C
    // during that check, the worker threads, or auxHooks.afterJoin() is all
    // handled uniformly.
    std::atomic<bool> stopRequested{false};
    SigintScope sigintScope(stopRequested);

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
        InterruptiblePredicate protoInterruptible(protoPredicate,
                                                  stopRequested);
        ConnectedInducedSubgraphEnumerator protoEnumerator(
            graph_.adjList.first, graph_.adjList.second, true,
            protoInterruptible);
        roots = protoEnumerator.getRoots();

        seedFoundCount = 1;
        if (protoEmbedding.isEmbedded() && protoEmbedding.satisfies(cond)) {
            seedSubgraphCount = 1;
            seedMaxFaces = static_cast<long long>(
                protoEmbedding.triangulation().size());
            seedFaceSum = seedMaxFaces;
            // graph_.graphToSkel[0] is the whole seed -- see
            // buildSeededGraph_.
            onSeedFound(graph_.graphToSkel[0]);
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
        InterruptiblePredicate interruptible(predicate, stopRequested);
        // ConnectedInducedSubgraphEnumerator holds a reference member, so
        // it isn't assignable -- construct it in place via optional rather
        // than choosing between two constructor calls via assignment.
        std::optional<ConnectedInducedSubgraphEnumerator> localEnumeratorOpt;
        if (isSeeded_)
            localEnumeratorOpt.emplace(graph_.adjList.first,
                                       graph_.adjList.second, true,
                                       interruptible);
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
        auto threadHook = makeThreadHook();

        while (true) {
            // Once interrupted, InterruptiblePredicate already prunes any
            // root already in progress down to nothing -- this just skips
            // starting fresh roots too, instead of doing so only to have
            // them immediately rejected.
            if (stopRequested.load(std::memory_order_relaxed))
                break;
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
                        auto faceCount = static_cast<long long>(
                            embedding.triangulation().size());
                        threadHook.onFound(embedding, U, faceCount);
                        localFaceSum += faceCount;
                        sinceFlushFaceSum += faceCount;
                        if (++sinceFlush >= FLUSH_EVERY_BDRY) {
                            globalSubgraphCount.fetch_add(
                                sinceFlush, std::memory_order_relaxed);
                            globalFaceSum.fetch_add(
                                sinceFlushFaceSum, std::memory_order_relaxed);
                            sinceFlush = 0;
                            sinceFlushFaceSum = 0;
                            threadHook.onFlush();
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
                interruptible);
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
        threadHook.onFlush();
        perThreadCount[tid] = localCount;
        perThreadFoundCount[tid] = localFoundCount;
        perThreadFaceSum[tid] = localFaceSum;
    };

    std::thread reporter([&]() {
        using namespace std::chrono_literals;
        size_t prevLines = 0;
        while (!workersFinished.load(std::memory_order_relaxed)) {
            std::this_thread::sleep_for(1s);

            // extraReportLines() (surface-type/boundary tallies) goes
            // first, status/elapsed/ETA last: the tallies can grow into a
            // long list, and keeping the status lines below it means
            // they're always the last thing printed -- right above the
            // cursor -- instead of getting pushed off-screen above a long
            // list.
            std::ostringstream report;
            report << extraReportLines();
            report << "[+] elapsed: "
                   << formatElapsed(std::chrono::steady_clock::now() -
                                    searchStart)
                   << "\n";
            report << "[+] roots completed: " << rootsCompleted.load() << "/"
                   << totalRoots
                   << "  | embedded submanifolds found so far" << ": "
                   << globalFoundCount.load() << "\n";
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

    std::thread aux = auxHooks.spawn(workersFinished);

    std::vector<std::thread> threads;
    threads.reserve(numThreads);
    for (unsigned t = 0; t < numThreads; ++t)
        threads.emplace_back(worker, t);
    for (auto &th : threads)
        th.join();

    workersFinished.store(true, std::memory_order_relaxed);
    reporter.join();

    // Printed only after reporter.join(), not before: the reporter thread
    // may already be past its own workersFinished check (asleep mid-tick)
    // when that flag flips, so it can still fire one more ANSI-erasing
    // redraw afterward -- printing this message any earlier meant that
    // redraw would immediately wipe it out.
    if (stopRequested.load(std::memory_order_relaxed))
        std::cerr << "\n[!] Interrupted -- moving on to whatever comes next "
                     "with what was found so far (Ctrl+C again to quit "
                     "immediately)\n";

    if (aux.joinable()) {
        aux.join();
        auxHooks.afterJoin();
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
    // extraReportLines() first, totals last -- same reasoning as the live
    // reporter above: keeps the totals visible right above the cursor
    // instead of pushed off-screen above a long tally list.
    std::string extra = extraReportLines();
    if (!extra.empty())
        std::cerr << "\n" << extra;
    std::cerr << "\nTotal elapsed: "
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
void EmbeddingSearch<dim, subdim>::search(const unsigned numThreads,
                                          BoundaryCondition cond) {
    struct NoopThreadHook {
        void onFound(EmbeddedSubmanifold<dim, subdim> &,
                    const std::vector<int> &, long long) {}
        void onFlush() {}
    };
    struct NoopAuxHooks {
        std::thread spawn(std::atomic<bool> &) { return {}; }
        void afterJoin() {}
    };
    runSearch_(
        numThreads, cond, [] { return NoopThreadHook{}; },
        [](const std::vector<int> &) {}, [] { return std::string(); },
        NoopAuxHooks{});
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

std::vector<std::vector<int>>
SurfaceSearch::PendingSurfaceBatch::popSome(size_t maxCount) {
    std::lock_guard<std::mutex> lock(mutex_);
    size_t count = std::min(maxCount, pending_.size());
    std::vector<std::vector<int>> out;
    out.reserve(count);
    for (size_t i = 0; i < count; ++i)
        out.push_back(std::move(pending_[pending_.size() - count + i]));
    pending_.resize(pending_.size() - count);
    return out;
}

void SurfaceSearch::LinkBoundaryTally::record(const std::string &descriptor,
                                              const SurfaceTypeKey &type) {
    std::lock_guard<std::mutex> lock(mutex_);
    ++descriptorSurfaceTypes_[descriptor][type];
}

std::string SurfaceSearch::LinkBoundaryTally::summary() const {
    std::lock_guard<std::mutex> lock(mutex_);
    std::ostringstream out;
    out << "Surface boundaries found:\n";
    if (descriptorSurfaceTypes_.empty()) {
        out << "  (none yet)\n";
    } else {
        std::vector<std::string> descriptors;
        descriptors.reserve(descriptorSurfaceTypes_.size());
        for (const auto &[descriptor, _] : descriptorSurfaceTypes_)
            descriptors.push_back(descriptor);
        std::sort(descriptors.begin(), descriptors.end());

        for (const std::string &descriptor : descriptors) {
            out << "  " << descriptor << " - ";
            bool first = true;
            for (const auto &[type, count] :
                 descriptorSurfaceTypes_.at(descriptor)) {
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

std::string SurfaceSearch::describeBoundary_(
    const std::vector<std::pair<size_t, Link>> &links) {
    std::ostringstream out;
    bool firstComponent = true;
    for (const auto &[component, link] : links) {
        if (!firstComponent)
            out << ", ";
        firstComponent = false;

        out << (component + 1) << ": ";
        bool firstCurve = true;
        for (const Knot &curve : link.comps_) {
            if (!firstCurve)
                out << ", ";
            firstCurve = false;
            out << curve.identify();
        }
        if (link.comps_.size() > 1)
            out << " (" << link.identify() << ")";
    }
    return out.str();
}

void SurfaceSearch::ThreadHook::onFound(EmbeddedSubmanifold<4, 2> &embedding,
                                        const std::vector<int> &U,
                                        long long) {
    ++localTypeCounts_[KnottedSurface::surfaceTypeKey(embedding.triangulation())];
    if (wantLinks_) {
        std::vector<int> faceIndices;
        for (int v : U)
            for (int f : owner_.graph_.graphToSkel[v - 1])
                faceIndices.push_back(f);
        localPending_.push_back(std::move(faceIndices));
    }
}

void SurfaceSearch::ThreadHook::onFlush() {
    tally_.merge(localTypeCounts_);
    if (wantLinks_)
        owner_.pendingSurfaces_.merge(localPending_);
}

std::thread SurfaceSearch::AuxHooks::spawn(std::atomic<bool> &workersFinished) {
    if (!wantLinks_)
        return {};
    return std::thread([this, &workersFinished]() {
        owner_.backgroundDrainLoop_(workersFinished);
    });
}

void SurfaceSearch::AuxHooks::afterJoin() {
    owner_.processRemainingSurfaceBoundaries(numThreads_);
}

void SurfaceSearch::backgroundDrainLoop_(
    const std::atomic<bool> &workersFinished) {
    using namespace std::chrono_literals;
    // Matches processBatchParallel_'s own CHUNK size -- not load-bearing
    // that they're equal, just a reasonable shared "small batch" constant.
    constexpr size_t POP_BATCH = 64;
    KnottedSurface embedding(skeleton_);
    while (!workersFinished.load(std::memory_order_relaxed)) {
        auto items = pendingSurfaces_.popSome(POP_BATCH);
        if (items.empty()) {
            std::this_thread::sleep_for(20ms);
            continue;
        }
        for (size_t i = 0; i < items.size(); ++i) {
            if (workersFinished.load(std::memory_order_relaxed)) {
                // The search ended partway through this small batch: hand
                // the unprocessed remainder back to the shared queue --
                // popSome() already removed it from pending_, so
                // processRemainingSurfaceBoundaries()'s drain() would
                // never otherwise see it.
                std::vector<std::vector<int>> remainder(
                    std::make_move_iterator(items.begin() +
                                            static_cast<ptrdiff_t>(i)),
                    std::make_move_iterator(items.end()));
                pendingSurfaces_.merge(remainder);
                return;
            }
            processEntry_(embedding, items[i]);
        }
    }
}

void SurfaceSearch::processRemainingSurfaceBoundaries(unsigned numThreads) {
    processBatchParallel_(pendingSurfaces_.drain(), numThreads);
}

void SurfaceSearch::processBatchParallel_(std::vector<std::vector<int>> batch,
                                          unsigned numThreads) {
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

            // linkTally_.summary() (the potentially long boundary list)
            // first, elapsed/processed/ETA last -- see runSearch_'s live
            // reporter for the same reasoning: keeps status visible right
            // above the cursor instead of pushed off-screen above a long
            // list.
            std::ostringstream report;
            report << linkTally_.summary();
            report << "[+] elapsed: " << formatElapsed(elapsed) << "\n";
            report << "[+] boundaries processed: " << processed << "/"
                   << total << " (" << std::fixed << std::setprecision(2)
                   << (total > 0 ? 100.0 * static_cast<double>(processed) /
                                       static_cast<double>(total)
                                 : 0.0)
                   << "%)\n";
            if (processed > 0 && processed < total && elapsedSeconds > 0) {
                long long etaSeconds = elapsedSeconds *
                                       static_cast<long long>(total - processed) /
                                       static_cast<long long>(processed);
                report << "[+] ETA: "
                       << formatElapsed(std::chrono::seconds(etaSeconds))
                       << "\n";
            }
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

void SurfaceSearch::processEntry_(KnottedSurface &embedding,
                                  const std::vector<int> &faceIndices) {
    for (int idx : faceIndices)
        embedding.addFace(idx);

    SurfaceTypeKey type = embedding.surfaceType();
    auto links = embedding.boundaryLinks();
    if (!links.empty())
        linkTally_.record(describeBoundary_(links), type);

    // Reverse order, mirroring how the DFS itself would back out --
    // resets embedding to empty for the next entry.
    for (auto it = faceIndices.rbegin(); it != faceIndices.rend(); ++it)
        embedding.removeFace(*it);
}

void SurfaceSearch::processBatchRange_(
    KnottedSurface &embedding, const std::vector<std::vector<int>> &batch,
    size_t begin, size_t end) {
    for (size_t i = begin; i < end; ++i)
        processEntry_(embedding, batch[i]);
}

void SurfaceSearch::search(unsigned numThreads, BoundaryCondition cond) {
    const bool wantLinks = cond == BoundaryCondition::proper ||
                           cond == BoundaryCondition::connected;
    SurfaceTypeTally surfaceTally;

    runSearch_(
        numThreads, cond,
        [this, &surfaceTally, wantLinks] {
            return ThreadHook(*this, surfaceTally, wantLinks);
        },
        [this, &surfaceTally, wantLinks](const std::vector<int> &seedFaces) {
            // One-off, not hot-path: rebuilds the seed as a KnottedSurface
            // (rather than reusing runSearch_'s generic proto embedding) to
            // get at surfaceType()/boundaryLinks(), which only KnottedSurface
            // exposes. seedFaces is added in the same order used to validate
            // it originally (see buildSeededGraph_), so this always embeds.
            KnottedSurface probe(skeleton_, seedFaces);
            SurfaceTypeKey type = probe.surfaceType();
            std::map<SurfaceTypeKey, long long> seedTypeCounts{{type, 1}};
            surfaceTally.merge(seedTypeCounts);
            if (wantLinks) {
                auto links = probe.boundaryLinks();
                if (!links.empty())
                    linkTally_.record(describeBoundary_(links), type);
            }
        },
        [this, &surfaceTally, wantLinks] {
            std::string out = surfaceTally.summary();
            if (wantLinks)
                out += linkTally_.summary();
            return out;
        },
        AuxHooks(*this, numThreads, wantLinks));
}
