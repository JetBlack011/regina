//
//  embeddingsearch.cpp
//
//  Created by John Teague on 07/15/2026.
//

#include "embeddingsearch.h"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <csignal>
#include <cstdlib>
#include <iomanip>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <thread>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

// satisfying-cond finds are rare relative to raw finds, so flush every one
// to keep the rolling progress report responsive.
#define FLUSH_EVERY_BDRY 1
// Like satisfying-cond finds, isEmbedded()-true finds are rare relative to
// raw finds (most incrementally-valid partial complexes still fail the
// final global embeddedness check) and can be exactly as common as
// satisfying-cond finds (e.g. under BoundaryCondition::all, satisfies() is
// trivially true whenever isEmbedded() is) -- so this is flushed every one,
// same as FLUSH_EVERY_BDRY, to keep the live embedded/sec rate accurate
// rather than lagging behind the (immediately-flushed) satisfying-count
// display.
#define FLUSH_EVERY_EMBEDDED 1
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

// The numeric counters shared across all worker threads in one runSearch_()
// call, updated concurrently as the search runs and periodically
// snapshotted (see runSearch_'s snapshotStats) into a caller-facing
// SearchStats. Grouped into one struct purely so the atomics it takes are
// declared/passed together instead of as four same-typed variables that are
// easy to mix up individually.
struct AtomicSearchStats {
    std::atomic<long long> foundCount;        // raw candidates visited, any cond
    std::atomic<long long> embeddedCount;     // isEmbedded()==true, any cond
    std::atomic<long long> satisfyingCount;   // satisfying the BoundaryCondition
    std::atomic<long long> largestSatisfying; // max faces among satisfying finds
    std::atomic<long long> satisfyingFaceSum;
    std::atomic<size_t> rootsCompleted{0};
};

// One worker thread's own running totals, plus the portion of each not yet
// folded into the shared AtomicSearchStats (flushed periodically -- see
// FLUSH_EVERY_BDRY/FLUSH_EVERY_FOUND -- rather than on every single find, to
// keep contention on the shared atomics down). One instance per thread,
// read back after all threads join to compute the exact final totals.
struct WorkerStats {
    long long foundCount = 0;
    long long embeddedCount = 0;
    long long satisfyingCount = 0;
    long long satisfyingFaceSum = 0;
    long long pendingFoundCount = 0;
    long long pendingEmbeddedCount = 0;
    long long pendingSatisfyingCount = 0;
    long long pendingFaceSum = 0;
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

// EmbeddednessPredicate's members are defined inline in embeddingsearch.h,
// not here: it is a member *template*, so the explicit instantiations of
// EmbeddingSearch<dim,subdim> at the bottom of this file do not force its
// specializations out, and its trivial constructor gets inlined away rather
// than emitted. Defining it here left any other translation unit that
// constructs one (e.g. tests/profile_driver.cpp) unable to link.

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
template <typename EmbeddingFactory, typename ThreadHookFactory,
          typename OnSeedFound, typename AuxHooks>
SearchStats EmbeddingSearch<dim, subdim>::runSearch_(
    unsigned numThreads, BoundaryCondition cond,
    EmbeddingFactory makeEmbedding, ThreadHookFactory makeThreadHook,
    OnSeedFound onSeedFound, const SearchCallbacks &callbacks,
    AuxHooks auxHooks) {
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
    long long seedEmbeddedCount = 0;
    long long seedSubgraphCount = 0;
    long long seedMaxFaces = 0;
    long long seedFaceSum = 0;
    if (isSeeded_) {
        auto protoEmbedding = makeEmbedding();
        EmbeddednessPredicate protoPredicate(protoEmbedding,
                                             graph_.graphToSkel);
        InterruptiblePredicate protoInterruptible(protoPredicate,
                                                  stopRequested);
        ConnectedInducedSubgraphEnumerator protoEnumerator(
            graph_.adjList.first, graph_.adjList.second, true,
            protoInterruptible);
        roots = protoEnumerator.getRoots();

        seedFoundCount = 1;
        if (protoEmbedding.isEmbedded()) {
            seedEmbeddedCount = 1;
            if (protoEmbedding.satisfies(cond)) {
                seedSubgraphCount = 1;
                seedMaxFaces = static_cast<long long>(
                    protoEmbedding.triangulation().size());
                seedFaceSum = seedMaxFaces;
                // graph_.graphToSkel[0] is the whole seed -- see
                // buildSeededGraph_.
                onSeedFound(graph_.graphToSkel[0]);
            }
        }
    } else {
        roots.resize(graph_.adjList.first);
        std::iota(roots.begin(), roots.end(), 1);
    }
    const size_t totalRoots = roots.size();

    std::atomic<size_t> nextRootIdx{0};
    std::atomic<bool> workersFinished{false};
    AtomicSearchStats stats{.foundCount = seedFoundCount,
                            .embeddedCount = seedEmbeddedCount,
                            .satisfyingCount = seedSubgraphCount,
                            .largestSatisfying = seedMaxFaces,
                            .satisfyingFaceSum = seedFaceSum};
    std::vector<WorkerStats> perThreadStats(numThreads);

    auto worker = [&](unsigned tid) {
        auto embedding = makeEmbedding();
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
        WorkerStats &local = perThreadStats[tid];
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
                    ++local.foundCount;
                    if (++local.pendingFoundCount >= FLUSH_EVERY_FOUND) {
                        stats.foundCount.fetch_add(local.pendingFoundCount,
                                                   std::memory_order_relaxed);
                        local.pendingFoundCount = 0;
                    }
                    if (embedding.isEmbedded()) {
                        ++local.embeddedCount;
                        if (++local.pendingEmbeddedCount >= FLUSH_EVERY_EMBEDDED) {
                            stats.embeddedCount.fetch_add(
                                local.pendingEmbeddedCount,
                                std::memory_order_relaxed);
                            local.pendingEmbeddedCount = 0;
                        }
                        if (embedding.satisfies(cond)) {
                            ++local.satisfyingCount;
                            auto faceCount = static_cast<long long>(
                                embedding.triangulation().size());
                            threadHook.onFound(embedding, U, faceCount);
                            local.satisfyingFaceSum += faceCount;
                            local.pendingFaceSum += faceCount;
                            if (++local.pendingSatisfyingCount >= FLUSH_EVERY_BDRY) {
                                stats.satisfyingCount.fetch_add(
                                    local.pendingSatisfyingCount,
                                    std::memory_order_relaxed);
                                stats.satisfyingFaceSum.fetch_add(
                                    local.pendingFaceSum, std::memory_order_relaxed);
                                local.pendingSatisfyingCount = 0;
                                local.pendingFaceSum = 0;
                                threadHook.onFlush();
                            }
                            auto prevMax = stats.largestSatisfying.load(
                                std::memory_order_relaxed);
                            while (faceCount > prevMax &&
                                   !stats.largestSatisfying.compare_exchange_weak(
                                       prevMax, faceCount,
                                       std::memory_order_relaxed))
                                ;
                        }
                    }
                },
                interruptible);
            stats.rootsCompleted.fetch_add(1, std::memory_order_relaxed);
        }
        // flush this thread's remainder so the global count ends up
        // exact
        if (local.pendingSatisfyingCount > 0) {
            stats.satisfyingCount.fetch_add(local.pendingSatisfyingCount,
                                            std::memory_order_relaxed);
            stats.satisfyingFaceSum.fetch_add(local.pendingFaceSum,
                                              std::memory_order_relaxed);
        }
        if (local.pendingEmbeddedCount > 0)
            stats.embeddedCount.fetch_add(local.pendingEmbeddedCount,
                                          std::memory_order_relaxed);
        if (local.pendingFoundCount > 0)
            stats.foundCount.fetch_add(local.pendingFoundCount,
                                       std::memory_order_relaxed);
        threadHook.onFlush();
    };

    auto snapshotStats = [&](long long foundCount, long long embeddedCount,
                             long long satisfyingCount, long long faceSum) {
        SearchStats result;
        result.elapsed = std::chrono::steady_clock::now() - searchStart;
        result.rootsCompleted = stats.rootsCompleted.load(std::memory_order_relaxed);
        result.totalRoots = totalRoots;
        result.foundCount = foundCount;
        result.embeddedCount = embeddedCount;
        result.satisfyingCount = satisfyingCount;
        result.satisfyingFaceSum = faceSum;
        result.largestSatisfying =
            stats.largestSatisfying.load(std::memory_order_relaxed);
        return result;
    };

    std::thread reporter([&]() {
        using namespace std::chrono_literals;
        while (!workersFinished.load(std::memory_order_relaxed)) {
            std::this_thread::sleep_for(1s);
            if (!callbacks.onProgress)
                continue;
            callbacks.onProgress(snapshotStats(
                stats.foundCount.load(std::memory_order_relaxed),
                stats.embeddedCount.load(std::memory_order_relaxed),
                stats.satisfyingCount.load(std::memory_order_relaxed) /
                    FLUSH_EVERY_BDRY,
                stats.satisfyingFaceSum.load(std::memory_order_relaxed)));
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

    // Fired only after reporter.join(), not before: the reporter thread may
    // already be past its own workersFinished check (asleep mid-tick) when
    // that flag flips, so it can still fire one more onProgress afterward --
    // firing this any earlier meant a caller redrawing in place on each
    // onProgress call would immediately overwrite it.
    if (stopRequested.load(std::memory_order_relaxed) && callbacks.onInterrupted)
        callbacks.onInterrupted();

    if (aux.joinable()) {
        aux.join();
        auxHooks.afterJoin();
    }

    long long total = seedSubgraphCount;
    long long totalFound = seedFoundCount;
    long long totalEmbedded = seedEmbeddedCount;
    long long totalFaceSum = seedFaceSum;
    for (const WorkerStats &local : perThreadStats) {
        total += local.satisfyingCount;
        totalFound += local.foundCount;
        totalEmbedded += local.embeddedCount;
        totalFaceSum += local.satisfyingFaceSum;
    }

    SearchStats finalStats =
        snapshotStats(totalFound, totalEmbedded, total, totalFaceSum);
    if (callbacks.onSearchComplete)
        callbacks.onSearchComplete(finalStats);

    assert(stats.foundCount.load() == totalFound);
    assert(stats.embeddedCount.load() == totalEmbedded);
    assert(stats.satisfyingCount.load() == total);
    assert(stats.satisfyingFaceSum.load() == totalFaceSum);

    return finalStats;
}

template <int dim, int subdim>
SearchStats EmbeddingSearch<dim, subdim>::search(const unsigned numThreads,
                                                 BoundaryCondition cond,
                                                 const SearchCallbacks &callbacks) {
    struct NoopThreadHook {
        void onFound(EmbeddedSubmanifold<dim, subdim> &,
                    const std::vector<int> &, long long) {}
        void onFlush() {}
    };
    struct NoopAuxHooks {
        std::thread spawn(std::atomic<bool> &) { return {}; }
        void afterJoin() {}
    };
    return runSearch_(
        numThreads, cond,
        [this] { return EmbeddedSubmanifold<dim, subdim>(skeleton_); },
        [] { return NoopThreadHook{}; },
        [](const std::vector<int> &) {}, callbacks, NoopAuxHooks{});
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

SurfaceSearch::SurfaceSearch(const regina::Triangulation<4> &tri,
                             const std::vector<int> &seedFaces)
    : EmbeddingSearch<4, 2>(tri, seedFaces) {
    // Additional validation beyond the base class's own (facet-level-only)
    // seeded-constructor check above: KnottedSurface's seeded constructor
    // runs the enhanced addFace()/addFaces() (transverse-self-intersection/
    // local-flatness checks included), throwing regina::InvalidArgument if
    // any face fails them -- this temporary exists purely for that side
    // effect.
    KnottedSurface(skeleton_, petalCache_, seedFaces);
}

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

BoundaryCondition SurfaceSearch::classifyByLinks_(
    const std::vector<std::pair<size_t, Link>> &links) {
    if (links.empty())
        return BoundaryCondition::closed;
    bool connected =
        std::all_of(links.begin(), links.end(), [](const auto &entry) {
            return entry.second.comps_.size() <= 1;
        });
    return connected ? BoundaryCondition::connected : BoundaryCondition::proper;
}

BoundaryCondition
SurfaceSearch::classifyCheaply_(const EmbeddedSubmanifold<4, 2> &embedding) {
    if (embedding.isClosed())
        return BoundaryCondition::closed;
    if (embedding.isProper())
        return BoundaryCondition::proper;
    return BoundaryCondition::all;
}

void SurfaceSearch::ThreadHook::onFound(EmbeddedSubmanifold<4, 2> &embedding,
                                        const std::vector<int> &U,
                                        long long faceCount) {
    auto type = KnottedSurface::surfaceTypeKey(embedding.triangulation());
    ++localTypeCounts_[type];
    if (wantLinks_) {
        std::vector<int> faceIndices;
        for (int v : U)
            for (int f : owner_.graph_.graphToSkel[v - 1])
                faceIndices.push_back(f);
        localPending_.push_back(std::move(faceIndices));
    } else if (callbacks_.onSurfaceFound) {
        auto [orientable, genus, punctures] = type;
        callbacks_.onSurfaceFound(SurfaceFoundInfo{
            .orientable = orientable,
            .genus = genus,
            .punctures = punctures,
            .triangleCount = faceCount,
            .mostRestrictive = classifyCheaply_(embedding)});
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
        owner_.backgroundDrainLoop_(workersFinished, callbacks_);
    });
}

void SurfaceSearch::AuxHooks::afterJoin() {
    owner_.processRemainingSurfaceBoundaries(numThreads_, callbacks_);
}

void SurfaceSearch::backgroundDrainLoop_(
    const std::atomic<bool> &workersFinished,
    const SurfaceSearchCallbacks &callbacks) {
    using namespace std::chrono_literals;
    // Matches processBatchParallel_'s own CHUNK size -- not load-bearing
    // that they're equal, just a reasonable shared "small batch" constant.
    constexpr size_t POP_BATCH = 64;
    KnottedSurface embedding(skeleton_, petalCache_);
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
            processEntry_(embedding, items[i], callbacks);
        }
    }
}

void SurfaceSearch::processRemainingSurfaceBoundaries(
    unsigned numThreads, const SurfaceSearchCallbacks &callbacks) {
    processBatchParallel_(pendingSurfaces_.drain(), numThreads, callbacks);
}

void SurfaceSearch::processBatchParallel_(
    std::vector<std::vector<int>> batch, unsigned numThreads,
    const SurfaceSearchCallbacks &callbacks) {
    if (batch.empty())
        return;

    const auto phaseStart = std::chrono::steady_clock::now();
    const size_t total = batch.size();
    const unsigned workerCount =
        static_cast<unsigned>(std::min<size_t>(numThreads, total));

    if (callbacks.onBoundaryProcessingStarted)
        callbacks.onBoundaryProcessingStarted(total, workerCount);

    constexpr size_t CHUNK = 64;
    std::atomic<size_t> nextIndex{0};
    std::atomic<size_t> processedCount{0};
    std::atomic<bool> done{false};

    std::thread reporter([&]() {
        using namespace std::chrono_literals;
        while (!done.load(std::memory_order_relaxed)) {
            std::this_thread::sleep_for(1s);
            if (!callbacks.onBoundaryProcessingProgress)
                continue;
            callbacks.onBoundaryProcessingProgress(
                processedCount.load(std::memory_order_relaxed), total,
                std::chrono::steady_clock::now() - phaseStart);
        }
    });

    auto worker = [&]() {
        KnottedSurface embedding(skeleton_, petalCache_);
        while (true) {
            size_t begin =
                nextIndex.fetch_add(CHUNK, std::memory_order_relaxed);
            if (begin >= total)
                break;
            size_t end = std::min(begin + CHUNK, total);
            processBatchRange_(embedding, batch, begin, end, callbacks);
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
    if (callbacks.onBoundaryProcessingComplete)
        callbacks.onBoundaryProcessingComplete(
            std::chrono::steady_clock::now() - phaseStart);
}

void SurfaceSearch::processEntry_(KnottedSurface &embedding,
                                  const std::vector<int> &faceIndices,
                                  const SurfaceSearchCallbacks &callbacks) {
    for (int idx : faceIndices)
        embedding.addFace(idx);

    SurfaceTypeKey type = embedding.surfaceType();
    auto links = embedding.boundaryLinks();
    std::string descriptor;
    if (!links.empty()) {
        descriptor = describeBoundary_(links);
        linkTally_.record(descriptor, type);
    }

    if (callbacks.onSurfaceBoundaryProcessed) {
        auto [orientable, genus, punctures] = type;
        callbacks.onSurfaceBoundaryProcessed(SurfaceBoundaryInfo{
            SurfaceFoundInfo{
                .orientable = orientable,
                .genus = genus,
                .punctures = punctures,
                .triangleCount = static_cast<long long>(faceIndices.size()),
                .mostRestrictive = classifyByLinks_(links)},
            descriptor});
    }

    // Reverse order, mirroring how the DFS itself would back out --
    // resets embedding to empty for the next entry.
    for (auto it = faceIndices.rbegin(); it != faceIndices.rend(); ++it)
        embedding.removeFace(*it);
}

void SurfaceSearch::processBatchRange_(
    KnottedSurface &embedding, const std::vector<std::vector<int>> &batch,
    size_t begin, size_t end, const SurfaceSearchCallbacks &callbacks) {
    for (size_t i = begin; i < end; ++i)
        processEntry_(embedding, batch[i], callbacks);
}

SearchStats SurfaceSearch::search(unsigned numThreads, BoundaryCondition cond,
                                  const SurfaceSearchCallbacks &callbacks) {
    const bool wantLinks = cond == BoundaryCondition::proper ||
                           cond == BoundaryCondition::connected;

    // Thread safety, not just performance: Vertex<4>::buildLink() caches
    // its result as a plain, unsynchronized lazily-constructed pointer
    // (see engine/triangulation/dim4/vertex4.cpp) -- since every worker
    // thread below shares the same ambient Triangulation<4> (and hence the
    // same Vertex<4> objects) even though each has its own KnottedSurface,
    // two workers concurrently reaching faces around the same ambient
    // vertex could otherwise race on that cache. Forcing every vertex's
    // link to build once, single-threaded, here -- before runSearch_
    // spawns any worker thread -- means every later call (even
    // concurrent) only ever takes the already-built, read-only path.
    for (auto v : skeleton_.triangulation().vertices())
        v->buildLink();

    return runSearch_(
        numThreads, cond,
        [this] { return KnottedSurface(skeleton_, petalCache_); },
        [this, wantLinks, &callbacks] {
            return ThreadHook(*this, surfaceTypeTally_, wantLinks, callbacks);
        },
        [this, wantLinks, &callbacks](const std::vector<int> &seedFaces) {
            // One-off, not hot-path: rebuilds the seed as a KnottedSurface
            // (rather than reusing runSearch_'s generic proto embedding) to
            // get at surfaceType()/boundaryLinks(), which only KnottedSurface
            // exposes. seedFaces is added in the same order used to validate
            // it originally (see buildSeededGraph_), so this always embeds.
            KnottedSurface probe(skeleton_, petalCache_, seedFaces);
            SurfaceTypeKey type = probe.surfaceType();
            std::map<SurfaceTypeKey, long long> seedTypeCounts{{type, 1}};
            surfaceTypeTally_.merge(seedTypeCounts);
            auto [orientable, genus, punctures] = type;
            auto triangleCount = static_cast<long long>(seedFaces.size());
            if (wantLinks) {
                auto links = probe.boundaryLinks();
                std::string descriptor;
                if (!links.empty()) {
                    descriptor = describeBoundary_(links);
                    linkTally_.record(descriptor, type);
                }
                if (callbacks.onSurfaceBoundaryProcessed)
                    callbacks.onSurfaceBoundaryProcessed(SurfaceBoundaryInfo{
                        SurfaceFoundInfo{.orientable = orientable,
                                        .genus = genus,
                                        .punctures = punctures,
                                        .triangleCount = triangleCount,
                                        .mostRestrictive =
                                            classifyByLinks_(links)},
                        descriptor});
            } else if (callbacks.onSurfaceFound) {
                callbacks.onSurfaceFound(SurfaceFoundInfo{
                    .orientable = orientable,
                    .genus = genus,
                    .punctures = punctures,
                    .triangleCount = triangleCount,
                    .mostRestrictive = classifyCheaply_(probe)});
            }
        },
        callbacks, AuxHooks(*this, numThreads, wantLinks, callbacks));
}
