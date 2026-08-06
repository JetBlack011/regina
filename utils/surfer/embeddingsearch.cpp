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
#include <memory>
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
    const regina::Triangulation<dim> &tri,
    std::optional<size_t> protectedBoundaryComponent)
    : skeleton_(tri),
      graph_(buildGraph_(skeleton_, protectedBoundaryComponent)) {}

template <int dim, int subdim>
EmbeddingSearch<dim, subdim>::EmbeddingSearch(
    const regina::Triangulation<dim> &tri, const std::vector<int> &seedFaces,
    std::optional<size_t> protectedBoundaryComponent)
    : skeleton_(tri),
      graph_(buildSeededGraph_(skeleton_, seedFaces,
                               protectedBoundaryComponent)),
      isSeeded_(true) {

    EmbeddedSubmanifold<dim, subdim>(skeleton_, seedFaces);
}

template <int dim, int subdim>
template <typename EmbeddingT>
SearchStats EmbeddingSearch<dim, subdim>::runSearch_(
    unsigned numThreads, BoundaryCondition cond,
    std::function<EmbeddingT()> makeEmbedding,
    std::function<std::unique_ptr<RunSearchThreadHook<dim, subdim>>()>
        makeThreadHook,
    std::function<void(const std::vector<int> &)> onSeedFound,
    const SearchCallbacks &callbacks, RunSearchAuxHooks &auxHooks,
    unsigned iddfsIterations, long long iddfsStep,
    std::optional<long long> iddfsStart,
    std::optional<unsigned> finalThreads, bool orientableOnly) {
    const auto searchStart = std::chrono::steady_clock::now();

    // See this method's doc comment (embeddingsearch.h) for the SIGINT
    // contract. stopRequested_/pauseRequested_ are class members (not
    // locals here) so SurfaceSearch's ThreadHook::onFlush() -- which only
    // has access to *this via its owner_ reference -- can also check/set
    // them; reset to false here since a class member could otherwise carry
    // a stale value from a previous search() call on the same instance.
    // SigintScope is constructed before the seeded proto-embedding check
    // below and destroyed only once this whole function returns, so a
    // Ctrl+C during that check, the worker threads, or auxHooks.afterJoin()
    // is all handled uniformly.
    stopRequested_.store(false, std::memory_order_relaxed);
    pauseRequested_.store(false, std::memory_order_relaxed);
    SigintScope sigintScope(stopRequested_);

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
        std::optional<OrientabilityPredicate<EmbeddingT>> protoOrientOpt;
        ConditionalPredicate *protoEmbeddingPredicate = &protoPredicate;
        if (orientableOnly) {
            protoOrientOpt.emplace(protoPredicate, protoEmbedding);
            protoEmbeddingPredicate = &*protoOrientOpt;
        }
        InterruptiblePredicate protoInterruptible(
            *protoEmbeddingPredicate, stopRequested_, pauseRequested_);
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

    // Shallow-first: descending vertex id. Unseeded roots are provably
    // bounded (extend()'s `w > s` rule means root s can only draw from
    // {s+1..n}, so its subtree depth is bounded by n - s) -- this ordering
    // guarantees the cheap, fast roots are claimed and drained first
    // regardless of thread count, which is what prevents starvation.
    // Seeded roots have no such bound (the seeded branch anchors at 1, not
    // s), so this is a harmless heuristic there.
    std::sort(roots.begin(), roots.end(), [](int a, int b) { return a > b; });

    const size_t totalRootsPerPass = roots.size();

    const unsigned resolvedFinalThreads = finalThreads.value_or(numThreads);

    std::atomic<size_t> nextRootIdx{0};
    std::atomic<bool> workersFinished{false};
    // Iterative-deepening progress, for SearchStats::iddfsRound/
    // iddfsTotalRounds/iddfsCapped/iddfsCap -- set right before each
    // round's threads are spawned (see the round loop below), so the
    // reporter thread (polling once a second) always sees which round is
    // currently in flight.
    std::atomic<unsigned> currentIddfsRound{1};
    std::atomic<bool> currentIddfsCapped{false};
    std::atomic<long long> currentIddfsCap{0};
    AtomicSearchStats stats{.foundCount = seedFoundCount,
                            .embeddedCount = seedEmbeddedCount,
                            .satisfyingCount = seedSubgraphCount,
                            .largestSatisfying = seedMaxFaces,
                            .satisfyingFaceSum = seedFaceSum};
    std::vector<WorkerStats> perThreadStats(
        std::max(numThreads, resolvedFinalThreads));

    auto worker = [&](unsigned tid, std::optional<long long> capFaces,
                      long long suppressBelow) {
        auto embedding = makeEmbedding();
        EmbeddednessPredicate predicate(embedding, graph_.graphToSkel);
        std::optional<OrientabilityPredicate<EmbeddingT>> orientOpt;
        ConditionalPredicate *embeddingPredicate = &predicate;
        if (orientableOnly) {
            orientOpt.emplace(predicate, embedding);
            embeddingPredicate = &*orientOpt;
        }
        InterruptiblePredicate interruptible(*embeddingPredicate, stopRequested_,
                                            pauseRequested_);
        // Only constructed for a capped (iterative-deepening) round; see
        // DepthCappedPredicate. Kept alive for the whole round exactly like
        // `interruptible` above, since it must see the seeded case's
        // one-time seed commit (inside the enumerator constructor below)
        // as well as every root's own tryAdd/undo calls.
        std::optional<DepthCappedPredicate> cappedOpt;
        ConditionalPredicate *activePredicate = &interruptible;
        if (capFaces) {
            long long maxDepth = iddfsMaxDepth(*capFaces, isSeeded_);
            cappedOpt.emplace(interruptible, static_cast<int>(maxDepth));
            activePredicate = &*cappedOpt;
        }
        // ConnectedInducedSubgraphEnumerator holds a reference member, so
        // it isn't assignable -- construct it in place via optional rather
        // than choosing between two constructor calls via assignment.
        std::optional<ConnectedInducedSubgraphEnumerator> localEnumeratorOpt;
        if (isSeeded_)
            localEnumeratorOpt.emplace(graph_.adjList.first,
                                       graph_.adjList.second, true,
                                       *activePredicate);
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
            if (stopRequested_.load(std::memory_order_relaxed))
                break;
            size_t idx = nextRootIdx.fetch_add(1);
            if (idx >= totalRootsPerPass)
                break;
            int s = roots[idx];
            localEnumerator.enumerateFromRootFiltered(
                s,
                [&](const std::vector<int> &U) {
                    // Already fully reported in an earlier (shallower)
                    // iterative-deepening round -- see the round loop
                    // below. Suppresses everything for this U, including
                    // the raw/embedded counters, since those were already
                    // counted in that earlier round too. Counts only faces
                    // *added* by the search (matching capFaces/suppressBelow's
                    // units -- see iddfsMaxDepth), so a seeded U's one
                    // contracted seed vertex doesn't count toward it.
                    if (suppressBelow > 0) {
                        long long addedFaceCount =
                            isSeeded_
                                ? static_cast<long long>(U.size()) - 1
                                : static_cast<long long>(U.size());
                        if (addedFaceCount <= suppressBelow)
                            return;
                    }
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
                            threadHook->onFound(embedding, U, faceCount);
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
                                threadHook->onFlush();
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
                *activePredicate);
            stats.rootsCompleted.fetch_add(1, std::memory_order_relaxed);
        }
        // Flush this thread's remainder so the global count ends up exact.
        // Must reset each pending counter back to 0 after flushing it, not
        // just when the in-loop threshold trips (above): perThreadStats[tid]
        // is reused across iterative-deepening rounds (a fresh std::thread
        // per round, but the same persistent WorkerStats slot for a given
        // tid), so a nonzero remainder left un-reset here would otherwise
        // get flushed a second time by a later round's own final flush,
        // double-counting it into stats.*.
        if (local.pendingSatisfyingCount > 0) {
            stats.satisfyingCount.fetch_add(local.pendingSatisfyingCount,
                                            std::memory_order_relaxed);
            stats.satisfyingFaceSum.fetch_add(local.pendingFaceSum,
                                              std::memory_order_relaxed);
            local.pendingSatisfyingCount = 0;
            local.pendingFaceSum = 0;
        }
        if (local.pendingEmbeddedCount > 0) {
            stats.embeddedCount.fetch_add(local.pendingEmbeddedCount,
                                          std::memory_order_relaxed);
            local.pendingEmbeddedCount = 0;
        }
        if (local.pendingFoundCount > 0) {
            stats.foundCount.fetch_add(local.pendingFoundCount,
                                       std::memory_order_relaxed);
            local.pendingFoundCount = 0;
        }
        threadHook->onFlush();
    };

    auto snapshotStats = [&](long long foundCount, long long embeddedCount,
                             long long satisfyingCount, long long faceSum) {
        SearchStats result;
        result.elapsed = std::chrono::steady_clock::now() - searchStart;
        result.rootsCompleted = stats.rootsCompleted.load(std::memory_order_relaxed);
        // Every iterative-deepening round (each capped round, plus the
        // final unbounded one) claims every root once, so the denominator
        // scales accordingly -- see the round loop below. With
        // iddfsIterations == 0 (the default), this is just
        // totalRootsPerPass, unchanged from before this feature existed.
        result.totalRoots = totalRootsPerPass * (iddfsIterations + 1);
        result.foundCount = foundCount;
        result.embeddedCount = embeddedCount;
        result.satisfyingCount = satisfyingCount;
        result.satisfyingFaceSum = faceSum;
        result.largestSatisfying =
            stats.largestSatisfying.load(std::memory_order_relaxed);
        result.iddfsRound = currentIddfsRound.load(std::memory_order_relaxed);
        result.iddfsTotalRounds = iddfsIterations + 1;
        result.iddfsCapped =
            currentIddfsCapped.load(std::memory_order_relaxed);
        result.iddfsCap = currentIddfsCap.load(std::memory_order_relaxed);
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

    // Iterative deepening: iddfsIterations capped passes (pass i's cap =
    // resolvedIddfsStart + (i - 1) * iddfsStep faces), each guaranteed
    // fast per-root since every root's subtree is bounded -- this is what
    // guarantees shallow results are found promptly regardless of how
    // deep/unbounded some roots' true subtrees are -- followed by one
    // final, fully unbounded pass that only reports what the capped
    // passes hadn't already reached. With iddfsIterations == 0 (the
    // default), the loop below never runs and this degenerates to exactly
    // one unbounded pass over every root, identical to this function's
    // behavior before this feature existed.
    const long long resolvedIddfsStart = iddfsStart.value_or(iddfsStep);
    long long prevCap = 0;
    for (unsigned iter = 1; iter <= iddfsIterations; ++iter) {
        long long cap = iddfsCapForRound(iter, resolvedIddfsStart, iddfsStep);
        currentIddfsRound.store(iter, std::memory_order_relaxed);
        currentIddfsCapped.store(true, std::memory_order_relaxed);
        currentIddfsCap.store(cap, std::memory_order_relaxed);
        nextRootIdx.store(0, std::memory_order_relaxed);
        std::vector<std::thread> roundThreads;
        roundThreads.reserve(numThreads);
        for (unsigned t = 0; t < numThreads; ++t)
            roundThreads.emplace_back(worker, t, cap, prevCap);
        for (auto &th : roundThreads)
            th.join();
        if (stopRequested_.load(std::memory_order_relaxed))
            break;
        prevCap = cap;
    }

    currentIddfsRound.store(iddfsIterations + 1, std::memory_order_relaxed);
    currentIddfsCapped.store(false, std::memory_order_relaxed);
    currentIddfsCap.store(0, std::memory_order_relaxed);
    nextRootIdx.store(0, std::memory_order_relaxed);
    std::vector<std::thread> threads;
    threads.reserve(resolvedFinalThreads);
    for (unsigned t = 0; t < resolvedFinalThreads; ++t)
        threads.emplace_back(worker, t, std::nullopt, prevCap);
    for (auto &th : threads)
        th.join();

    workersFinished.store(true, std::memory_order_relaxed);
    reporter.join();

    // Fired only after reporter.join(), not before: the reporter thread may
    // already be past its own workersFinished check (asleep mid-tick) when
    // that flag flips, so it can still fire one more onProgress afterward --
    // firing this any earlier meant a caller redrawing in place on each
    // onProgress call would immediately overwrite it.
    if (stopRequested_.load(std::memory_order_relaxed) && callbacks.onInterrupted)
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
SearchStats EmbeddingSearch<dim, subdim>::search(
    const unsigned numThreads, BoundaryCondition cond,
    const SearchCallbacks &callbacks, unsigned iddfsIterations,
    long long iddfsStep, std::optional<long long> iddfsStart,
    std::optional<unsigned> finalThreads, bool orientableOnly) {
    struct NoopThreadHook : RunSearchThreadHook<dim, subdim> {
        void onFound(EmbeddedSubmanifold<dim, subdim> &,
                    const std::vector<int> &, long long) override {}
        void onFlush() override {}
    };
    struct NoopAuxHooks : RunSearchAuxHooks {
        std::thread spawn(std::atomic<bool> &) override { return {}; }
        void afterJoin() override {}
    };
    NoopAuxHooks noopAuxHooks;
    return runSearch_<EmbeddedSubmanifold<dim, subdim>>(
        numThreads, cond,
        [this] { return EmbeddedSubmanifold<dim, subdim>(skeleton_); },
        [] { return std::make_unique<NoopThreadHook>(); },
        [](const std::vector<int> &) {}, callbacks, noopAuxHooks,
        iddfsIterations, iddfsStep, iddfsStart, finalThreads, orientableOnly);
}

template <int dim, int subdim>
typename EmbeddingSearch<dim, subdim>::Graph
EmbeddingSearch<dim, subdim>::buildGraph_(
    const Skeleton<dim, subdim> &skeleton,
    std::optional<size_t> protectedBoundaryComponent,
    const std::vector<int> &exemptSkeletonIndices) {
    const auto &nodes = skeleton.getNodes();

    std::set<int> exempt(exemptSkeletonIndices.begin(),
                         exemptSkeletonIndices.end());

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

        // Reject any face with an edge on protectedBoundaryComponent,
        // unless it's explicitly exempted (a seed face, which legitimately
        // touches that boundary component by construction -- see
        // buildSeededGraph_()'s doc comment). This is edge-level, not
        // face-level: a face with just one edge on the protected
        // component, while the face itself is interior, must be rejected
        // too -- including it could turn one of that boundary's own edges
        // into an interior edge of the constructed surface, silently
        // changing which edges end up exposed as the surface's own
        // boundary. A face that's itself entirely one of the protected
        // component's own 2D faces trivially has all of its edges on that
        // boundary too, so this one check subsumes both cases.
        if (protectedBoundaryComponent && !exempt.contains(static_cast<int>(i))) {
            bool touchesProtected = false;
            for (int e = 0; e < 3; ++e) {
                auto *bc = nodes[i].face->edge(e)->boundaryComponent();
                if (bc && bc->index() == *protectedBoundaryComponent) {
                    touchesProtected = true;
                    break;
                }
            }
            if (touchesProtected)
                continue;
        }

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
    const Skeleton<dim, subdim> &skeleton, const std::vector<int> &seedFaces,
    std::optional<size_t> protectedBoundaryComponent) {
    // seedFaces is forwarded as the exemption list -- see this function's
    // own doc comment in embeddingsearch.h for why the seed must never be
    // excluded by protectedBoundaryComponent even though it legitimately
    // touches it.
    Graph base =
        buildGraph_(skeleton, protectedBoundaryComponent, seedFaces);

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

// Explicit instantiations of runSearch_(), one per concrete EmbeddingT
// actually used -- see the extern template declarations in
// embeddingsearch.h (and runSearch_()'s own doc comment there) for why
// this works despite runSearch_ being a member template: EmbeddingT is a
// small, fixed, nameable set of concrete types (unlike the lambda-closure
// types the other hook parameters used to be), so ordinary explicit
// instantiation applies. KnottedSurface's full definition is visible here
// via embeddedsubmanifold.h (included transitively through
// embeddingsearch.h), so this instantiation -- used only by
// SurfaceSearch::search(), in surfacesearch.cpp -- can live in this file
// rather than surfacesearch.cpp itself.
template SearchStats
EmbeddingSearch<3, 2>::runSearch_<EmbeddedSubmanifold<3, 2>>(
    unsigned, BoundaryCondition, std::function<EmbeddedSubmanifold<3, 2>()>,
    std::function<std::unique_ptr<RunSearchThreadHook<3, 2>>()>,
    std::function<void(const std::vector<int> &)>, const SearchCallbacks &,
    RunSearchAuxHooks &, unsigned, long long, std::optional<long long>,
    std::optional<unsigned>, bool);
template SearchStats
EmbeddingSearch<4, 2>::runSearch_<EmbeddedSubmanifold<4, 2>>(
    unsigned, BoundaryCondition, std::function<EmbeddedSubmanifold<4, 2>()>,
    std::function<std::unique_ptr<RunSearchThreadHook<4, 2>>()>,
    std::function<void(const std::vector<int> &)>, const SearchCallbacks &,
    RunSearchAuxHooks &, unsigned, long long, std::optional<long long>,
    std::optional<unsigned>, bool);
template SearchStats
EmbeddingSearch<4, 2>::runSearch_<KnottedSurface>(
    unsigned, BoundaryCondition, std::function<KnottedSurface()>,
    std::function<std::unique_ptr<RunSearchThreadHook<4, 2>>()>,
    std::function<void(const std::vector<int> &)>, const SearchCallbacks &,
    RunSearchAuxHooks &, unsigned, long long, std::optional<long long>,
    std::optional<unsigned>, bool);

