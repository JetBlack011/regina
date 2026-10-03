//
//  embeddingsearch.cpp
//
//  Created by John Teague on 07/15/2026.
//

#include "surfer/enumeration/submanifoldsearch.h"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <deque>
#include <limits>
#include <cassert>
#include <csignal>
#include <cstdlib>
#include <iomanip>
#include <memory>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <type_traits>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include "surfer/pairsig/sha1.h"

#define FLUSH_EVERY_BDRY 1
#define FLUSH_EVERY_EMBEDDED 1
#define FLUSH_EVERY_FOUND 10'000

namespace {
// Ctrl+C handling for runSearch_(). g_stopRequested/g_sigintCount are
// process-wide, not per-search.
std::atomic<bool> *g_stopRequested = nullptr;
std::atomic<int> g_sigintCount{0};

void handleSigint(int) {
    if (g_sigintCount.fetch_add(1, std::memory_order_relaxed) > 0)
        std::_Exit(130); // 128 + SIGINT
    if (g_stopRequested)
        g_stopRequested->store(true, std::memory_order_relaxed);
}

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
// SearchStats.
struct AtomicSearchStats {
    std::atomic<long long> foundCount;      // raw candidates visited, any cond
    std::atomic<long long> embeddedCount;   // isEmbedded()==true, any cond
    std::atomic<long long> satisfyingCount; // satisfying the BoundaryCondition
    std::atomic<long long>
        largestSatisfying; // max faces among satisfying finds
    std::atomic<long long> satisfyingFaceSum;
    std::atomic<size_t> rootsCompleted{0};
    std::atomic<long long> resolvedCount{0}; // satisfying, but not embedded
};

// One worker thread's own running totals, plus the portion of each not yet
// folded into the shared AtomicSearchStats.
struct WorkerStats {
    long long foundCount = 0;
    long long embeddedCount = 0;
    long long satisfyingCount = 0;
    long long satisfyingFaceSum = 0;
    long long pendingFoundCount = 0;
    long long pendingEmbeddedCount = 0;
    long long pendingSatisfyingCount = 0;
    long long pendingFaceSum = 0;
    long long pendingResolvedCount = 0;
    // For SearchStats::Profile; see its fields.
    long long nodes = 0;
    long long attempts = 0;
    long long evaluated = 0;
    long long charged = 0;
    long long replayed = 0;
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
    : skeleton_(tri), graph_(buildSeededGraph_(skeleton_, seedFaces,
                                               protectedBoundaryComponent)),
      isSeeded_(true) {

    EmbeddedSubmanifold<dim, subdim>(skeleton_, seedFaces);
}

template <int dim, int subdim>
std::string EmbeddingSearch<dim, subdim>::frontierFingerprint_(
    const std::vector<int> &roots, BoundaryCondition cond,
    unsigned iddfsIterations, long long iddfsStart, long long iddfsStep,
    std::optional<long long> hardFaceCap, bool orientableOnly,
    long long rootBudgetStart, long long rootBudgetGrowth) const {
    std::ostringstream s;
    s << "traversal " << SearchFrontier::kTraversalVersion << " dim " << dim
      << " subdim " << subdim << "\ntri";
    // The labelled triangulation: the graph's faces and the predicate's
    // checks are both read off it.
    const regina::Triangulation<dim> &tri = skeleton_.triangulation();
    s << ' ' << tri.size();
    for (size_t i = 0; i < tri.size(); ++i) {
        const auto *simplex = tri.simplex(i);
        for (int f = 0; f <= dim; ++f) {
            const auto *adj = simplex->adjacentSimplex(f);
            s << ' ';
            if (adj)
                s << adj->index() << ':' << simplex->adjacentGluing(f).str();
            else
                s << '-';
        }
    }
    s << "\ngraph " << graph_.adjList.first;
    for (const auto &neighbours : graph_.adjList.second) {
        s << " |";
        for (int v : neighbours)
            s << ' ' << v;
    }
    // Each vertex's faces as a set: the seed's arrive in no fixed order
    // (measured: two hops on one row), and it is committed whole, so only
    // which faces it holds can matter to the search.
    s << "\nskeleton";
    for (const auto &faces : graph_.graphToSkel) {
        std::vector<int> sorted(faces);
        std::sort(sorted.begin(), sorted.end());
        s << " |";
        for (int f : sorted)
            s << ' ' << f;
    }
    s << "\nseeded " << (isSeeded_ ? 1 : 0) << "\nroots";
    for (int r : roots)
        s << ' ' << r;
    s << "\niddfs " << iddfsIterations << ' ' << iddfsStart << ' ' << iddfsStep
      << " hard " << (hardFaceCap ? std::to_string(*hardFaceCap) : "-")
      << "\nbudget " << rootBudgetStart << ' ' << rootBudgetGrowth
      << "\ncondition " << boundaryConditionName(cond) << " orientable "
      << (orientableOnly ? 1 : 0) << "\ncontext " << frontierContext_ << '\n';
    return pairsig::sha1Hex(s.str());
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
    std::optional<long long> iddfsStart, std::optional<unsigned> finalThreads,
    bool orientableOnly, std::optional<long long> hardFaceCap,
    long long rootBudgetStart, long long rootBudgetGrowth) {
    const auto searchStart = std::chrono::steady_clock::now();

    stopRequested_.store(false, std::memory_order_relaxed);
    pauseRequested_.store(false, std::memory_order_relaxed);
    std::optional<SigintScope> sigintScope;
    if (handleSigint_)
        sigintScope.emplace(stopRequested_);

    // Shared dynamic work queue over roots. Unseeded: every graph vertex,
    // unconditionally (matches the old s = 1..n sweep exactly). Seeded:
    // every sibling of the seed surviving the predicate, via one prototype
    // embedding/enumerator built once here on the calling thread. 
    std::vector<int> roots;
    long long seedFoundCount = 0;
    long long seedEmbeddedCount = 0;
    long long seedSubgraphCount = 0;
    long long seedMaxFaces = 0;
    long long seedFaceSum = 0;
    long long seedFaceCount = 0;
    long long seedResolvedCount = 0;
    bool seedSatisfies = false;
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
        seedFaceCount =
            static_cast<long long>(protoEmbedding.triangulation().size());

        seedFoundCount = 1;
        if (protoEmbedding.isEmbedded())
            seedEmbeddedCount = 1;
        // Same acceptance rule as the workers' visit callback below.
        if (protoEmbedding.satisfies(cond) && protoEmbedding.isAcceptable()) {
            seedSubgraphCount = 1;
            if (!protoEmbedding.isEmbedded())
                seedResolvedCount = 1;
            seedMaxFaces = static_cast<long long>(
                protoEmbedding.triangulation().size());
            seedFaceSum = seedMaxFaces;
            seedSatisfies = true; // reported below, unless resuming
        }
    } else {
        roots.resize(graph_.adjList.first);
        std::iota(roots.begin(), roots.end(), 1);
    }

    // Shallow-first: descending vertex id.
    std::sort(roots.begin(), roots.end(), [](int a, int b) { return a > b; });

    const unsigned resolvedFinalThreads = finalThreads.value_or(numThreads);
    const long long resolvedIddfsStart = iddfsStart.value_or(iddfsStep);

    // The frontier this search carries on from, if it was given one that is
    // its own, and the one it records; see SearchFrontier.
    frontier_.reset();
    resumed_ = false;
    resumeRefusal_.clear();
    frontierSeconds_ = 0;
    auto frontierClock = [this](std::chrono::steady_clock::time_point since) {
        frontierSeconds_ += std::chrono::duration<double>(
                                std::chrono::steady_clock::now() - since)
                                .count();
    };
    std::string fingerprint;
    if (recordFrontier_ || resumeFrontier_) {
        const auto t = std::chrono::steady_clock::now();
        fingerprint = frontierFingerprint_(
            roots, cond, iddfsIterations, resolvedIddfsStart, iddfsStep,
            hardFaceCap, orientableOnly, rootBudgetStart, rootBudgetGrowth);
        frontierClock(t);
    }
    const SearchFrontier *resume = nullptr;
    if (resumeFrontier_) {
        if (resumeFrontier_->fingerprint != fingerprint)
            resumeRefusal_ = "its fingerprint is not this search's";
        else if (resumeFrontier_->roots.size() != roots.size() ||
                 resumeFrontier_->rounds != iddfsIterations + 1 ||
                 resumeFrontier_->round < 1 ||
                 resumeFrontier_->round > iddfsIterations + 1)
            resumeRefusal_ = "its roots or rounds are not this search's";
        else {
            resume = resumeFrontier_;
            resumed_ = true;
        }
    }

    // The seed alone: reported by the run a frontier started with, so a
    // resumed run neither reports nor counts it again.
    if (seedSatisfies && resume) {
        seedFoundCount = seedEmbeddedCount = seedSubgraphCount = 0;
        seedMaxFaces = seedFaceSum = seedResolvedCount = 0;
    } else if (seedSatisfies) {
        // graph_.graphToSkel[0] is the whole seed -- see buildSeededGraph_.
        onSeedFound(graph_.graphToSkel[0]);
    } else if (resume) {
        seedFoundCount = seedEmbeddedCount = 0;
    }

    const auto prototypeTime = std::chrono::steady_clock::now() - searchStart;
    if (callbacks.onRootsReady)
        callbacks.onRootsReady();

    const size_t totalRootsPerPass = roots.size();

    std::atomic<size_t> nextRootIdx{0};
    std::atomic<bool> workersFinished{false};
    std::atomic<bool> surfaceTargetReached{false}; // see SearchCallbacks
    // Iterative-deepening progress, for SearchStats::iddfsRound/
    // iddfsTotalRounds/iddfsCapped/iddfsCap.
    std::atomic<unsigned> currentIddfsRound{1};
    std::atomic<bool> currentIddfsCapped{false};
    std::atomic<long long> currentIddfsCap{0};
    // Per-root budget currently in force (0 = unbudgeted), for progress.
    std::atomic<long long> currentRootBudget{0};
    // The largest face cap whose round finished with EVERY root enumerated
    // to completion. -1 means no round did, so nothing exhaustive can be
    // claimed. This is what turns "we found nothing" into "nothing exists
    // with at most this many added faces".
    std::atomic<long long> deepestExhausted{-1};
    AtomicSearchStats stats{.foundCount = seedFoundCount,
                            .embeddedCount = seedEmbeddedCount,
                            .satisfyingCount = seedSubgraphCount,
                            .largestSatisfying = seedMaxFaces,
                            .satisfyingFaceSum = seedFaceSum,
                            .resolvedCount = seedResolvedCount};
    std::vector<WorkerStats> perThreadStats(
        std::max(numThreads, resolvedFinalThreads));

    // The surface target is the search's breadth (see
    // SearchCallbacks::surfaceTarget): a resumed search has already covered
    // its frontier's surfaces, and one already that broad stops at once.
    long long surfaceTarget = callbacks.surfaceTarget;
    if (surfaceTarget > 0 && resume) {
        surfaceTarget -= resume->satisfying;
        if (surfaceTarget <= 0) {
            surfaceTarget = 0;
            surfaceTargetReached.store(true);
            if (callbacks.onSurfaceTarget)
                callbacks.onSurfaceTarget();
            requestStop();
        }
    }

    // Per-root scheduling state, shared across budget passes within one
    // depth round. Plain values rather than atomics: within a pass each root
    // index is claimed by exactly one thread (via nextRootIdx), and passes
    // are separated by thread joins, which supply the happens-before edge.
    // uint8_t rather than bool because std::vector<bool> packs bits, making
    // concurrent writes to distinct indices unsafe.
    std::vector<uint8_t> rootDone(roots.size(), 0);
    // positions[i]: where root i's previous budget pass stopped. Its next
    // pass carries on from exactly there, on whichever thread takes it (see
    // ConnectedInducedSubgraphEnumerator::Position), so the passes together
    // walk the root once: nothing is retraced and nothing reported twice.
    // Until 2026-09-28 each pass instead re-walked its root from the start
    // and skipped the visits already reported, which was about half of all
    // the search's attempts.
    std::vector<ConnectedInducedSubgraphEnumerator::Position> positions(
        roots.size());
    // rootLevel[i]: how many budget passes root i has had. Its ration is
    // rootBudgetStart * growth^rootLevel[i]. Levels are sequential PER ROOT
    // (pass k+1 needs pass k's visit count to know what to skip) but wholly
    // independent ACROSS roots -- which is what lets the schedule below run
    // without a barrier.
    std::vector<unsigned> rootLevel(roots.size(), 0);

    // Barrier-free work queue over roots. A worker pops a root, runs it at
    // its own next level, and pushes it back if it still has work; finished
    // roots simply drop out.
    //
    // The previous design ran one level at a time and joined between them.
    // That starved the machine badly: root costs here vary by orders of
    // magnitude, so all roots would be claimed early and 10 of 12 threads sat
    // at the join while two ground out the expensive stragglers -- measured
    // at ~40% utilisation, with candidate throughput dropping to zero for
    // stretches. Recycling roots through a queue keeps every thread busy for
    // as long as ANY root has work left.
    // Ration for a root at pass `level`: start * growth^level, saturating
    // (a budget that large is unlimited in practice anyway). A
    // rootBudgetStart of <= 0 means unbudgeted, signalled as -1 so
    // BudgetedPredicate stays transparent.
    auto budgetForLevel = [&](unsigned level) -> long long {
        if (rootBudgetStart <= 0)
            return -1;
        long long b = rootBudgetStart;
        for (unsigned i = 0; i < level; ++i) {
            if (b > std::numeric_limits<long long>::max() / rootBudgetGrowth)
                return std::numeric_limits<long long>::max();
            b *= rootBudgetGrowth;
        }
        return b;
    };

    std::mutex queueMutex;
    std::condition_variable queueCv;
    std::deque<size_t> rootQueue;
    size_t rootsInFlight = 0;

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
        // Before the predicate, which calls it while a pause holds this
        // worker (RunSearchThreadHook::onPaused()).
        auto threadHook = makeThreadHook();
        InterruptiblePredicate interruptible(
            *embeddingPredicate, stopRequested_, pauseRequested_,
            [&threadHook] { return threadHook->onPaused(); });
        // Outermost, so its counter sees every attempt. The depth cap is the
        // enumerator's own (setMaxSize() below), so every attempt it charges
        // reaches the embedding checks -- the cost it exists to ration.
        // Constructed even when unbudgeted, where it is transparent -- that
        // keeps one code path for both modes.
        //
        // Constructed UNLIMITED (-1) whatever the per-root ration will be:
        // a seeded enumerator commits the seed via one tryAdd() during
        // construction below, before any root is reached, and that commit
        // must succeed or the enumerator's invariants break (it later
        // removes faces it believes it added). The real ration is installed
        // per root by reset().
        BudgetedPredicate budgeted(interruptible, -1);
        ConditionalPredicate *activePredicate = &budgeted;
        std::optional<ConnectedInducedSubgraphEnumerator> localEnumeratorOpt;
        if (isSeeded_)
            localEnumeratorOpt.emplace(graph_.adjList.first,
                                       graph_.adjList.second, true,
                                       *activePredicate);
        else
            localEnumeratorOpt.emplace(graph_.adjList.first,
                                       graph_.adjList.second);
        auto &localEnumerator = *localEnumeratorOpt;
        // Recording a frontier: a stop suspends each root where it stands.
        // Only now, after the seed's commit and root probes above, which
        // must never be refused by it.
        if (recordFrontier_)
            budgeted.suspendOnStop(&stopRequested_);
        // At least 1, so a seeded search can always commit its seed (as
        // DepthCappedPredicate's own clamp did).
        localEnumerator.setMaxSize(
            capFaces ? static_cast<size_t>(std::max<long long>(
                           1, iddfsMaxDepth(*capFaces, isSeeded_)))
                     : 0);
        WorkerStats &local = perThreadStats[tid];
        // The seeded constructor's own seed commit and root probes are not
        // root enumeration; the profile counts from here on.
        const long long attemptsAtStart = budgeted.attempts();

        // Dispatch is safe only because resetCandidateOrder() below makes a
        // root's traversal a function of that root alone, not of which roots
        // this thread handled first.
        while (true) {
            size_t idx;
            long long rootBudget;
            {
                std::unique_lock<std::mutex> lock(queueMutex);
                queueCv.wait(lock, [&] {
                    return !rootQueue.empty() || rootsInFlight == 0 ||
                           stopRequested_.load(std::memory_order_relaxed);
                });
                if (stopRequested_.load(std::memory_order_relaxed))
                    break;
                if (rootQueue.empty()) {
                    if (rootsInFlight == 0)
                        break; // every root finished; the round is done
                    continue;  // someone may yet push a root back
                }
                idx = rootQueue.front();
                rootQueue.pop_front();
                ++rootsInFlight;
                rootBudget = budgetForLevel(rootLevel[idx]);
                currentRootBudget.store(rootBudget, std::memory_order_relaxed);
            }
            int s = roots[idx];
            // A finished root never comes back (it leaves the queue below),
            // and the enumerator leaves its candidate list exactly as it found
            // it, so a root's traversal is a function of that root alone.
            // resetCandidateOrder() is kept as a cheap guarantee of that.
            localEnumerator.resetCandidateOrder();
            budgeted.reset(rootBudget <= 0 ? -1 : rootBudget);
            // Carry on from where this root's previous pass stopped; the
            // re-adds back down to there were paid for then.
            ConnectedInducedSubgraphEnumerator::Position &position =
                positions[idx];
            const long long rebuild = position.rebuildAttempts();
            budgeted.freeAttempts(rebuild);
            long long visits = 0;
            const auto outcome = localEnumerator.enumerateFromRootFiltered(
                s,
                [&](const std::vector<int> &U) {
                    ++visits;
                    if (suppressBelow > 0) {
                        long long addedFaceCount =
                            isSeeded_ ? static_cast<long long>(U.size()) - 1
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
                    const bool embedded = embedding.isEmbedded();
                    if (embedded) {
                        ++local.embeddedCount;
                        if (++local.pendingEmbeddedCount >=
                            FLUSH_EVERY_EMBEDDED) {
                            stats.embeddedCount.fetch_add(
                                local.pendingEmbeddedCount,
                                std::memory_order_relaxed);
                            local.pendingEmbeddedCount = 0;
                        }
                    }
                    // Acceptance. An embedded candidate is accepted when it
                    // satisfies `cond` and isAcceptable() -- for
                    // KnottedSurface, when it is also smooth at the boundary,
                    // a post hoc check because it is not hereditary. A
                    // non-embedded one is examined further only when the
                    // embedding type can say more (mayResolve(): KnottedSurface
                    // with --resolve-unlinked, or measuring). `cond` always
                    // comes first because it is O(1) and the rest is not.
                    bool accepted = false;
                    bool resolved = false;
                    if (embedded) {
                        accepted = embedding.satisfies(cond) &&
                                   embedding.isAcceptable();
                    } else if (embedding.mayResolve() &&
                               embedding.satisfies(cond)) {
                        if constexpr (std::is_same_v<EmbeddingT,
                                                     KnottedSurface>)
                            embedding.tallySelfIntersection();
                        resolved = accepted = embedding.isAcceptable();
                    }
                    if (resolved &&
                        ++local.pendingResolvedCount >= FLUSH_EVERY_EMBEDDED) {
                        stats.resolvedCount.fetch_add(
                            local.pendingResolvedCount,
                            std::memory_order_relaxed);
                        local.pendingResolvedCount = 0;
                    }
                    if (accepted) {
                        ++local.satisfyingCount;
                        auto faceCount = static_cast<long long>(
                            embedding.triangulation().size());
                        threadHook->onFound(embedding, U, faceCount);
                        local.satisfyingFaceSum += faceCount;
                        local.pendingFaceSum += faceCount;
                        if (++local.pendingSatisfyingCount >=
                            FLUSH_EVERY_BDRY) {
                            const long long total =
                                stats.satisfyingCount.fetch_add(
                                    local.pendingSatisfyingCount,
                                    std::memory_order_relaxed) +
                                local.pendingSatisfyingCount;
                            stats.satisfyingFaceSum.fetch_add(
                                local.pendingFaceSum,
                                std::memory_order_relaxed);
                            local.pendingSatisfyingCount = 0;
                            local.pendingFaceSum = 0;
                            threadHook->onFlush();
                            // The surface target, exactly; see
                            // SearchCallbacks::surfaceTarget.
                            if (surfaceTarget > 0 && total >= surfaceTarget &&
                                !surfaceTargetReached.exchange(true)) {
                                if (callbacks.onSurfaceTarget)
                                    callbacks.onSurfaceTarget();
                                requestStop();
                            }
                        }
                        auto prevMax = stats.largestSatisfying.load(
                            std::memory_order_relaxed);
                        while (
                            faceCount > prevMax &&
                            !stats.largestSatisfying.compare_exchange_weak(
                                prevMax, faceCount,
                                std::memory_order_relaxed))
                            ;
                    }
                },
                *activePredicate, &position);
            // A path that re-added fine when it was recorded can only be
            // refused again by a stop.
            if (outcome ==
                    ConnectedInducedSubgraphEnumerator::Outcome::stopped &&
                !stopRequested_.load(std::memory_order_relaxed))
                throw std::logic_error(
                    "EmbeddingSearch: a budget pass could not resume where "
                    "the previous one stopped");
            local.nodes += visits;
            local.charged += budgeted.spent();
            local.replayed += rebuild; // the uncharged re-adds, nothing more
            // A root whose pass ran to the end was enumerated exhaustively
            // (within this round's depth cap), so it needs no further pass.
            // An interrupted run must not claim this: stopRequested_ prunes
            // via InterruptiblePredicate without tripping the budget, which
            // would otherwise look like completion. Recording a frontier, a
            // stop always suspends instead (suspendOnStop()), so completion
            // is real even with a stop raised -- and a frontier must say so,
            // or its root would be walked again from scratch.
            const bool finished =
                outcome ==
                    ConnectedInducedSubgraphEnumerator::Outcome::completed &&
                (recordFrontier_ ||
                 !stopRequested_.load(std::memory_order_relaxed));
            stats.rootsCompleted.fetch_add(1, std::memory_order_relaxed);
            {
                std::lock_guard<std::mutex> lock(queueMutex);
                --rootsInFlight;
                if (finished)
                    rootDone[idx] = 1;
                else if (!stopRequested_.load(std::memory_order_relaxed)) {
                    ++rootLevel[idx]; // more ration next time round
                    rootQueue.push_back(idx);
                }
                queueCv.notify_all();
            }
        }
        const long long rootAttempts = budgeted.attempts() - attemptsAtStart;
        local.attempts += rootAttempts;
        // The cap is the enumerator's (setMaxSize()), so nothing between the
        // budget and the embedding checks refuses: every attempt is
        // evaluated.
        local.evaluated += rootAttempts;
        // Flush this thread's remainder so the global count ends up exact.
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
        if (local.pendingResolvedCount > 0) {
            stats.resolvedCount.fetch_add(local.pendingResolvedCount,
                                          std::memory_order_relaxed);
            local.pendingResolvedCount = 0;
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
        result.rootsCompleted =
            stats.rootsCompleted.load(std::memory_order_relaxed);
        // Every iterative-deepening round (each capped round, plus the
        // final unbounded one) claims every root once, so the denominator
        // scales accordingly
        result.totalRoots = totalRootsPerPass * (iddfsIterations + 1);
        result.seedFaces = seedFaceCount;
        result.rootBudget = currentRootBudget.load(std::memory_order_relaxed);
        result.rootsPerPass = totalRootsPerPass;
        // Read without synchronisation while workers may be writing: these
        // are byte-sized flags and this is a progress counter, so a torn or
        // stale read costs at most a slightly-off display.
        result.rootsExhausted =
            static_cast<size_t>(std::count(rootDone.begin(), rootDone.end(),
                                           static_cast<uint8_t>(1)));
        long long exhausted = deepestExhausted.load(std::memory_order_relaxed);
        result.deepestExhaustedCap =
            exhausted >= 0 ? std::optional<long long>(exhausted) : std::nullopt;
        result.foundCount = foundCount;
        result.embeddedCount = embeddedCount;
        result.satisfyingCount = satisfyingCount;
        result.resolvedCount =
            stats.resolvedCount.load(std::memory_order_relaxed);
        result.satisfyingFaceSum = faceSum;
        result.largestSatisfying =
            stats.largestSatisfying.load(std::memory_order_relaxed);
        result.iddfsRound = currentIddfsRound.load(std::memory_order_relaxed);
        result.iddfsTotalRounds = iddfsIterations + 1;
        result.iddfsCapped = currentIddfsCapped.load(std::memory_order_relaxed);
        result.iddfsCap = currentIddfsCap.load(std::memory_order_relaxed);
        return result;
    };

    // The reporter reports once a second, and is woken the moment the
    // workers finish (it used to sleep out its second, which a short search
    // -- a cascade hop -- paid as idle wall time on every call).
    std::mutex reporterMutex;
    std::condition_variable reporterWake;
    std::thread reporter([&]() {
        using namespace std::chrono_literals;
        while (!workersFinished.load(std::memory_order_relaxed)) {
            {
                std::unique_lock<std::mutex> lock(reporterMutex);
                reporterWake.wait_for(lock, 1s, [&] {
                    return workersFinished.load(std::memory_order_relaxed);
                });
            }
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

    // Runs one depth round to completion, as a sequence of budget passes.
    //
    // Unbudgeted (rootBudgetStart <= 0) this is a single pass over every
    // root. Budgeted, each pass gives every not-yet-finished root a ration,
    // multiplied pass by pass, and carries on where that root's previous
    // pass stopped (positions), so the roots that need more effort get it,
    // the cheap ones are not revisited (rootDone), and nothing is retraced.
    //
    // Runs one depth round to completion. Workers are spawned ONCE and drain
    // the recycling root queue, so there is no barrier between budget levels
    // and a thread that finishes a cheap root immediately takes another.
    //
    // Returns whether every root finished, i.e. whether this round enumerated
    // its depth exhaustively -- a statement about the object, not about how
    // long we ran.
    //
    // `from`, if set, is a frontier stopped in this round: its finished roots
    // stay finished, and the rest carry on at their recorded level and
    // position, exactly as its own run would have.
    std::vector<std::chrono::steady_clock::duration> roundTimes;
    auto runRound = [&](std::optional<long long> capFaces,
                        long long suppressBelow, unsigned threadCount,
                        const SearchFrontier *from) -> bool {
        const auto roundStart = std::chrono::steady_clock::now();
        std::fill(rootDone.begin(), rootDone.end(), 0);
        for (auto &position : positions)
            position = {}; // each round starts every root from scratch
        std::fill(rootLevel.begin(), rootLevel.end(), 0u);
        if (from)
            for (size_t i = 0; i < roots.size(); ++i) {
                rootDone[i] = from->roots[i].done ? 1 : 0;
                rootLevel[i] = from->roots[i].level;
                positions[i] = from->roots[i].position;
            }
        {
            std::lock_guard<std::mutex> lock(queueMutex);
            rootQueue.clear();
            for (size_t i = 0; i < roots.size(); ++i)
                if (!rootDone[i])
                    rootQueue.push_back(i);
            rootsInFlight = 0;
        }

        std::vector<std::thread> roundThreads;
        roundThreads.reserve(threadCount);
        for (unsigned t = 0; t < threadCount; ++t)
            roundThreads.emplace_back(worker, t, capFaces, suppressBelow);
        for (auto &th : roundThreads)
            th.join();
        roundTimes.push_back(std::chrono::steady_clock::now() - roundStart);

        if (stopRequested_.load(std::memory_order_relaxed))
            return false;
        return std::all_of(rootDone.begin(), rootDone.end(),
                           [](uint8_t d) { return d != 0; });
    };

    // The frontier, taken once: when a round ends on a stop (its threads
    // joined, so every root's state is final), or when every round ran.
    auto takeFrontier = [&](unsigned round, std::optional<long long> cap,
                            long long suppressBelow, bool complete) {
        if (!recordFrontier_ || frontier_)
            return;
        const auto t = std::chrono::steady_clock::now();
        SearchFrontier f;
        f.fingerprint = fingerprint;
        f.round = round;
        f.rounds = iddfsIterations + 1;
        f.cap = cap;
        f.suppressBelow = suppressBelow;
        f.complete = complete;
        f.roots.resize(roots.size());
        for (size_t i = 0; i < roots.size(); ++i) {
            SearchFrontier::Root &r = f.roots[i];
            r.done = complete || rootDone[i] != 0;
            if (!r.done) {
                r.level = rootLevel[i];
                r.position = positions[i];
            }
        }
        frontier_ = std::move(f);
        frontierClock(t);
    };

    // A resumed search starts in its frontier's round, the rounds before
    // having run to the end; a complete frontier leaves nothing to run.
    unsigned firstRound = 1;
    long long prevCap = 0;
    if (resume) {
        firstRound = resume->complete ? iddfsIterations + 2 : resume->round;
        prevCap = resume->suppressBelow;
        if (resume->deepestExhausted)
            deepestExhausted.store(*resume->deepestExhausted,
                                   std::memory_order_relaxed);
    }
    for (unsigned iter = firstRound; iter <= iddfsIterations; ++iter) {
        long long cap = iddfsCapForRound(iter, resolvedIddfsStart, iddfsStep);
        currentIddfsRound.store(iter, std::memory_order_relaxed);
        currentIddfsCapped.store(true, std::memory_order_relaxed);
        currentIddfsCap.store(cap, std::memory_order_relaxed);
        if (runRound(cap, prevCap, numThreads,
                     resume && iter == firstRound ? resume : nullptr))
            deepestExhausted.store(cap, std::memory_order_relaxed);
        if (stopRequested_.load(std::memory_order_relaxed)) {
            takeFrontier(iter, cap, prevCap, false);
            break;
        }
        prevCap = cap;
    }

    const bool finalPassRedundant = hardFaceCap && prevCap >= *hardFaceCap;
    currentIddfsRound.store(iddfsIterations + 1, std::memory_order_relaxed);
    currentIddfsCapped.store(hardFaceCap.has_value(),
                             std::memory_order_relaxed);
    currentIddfsCap.store(hardFaceCap.value_or(0), std::memory_order_relaxed);
    if (!finalPassRedundant && firstRound <= iddfsIterations + 1) {
        const SearchFrontier *from =
            resume && firstRound == iddfsIterations + 1 ? resume : nullptr;
        if (runRound(hardFaceCap, prevCap, resolvedFinalThreads, from) &&
            hardFaceCap)
            deepestExhausted.store(*hardFaceCap, std::memory_order_relaxed);
        if (stopRequested_.load(std::memory_order_relaxed))
            takeFrontier(iddfsIterations + 1, hardFaceCap, prevCap, false);
    }
    if (!stopRequested_.load(std::memory_order_relaxed))
        takeFrontier(iddfsIterations + 1, hardFaceCap, prevCap, true);

    {
        std::lock_guard<std::mutex> lock(reporterMutex);
        workersFinished.store(true, std::memory_order_relaxed);
    }
    reporterWake.notify_all();
    reporter.join();

    if (stopRequested_.load(std::memory_order_relaxed) &&
        callbacks.onInterrupted)
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
    finalStats.profile.prototype = prototypeTime;
    finalStats.profile.rounds = roundTimes;
    for (const WorkerStats &local : perThreadStats) {
        finalStats.profile.nodes += local.nodes;
        finalStats.profile.attempts += local.attempts;
        finalStats.profile.evaluated += local.evaluated;
        finalStats.profile.charged += local.charged;
        finalStats.profile.replayed += local.replayed;
    }
    if (frontier_) {
        // Cumulative over the frontier this run carried on from.
        frontier_->deepestExhausted = finalStats.deepestExhaustedCap;
        frontier_->runs = (resume ? resume->runs : 0) + 1;
        frontier_->satisfying =
            (resume ? resume->satisfying : 0) + finalStats.satisfyingCount;
        frontier_->found =
            (resume ? resume->found : 0) + finalStats.foundCount;
        frontier_->attempts =
            (resume ? resume->attempts : 0) + finalStats.profile.attempts;
    }
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
    std::optional<unsigned> finalThreads, bool orientableOnly,
    std::optional<long long> hardFaceCap, long long rootBudgetStart,
    long long rootBudgetGrowth) {
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
        iddfsIterations, iddfsStep, iddfsStart, finalThreads, orientableOnly,
        hardFaceCap, rootBudgetStart, rootBudgetGrowth);
}

template <int dim, int subdim>
size_t EmbeddingSearch<dim, subdim>::countSearchableFacesTouching(
    size_t boundaryComponent) const {
    const auto &nodes = skeleton_.getNodes();
    size_t count = 0;
    // Graph vertex 1 of a seeded graph is the seed itself; every other vertex
    // is a single searchable face.
    for (size_t v = isSeeded_ ? 1 : 0; v < graph_.graphToSkel.size(); ++v) {
        for (int f : graph_.graphToSkel[v]) {
            const auto *face = nodes[f].face;
            bool touches = false;
            for (int e = 0; e < 3 && !touches; ++e) {
                const auto *bc = face->edge(e)->boundaryComponent();
                touches = bc && bc->index() == boundaryComponent;
            }
            count += touches;
        }
    }
    return count;
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

        if (protectedBoundaryComponent &&
            !exempt.contains(static_cast<int>(i))) {
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
    Graph base = buildGraph_(skeleton, protectedBoundaryComponent, seedFaces);

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
            // Phase 2 in addFace(); see buildGraph_() above.
            //
            // else if (EmbeddedSubmanifold<dim, subdim>::
            //              hasUnexplainedSelfCollision(node.face,
            //              node.gluings))
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

template SearchStats
EmbeddingSearch<3, 2>::runSearch_<EmbeddedSubmanifold<3, 2>>(
    unsigned, BoundaryCondition, std::function<EmbeddedSubmanifold<3, 2>()>,
    std::function<std::unique_ptr<RunSearchThreadHook<3, 2>>()>,
    std::function<void(const std::vector<int> &)>, const SearchCallbacks &,
    RunSearchAuxHooks &, unsigned, long long, std::optional<long long>,
    std::optional<unsigned>, bool, std::optional<long long>, long long,
    long long);
template SearchStats
EmbeddingSearch<4, 2>::runSearch_<EmbeddedSubmanifold<4, 2>>(
    unsigned, BoundaryCondition, std::function<EmbeddedSubmanifold<4, 2>()>,
    std::function<std::unique_ptr<RunSearchThreadHook<4, 2>>()>,
    std::function<void(const std::vector<int> &)>, const SearchCallbacks &,
    RunSearchAuxHooks &, unsigned, long long, std::optional<long long>,
    std::optional<unsigned>, bool, std::optional<long long>, long long,
    long long);
template SearchStats EmbeddingSearch<4, 2>::runSearch_<KnottedSurface>(
    unsigned, BoundaryCondition, std::function<KnottedSurface()>,
    std::function<std::unique_ptr<RunSearchThreadHook<4, 2>>()>,
    std::function<void(const std::vector<int> &)>, const SearchCallbacks &,
    RunSearchAuxHooks &, unsigned, long long, std::optional<long long>,
    std::optional<unsigned>, bool, std::optional<long long>, long long,
    long long);
