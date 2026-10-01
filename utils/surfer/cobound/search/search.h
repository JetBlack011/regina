//
//  search.h
//
//  One search from a link, as every driver runs it.
//

#ifndef SURFER_COBOUND_SEARCH_H
#define SURFER_COBOUND_SEARCH_H

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <optional>
#include <thread>

#include "surfer/submanifold/submanifold.h"

/*! \file utils/surfer/cobound/search/search.h
 *  \brief One search from a link: which BoundaryCondition it runs under
 *  (conditionFor()), and the watchdog that ends it at its surface target or
 *  deadline (RowWatchdog).
 */

namespace rowsearch {

/** Which BoundaryCondition to search a row under; see conditionFor(). */
enum class BoundaryConditionMode { automatic, connected, proper };

/**
 * The BoundaryCondition for a row of `componentCount` components.
 *
 * `connected` (one curve per ambient boundary component) prunes far harder,
 * but it caps the curve count on the far side too, so a knot row searched
 * under it can never find a knot-to-link cobordism. It is honoured only for
 * a knot: a multi-component link can never meet it on its own search side.
 * `automatic` is `connected` for a knot and `proper` otherwise.
 */
BoundaryCondition conditionFor(BoundaryConditionMode mode, int componentCount);

/** When RowWatchdog ends a row's search. Unset limits never fire. */
struct WatchdogLimits {
    std::optional<long long> surfaceTarget;
    std::optional<double> rowSeconds;
    std::optional<double> sweepSeconds;
    std::chrono::steady_clock::time_point sweepStart{};
    std::optional<double> quiescenceSeconds;
    /** Milliseconds since the row last taught us something new; read only
     *  with `quiescenceSeconds`. */
    std::function<long long()> idleMillis;

    bool any() const {
        return surfaceTarget || rowSeconds || sweepSeconds || quiescenceSeconds;
    }
};

/**
 * Polls every 200 ms, while the search and its drain run, and calls
 * `endRow(why)` at most once: "surface-target", "timeout" or "quiescent".
 * stop() wakes it at once rather than waiting out its poll.
 * The surface target is checked before the clocks, so a row reaching it in
 * the same tick as a deadline records "surface-target": the two mean
 * different things to anyone later reading the negative.
 *
 * No thread is started when no limit is set.
 */
class RowWatchdog {
public:
    RowWatchdog(WatchdogLimits limits, std::function<void(const char *)> endRow);
    ~RowWatchdog();
    RowWatchdog(const RowWatchdog &) = delete;
    RowWatchdog &operator=(const RowWatchdog &) = delete;

    /** The live count of surfaces satisfying the boundary condition, from
     *  SearchCallbacks::onProgress (the only place it is handed out). */
    void publishSatisfying(long long count) {
        satisfying_.store(count, std::memory_order_relaxed);
    }
    /** Stops polling and joins; idempotent. Call once the search returns. */
    void stop();

private:
    WatchdogLimits limits_;
    std::function<void(const char *)> endRow_;
    std::atomic<long long> satisfying_{0};
    std::atomic<bool> done_{false};
    std::mutex wakeMutex_;
    std::condition_variable wake_;
    std::thread thread_;
};

} // namespace rowsearch

#endif // SURFER_COBOUND_SEARCH_H
