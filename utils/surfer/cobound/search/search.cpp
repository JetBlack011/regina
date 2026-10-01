//
//  search.cpp
//

#include "cobound/search/search.h"

namespace rowsearch {

BoundaryCondition conditionFor(BoundaryConditionMode mode, int componentCount) {
    switch (mode) {
    case BoundaryConditionMode::proper:
        return BoundaryCondition::proper;
    case BoundaryConditionMode::connected:
    case BoundaryConditionMode::automatic:
    default:
        return componentCount == 1 ? BoundaryCondition::connected
                                   : BoundaryCondition::proper;
    }
}

RowWatchdog::RowWatchdog(WatchdogLimits limits,
                         std::function<void(const char *)> endRow)
    : limits_(std::move(limits)), endRow_(std::move(endRow)) {
    if (!limits_.any())
        return;
    thread_ = std::thread([this] {
        const auto rowDeadline =
            std::chrono::steady_clock::now() +
            std::chrono::duration<double>(limits_.rowSeconds.value_or(0));
        while (!done_.load(std::memory_order_relaxed)) {
            {
                std::unique_lock<std::mutex> lock(wakeMutex_);
                wake_.wait_for(lock, std::chrono::milliseconds(200), [this] {
                    return done_.load(std::memory_order_relaxed);
                });
            }
            if (done_.load(std::memory_order_relaxed))
                break;
            // Checked before the clocks: see the class comment.
            if (limits_.surfaceTarget &&
                satisfying_.load(std::memory_order_relaxed) >=
                    *limits_.surfaceTarget) {
                endRow_("surface-target");
                break;
            }
            const auto now = std::chrono::steady_clock::now();
            if (limits_.rowSeconds && now >= rowDeadline) {
                endRow_("timeout");
                break;
            }
            if (limits_.sweepSeconds &&
                now - limits_.sweepStart >
                    std::chrono::duration<double>(*limits_.sweepSeconds)) {
                endRow_("timeout");
                break;
            }
            // Quiescence: this row has stopped teaching us anything new.
            // Meaningful during the drain too, which is where witnesses are
            // actually identified.
            if (limits_.quiescenceSeconds &&
                limits_.idleMillis() >
                    static_cast<long long>(*limits_.quiescenceSeconds * 1000)) {
                endRow_("quiescent");
                break;
            }
        }
    });
}

void RowWatchdog::stop() {
    {
        std::lock_guard<std::mutex> lock(wakeMutex_);
        done_.store(true, std::memory_order_relaxed);
    }
    wake_.notify_all();
    if (thread_.joinable())
        thread_.join();
}

RowWatchdog::~RowWatchdog() { stop(); }

} // namespace rowsearch
