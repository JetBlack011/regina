//
//  timers.h
//
//  Wall and CPU clocks, one of each.
//

#ifndef SURFER_COBOUND_TIMERS_H
#define SURFER_COBOUND_TIMERS_H

#include <chrono>

#include <sys/resource.h>

/*! \file utils/surfer/cobound/driver/timers.h
 *  \brief The clocks cobound reads: wall time since a point, and the
 *  process's CPU time. (linknaming, below cobound, keeps its own
 *  microsSince() for the namer's statistics: the dependency rule lets it
 *  see nothing of cobound's.)
 */

namespace timers {

using Clock = std::chrono::steady_clock;

/** Wall seconds since `t`. */
inline double secondsSince(Clock::time_point t) {
    return std::chrono::duration<double>(Clock::now() - t).count();
}

/** Wall microseconds since `t`. */
inline long long microsSince(Clock::time_point t) {
    return std::chrono::duration_cast<std::chrono::microseconds>(Clock::now() - t).count();
}

/** This process's CPU seconds so far, user and system, every thread. */
inline double processCpuSeconds() {
    rusage u{};
    getrusage(RUSAGE_SELF, &u);
    return static_cast<double>(u.ru_utime.tv_sec + u.ru_stime.tv_sec) +
           static_cast<double>(u.ru_utime.tv_usec + u.ru_stime.tv_usec) / 1e6;
}

} // namespace timers

#endif // SURFER_COBOUND_TIMERS_H
