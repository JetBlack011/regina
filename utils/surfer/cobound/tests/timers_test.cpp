// timers_test.cpp: the clocks cobound reads (../driver/timers.h). Wall
// time since a point, in seconds and microseconds, agree with each other and
// with a known pause; the process's CPU time counts every thread's.

#include <chrono>
#include <ctime>
#include <thread>
#include <vector>

#include "cobound/driver/timers.h"
#include "linknaming/tests/check.h"

namespace {

// This thread's own CPU seconds (not the process's), to spin on.
double threadCpuSeconds() {
  timespec ts{};
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts);
  return static_cast<double>(ts.tv_sec) + static_cast<double>(ts.tv_nsec) / 1e9;
}

// Burns `seconds` of this thread's CPU, however busy the machine is.
void spin(double seconds) {
  const double start = threadCpuSeconds();
  volatile unsigned long sink = 0;
  while (threadCpuSeconds() - start < seconds)
    for (int i = 0; i < 10000; ++i) sink = sink + static_cast<unsigned long>(i);
}

void testWallClocks() {
  const timers::Clock::time_point t = timers::Clock::now();
  CHECK(timers::secondsSince(t) >= 0.0, "no time has gone backwards");
  std::this_thread::sleep_for(std::chrono::milliseconds(30));
  const long long us = timers::microsSince(t);
  const double s = timers::secondsSince(t);
  CHECK(us >= 30000, "a 30 ms pause is at least 30,000 microseconds");
  CHECK(s >= 0.030, "...and at least 0.030 seconds");
  CHECK(s >= us / 1e6, "the later reading is no earlier than the first");
  CHECK(s - us / 1e6 < 0.5, "the two clocks agree");
}

void testProcessCpu() {
  const double before = timers::processCpuSeconds();
  CHECK(before >= 0.0, "CPU time is non-negative");
  std::vector<std::thread> threads;
  for (int k = 0; k < 2; ++k) threads.emplace_back([] { spin(0.05); });
  for (std::thread &th : threads) th.join();
  const double after = timers::processCpuSeconds();
  CHECK(after - before >= 0.09,
        "two threads' 0.05 s each are both counted (every thread, not the caller's)");
}

} // namespace

int main() {
  testWallClocks();
  testProcessCpu();
  return checks::finish("timers_test");
}
