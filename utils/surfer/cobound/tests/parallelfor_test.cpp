// parallelfor_test.cpp: parallelFor() (../parallelfor.h) calls work(i)
// exactly once for every i, whatever the thread count, and none for an
// empty range.

#include <atomic>
#include <string>
#include <thread>
#include <vector>

#include "cobound/parallelfor.h"
#include "linknaming/tests/check.h"

int main() {
  for (unsigned threads : {0u, 1u, 2u, 7u, 64u})
    for (size_t n : {size_t(0), size_t(1), size_t(5), size_t(1000)}) {
      std::vector<std::atomic<int>> calls(n);
      parallelFor(n, threads, [&](size_t i) { calls[i].fetch_add(1); });
      bool once = true;
      for (const auto &c : calls) once = once && c.load() == 1;
      CHECK(once, "n " + std::to_string(n) + ", threads " + std::to_string(threads) +
                      ": every index exactly once");
    }
  // The calling thread works too: with one thread, everything runs on it.
  const std::thread::id self = std::this_thread::get_id();
  bool onCaller = true;
  parallelFor(10, 1, [&](size_t) { onCaller = onCaller && std::this_thread::get_id() == self; });
  CHECK(onCaller, "one thread is the calling thread");
  return cascadetest::finish("parallelfor_test");
}
