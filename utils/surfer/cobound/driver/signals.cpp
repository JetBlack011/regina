//
//  signals.cpp
//

#include "cobound/driver/signals.h"

#include <atomic>
#include <csignal>
#include <cstdlib>

namespace {

// Lock-free atomics only: a signal handler may touch nothing else.
std::atomic<int> g_signal{0};
std::atomic<int> g_count{0};
static_assert(std::atomic<int>::is_always_lock_free);

extern "C" void onSignal(int sig) {
  if (g_count.fetch_add(1, std::memory_order_relaxed) > 0)
    std::_Exit(128 + sig);
  g_signal.store(sig, std::memory_order_relaxed);
}

} // namespace

namespace runsignals {

void install() {
  struct sigaction sa {};
  sa.sa_handler = onSignal;
  sigemptyset(&sa.sa_mask);
  sa.sa_flags = SA_RESTART;
  sigaction(SIGINT, &sa, nullptr);
  sigaction(SIGTERM, &sa, nullptr);
}

bool interrupted() { return g_signal.load(std::memory_order_relaxed) != 0; }

int which() { return g_signal.load(std::memory_order_relaxed); }

const char *name() {
  switch (which()) {
  case SIGINT:
    return "SIGINT";
  case SIGTERM:
    return "SIGTERM";
  default:
    return "";
  }
}

} // namespace runsignals
