//
//  fatal.cpp
//

#include "cobound/driver/fatal.h"

#include <atomic>
#include <cstdlib>
#include <iostream>
#include <mutex>

namespace fatal {

namespace {

// Set (once) from any thread via flag(); only ever read/acted on back on the
// main thread, after the offending search has returned and any watchdog
// thread has joined -- never torn down the process from within a callback,
// which may be running on one of several concurrent worker threads.
std::atomic<bool> detected{false};
std::mutex messageMutex;
std::string message_;
std::function<void()> beforeHalt_;

} // namespace

void flag(std::string message) {
  bool expected = false;
  if (!detected.compare_exchange_strong(expected, true))
    return;
  std::lock_guard<std::mutex> lock(messageMutex);
  message_ = std::move(message);
}

bool flagged() { return detected.load(); }

void beforeHalt(std::function<void()> write) { beforeHalt_ = std::move(write); }

void haltIfFlagged() {
  if (!detected.load())
    return;
  if (beforeHalt_)
    beforeHalt_();
  std::string message;
  {
    std::lock_guard<std::mutex> lock(messageMutex);
    message = message_;
  }
  std::cerr
      << "\n\x1b[1;31m"
      << "################################################################"
         "################\n"
         "#####  FATAL: SEARCH LIBRARY BUG DETECTED -- HALTING NOW  #####\n"
         "################################################################"
         "################\x1b[0m\n\n"
      << message << "\n\n"
      << "This is a mathematical impossibility, not a data or literature "
         "issue -- it means\n"
         "surfacesearch.h/embeddedsubmanifold.h computed something wrong "
         "(a bad genus, a\n"
         "false claim of embeddedness, etc.). Nothing else this run could "
         "report from here\n"
         "on can be trusted, so the program is stopping immediately rather "
         "than continuing\n"
         "to the next knot/link. Whatever was already durably written to "
         "--output before\n"
         "this point is unaffected and safe to keep.\n";
  std::exit(2);
}

} // namespace fatal
