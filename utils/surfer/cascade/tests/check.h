// check.h: the minimal assertion helpers the cascade tests share.

#pragma once

#include <iostream>
#include <sstream>
#include <string>

namespace cascadetest {

inline int passed = 0, failed = 0;

inline void report(bool ok, const std::string &what, const std::string &detail) {
  if (ok) {
    ++passed;
  } else {
    ++failed;
    std::cout << "FAIL: " << what;
    if (!detail.empty())
      std::cout << " -- " << detail;
    std::cout << "\n";
  }
}

inline int finish(const char *name) {
  std::cout << name << ": " << passed << " passed, " << failed << " failed\n";
  return failed == 0 ? 0 : 1;
}

} // namespace cascadetest

#define CHECK(cond, what)                                                      \
  cascadetest::report(static_cast<bool>(cond), what, #cond)

#define CHECK_EQ(actual, expected, what)                                       \
  do {                                                                         \
    auto _a = (actual);                                                        \
    auto _e = (expected);                                                      \
    std::ostringstream _o;                                                     \
    if (!(_a == _e))                                                           \
      _o << "got " << _a << ", expected " << _e;                               \
    cascadetest::report(_a == _e, what, _o.str());                             \
  } while (0)
