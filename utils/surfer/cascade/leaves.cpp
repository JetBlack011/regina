// leaves.cpp

#include "leaves.h"

#include <cctype>

namespace cascade {

namespace {
std::optional<int> parseInt(const std::string &s) {
  if (s.empty()) return std::nullopt;
  for (char c : s)
    if (!std::isdigit(static_cast<unsigned char>(c))) return std::nullopt;
  return std::stoi(s);
}
} // namespace

std::optional<std::pair<int, int>> parseTableG4(const std::string &s) {
  if (s.size() >= 5 && s.front() == '[' && s.back() == ']') {
    const auto semi = s.find(';');
    if (semi == std::string::npos) return std::nullopt;
    auto lo = parseInt(s.substr(1, semi - 1));
    auto hi = parseInt(s.substr(semi + 1, s.size() - semi - 2));
    if (!lo || !hi || *lo > *hi) return std::nullopt;
    return std::make_pair(*lo, *hi);
  }
  if (auto v = parseInt(s)) return std::make_pair(*v, *v);
  return std::nullopt;
}

bool mayUseLiteratureUpperBound(const std::string &nodeClass,
                                const std::string &targetClass,
                                bool literatureAllowed) {
  if (!literatureAllowed || nodeClass.empty()) return false;
  return targetClass.empty() || nodeClass != targetClass;
}

} // namespace cascade
