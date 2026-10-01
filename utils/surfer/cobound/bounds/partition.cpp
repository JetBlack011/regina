// partition.cpp

#include "partition.h"

#include <functional>
#include <map>
#include <mutex>
#include <sstream>
#include <stdexcept>

namespace cascade {

Partition Partition::singletons(int n) {
  std::vector<int> l(n);
  for (int i = 0; i < n; ++i)
    l[i] = i;
  return fromLabels(l);
}

Partition Partition::coarsest(int n) {
  return fromLabels(std::vector<int>(n, 0));
}

Partition Partition::fromLabels(const std::vector<int> &labels) {
  Partition p;
  p.label_.resize(labels.size());
  std::map<int, int> renumber;
  for (size_t i = 0; i < labels.size(); ++i) {
    auto [it, inserted] =
        renumber.try_emplace(labels[i], static_cast<int>(renumber.size()));
    p.label_[i] = it->second;
  }
  p.blocks_ = static_cast<int>(renumber.size());
  return p;
}

bool Partition::refines(const Partition &coarser) const {
  if (coarser.size() != size())
    throw std::invalid_argument("Partition::refines: sizes differ");
  // Each of our blocks must map into a single block of `coarser`.
  std::vector<int> image(blocks_, -1);
  for (int i = 0; i < size(); ++i) {
    int &im = image[label_[i]];
    if (im < 0)
      im = coarser.label_[i];
    else if (im != coarser.label_[i])
      return false;
  }
  return true;
}

std::vector<std::vector<int>> Partition::blockList() const {
  std::vector<std::vector<int>> ans(blocks_);
  for (int i = 0; i < size(); ++i)
    ans[label_[i]].push_back(i);
  return ans;
}

std::string Partition::str() const {
  std::ostringstream o;
  for (const auto &b : blockList()) {
    o << '{';
    for (size_t k = 0; k < b.size(); ++k)
      o << (k ? "," : "") << b[k];
    o << '}';
  }
  if (label_.empty())
    o << "{}";
  return o.str();
}

namespace {
std::vector<Partition> buildAllPartitions(int n) {
  // Restricted growth strings: a[0] = 0, a[i] <= 1 + max(a[0..i-1]).
  std::vector<Partition> ans;
  std::vector<int> a(n, 0);
  std::function<void(int, int)> rec = [&](int i, int mx) {
    if (i == n) {
      ans.push_back(Partition::fromLabels(a));
      return;
    }
    for (int v = 0; v <= mx + 1; ++v) {
      a[i] = v;
      rec(i + 1, std::max(mx, v));
    }
  };
  if (n == 0)
    ans.push_back(Partition::fromLabels({}));
  else {
    a[0] = 0;
    rec(1, 0);
  }
  return ans;
}
} // namespace

const std::vector<Partition> &allPartitions(int n) {
  // Fixed for each n, so built once and shared: propagateLower() and the
  // lower report ask for the same small n millions of times, and rebuilding
  // them was half of a run's single-threaded time (2026-09-30 profile).
  // call_once keeps it lock-free after the first build, for callers on many
  // threads at once.
  constexpr int kShared = 16;
  static std::once_flag once[kShared];
  static std::vector<Partition> built[kShared];
  if (n >= 0 && n < kShared) {
    std::call_once(once[n], [n] { built[n] = buildAllPartitions(n); });
    return built[n];
  }
  thread_local std::map<int, std::vector<Partition>> large;
  auto it = large.find(n);
  if (it == large.end()) it = large.emplace(n, buildAllPartitions(n)).first;
  return it->second;
}

} // namespace cascade
