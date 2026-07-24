//
//  rollbackunionfind.cpp
//

#include "rollbackunionfind.h"

#include <utility>

RollbackUnionFind::RollbackUnionFind(size_t n) : parent_(n), size_(n, 1) {
  for (size_t i = 0; i < n; ++i)
    parent_[i] = static_cast<int>(i);
}

int RollbackUnionFind::find(int x) const {
  while (parent_[x] != x)
    x = parent_[x];
  return x;
}

bool RollbackUnionFind::unite(int x, int y) {
  int rx = find(x), ry = find(y);
  if (rx == ry)
    return false;

  if (size_[rx] < size_[ry])
    std::swap(rx, ry);

  int survivorOldSize = size_[rx];
  parent_[ry] = rx;
  size_[rx] += size_[ry];
  log_.push_back({.child = ry, .survivorOldSize = survivorOldSize});
  return true;
}

void RollbackUnionFind::rollbackTo(size_t mark) {
  while (log_.size() > mark) {
    UndoEntry e = log_.back();
    log_.pop_back();
    int rx = parent_[e.child];
    parent_[e.child] = e.child;
    size_[rx] = e.survivorOldSize;
  }
}
