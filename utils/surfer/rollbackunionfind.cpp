//
//  rollbackunionfind.cpp
//
//  Created by John Teague on 07/23/2026.
//

#include "rollbackunionfind.h"

#include <utility>

RollbackUnionFind::RollbackUnionFind(size_t n)
    : parent_(n), size_(n, 1), parityToParent_(n, false) {
  for (size_t i = 0; i < n; ++i)
    parent_[i] = static_cast<int>(i);
}

int RollbackUnionFind::find(int x) const {
  while (parent_[x] != x)
    x = parent_[x];
  return x;
}

std::pair<int, bool> RollbackUnionFind::findWithParity_(int x) const {
  bool parity = false;
  while (parent_[x] != x) {
    parity ^= parityToParent_[x];
    x = parent_[x];
  }
  return {x, parity};
}

bool RollbackUnionFind::unite(int x, int y) { return unite(x, y, true); }

bool RollbackUnionFind::unite(int x, int y, bool sameOrientation) {
  auto [rx, px] = findWithParity_(x);
  auto [ry, py] = findWithParity_(y);
  if (rx == ry)
    return false;

  bool mismatch = !sameOrientation;

  if (size_[rx] < size_[ry]) {
    std::swap(rx, ry);
    std::swap(px, py);
  }

  int survivorOldSize = size_[rx];
  // Chosen so that, after this merge, walking from y up through ry to rx
  // accumulates a total parity of (px ^ mismatch) relative to rx -- i.e.
  // parity(x) ^ parity(y) == mismatch, matching the requested relationship
  // (see the class documentation's derivation).
  parityToParent_[ry] = px ^ py ^ mismatch;
  parent_[ry] = rx;
  size_[rx] += size_[ry];
  log_.push_back({.child = ry, .survivorOldSize = survivorOldSize});
  return true;
}

bool RollbackUnionFind::sameOrientation(int x, int y) const {
  auto [rx, px] = findWithParity_(x);
  auto [ry, py] = findWithParity_(y);
  return px == py; // \pre rx == ry, i.e. x and y are already in the same set
}

void RollbackUnionFind::rollbackTo(size_t mark) {
  while (log_.size() > mark) {
    UndoEntry e = log_.back();
    log_.pop_back();
    int rx = parent_[e.child];
    parent_[e.child] = e.child;
    parityToParent_[e.child] = false; // a root's parity relative to itself is trivially 0
    size_[rx] = e.survivorOldSize;
  }
}
