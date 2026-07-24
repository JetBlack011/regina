#ifndef ROLLBACKUNIONFIND_H

#define ROLLBACKUNIONFIND_H

#include <cstddef>
#include <vector>

// A union-find over a fixed n-element universe supporting exact rollback of
// the most recently performed unite() calls, in LIFO order. Deliberately
// omits path compression (which would make rollback unsound) and instead
// unions by size, keeping find() at O(log(component size)) while every
// unite()/rollback step stays O(1).
class RollbackUnionFind {
public:
  explicit RollbackUnionFind(size_t n);

  int find(int x) const;

  // Merges the sets containing x and y. Returns false (no-op, nothing
  // logged) if they were already in the same set.
  bool unite(int x, int y);

  size_t checkpoint() const { return log_.size(); }

  // Reverses every unite() call back to (but not including) the state at
  // `mark` (a value previously returned by checkpoint()), in strict LIFO
  // order.
  void rollbackTo(size_t mark);

private:
  std::vector<int> parent_;
  std::vector<int> size_;

  struct UndoEntry {
    int child;
    int survivorOldSize;
  };
  std::vector<UndoEntry> log_;
};

#endif // ROLLBACKUNIONFIND_H
