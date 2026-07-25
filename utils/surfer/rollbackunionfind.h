//
//  rollbackunionfind.h
//
//  Created by John Teague on 07/23/2026.
//

#ifndef ROLLBACKUNIONFIND_H

#define ROLLBACKUNIONFIND_H

#include <cstddef>
#include <vector>

/*! \file utils/surfer/rollbackunionfind.h
 *  \brief A union-find data structure supporting rollback of recent merges.
 */

/**
 * A union-find over a fixed n-element universe supporting exact rollback of
 * the most recently performed unite() calls, in LIFO order.
 *
 * This deliberately omits path compression, which would make rollback
 * unsound, and instead unions by size: find() stays O(log(component size))
 * while every unite()/rollback step stays O(1).
 */
class RollbackUnionFind {
public:
  /**
   * Creates a union-find over the universe {0, ..., n-1}, with every
   * element initially in its own singleton set.
   */
  explicit RollbackUnionFind(size_t n);

  /**
   * Returns the representative of the set currently containing x.
   */
  int find(int x) const;

  /**
   * Merges the sets containing x and y.
   *
   * \note Named unite() rather than union(), since \c union is a C++
   * keyword.
   *
   * \return \c false if x and y were already in the same set (in which
   * case nothing changes and nothing is logged), or \c true otherwise.
   */
  bool unite(int x, int y);

  /**
   * Returns an opaque marker for the current state, to be passed to a
   * later rollbackTo() call.
   */
  size_t checkpoint() const { return log_.size(); }

  /**
   * Reverses every unite() call performed since `mark` was returned by
   * checkpoint(), in strict LIFO order.
   */
  void rollbackTo(size_t mark);

private:
  std::vector<int> parent_;
      /**< parent_[x] is x's parent in its tree, or x itself if x is a
           root. */
  std::vector<int> size_;
      /**< size_[x] is the size of x's set; only meaningful when x is a
           root. */

  /**
   * A single logged unite() call, sufficient to undo it exactly.
   */
  struct UndoEntry {
    int child;
        /**< The root that was merged into another (the non-survivor). */
    int survivorOldSize;
        /**< The surviving root's size() before this merge. */
  };
  std::vector<UndoEntry> log_; /**< The undo log for unite(), in call order. */
};

#endif // ROLLBACKUNIONFIND_H
