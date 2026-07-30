//
//  rollbackunionfind.h
//
//  Created by John Teague on 07/23/2026.
//

#ifndef ROLLBACKUNIONFIND_H

#define ROLLBACKUNIONFIND_H

#include <cstddef>
#include <utility>
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
   * As unite(x, y), but additionally records whether x and y should be
   * considered in the "same" or "different" relative state (e.g. the same
   * vs. opposite orientation) -- a weighted/signed union-find, the
   * standard incremental technique for problems like bipartiteness or
   * orientability checking. The relationship is tracked as a parity bit
   * per element, accumulated along the union-find tree; see
   * sameOrientation().
   *
   * Does not itself detect or report a contradiction (x and y already in
   * the same set, united here with a relationship inconsistent with their
   * existing accumulated parity) -- callers that care (e.g. orientability
   * tracking, where a contradiction means "not orientable") must check
   * sameOrientation(x, y) themselves *before* calling this, while x and y
   * are still known to be in the same set.
   *
   * \return \c false if x and y were already in the same set (in which
   * case nothing changes and nothing is logged, regardless of
   * `sameOrientation`), or \c true otherwise.
   */
  bool unite(int x, int y, bool sameOrientation);

  /**
   * Returns whether x and y are in the "same" relative state, per the
   * accumulated chain of unite(x, y, sameOrientation) calls connecting
   * them (unite(x, y) with no explicit relationship is equivalent to
   * unite(x, y, true), so this is well-defined for sets built with a mix
   * of both overloads).
   *
   * \pre x and y are in the same set (find(x) == find(y)).
   */
  bool sameOrientation(int x, int y) const;

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
  std::vector<bool> parityToParent_;
      /**< parityToParent_[x]: x's accumulated parity relative to
           parent_[x] -- only meaningful (and only ever read) when x is
           not currently a root; see findWithParity_(). Maintained
           alongside parent_/size_ by unite(x, y, bool). */

  /**
   * Walks from `x` up to its current root, accumulating the XOR of every
   * parityToParent_ entry along the way -- the shared traversal both
   * find() and the parity-aware operations above are built on. No path
   * compression (see the class documentation), so this is the same
   * O(log(component size)) walk find() already does, just also tracking
   * parity.
   */
  std::pair<int, bool> findWithParity_(int x) const;

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
