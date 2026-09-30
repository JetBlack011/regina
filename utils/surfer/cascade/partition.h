// partition.h
//
// A partition of the components {0, ..., n-1} of a link: which components of
// the link bound the same piece of a surface. See README.md, "Profiles".

#pragma once

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace cascade {

/**
 * A partition of {0, ..., n-1}, stored as a restricted growth string: element
 * i lies in block label()[i], and blocks are numbered 0, 1, 2, ... in order of
 * their first element. That normal form makes equality of partitions equality
 * of the label vectors.
 */
class Partition {
public:
  Partition() = default;

  /// Every element in its own block.
  static Partition singletons(int n);
  /// One block holding every element.
  static Partition coarsest(int n);
  /// Any labelling; elements with equal labels share a block. Labels may be
  /// arbitrary integers, and are renumbered into normal form.
  static Partition fromLabels(const std::vector<int> &labels);

  int size() const { return static_cast<int>(label_.size()); }
  int blocks() const { return blocks_; }
  int blockOf(int i) const { return label_[i]; }
  const std::vector<int> &labels() const { return label_; }

  /// Whether every block of *this lies inside a block of `coarser`, i.e.
  /// `coarser` can be reached from *this by merging blocks. Reflexive. Both
  /// must have the same size.
  bool refines(const Partition &coarser) const;

  /// The elements of each block, blocks in label order.
  std::vector<std::vector<int>> blockList() const;

  /// E.g. "{0,2}{1}".
  std::string str() const;

  bool operator==(const Partition &o) const { return label_ == o.label_; }
  bool operator<(const Partition &o) const { return label_ < o.label_; }

private:
  std::vector<int> label_;
  int blocks_ = 0;
};

/// Every partition of {0, ..., n-1}, in no particular order (Bell(n) of
/// them), built once per n and shared (safe from many threads).
const std::vector<Partition> &allPartitions(int n);

} // namespace cascade

template <> struct std::hash<cascade::Partition> {
  std::size_t operator()(const cascade::Partition &p) const noexcept {
    std::size_t h = 0x9e3779b97f4a7c15ULL;
    for (int l : p.labels())
      h = (h ^ static_cast<std::size_t>(l + 1)) * 0x100000001b3ULL;
    return h;
  }
};
