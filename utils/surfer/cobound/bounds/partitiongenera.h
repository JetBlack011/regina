// profile.h
//
// What the cascade knows about how a link bounds surfaces in B^4, and how a
// cobordism transports that knowledge. See README.md, "Profiles".
//
// Conventions. A link's components are numbered 0, ..., n-1. A surface F in
// B^4 bounded by it (smooth or locally flat, oriented, proper, with no closed
// components) determines
//   - a partition of {0, ..., n-1}: which components bound the same piece of F;
//   - its total genus: the sum of the genera of its pieces.
// A profile records, for a link, the pairs (partition, total genus) that
// surfaces we can exhibit achieve.
//
// Two facts make a profile a Pareto set (README.md, "Profiles"):
//   - tubing two pieces together keeps the total genus and merges their blocks
//     (paper lem:tubing), so (P, g) achievable implies (Q, g) achievable for
//     every Q coarser than P;
//   - so an entry (P, g) makes every (Q, h) with P refining Q and g <= h
//     redundant.
// Connected slice genus is the least genus over entries (every P refines the
// one-block partition); "strongly slice" (bounds disjoint discs) is an entry
// (singletons, 0).

#pragma once

#include <optional>
#include <string>
#include <vector>

#include "cobound/bounds/partition.h"

namespace bounds {

/**
 * The part of a cobordism C in S^3 x [0,1] that bounds depend on.
 *
 * C runs from its incoming link (at 0; for a witness, the searched link) to
 * its outgoing link (at 1; the far side). Its components are numbered
 * 0, ..., components-1; each boundary curve lies on exactly one of them.
 */
struct CobordismShape {
  int components = 0;
  /// The sum of the genera of C's components. A witness's `genus` column is
  /// exactly this (its tubed genus: see cobordismgraph.h).
  int genus = 0;
  /// For each incoming curve, in the incoming link's component order, the
  /// component of C it lies on.
  std::vector<int> inComponent;
  /// Likewise for the outgoing curves.
  std::vector<int> outComponent;

  /// Throws std::invalid_argument unless every curve's component is in range
  /// and every component carries at least one curve (a component with no
  /// boundary would be closed, which a witness never has).
  void validate() const;

  /// The same cobordism read in the other direction: incoming and outgoing
  /// swap. Genus and components are unchanged.
  CobordismShape reversed() const;

  /// The shape of the product cobordism L x [0,1] of an n-component link:
  /// n annuli, curve i on component i at both ends.
  static CobordismShape product(int n);
};

enum class Side { incoming, outgoing };

/// The result of glue().
struct Glued {
  /// How the resulting surface's pieces partition the free side's curves.
  Partition partition;
  /// An upper bound on its total genus (exact when no piece closes up).
  int genus = 0;
};

/**
 * Caps one side of a cobordism with a surface in B^4 and describes the
 * result as a surface bounded by the other ("free") side.
 *
 * `gluedPartition` partitions the glued side's curves (in that side's curve
 * order) into the pieces of a surface G bounded by them, of total genus
 * `gluedGenus`. The union S = C u G is then a surface in B^4 (after the
 * collar identification of README.md, "Composing hops") bounded by the free
 * side. Its pieces are the connected components of the graph whose vertices
 * are C's components and G's pieces and whose edges are the glued curves; a
 * graph component K has genus
 *     sum of its vertices' genera + b1(K),   b1(K) = E(K) - V(K) + 1,
 * because gluing along circles adds nothing to the Euler characteristic
 * (README.md derives this). Pieces carrying no free-side curve are closed and
 * are discarded (lem:tubing), which only lowers the genus: the returned genus
 * counts them anyway, so it is an upper bound, exact when nothing closes up.
 *
 * Returns nullopt if the free side has no curves.
 */
std::optional<Glued> glue(const CobordismShape &c, Side glued,
                          const Partition &gluedPartition, int gluedGenus);

/**
 * Whether a surface in B^4 whose pieces partition a link's components as `p`
 * can exist, as far as linking numbers decide it. Distinct pieces F_a, F_b are
 * disjoint, so lk(dF_a, dF_b) = F_a . F_b = 0: the total linking number between
 * the components of any two different blocks must vanish. (Not each pair
 * separately: blocks {0} and {1,2} with lk(0,1) = 1, lk(0,2) = -1 are fine.)
 *
 * `lk` is the symmetric linking matrix (diagonal ignored).
 */
bool linkingAllows(const Partition &p, const std::vector<std::vector<int>> &lk);

/// One achievable (partition, total genus) and the proof record behind it.
struct ProfileEntry {
  Partition partition;
  int genus = 0;
  long derivation = -1;
};

/**
 * The Pareto set of achievable (partition, genus) pairs for one link.
 */
class Profile {
public:
  Profile() = default;
  explicit Profile(int components) : n_(components) {}

  int components() const { return n_; }
  const std::vector<ProfileEntry> &entries() const { return entries_; }
  bool empty() const { return entries_.empty(); }

  /// Whether some entry already implies (p, g): its partition refines p and
  /// its genus is at most g.
  bool implies(const Partition &p, int g) const;

  /// Adds (p, g) unless implied, removing the entries it implies. Returns
  /// whether it was added.
  bool insert(const Partition &p, int g, long derivation);

  /// The least genus known for a surface whose partition refines `target`
  /// (and the record proving it), or nullopt. For the one-block partition
  /// this is the connected slice genus bound.
  std::optional<ProfileEntry> best(const Partition &target) const;

private:
  int n_ = 0;
  std::vector<ProfileEntry> entries_;
};

} // namespace bounds
