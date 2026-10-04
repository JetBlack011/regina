// leaves.h
//
// Which outside facts may become proof-graph leaves, and the links that
// receive them. See README.md, "Leaf facts".

#pragma once

#include <map>
#include <ostream>
#include <string>
#include <vector>

#include "cobound/bounds/cobordismgraph.h"
#include "cobound/bounds/links.h"
#include "linknaming/linknamer.h"
#include "linknaming/tables.h"

namespace bounds {

// A table's literature 4-genus is parsed by linknaming::parseTableG4()
// (linknaming/tables.h), the one parser of that field.

/**
 * Whether a node whose exact table class is `nodeClass` may receive that
 * class's literature UPPER bound as a leaf. Never the target's own class:
 * the target's literature value proving the target would be circular, and a
 * duplicate node of the target (a diagram the registry did not recognise)
 * would otherwise do exactly that. Lower bounds are always kept: they only
 * gate contradictions.
 */
bool mayUseLiteratureUpperBound(const std::string &linkClass,
                                const std::string &targetClass,
                                bool literatureAllowed);

/**
 * The outside facts a cobordism graph's links rest on, attached as each link
 * joins the graph: the link named exactly (ExactNamer::identify(); a knot
 * the tables hold only summand by summand, by the whole-diagram namer), and
 * given what its name proves -- its literature lower bound always (the
 * contradiction gates read it), its literature upper bound as a leaf unless
 * it is the target's own class, and the unknot's disc when it is a slice
 * composite (an anchor). An untabulated link with visible sum spheres is cut
 * into its prime summands, each a link of its own, joined to it by a sum
 * edge and named in turn.
 *
 * A goal run (cascadesearch) and a depth-0 search's own graph (SearchJudge)
 * name their links this one way.
 */
class LinkAxioms {
public:
  struct Options {
    /// Literature upper bounds may be leaves (cascadesearch --literature).
    bool literature = true;
    /// Compare table names by their link class (ExactNamer::canonicalName(),
    /// computed lazily per base); false, by the name itself, which costs
    /// nothing at depth 0 (plan, "Startup per process").
    bool classes = true;
    unsigned threads = 1;
    /// Where sums found and slice composites anchored are announced; none
    /// for a graph that is only a search's judge.
    std::ostream *log = nullptr;
  };

  /// The limits every graph's link naming runs under: the cheap ones
  /// (no Reidemeister searches; a few simplification tries).
  static linknaming::NamerLimits namerLimits();

  /// All references must outlive this. `symmetries` may change content
  /// later (cascadesearch fills its NameTable after constructing this).
  LinkAxioms(CobordismGraph &graph, LinkRegistry &links, const linknaming::ExactTables &tables,
             const linknaming::ExactNamer &namer,
             const linknaming::SymmetryTable &symmetries, Options options);
  LinkAxioms(const LinkAxioms &) = delete;
  LinkAxioms &operator=(const LinkAxioms &) = delete;

  /// Names every node of `ns` (on the threads: names are independent of
  /// each other), then records each at `depth` with its table name and
  /// outside facts, in node order, as one at a time would; summands of a
  /// sum are named at depth + 1.
  void name(const std::vector<LinkId> &ns, int depth);

  /// The class a table name stands for (Options::classes), else the name.
  std::string classOf(const std::string &tableName) const;

  /// The target: never anchored as a composite, and its class's literature
  /// upper bound is no leaf, even of a duplicate node of it.
  LinkId target = -1;
  std::string targetClass;

  // What naming found.
  std::map<LinkId, std::string> tableName; ///< a table (or composite) name
  std::map<LinkId, int> depth;             ///< hops from the target when met
  std::map<LinkId, std::vector<LinkId>> sumOf; ///< each sum node's summands
  int anchors = 0;

private:
  void applyName(LinkId n, const linknaming::PieceName &pn);
  void applyComposite(LinkId n, const linknaming::LinkName &fs);
  void applySum(LinkId n, const std::vector<linknaming::GaussDiagram> &primes, int depth);
  std::vector<LinkId> linksSince(size_t first) const;

  CobordismGraph &g_;
  LinkRegistry &reg_;
  const linknaming::ExactTables &tables_;
  const linknaming::ExactNamer &namer_;
  const linknaming::SymmetryTable &symmetries_;
  Options options_;
};

} // namespace bounds
