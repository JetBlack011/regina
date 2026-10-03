// leaves.cpp

#include "cobound/bounds/axioms.h"

#include <algorithm>
#include <numeric>
#include <optional>
#include <sstream>

#include "cobound/parallelfor.h"

using exactnaming::GaussDiagram;

namespace cascade {

bool mayUseLiteratureUpperBound(const std::string &nodeClass,
                                const std::string &targetClass,
                                bool literatureAllowed) {
  if (!literatureAllowed || nodeClass.empty()) return false;
  return targetClass.empty() || nodeClass != targetClass;
}

exactnaming::NamerLimits NodeAxioms::namerLimits() {
  exactnaming::NamerLimits l;
  l.simplifyTries = 4;
  l.exhaustiveHeight = 0;
  l.searchHeight = -1;
  l.maxDeepCrossings = 0;
  l.tableSideHeight = -1;
  return l;
}

NodeAxioms::NodeAxioms(ProofGraph &graph, NodeRegistry &nodes,
                       const exactnaming::ExactTables &tables,
                       const exactnaming::ExactNamer &namer,
                       const exactnaming::SymmetryTable &symmetries, Options options)
    : g_(graph), reg_(nodes), tables_(tables), namer_(namer), symmetries_(symmetries),
      options_(options) {}

std::string NodeAxioms::classOf(const std::string &name) const {
  if (!options_.classes) return name;
  const exactnaming::TableEntry *e = tables_.entry(name);
  return e ? namer_.canonicalName(*e) : name;
}

std::vector<NodeId> NodeAxioms::nodesSince(size_t first) const {
  std::vector<NodeId> ns;
  for (size_t m = first; m < g_.nodeCount(); ++m) ns.push_back(static_cast<NodeId>(m));
  return ns;
}

void NodeAxioms::name(const std::vector<NodeId> &ns, int atDepth) {
  // Naming is most of what a hop does outside its search, and each node's
  // name is independent of the others', so the names are found on a pool
  // (ExactNamer is safe to share: its caches are locked, and its SnapPea
  // calls serialised) and applied here in node order, as one at a time would.
  std::vector<std::optional<exactnaming::PieceName>> names(ns.size());
  // A knot identify() leaves untabulated may be a connected sum the tables
  // hold only summand by summand: the whole-diagram namer cuts it at its
  // visible sum spheres and composes `A#mB` (exactnamer.h, step 3).
  std::vector<std::optional<exactnaming::FarSideName>> composites(ns.size());
  // Any untabulated node (knot or link) with visible sum spheres is cut into
  // its prime summands, which become nodes joined to it by a sum edge
  // (applySum; paper lem:sum-partitions): the tables list prime links only,
  // and most far sides of small links are sums of smaller ones.
  std::vector<std::vector<GaussDiagram>> summands(ns.size());
  parallelFor(ns.size(), std::max(options_.threads, 1u), [&](size_t i) {
    const NodeId n = ns[i];
    if (!reg_.known(n) || n == reg_.unknot()) return;
    try {
      const GaussDiagram &d = reg_.info(n).diagram;
      names[i] = namer_.identify(d);
      if (names[i]->by == exactnaming::PieceName::By::untabulated && d.crossings() > 0) {
        if (d.components() == 1) {
          exactnaming::FarSideName fs = namer_.name(d.link());
          if (fs.exact && fs.pinned && fs.pieces.size() >= 2 &&
              fs.name.find('#') != std::string::npos)
            composites[i] = std::move(fs);
        }
        GaussDiagram own = d;
        own.origin.resize(own.components());
        std::iota(own.origin.begin(), own.origin.end(), 0);
        std::vector<GaussDiagram> primes;
        namer_.decompose(own, primes);
        if (primes.size() >= 2) summands[i] = std::move(primes);
      }
    } catch (const std::exception &) {
    }
  });
  for (size_t i = 0; i < ns.size(); ++i) {
    depth.emplace(ns[i], atDepth);
    if (composites[i]) applyComposite(ns[i], *composites[i]);
    else if (names[i]) applyName(ns[i], *names[i]);
  }
  for (size_t i = 0; i < ns.size(); ++i)
    if (!summands[i].empty()) applySum(ns[i], summands[i], atDepth);
}

void NodeAxioms::applySum(NodeId n, const std::vector<GaussDiagram> &primes, int atDepth) {
  // Each prime is interned as a node (a duplicate costs search, never
  // soundness), and its components are mapped to the whole's through the
  // registry's component map and the prime's origins. The new summand
  // nodes are then named like any other (literature leaves included), so
  // the sum rule can combine their bounds.
  const size_t before = g_.nodeCount();
  std::vector<NodeId> pieces;
  std::vector<std::vector<int>> maps;
  for (size_t k = 0; k < primes.size(); ++k) {
    const GaussDiagram &p = primes[k];
    NodeMatch m = reg_.intern(p, "summand " + std::to_string(k) + " of node " + std::to_string(n));
    std::vector<int> map(p.components(), -1);
    for (size_t c = 0; c < p.components(); ++c)
      map[static_cast<size_t>(m.componentMap[c])] = static_cast<int>(p.origin[c]);
    pieces.push_back(m.node);
    maps.push_back(std::move(map));
  }
  try {
    g_.addSum(n, pieces, maps);
  } catch (const std::exception &e) {
    if (options_.log)
      *options_.log << "[!] node " << n << ": summands not recorded: " << e.what() << "\n";
    return;
  }
  sumOf[n] = pieces;
  name(nodesSince(before), atDepth + 1);
  if (!options_.log) return;
  std::ostringstream o;
  o << "[+] node " << n << " is a sum along components of";
  for (NodeId p : pieces) {
    auto t = tableName.find(p);
    o << " node " << p << (t == tableName.end() ? "" : " (" + t->second + ")");
  }
  *options_.log << o.str() << "\n";
}

void NodeAxioms::applyComposite(NodeId n, const exactnaming::FarSideName &fs) {
  // The composite's name is recorded (certificates, node bounds, the
  // subject name stays cascade:, since no table row holds it). It is an
  // ANCHOR when its summands cancel in concordance (cobordismgraph.h
  // isElementarySlice: the explicit allowlist, or symmetry types from
  // --knot-symmetry): K # m(K^r) bounds a ribbon disc, so the node gets the
  // unknot's leaf, constructive like the unknot's, never for the target.
  tableName[n] = fs.name;
  if (n == target) return;
  if (!exactnaming::isElementarySlice(fs.name, symmetries_)) return;
  g_.addLeaf(n, Partition::coarsest(1), 0, "anchor " + fs.name);
  ++anchors;
  if (options_.log)
    *options_.log << "[+] node " << n << " is " << fs.name << " (" << fs.proof()
                  << "): a slice composite, anchored\n";
}

void NodeAxioms::applyName(NodeId n, const exactnaming::PieceName &pn) {
  if (pn.by == exactnaming::PieceName::By::untabulated || !pn.pinned() || pn.names.size() != 1)
    return;
  const std::string &name = pn.names.front();
  tableName[n] = name;
  const exactnaming::TableEntry *e = tables_.entry(name);
  if (!e) return;
  auto g4 = exactnaming::parseTableG4(e->g4);
  if (!g4) return;
  const auto [lo, hi] = *g4;
  g_.setGenusLowerBound(n, lo, "literature " + name + " " + e->g4);
  // Never let the target's own literature value prove the target, even
  // through a duplicate node of it (README.md, "Leaf facts").
  if (!mayUseLiteratureUpperBound(classOf(name), targetClass, options_.literature))
    return;
  g_.addLeaf(n, Partition::coarsest(g_.node(n).components), hi,
             "literature " + name + " " + e->g4);
}

} // namespace cascade
