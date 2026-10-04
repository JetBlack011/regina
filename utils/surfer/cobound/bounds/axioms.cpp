// axioms.cpp

#include "cobound/bounds/axioms.h"

#include <algorithm>
#include <numeric>
#include <optional>
#include <sstream>

#include "cobound/parallelfor.h"
#include "cobound/frozen.h"

using linknaming::GaussDiagram;

namespace bounds {

bool mayUseLiteratureUpperBound(const std::string &linkClass,
                                const std::string &targetClass,
                                bool literatureAllowed) {
  if (!literatureAllowed || linkClass.empty()) return false;
  return targetClass.empty() || linkClass != targetClass;
}

linknaming::NamerLimits LinkAxioms::namerLimits() {
  linknaming::NamerLimits l;
  l.simplifyTries = 4;
  l.exhaustiveHeight = 0;
  l.searchHeight = -1;
  l.maxDeepCrossings = 0;
  l.tableSideHeight = -1;
  return l;
}

LinkAxioms::LinkAxioms(CobordismGraph &graph, LinkRegistry &links,
                       const linknaming::Tables &tables,
                       const linknaming::LinkNamer &namer,
                       const linknaming::SymmetryTable &symmetries, Options options)
    : g_(graph), reg_(links), tables_(tables), namer_(namer), symmetries_(symmetries),
      options_(options) {}

std::string LinkAxioms::classOf(const std::string &name) const {
  if (!options_.classes) return name;
  const linknaming::TableEntry *e = tables_.entry(name);
  return e ? namer_.canonicalName(*e) : name;
}

std::vector<LinkId> LinkAxioms::linksSince(size_t first) const {
  std::vector<LinkId> ns;
  for (size_t m = first; m < g_.linkCount(); ++m) ns.push_back(static_cast<LinkId>(m));
  return ns;
}

void LinkAxioms::name(const std::vector<LinkId> &ns, int atDepth) {
  // Naming is most of what a goal run does outside its searches, and each link's
  // name is independent of the others', so the names are found on a pool
  // (LinkNamer is safe to share: its caches are locked, and its SnapPea
  // calls serialised) and applied here in link order, as one at a time would.
  std::vector<std::optional<linknaming::PieceName>> names(ns.size());
  // A knot namePiece() leaves untabulated may be a connected sum the tables
  // hold only summand by summand: the whole-diagram namer cuts it at its
  // visible sum spheres and composes `A#mB` (linknamer.h, step 3).
  std::vector<std::optional<linknaming::LinkName>> composites(ns.size());
  // Any untabulated link (a knot included) with visible sum spheres is cut into
  // its prime summands, which become links joined to it by a sum
  // (applySum; paper lem:sum-partitions): the tables list prime links only,
  // and most outgoing links of small links are sums of smaller ones.
  std::vector<std::vector<GaussDiagram>> summands(ns.size());
  parallelFor(ns.size(), std::max(options_.threads, 1u), [&](size_t i) {
    const LinkId n = ns[i];
    if (!reg_.known(n) || n == reg_.unknot()) return;
    try {
      const GaussDiagram &d = reg_.info(n).diagram;
      names[i] = namer_.namePiece(d);
      if (names[i]->by == linknaming::PieceName::By::untabulated && d.crossings() > 0) {
        if (d.components() == 1) {
          linknaming::LinkName fs = namer_.name(d.link());
          if (fs.isName && fs.pinned && fs.pieces.size() >= 2 &&
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

void LinkAxioms::applySum(LinkId n, const std::vector<GaussDiagram> &primes, int atDepth) {
  // Each prime is interned as a link (a duplicate costs search, never
  // soundness), and its components are mapped to the whole's through the
  // registry's component map and the prime's origins. The new summand
  // links are then named like any other (literature leaves included), so
  // the sum rule can combine their bounds.
  const size_t before = g_.linkCount();
  std::vector<LinkId> pieces;
  std::vector<std::vector<int>> maps;
  for (size_t k = 0; k < primes.size(); ++k) {
    const GaussDiagram &p = primes[k];
    LinkMatch m = reg_.intern(p, "summand " + std::to_string(k) + kFrozenSummandOfNodeLabel +
                                        std::to_string(n));
    std::vector<int> map(p.components(), -1);
    for (size_t c = 0; c < p.components(); ++c)
      map[static_cast<size_t>(m.componentMap[c])] = static_cast<int>(p.origin[c]);
    pieces.push_back(m.link);
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
  name(linksSince(before), atDepth + 1);
  if (!options_.log) return;
  std::ostringstream o;
  o << "[+] node " << n << " is a sum along components of";
  for (LinkId p : pieces) {
    auto t = tableName.find(p);
    o << " node " << p << (t == tableName.end() ? "" : " (" + t->second + ")");
  }
  *options_.log << o.str() << "\n";
}

void LinkAxioms::applyComposite(LinkId n, const linknaming::LinkName &fs) {
  // The composite's name is recorded (certificates, link bounds, the
  // subject name stays cascade:, since no table entry holds it). It is an
  // ANCHOR when its summands cancel in concordance (cobordismgraph.h
  // isElementarySlice: the explicit allowlist, or symmetry types from
  // --knot-symmetry): K # m(K^r) bounds a ribbon disc, so the link gets the
  // unknot's leaf, constructive like the unknot's, never for the target.
  tableName[n] = fs.name;
  if (n == target) return;
  if (!linknaming::isElementarySlice(fs.name, symmetries_)) return;
  g_.addLeaf(n, Partition::coarsest(1), 0, "anchor " + fs.name);
  ++anchors;
  if (options_.log)
    *options_.log << "[+] node " << n << " is " << fs.name << " (" << fs.proof()
                  << "): a slice composite, anchored\n";
}

void LinkAxioms::applyName(LinkId n, const linknaming::PieceName &pn) {
  if (pn.by == linknaming::PieceName::By::untabulated || !pn.pinned() || pn.names.size() != 1)
    return;
  const std::string &name = pn.names.front();
  tableName[n] = name;
  const linknaming::TableEntry *e = tables_.entry(name);
  if (!e) return;
  auto g4 = linknaming::parseTableG4(e->g4);
  if (!g4) return;
  const auto [lo, hi] = *g4;
  g_.setGenusLowerBound(n, lo, "literature " + name + " " + e->g4);
  // Never let the target's own literature value prove the target, even
  // through a duplicate link of it (README.md, "Leaf facts").
  if (!mayUseLiteratureUpperBound(classOf(name), targetClass, options_.literature))
    return;
  g_.addLeaf(n, Partition::coarsest(g_.link(n).components), hi,
             "literature " + name + " " + e->g4);
}

} // namespace bounds
