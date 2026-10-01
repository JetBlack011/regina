// nodes.cpp

#include "cobound/bounds/nodes.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

#include <link/link.h>

#include "linknaming/diagrams/diagramiso.h"

using exactnaming::GaussDiagram;
using exactnaming::KernelLink;

namespace cascade {

std::vector<std::vector<int>> linkingMatrix(const GaussDiagram &d) {
  const size_t m = d.components(), n = d.crossings();
  // For each crossing, the components of its over and under strands.
  std::vector<long> over(n, -1), under(n, -1);
  for (size_t c = 0; c < m; ++c)
    for (long x : d.comps[c])
      (x > 0 ? over : under)[std::labs(x) - 1] = static_cast<long>(c);
  std::vector<std::vector<int>> lk(m, std::vector<int>(m, 0));
  for (size_t k = 0; k < n; ++k)
    if (over[k] != under[k] && over[k] >= 0 && under[k] >= 0) {
      lk[over[k]][under[k]] += d.signs[k];
      lk[under[k]][over[k]] += d.signs[k];
    }
  for (auto &row : lk)
    for (int &v : row) {
      if (v % 2 != 0)
        throw std::logic_error("linkingMatrix: odd crossing sum");
      v /= 2;
    }
  return lk;
}

namespace {

// The crossings on the far side of nugatory-candidate crossing `c` from the
// rest of the diagram: everything reachable from the arc c's component runs
// between its two visits to c, without passing through c. `reachesBack` is
// set when that side reaches the component's other arc (so c is not
// nugatory). `c` must be a self-crossing.
std::vector<char> sideOf(const GaussDiagram &d, size_t c, bool &reachesBack) {
  const size_t n = d.crossings();
  // Neighbours along every component, c excluded.
  std::vector<std::vector<size_t>> adj(n);
  for (const auto &w : d.comps)
    for (size_t i = 0; i < w.size(); ++i) {
      const size_t a = std::labs(w[i]) - 1, b = std::labs(w[(i + 1) % w.size()]) - 1;
      if (a == c || b == c || a == b) continue;
      adj[a].push_back(b);
      adj[b].push_back(a);
    }
  // c's component and its two visits.
  size_t comp = 0, i = 0, j = 0;
  bool found = false;
  for (size_t k = 0; k < d.components() && !found; ++k) {
    std::vector<size_t> at;
    for (size_t p = 0; p < d.comps[k].size(); ++p)
      if (static_cast<size_t>(std::labs(d.comps[k][p]) - 1) == c) at.push_back(p);
    if (at.size() == 2) { comp = k; i = at[0]; j = at[1]; found = true; }
    else if (!at.empty()) throw std::logic_error("sideOf: not a self-crossing");
  }
  if (!found) throw std::logic_error("sideOf: crossing on no component");
  const auto &w = d.comps[comp];
  std::vector<char> side(n, 0), other(n, 0);
  std::vector<size_t> stack;
  for (size_t p = i + 1; p < j; ++p) {
    const size_t x = std::labs(w[p]) - 1;
    if (!side[x]) { side[x] = 1; stack.push_back(x); }
  }
  for (size_t p = j + 1; p < w.size() + i; ++p)
    other[std::labs(w[p % w.size()]) - 1] = 1;
  reachesBack = false;
  while (!stack.empty()) {
    const size_t x = stack.back();
    stack.pop_back();
    if (other[x]) reachesBack = true;
    for (size_t y : adj[x])
      if (!side[y]) { side[y] = 1; stack.push_back(y); }
  }
  return side;
}

bool isSelfCrossing(const GaussDiagram &d, size_t c) {
  for (const auto &w : d.comps) {
    int hits = 0;
    for (long x : w) if (static_cast<size_t>(std::labs(x) - 1) == c) ++hits;
    if (hits == 2) return true;
    if (hits == 1) return false;
  }
  return false;
}

} // namespace

std::optional<size_t> nugatoryCrossing(const GaussDiagram &d) {
  for (size_t c = 0; c < d.crossings(); ++c) {
    if (!isSelfCrossing(d, c)) continue;
    bool reachesBack = false;
    sideOf(d, c, reachesBack);
    if (!reachesBack) return c;
  }
  return std::nullopt;
}

GaussDiagram removeNugatoryCrossings(GaussDiagram d) {
  while (auto c = nugatoryCrossing(d)) {
    bool reachesBack = false;
    const std::vector<char> side = sideOf(d, *c, reachesBack);
    // Turn the side c cuts off over (a half turn about an axis in the plane
    // through c): c goes, every crossing on that side swaps over and under,
    // and every sign stays (a rigid motion keeps each crossing's handedness).
    GaussDiagram r;
    r.origin = d.origin;
    for (size_t k = 0; k < d.crossings(); ++k)
      if (k != *c) r.signs.push_back(d.signs[k]);
    for (const auto &w : d.comps) {
      std::vector<long> nw;
      for (long x : w) {
        const size_t k = std::labs(x) - 1;
        if (k == *c) continue;
        const long renumbered = static_cast<long>(k < *c ? k : k - 1) + 1;
        const bool over = (x > 0) != (side[k] != 0);
        nw.push_back(over ? renumbered : -renumbered);
      }
      r.comps.push_back(std::move(nw));
    }
    d = std::move(r);
  }
  return d;
}

GaussDiagram liftSplitComponents(GaussDiagram d) {
  for (;;) {
    // A component over (or under) at every crossing it meets has no
    // self-crossing, so it is an unknot lying above (below) everything else:
    // lifting it off is an isotopy, and it becomes a split unknot.
    std::optional<size_t> lift;
    for (size_t c = 0; c < d.components() && !lift; ++c) {
      const auto &w = d.comps[c];
      if (w.empty()) continue;
      const bool over = std::all_of(w.begin(), w.end(), [](long x) { return x > 0; });
      const bool under = std::all_of(w.begin(), w.end(), [](long x) { return x < 0; });
      if (over || under) lift = c;
    }
    if (!lift) return d;
    std::vector<char> drop(d.crossings(), 0);
    for (long x : d.comps[*lift]) drop[std::labs(x) - 1] = 1;
    std::vector<long> renumbered(d.crossings(), 0);
    GaussDiagram r;
    r.origin = d.origin;
    for (size_t k = 0; k < d.crossings(); ++k)
      if (!drop[k]) {
        r.signs.push_back(d.signs[k]);
        renumbered[k] = static_cast<long>(r.signs.size());
      }
    for (const auto &w : d.comps) {
      std::vector<long> nw;
      for (long x : w) {
        const size_t k = std::labs(x) - 1;
        if (!drop[k]) nw.push_back(x > 0 ? renumbered[k] : -renumbered[k]);
      }
      r.comps.push_back(std::move(nw));
    }
    d = std::move(r);
  }
}

GaussDiagram simplifyKeepingComponents(const GaussDiagram &d) {
  regina::Link l = d.link();
  l.simplify();
  // Reduced as well: knotbuilder's drawer cannot draw a diagram with a
  // nugatory crossing back (a block's corners meet), so a hop on one could
  // not be certified. Regina's simplify() removes kinks but not every
  // nugatory crossing. And with every component that lies above (or below)
  // everything lifted off: its PD code cannot carry its orientation, so its
  // row could not be certified either; lifted, it is a split unknot.
  GaussDiagram s = liftSplitComponents(
      removeNugatoryCrossings(liftSplitComponents(GaussDiagram::of(l, d.origin))));
  if (s.components() != d.components())
    throw std::logic_error("simplify changed the number of components");
  if (linkingMatrix(s) != linkingMatrix(d))
    throw std::logic_error(
        "simplify changed a pairwise linking number: component order or "
        "orientation was not preserved");
  return s;
}

std::string NodeRegistry::diagramKey(const GaussDiagram &d) {
  std::vector<size_t> lengths;
  for (const auto &w : d.comps) lengths.push_back(w.size());
  std::sort(lengths.begin(), lengths.end());
  long writhe = 0;
  for (int s : d.signs) writhe += s;
  std::ostringstream o;
  o << d.components() << '/' << d.crossings() << '/';
  // Writhe up to mirror (the isomorphism test allows mirrors).
  o << std::labs(writhe) << '/';
  for (size_t l : lengths) o << l << ',';
  return o.str();
}

NodeRegistry::NodeRegistry(ProofGraph &graph) : g_(graph) {}

NodeId NodeRegistry::unknot() {
  if (unknot_ < 0) {
    unknot_ = g_.addNode(1, "unknot", std::vector<std::vector<int>>{{0}});
    NodeInfo ni;
    ni.diagram.comps = {{}};
    ni.diagram.origin = {0};
    ni.linking = {{0}};
    info_[unknot_] = ni;
    g_.addLeaf(unknot_, Partition::coarsest(1), 0, "unknot: bounds a disc");
  }
  return unknot_;
}

NodeMatch NodeRegistry::intern(const GaussDiagram &piece, const std::string &label) {
  ++stats_.lookups;
  if (piece.components() == 0)
    throw std::invalid_argument("intern: empty diagram");
  if (exactnaming::splitPieces(piece).size() != 1)
    throw std::invalid_argument("intern: not one split piece");

  // Crossingless: the unknot.
  if (piece.crossings() == 0) {
    if (piece.components() != 1)
      throw std::logic_error("intern: a crossingless piece has one component");
    ++stats_.unknots;
    return NodeMatch{unknot(), {0}, false, false, false, "unknot"};
  }

  // 1. The same diagram.
  const std::string key = diagramKey(piece);
  for (auto [it, end] = byDiagramKey_.equal_range(key); it != end; ++it) {
    const NodeId n = it->second;
    if (auto iso = findDiagramIsomorphism(piece, info_.at(n).diagram,
                                          /*allowMirror=*/true, /*allowReverse=*/true)) {
      ++stats_.diagramHits;
      return NodeMatch{n, iso->componentMap, iso->mirrored, iso->reversed, false,
                       "diagram"};
    }
  }

  // 2. For hyperbolic pieces: the same link by an isometry carrying meridians.
  auto kl = std::make_unique<KernelLink>(piece.link());
  const size_t m = piece.components();
  auto tryIsometry = [&](const KernelLink &k) -> std::optional<NodeMatch> {
    const long long bucket = std::llround(k.volume() * 1e6);
    for (long long b = bucket - 2; b <= bucket + 2; ++b)
      for (auto [it, end] = byVolume_.equal_range({m, b}); it != end; ++it) {
        const NodeId n = it->second;
        for (const KernelLink::Meridional &mi : k.meridionalIsometriesTo(*kernel_.at(n)))
          if (mi.uniform()) {
            NodeMatch nm;
            nm.node = n;
            nm.componentMap = mi.image;
            nm.mirrored = mi.reflects;
            nm.reversed = !mi.sign.empty() && mi.sign.front() < 0;
            nm.method = "isometry";
            return nm;
          }
      }
    return std::nullopt;
  };
  auto hasVolumeCandidate = [&](const KernelLink &k) {
    const long long bucket = std::llround(k.volume() * 1e6);
    for (long long b = bucket - 2; b <= bucket + 2; ++b)
      if (byVolume_.count({m, b})) return true;
    return false;
  };
  if (kl->hyperbolic()) {
    if (auto nm = tryIsometry(*kl)) {
      ++stats_.isometryHits;
      return *nm;
    }
    // A miss is not an answer when a candidate of that volume exists: the
    // canonisation can go astray. Retry once from a randomised triangulation.
    if (hasVolumeCandidate(*kl)) {
      ++stats_.isometryRetries;
      KernelLink retry(piece.link(), /*randomisations=*/1);
      if (retry.hyperbolic())
        if (auto nm = tryIsometry(retry)) {
          ++stats_.isometryHits;
          return *nm;
        }
    }
  }

  // 3. A new node.
  NodeInfo ni;
  ni.diagram = piece;
  ni.linking = linkingMatrix(piece);
  ni.hyperbolic = kl->hyperbolic();
  ni.volume = kl->volume();
  const NodeId n = g_.addNode(static_cast<int>(m), label, ni.linking);
  info_[n] = ni;
  byDiagramKey_.emplace(key, n);
  if (ni.hyperbolic)
    byVolume_.emplace(std::make_pair(m, std::llround(ni.volume * 1e6)), n);
  kernel_[n] = std::move(kl);
  ++stats_.created;
  NodeMatch nm;
  nm.node = n;
  nm.componentMap.resize(m);
  for (size_t i = 0; i < m; ++i) nm.componentMap[i] = static_cast<int>(i);
  nm.created = true;
  nm.method = "new";
  return nm;
}

} // namespace cascade
