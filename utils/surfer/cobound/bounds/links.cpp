// links.cpp

#include "cobound/bounds/links.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

#include <link/link.h>

#include "linknaming/diagrams/diagramiso.h"

using linknaming::GaussDiagram;
using linknaming::KernelLink;

namespace bounds {

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
  if (linknaming::splitPieces(piece).size() != 1)
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
    if (auto iso = linknaming::findDiagramIsomorphism(piece, info_.at(n).diagram,
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
  ni.linking = linknaming::linkingMatrix(piece);
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

} // namespace bounds
