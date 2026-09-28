// hopedges.cpp

#include "hopedges.h"

#include <map>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

#include <link/link.h>

#include "diagramiso.h"
#include "embeddedsubmanifold.h"
#include "farsidecurves.h"

using exactnaming::GaussDiagram;

namespace cascade {

GaussDiagram gaussOf(const knotbuilder::Diagram &d) {
  GaussDiagram g;
  for (const auto &c : d.crossings) g.signs.push_back(c.sign);
  g.comps = d.gauss;
  g.origin.resize(d.components);
  for (size_t i = 0; i < d.components; ++i) g.origin[i] = i;
  if (g.comps.size() != d.components)
    throw std::logic_error("gaussOf: gauss codes do not cover every component");
  return g;
}

std::string rowPD(const GaussDiagram &d) {
  regina::Link l = d.link();
  if (l.pdAmbiguous()) {
    // Only harmless when each ambiguous component (over at every crossing it
    // meets, no self-crossings) lies above everything: then it is an unknot
    // split from the rest, and reversing it changes nothing. Such a
    // component has linking number 0 with every other; anything else is a
    // PD we cannot trust to carry orientation.
    auto lk = linkingMatrix(d);
    for (size_t c = 0; c < d.components(); ++c) {
      bool allOver = !d.comps[c].empty();
      for (long x : d.comps[c]) if (x < 0) allOver = false;
      if (!allOver) continue;
      for (size_t o = 0; o < d.components(); ++o)
        if (o != c && lk[c][o] != 0)
          throw std::runtime_error("rowPD: an over-everywhere component is linked");
    }
  }
  std::ostringstream o;
  o << '[';
  bool first = true;
  for (const auto &x : l.pdData()) {
    o << (first ? "" : ";") << '[' << x[0] << ';' << x[1] << ';' << x[2] << ';' << x[3] << ']';
    first = false;
  }
  o << ']';
  return o.str();
}

HopAssembler::HopAssembler(ProofGraph &graph, NodeRegistry &nodes, HopRow row,
                           Read read)
    : read_(read), g_(graph), nodes_(nodes), row_(std::move(row)) {
  redraw_ = std::make_unique<farside::WitnessRedrawer>(row_.pd, row_.layers);
  const auto &built = redraw_->built();
  const auto &cycles = redraw_->rowCycles();
  {
    std::unordered_map<size_t, size_t> compOfT;
    for (size_t c = 0; c < cycles.size(); ++c)
      for (const auto &de : cycles[c]) compOfT[de.edge] = c;
    componentOfRowEdge_.resize(built.edges.size());
    for (size_t i = 0; i < built.edges.size(); ++i)
      componentOfRowEdge_[i] = compOfT.at(built.edges[i]->index());
  }
  // Certify the row: knotbuilder's link, drawn back, is row.diagram.
  const GaussDiagram drawn = gaussOf(redraw_->drawer().draw(cycles));
  auto iso = findDiagramIsomorphism(drawn, row_.diagram, /*allowMirror=*/false,
                                    /*allowReverse=*/false);
  if (!iso)
    throw std::runtime_error(
        "HopAssembler: the row's triangulated link does not redraw as its "
        "diagram (no orientation-preserving isomorphism); refusing the row");
  if (row_.nodeMap.size() != row_.diagram.components())
    throw std::invalid_argument("HopAssembler: nodeMap size");
  rowToNode_.resize(drawn.components());
  for (size_t i = 0; i < drawn.components(); ++i)
    rowToNode_[i] = row_.nodeMap[iso->componentMap[i]];
}

std::optional<HopAssembler::ReadBack>
HopAssembler::readBack(const std::string &pairsig, std::string &why) const {
  const size_t n = redraw_->rowCycles().size();
  ReadBack rb;
  rb.surfaceOfRowComponent.assign(n, static_cast<size_t>(-1));
  auto place = [&](size_t firstEdgeIndex, size_t sc) -> bool {
    const size_t rowEdge = redraw_->row().rowIndexOf.at(firstEdgeIndex);
    const size_t rc = componentOfRowEdge_[rowEdge];
    if (rb.surfaceOfRowComponent[rc] != static_cast<size_t>(-1)) {
      why = "a row component met twice";
      return false;
    }
    rb.surfaceOfRowComponent[rc] = sc;
    return true;
  };
  if (read_ == Read::fast) {
    auto link = redraw_->outgoingLinkFast(pairsig, why);
    if (!link) return std::nullopt;
    for (size_t i = 0; i < link->incomingFirstEdge.size(); ++i)
      if (!place(link->incomingFirstEdge[i], link->incomingSurfaceComponent[i]))
        return std::nullopt;
    rb.link = std::move(*link);
  } else {
    auto carried = redraw_->carry(pairsig, why);
    if (!carried) {
      why = "carry: " + why;
      return std::nullopt;
    }
    KnottedSurface surface(redraw_->skeleton(), *carried);
    auto link = farside::orientedOutgoingLink(surface, redraw_->outgoing(),
                                              redraw_->row(), redraw_->incomingBC());
    if (!link) {
      why = "incoming orientation inconsistent";
      return std::nullopt;
    }
    const auto surfaceOf = surface.boundaryEdgeSurfaceComponent();
    for (const auto &[bc, curves] : surface.orientedBoundaryLinks()) {
      if (bc != redraw_->incomingBC()) continue;
      for (const OrientedCurve &curve : curves)
        if (!curve.empty() &&
            !place(curve.front().edge->index(), surfaceOf.at(curve.front().edge)))
          return std::nullopt;
    }
    rb.link = std::move(*link);
  }
  for (size_t sc : rb.surfaceOfRowComponent)
    if (sc == static_cast<size_t>(-1)) {
      why = "a row component is on no surface component";
      return std::nullopt;
    }
  return rb;
}

HopEdge HopAssembler::add(const HopWitness &w) {
  HopEdge out;
  std::string why;
  auto rb = readBack(w.pairsig, why);
  if (!rb) {
    out.why = why;
    return out;
  }
  const farside::OutgoingLink *link = &rb->link;

  // Surface components, renumbered 0..c-1 in order of first appearance.
  std::map<size_t, int> compIndex;
  auto idx = [&](size_t sc) {
    auto [it, ins] = compIndex.try_emplace(sc, static_cast<int>(compIndex.size()));
    return it->second;
  };
  const size_t n = redraw_->rowCycles().size();
  std::vector<int> inComp(n, -1);
  for (size_t rc = 0; rc < n; ++rc)
    inComp[rc] = idx(rb->surfaceOfRowComponent[rc]);
  CobordismShape shape;
  shape.genus = w.genus;
  shape.inComponent.resize(n);
  // Incoming curves in NODE order: node component rowToNode_[rc] <- rc.
  for (size_t rc = 0; rc < n; ++rc)
    shape.inComponent[rowToNode_[rc]] = inComp[rc];
  for (size_t sc : link->surfaceComponent) shape.outComponent.push_back(idx(sc));
  for (const auto &cyc : link->curves) {
    std::vector<size_t> es;
    for (const auto &de : cyc) es.push_back(de.edge);
    std::sort(es.begin(), es.end());
    out.farCurveEdges.push_back(std::move(es));
  }
  shape.components = static_cast<int>(compIndex.size());
  shape.validate();
  out.shape = shape;

  std::vector<int> inMap(n);
  for (size_t i = 0; i < n; ++i) inMap[i] = static_cast<int>(i); // already node order

  if (link->curves.empty()) {
    // A surface bounding the row alone: a leaf for the row's node.
    std::vector<int> labels(n);
    for (size_t i = 0; i < n; ++i) labels[i] = shape.inComponent[i];
    g_.addLeaf(row_.node, Partition::fromLabels(labels), w.genus,
               "direct witness " + w.key);
    out.ok = true;
    out.direct = true;
    return out;
  }

  // The far side, drawn, split into pieces, each simplified and interned.
  const knotbuilder::Diagram d = redraw_->drawer().draw(link->curves);
  const GaussDiagram whole = gaussOf(d);
  const size_t m = whole.components();
  // One pass: simplify() is randomised, so each simplified piece must stay
  // together with its own match and origins.
  std::vector<std::vector<size_t>> pieceOrigins;
  for (const GaussDiagram &p : exactnaming::splitPieces(whole)) {
    const GaussDiagram s = simplifyKeepingComponents(p);
    // simplify() can make a piece split further (a component unlinked by
    // Reidemeister moves): intern each resulting piece separately.
    for (const GaussDiagram &q : exactnaming::splitPieces(s)) {
      out.pieces.push_back(nodes_.intern(q, "far side of " + w.key));
      pieceOrigins.push_back(q.origin);
    }
  }

  std::vector<int> outMap(m, -1);
  if (out.pieces.size() == 1) {
    out.farNode = out.pieces[0].node;
    for (size_t i = 0; i < pieceOrigins[0].size(); ++i)
      outMap[pieceOrigins[0][i]] = out.pieces[0].componentMap[i];
  } else {
    // A split far side: a fresh whole node (never merged: mirroring or
    // reversing ONE piece changes a split link), joined to its pieces.
    out.farNode = g_.addNode(static_cast<int>(m), "split far side of " + w.key,
                             linkingMatrix(whole));
    std::vector<NodeId> pn;
    std::vector<std::vector<int>> pmap;
    for (size_t k = 0; k < out.pieces.size(); ++k) {
      pn.push_back(out.pieces[k].node);
      std::vector<int> mapK(g_.node(out.pieces[k].node).components, -1);
      for (size_t i = 0; i < pieceOrigins[k].size(); ++i)
        mapK[out.pieces[k].componentMap[i]] = static_cast<int>(pieceOrigins[k][i]);
      pmap.push_back(mapK);
    }
    out.splitEdge = g_.addSplit(out.farNode, pn, pmap);
    for (size_t j = 0; j < m; ++j) outMap[j] = static_cast<int>(j);
  }
  out.pieceOrigins = pieceOrigins;
  for (int v : outMap)
    if (v < 0) throw std::logic_error("hop: a far-side curve is in no piece");
  out.edge = g_.addWitness(row_.node, out.farNode, shape, inMap, outMap, w.key);
  out.ok = true;
  return out;
}

} // namespace cascade
