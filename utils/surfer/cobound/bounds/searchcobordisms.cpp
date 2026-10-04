// hopedges.cpp

#include "cobound/bounds/searchcobordisms.h"

#include <map>
#include <sstream>
#include <stdexcept>

#include <link/link.h>

#include "diagramtriangulation/pdcode.h"
#include "linknaming/diagrams/diagramiso.h"
#include "surfer/submanifold/submanifold.h"
#include "cobound/outgoing/outgoinglink.h"
#include "cobound/frozen.h"

using linknaming::GaussDiagram;

namespace bounds {

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

std::string diagramPD(const GaussDiagram &d) {
  regina::Link l = d.link();
  if (l.pdAmbiguous()) {
    // Only harmless when each ambiguous component (over at every crossing it
    // meets, no self-crossings) lies above everything: then it is an unknot
    // split from the rest, and reversing it changes nothing. Such a
    // component has linking number 0 with every other; anything else is a
    // PD we cannot trust to carry orientation.
    auto lk = linknaming::linkingMatrix(d);
    for (size_t c = 0; c < d.components(); ++c) {
      bool allOver = !d.comps[c].empty();
      for (long x : d.comps[c]) if (x < 0) allOver = false;
      if (!allOver) continue;
      for (size_t o = 0; o < d.components(); ++o)
        if (o != c && lk[c][o] != 0)
          throw std::runtime_error("rowPD: an over-everywhere component is linked");
    }
  }
  return knotbuilder::formatPDCode(l.pdData(), knotbuilder::PDSpelling::semicolons);
}

CobordismAssembler::CobordismAssembler(ProofGraph &graph, LinkRegistry &links, SearchedLink searched,
                           Read read)
    : read_(read), g_(graph), links_(links), searched_(std::move(searched)) {
  redraw_ = std::make_unique<outgoing::OutgoingReader>(searched_.pd, searched_.layers);
  certifyIncoming_();
}

CobordismAssembler::CobordismAssembler(ProofGraph &graph, LinkRegistry &links, SearchedLink searched,
                           std::unique_ptr<outgoing::OutgoingReader> built, Read read)
    : read_(read), g_(graph), links_(links), searched_(std::move(searched)), redraw_(std::move(built)) {
  if (!redraw_) throw std::invalid_argument("HopAssembler: no redrawer");
  certifyIncoming_();
}

void CobordismAssembler::certifyIncoming_() {
  const auto &cycles = redraw_->incomingCycles();
  // Certify the row: knotbuilder's link, drawn back, is row.diagram.
  const GaussDiagram drawn = gaussOf(redraw_->drawer().draw(cycles));
  auto iso = linknaming::findDiagramIsomorphism(drawn, searched_.diagram, /*allowMirror=*/false,
                                    /*allowReverse=*/false);
  if (!iso)
    throw std::runtime_error(
        "HopAssembler: the row's triangulated link does not redraw as its "
        "diagram (no orientation-preserving isomorphism); refusing the row");
  if (searched_.linkMap.size() != searched_.diagram.components())
    throw std::invalid_argument("HopAssembler: nodeMap size");
  incomingToLink_.resize(drawn.components());
  for (size_t i = 0; i < drawn.components(); ++i)
    incomingToLink_[i] = searched_.linkMap[iso->componentMap[i]];
}

std::optional<outgoing::OutgoingLink>
CobordismAssembler::readBack(const std::string &pairsig, std::string &why) const {
  return read_ == Read::fast ? redraw_->outgoingLinkFast(pairsig, why)
                            : redraw_->outgoingLink(pairsig, why);
}

std::optional<std::vector<size_t>>
CobordismAssembler::surfaceOfIncomingComponents(const outgoing::OutgoingLink &link,
                                     std::string &why) const {
  std::vector<size_t> of(redraw_->incomingCycles().size(), static_cast<size_t>(-1));
  for (size_t i = 0; i < link.incomingFirstEdge.size(); ++i) {
    const size_t rc = redraw_->incomingComponentOf(link.incomingFirstEdge[i]);
    if (of[rc] != static_cast<size_t>(-1)) {
      why = "a row component met twice";
      return std::nullopt;
    }
    of[rc] = link.incomingSurfaceComponent[i];
  }
  for (size_t sc : of)
    if (sc == static_cast<size_t>(-1)) {
      why = "a row component is on no surface component";
      return std::nullopt;
    }
  return of;
}

AddedCobordism CobordismAssembler::add(const SignedCobordism &w) {
  std::string why;
  auto link = readBack(w.pairsig, why);
  if (!link) {
    AddedCobordism out;
    out.why = why;
    return out;
  }
  return addRead(*link, w.genus, w.key);
}

AddedCobordism CobordismAssembler::addRead(const outgoing::OutgoingLink &read, int genus,
                              const std::string &key) {
  AddedCobordism out;
  std::string why;
  auto surfaceOfIncomingComponent = surfaceOfIncomingComponents(read, why);
  if (!surfaceOfIncomingComponent) {
    out.why = why;
    return out;
  }
  const outgoing::OutgoingLink *link = &read;

  // Surface components, renumbered 0..c-1 in order of first appearance.
  std::map<size_t, int> compIndex;
  auto idx = [&](size_t sc) {
    auto [it, ins] = compIndex.try_emplace(sc, static_cast<int>(compIndex.size()));
    return it->second;
  };
  const size_t n = redraw_->incomingCycles().size();
  std::vector<int> inComp(n, -1);
  for (size_t rc = 0; rc < n; ++rc)
    inComp[rc] = idx((*surfaceOfIncomingComponent)[rc]);
  CobordismShape shape;
  shape.genus = genus;
  shape.inComponent.resize(n);
  // Incoming curves in NODE order: node component rowToNode_[rc] <- rc.
  for (size_t rc = 0; rc < n; ++rc)
    shape.inComponent[incomingToLink_[rc]] = inComp[rc];
  for (size_t sc : link->surfaceComponent) shape.outComponent.push_back(idx(sc));
  for (const auto &cyc : link->curves) {
    std::vector<size_t> es;
    for (const auto &de : cyc) es.push_back(de.edge);
    std::sort(es.begin(), es.end());
    out.outgoingCurveEdges.push_back(std::move(es));
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
    g_.addLeaf(searched_.link, Partition::fromLabels(labels), genus,
               kFrozenDirectWitnessSource + key);
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
  for (const GaussDiagram &p : linknaming::splitPieces(whole)) {
    const GaussDiagram s = linknaming::simplifyKeepingComponents(p);
    // simplify() can make a piece split further (a component unlinked by
    // Reidemeister moves): intern each resulting piece separately.
    for (const GaussDiagram &q : linknaming::splitPieces(s)) {
      out.pieces.push_back(links_.intern(q, kFrozenFarSideLabel + key));
      pieceOrigins.push_back(q.origin);
    }
  }

  std::vector<int> outMap(m, -1);
  if (out.pieces.size() == 1) {
    out.outgoing = out.pieces[0].link;
    for (size_t i = 0; i < pieceOrigins[0].size(); ++i)
      outMap[pieceOrigins[0][i]] = out.pieces[0].componentMap[i];
  } else {
    // A split far side: a fresh whole node (never merged: mirroring or
    // reversing ONE piece changes a split link), joined to its pieces.
    out.outgoing = g_.addLink(static_cast<int>(m), kFrozenSplitFarSideLabel + key,
                             linknaming::linkingMatrix(whole));
    std::vector<LinkId> pn;
    std::vector<std::vector<int>> pmap;
    for (size_t k = 0; k < out.pieces.size(); ++k) {
      pn.push_back(out.pieces[k].link);
      std::vector<int> mapK(g_.link(out.pieces[k].link).components, -1);
      for (size_t i = 0; i < pieceOrigins[k].size(); ++i)
        mapK[out.pieces[k].componentMap[i]] = static_cast<int>(pieceOrigins[k][i]);
      pmap.push_back(mapK);
    }
    out.splitEdge = g_.addSplit(out.outgoing, pn, pmap);
    for (size_t j = 0; j < m; ++j) outMap[j] = static_cast<int>(j);
  }
  out.pieceOrigins = pieceOrigins;
  for (int v : outMap)
    if (v < 0) throw std::logic_error("hop: a far-side curve is in no piece");
  out.edge = g_.addCobordism(searched_.link, out.outgoing, shape, inMap, outMap, key);
  out.ok = true;
  return out;
}

} // namespace bounds
