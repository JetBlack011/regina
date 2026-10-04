// hopedges.h
//
// One hop's witnesses, turned into proof-graph edges. See README.md,
// "Composing hops".

#pragma once

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "linknaming/diagrams/gaussdiagram.h"
#include "cobound/outgoing/fromdatabase.h"
#include "cobound/bounds/links.h"
#include "cobound/bounds/cobordismgraph.h"

namespace bounds {

/// What a hop searched: a node, as the diagram written into the row's PD.
struct SearchedLink {
  LinkId link = -1;
  /// The diagram the row was built from (its PD), components in its order.
  linknaming::GaussDiagram diagram;
  /// diagram component i is the node's component nodeMap[i].
  std::vector<int> linkMap;
  std::string pd;   ///< as written to the row, `;`-separated
  int layers = 2;   ///< thicken_layers of the witnesses
};

struct SignedCobordism {
  std::string pairsig;
  int genus = 0; ///< the witness's tubed genus (cobordisms.csv `genus`)
  std::string key; ///< provenance (witness key)
};

struct AddedCobordism {
  bool ok = false;
  std::string why;          ///< when !ok
  bool direct = false;      ///< no far side: the surface bounds the row alone
  EdgeId edge = -1;         ///< the witness edge (when not direct)
  LinkId outgoing = -1;      ///< the far side's node (a split whole or a piece)
  std::vector<LinkMatch> pieces; ///< its split pieces, interned
  /// Per piece: which far-side curves (drawn order) its components are.
  std::vector<std::vector<size_t>> pieceOrigins;
  EdgeId splitEdge = -1;    ///< the split edge, when the far side is split
  /// Per far-side curve (in this read's order): its edges of T, sorted.
  /// Curve order is a property of the read, not of the cobordism; tests use
  /// these to compare two reads up to relabelling the curves.
  std::vector<std::vector<size_t>> outgoingCurveEdges;
  CobordismShape shape;
};

/**
 * Turns witnesses of one hop row into edges of the proof graph.
 *
 * Construction certifies the row: the row's own link, drawn from knotbuilder's
 * triangulation of `row.pd`, must be isomorphic as a diagram to `row.diagram`
 * (orientation kept, no mirror). That isomorphism is the map from the
 * witnesses' incoming curves (knotbuilder's component order) to the node's
 * components, and it is the per-node certificate that the triangulation
 * searched carries the node's oriented link. Throws std::runtime_error if it
 * does not exist: a hop whose row cannot be certified contributes nothing.
 */
class CobordismAssembler {
public:
  /// How a witness is read back from its pair signature. `fast` is
  /// WitnessRedrawer::outgoingLinkFast() (no re-run of the embeddedness
  /// checks a stored witness already passed; milliseconds); `reference`
  /// rebuilds a KnottedSurface (~1 s). hopedges_test checks they agree.
  enum class Read { fast, reference };

  CobordismAssembler(ProofGraph &graph, LinkRegistry &links, SearchedLink searched,
               Read read = Read::fast);
  /// As above, with the row's redrawer already built (for row.pd and
  /// row.layers): master rows are built on worker threads, then assembled.
  CobordismAssembler(ProofGraph &graph, LinkRegistry &links, SearchedLink searched,
               std::unique_ptr<outgoing::OutgoingReader> built, Read read = Read::fast);

  /// A stored witness: read back from its pair signature, then addRead().
  AddedCobordism add(const SignedCobordism &w);

  /// A surface already read: its oriented far side and incoming side, as
  /// outgoing::orientedOutgoingLink() gives them for a surface in
  /// redrawer().thickening() (an in-process hop's search; see hoprunner.h).
  AddedCobordism addRead(const outgoing::OutgoingLink &link, int genus,
                  const std::string &key);

  /// knotbuilder's row component i is the node's component rowToNode()[i].
  const std::vector<int> &incomingToLink() const { return incomingToLink_; }

  /// The row's thickening and everything read from it.
  const outgoing::OutgoingReader &redrawer() const { return *redraw_; }

private:
  std::optional<outgoing::OutgoingLink> readBack(const std::string &pairsig,
                                                std::string &why) const;
  /// Per knotbuilder row component, the surface component it lies on.
  std::optional<std::vector<size_t>>
  surfaceOfIncomingComponents(const outgoing::OutgoingLink &link, std::string &why) const;
  void certifyIncoming_();
  Read read_;
  ProofGraph &g_;
  LinkRegistry &links_;
  SearchedLink searched_;
  std::unique_ptr<outgoing::OutgoingReader> redraw_;
  std::vector<int> incomingToLink_;
};

/// The GaussDiagram of a drawn diagram, `origin` = drawn component index.
linknaming::GaussDiagram gaussOf(const knotbuilder::Diagram &d);

/// A PD code for a connected diagram, as a row's `PD Notation` (`;`-separated,
/// labels from 1). Throws if the PD would not fix every orientation
/// (regina::Link::pdAmbiguous()), unless every ambiguous component is split
/// from the rest, when orientation there cannot matter.
std::string diagramPD(const linknaming::GaussDiagram &d);

} // namespace bounds
