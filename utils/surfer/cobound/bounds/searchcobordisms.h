// searchcobordisms.h
//
// One search's cobordisms, turned into cobordisms of the graph. See README.md,
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

/// The link a search searched, as the diagram its PD was written from.
struct SearchedLink {
  LinkId link = -1;
  /// The diagram the thickening was built from (its PD), components in its order.
  linknaming::GaussDiagram diagram;
  /// diagram component i is the link's component linkMap[i].
  std::vector<int> linkMap;
  std::string pd;   ///< as the search was given it, `;`-separated
  int layers = 2;   ///< thicken_layers of the cobordisms
};

struct SignedCobordism {
  std::string pairsig;
  int genus = 0; ///< the cobordism's tubed genus (cobordisms.csv `genus`)
  std::string key; ///< provenance (cobordism key)
};

struct AddedCobordism {
  bool ok = false;
  std::string why;          ///< when !ok
  bool direct = false;      ///< no outgoing link: the surface bounds the searched link alone
  RelationId cobordism = -1;         ///< the cobordism (when not direct)
  LinkId outgoing = -1;      ///< the outgoing link as a graph link (a split whole or a piece)
  std::vector<LinkMatch> pieces; ///< its split pieces, interned
  /// Per piece: which outgoing curves (drawn order) its components are.
  std::vector<std::vector<size_t>> pieceOrigins;
  RelationId split = -1;    ///< the split, when the outgoing link is split
  /// Per outgoing curve (in this read's order): its edges of T, sorted.
  /// Curve order is a property of the read, not of the cobordism; tests use
  /// these to compare two reads up to relabelling the curves.
  std::vector<std::vector<size_t>> outgoingCurveEdges;
  CobordismShape shape;
};

/**
 * Turns cobordisms of one search into cobordisms of the graph.
 *
 * Construction certifies the incoming link: the link drawn from knotbuilder's
 * triangulation of `searched.pd` must be isomorphic as a diagram to
 * `searched.diagram`
 * (orientation kept, no mirror). That isomorphism is the map from the
 * cobordisms' incoming curves (knotbuilder's component order) to the link's
 * components, and it is the per-link certificate that the triangulation
 * searched carries the link's oriented link. Throws std::runtime_error if it
 * does not exist: a search whose incoming link cannot be certified contributes nothing.
 */
class CobordismAssembler {
public:
  /// How a cobordism is read back from its pair signature. `fast` is
  /// OutgoingReader::outgoingLinkFast() (no re-run of the embeddedness
  /// checks a stored cobordism already passed; milliseconds); `reference`
  /// rebuilds a KnottedSurface (~1 s). searchcobordisms_test checks they agree.
  enum class Read { fast, reference };

  CobordismAssembler(CobordismGraph &graph, LinkRegistry &links, SearchedLink searched,
               Read read = Read::fast);
  /// As above, with the redrawer already built (for searched.pd and
  /// searched.layers): master subjects are built on worker threads, then assembled.
  CobordismAssembler(CobordismGraph &graph, LinkRegistry &links, SearchedLink searched,
               std::unique_ptr<outgoing::OutgoingReader> built, Read read = Read::fast);

  /// A stored cobordism: read back from its pair signature, then addRead().
  AddedCobordism add(const SignedCobordism &w);

  /// A surface already read: its oriented outgoing link and incoming side, as
  /// outgoing::orientedOutgoingLink() gives them for a surface in
  /// redrawer().thickening() (an in-process search; see search/search.h).
  AddedCobordism addRead(const outgoing::OutgoingLink &link, int genus,
                  const std::string &key);

  /// knotbuilder's incoming component i is the link's component incomingToLink()[i].
  const std::vector<int> &incomingToLink() const { return incomingToLink_; }

  /// The searched link's thickening and everything read from it.
  const outgoing::OutgoingReader &redrawer() const { return *redraw_; }

private:
  std::optional<outgoing::OutgoingLink> readBack(const std::string &pairsig,
                                                std::string &why) const;
  /// Per knotbuilder incoming component, the surface component it lies on.
  std::optional<std::vector<size_t>>
  surfaceOfIncomingComponents(const outgoing::OutgoingLink &link, std::string &why) const;
  void certifyIncoming_();
  Read read_;
  CobordismGraph &g_;
  LinkRegistry &links_;
  SearchedLink searched_;
  std::unique_ptr<outgoing::OutgoingReader> redraw_;
  std::vector<int> incomingToLink_;
};

/// The GaussDiagram of a drawn diagram, `origin` = drawn component index.
linknaming::GaussDiagram gaussOf(const knotbuilder::Diagram &d);

/// A PD code for a connected diagram, as a table row's `PD Notation` (`;`-separated,
/// labels from 1). Throws if the PD would not fix every orientation
/// (regina::Link::pdAmbiguous()), unless every ambiguous component is split
/// from the rest, when orientation there cannot matter.
std::string diagramPD(const linknaming::GaussDiagram &d);

} // namespace bounds
