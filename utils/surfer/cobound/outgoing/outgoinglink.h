//
//  outgoinglink.h
//
//  A cobordism's outgoing link, as oriented curves in the diagram's own
//  triangulation T of the incoming link.
//

/*! \file utils/surfer/cobound/outgoing/outgoinglink.h
 *  \brief Carries a surface's outgoing boundary curves onto the diagram's
 *  triangulation T of the incoming link (OutgoingMap, diagramtriangulation/thickening),
 *  oriented as a cobordism from the oriented incoming link -- ready for
 *  diagramtriangulation::DiagramDrawer.
 *
 *  Orientation. A surface's boundary is oriented per connected component of
 *  the surface, each component independently
 *  (KnottedSurface::orientedBoundaryLinks()). Every component meets the
 *  incoming boundary (it contains a seed annulus), and the search accepts a
 *  surface only if each component's incoming curves all run with the incoming link's
 *  orientation or all against it (search::classifyIncomingOrientation()).
 *  Reversing the components that run against it orients the surface as an
 *  oriented cobordism from the oriented incoming link L; its outgoing curves
 *  are then the oriented outgoing link, up to reversing every component at
 *  once -- which neither a table name nor a slice genus can see.
 */

#ifndef SURFER_OUTGOINGLINK_H
#define SURFER_OUTGOINGLINK_H

#include <map>
#include <string>
#include <optional>
#include <vector>

#include <triangulation/dim3.h>

#include "diagramtriangulation/thickening/thickening.h"
#include "cobound/search/incoming.h"
#include "surfer/submanifold/submanifold.h"
#include "diagramtriangulation/todiagram.h"

namespace outgoing {

/** A boundary curve as a search hands it out (KnottedSurface's
 *  OrientedCurve), as the OutgoingMap reads curves. */
OutgoingCurve outgoingCurve(const OrientedCurve &curve);

/**
 * Per surface component, +1 to keep its orientation or -1 to reverse it,
 * so that its incoming curves run as the incoming link does. nullopt when some
 * component's incoming curves disagree among themselves, a curve has an
 * edge off the incoming link, or a curve's component is unknown -- exactly the
 * surfaces classifyIncomingOrientation() rejects; an empty map for no curves.
 * (search::judgeIncomingOrientation()'s consistentFlips(): one walk.)
 */
std::optional<std::map<size_t, int>> incomingFlips(
    const search::IncomingOrientation &incoming,
    const std::vector<OrientedCurve> &incomingCurves,
    const std::map<const regina::Edge<3> *, size_t> &surfaceComponentOf);

/** A surface's outgoing link, oriented as a cobordism from the incoming link. */
struct OutgoingLink {
    std::vector<diagramtriangulation::EdgeCycle> curves; /**< In the diagram's triangulation T. */
    std::vector<size_t> surfaceComponent;       /**< Per curve. */
    /** The incoming side, per incoming curve: the index of its first edge in
     *  the incoming boundary component's built triangulation, and the surface
     *  component it lies on (goal runs need both ends). */
    std::vector<size_t> incomingFirstEdge;
    std::vector<size_t> incomingSurfaceComponent;
};

/**
 * `surface`'s outgoing curves, carried onto T by `map` and oriented per
 * incomingFlips() against the incoming link (`incoming`, built on boundary component
 * `incomingBC`). nullopt when incomingFlips() is.
 */
std::optional<OutgoingLink> orientedOutgoingLink(
    const KnottedSurface &surface, const OutgoingMap &map,
    const search::IncomingOrientation &incoming, size_t incomingBC);

/**
 * As above, from a surface's boundary as a search hands it out
 * (SurfaceBoundaryInfo::captureOrientedBoundaryLinks() and
 * captureBoundaryEdgeSurfaceComponent()) rather than a KnottedSurface.
 */
std::optional<OutgoingLink> orientedOutgoingLink(
    const std::vector<std::pair<size_t, std::vector<OrientedCurve>>> &oriented,
    const std::map<const regina::Edge<3> *, size_t> &surfaceOf,
    const OutgoingMap &map, const search::IncomingOrientation &incoming,
    size_t incomingBC, std::string *why = nullptr);

/**
 * As above, with the incoming curves' flips already judged (incomingFlips(),
 * or a search's gate, search::GatedSurface::flips). nullopt ("a surface
 * component misses the row") when an outgoing curve's component has none.
 */
std::optional<OutgoingLink> orientedOutgoingLink(
    const std::vector<std::pair<size_t, std::vector<OrientedCurve>>> &oriented,
    const std::map<const regina::Edge<3> *, size_t> &surfaceOf,
    const OutgoingMap &map, const std::map<size_t, int> &flips, size_t incomingBC,
    std::string *why = nullptr);

} // namespace outgoing

#endif
