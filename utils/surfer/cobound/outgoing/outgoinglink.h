//
//  outgoinglink.h
//
//  A cobordism's outgoing link, as oriented curves in knotbuilder's own
//  triangulation of the row.
//

/*! \file utils/surfer/cobound/outgoing/outgoinglink.h
 *  \brief Carries a surface's outgoing boundary curves onto knotbuilder's
 *  triangulation T of the row (OutgoingMap, diagramtriangulation/thickening),
 *  oriented as a cobordism from the row's own oriented link -- ready for
 *  knotbuilder::DiagramDrawer.
 *
 *  Orientation. A surface's boundary is oriented per connected component of
 *  the surface, each component independently
 *  (KnottedSurface::orientedBoundaryLinks()). Every component meets the
 *  incoming boundary (it contains a seed annulus), and the search accepts a
 *  surface only if each component's incoming curves all run with the row's
 *  orientation or all against it (cobordismgraph::classifyRowOrientation()).
 *  Reversing the components that run against it orients the surface as an
 *  oriented cobordism from the row's oriented link L; its outgoing curves
 *  are then the oriented outgoing link, up to reversing every component at
 *  once -- which neither a table name nor a slice genus can see.
 */

#ifndef SURFER_OUTGOINGLINK_H
#define SURFER_OUTGOINGLINK_H

#include <map>
#include <optional>
#include <vector>

#include <triangulation/dim3.h>

#include "diagramtriangulation/thickening/thickening.h"
#include "cobound/search/incoming.h"
#include "surfer/submanifold/submanifold.h"
#include "diagramtriangulation/todiagram.h"

namespace farside {

/** A boundary curve as a search hands it out (KnottedSurface's
 *  OrientedCurve), as the OutgoingMap reads curves. */
OutgoingCurve outgoingCurve(const OrientedCurve &curve);

/**
 * Per surface component, +1 to keep its orientation or -1 to reverse it,
 * so that its incoming curves run as the row's link does. nullopt when some
 * component's incoming curves disagree among themselves, a curve has an
 * edge off the row's link, or a curve's component is unknown -- exactly the
 * surfaces classifyRowOrientation() rejects.
 */
std::optional<std::map<size_t, int>> incomingFlips(
    const cobordismgraph::RowOrientation &row,
    const std::vector<OrientedCurve> &incomingCurves,
    const std::map<const regina::Edge<3> *, size_t> &surfaceComponentOf);

/** A surface's outgoing link, oriented as a cobordism from the row's link. */
struct OutgoingLink {
    std::vector<knotbuilder::EdgeCycle> curves; /**< In knotbuilder's T. */
    std::vector<size_t> surfaceComponent;       /**< Per curve. */
    /** The incoming side, per incoming curve: the index of its first edge in
     *  the incoming boundary component's built triangulation, and the surface
     *  component it lies on (the cascade needs both ends). */
    std::vector<size_t> incomingFirstEdge;
    std::vector<size_t> incomingSurfaceComponent;
};

/**
 * `surface`'s outgoing curves, carried onto T by `map` and oriented per
 * incomingFlips() against the row (`row`, built on boundary component
 * `incomingBC`). nullopt when incomingFlips() is.
 */
std::optional<OutgoingLink> orientedOutgoingLink(
    const KnottedSurface &surface, const OutgoingMap &map,
    const cobordismgraph::RowOrientation &row, size_t incomingBC);

/**
 * As above, from a surface's boundary as a search hands it out
 * (SurfaceBoundaryInfo::captureOrientedBoundaryLinks() and
 * captureBoundaryEdgeSurfaceComponent()) rather than a KnottedSurface.
 */
std::optional<OutgoingLink> orientedOutgoingLink(
    const std::vector<std::pair<size_t, std::vector<OrientedCurve>>> &oriented,
    const std::map<const regina::Edge<3> *, size_t> &surfaceOf,
    const OutgoingMap &map, const cobordismgraph::RowOrientation &row,
    size_t incomingBC);

} // namespace farside

#endif
