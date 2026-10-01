//
//  farsidecurves.h
//
//  A cobordism's outgoing link, as oriented curves in knotbuilder's own
//  triangulation of the row.
//

/*! \file utils/surfer/farsidecurves.h
 *  \brief Carries a surface's outgoing boundary curves onto knotbuilder's
 *  triangulation T of the row, oriented as a cobordism from the row's own
 *  oriented link -- ready for knotbuilder::DiagramDrawer.
 *
 *  verifyslicegenus searches S^3 x [0,2] = T x [0,2], thickened from
 *  knotbuilder's triangulation T by CobordismBuilder. Its outgoing boundary
 *  is literally T x {2}: over each base tetrahedron sigma, the top prism
 *  piece P_3(sigma) has sigma x {2} as a facet, with base vertex v's top copy
 *  at SimplicialPrism::localVertex(v, true). OutgoingMap uses exactly that
 *  to carry each edge of the outgoing boundary onto an edge of T. There is
 *  no isomorphism search anywhere on the outgoing side, so none of T's
 *  automorphisms -- some of which reverse orientation -- can relabel or
 *  mirror the far side.
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

#ifndef SURFER_FARSIDECURVES_H
#define SURFER_FARSIDECURVES_H

#include <map>
#include <optional>
#include <unordered_map>
#include <vector>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include "cobordismbuilder.h"
#include "cobordismgraph.h"
#include "embeddedsubmanifold.h"
#include "knotbuilder/diagramdrawer.h"

namespace farside {

/**
 * The outgoing boundary of a thickening, edge by edge, as knotbuilder's T.
 */
class OutgoingMap {
  public:
    /**
     * \param knotT knotbuilder::buildLink()'s triangulation, unmodified.
     * \param cob built from `knotT` (CobordismBuilder takes an ordered copy,
     *        relabelling vertices within tetrahedra), after its last
     *        thicken() and with no cone().
     */
    OutgoingMap(const regina::Triangulation<3> &knotT,
                const CobordismBuilder<3> &cob);

    /** The outgoing boundary component of cob.getCobordism(). */
    size_t boundaryComponent() const { return bc_; }

    /**
     * A curve on the outgoing boundary -- oriented edges of that boundary
     * component's built triangulation, as orientedBoundaryLinks() gives
     * them -- as a directed edge cycle of `knotT`.
     */
    knotbuilder::EdgeCycle carry(const OrientedCurve &curve) const;

    /**
     * A closed curve given as its edges in any order (a boundary link's
     * component, which carries no orientation) as a directed edge cycle of
     * `knotT`, traversed in the direction of `edges.front()`.
     *
     * \exception regina::InvalidArgument the edges do not form one closed
     * curve.
     */
    knotbuilder::EdgeCycle carryCycle(const std::vector<const regina::Edge<3> *> &edges) const;

  private:
    size_t bc_ = 0;
    std::unordered_map<size_t, size_t> edgeToT_;   /**< boundary edge -> T edge */
    std::unordered_map<size_t, size_t> vertexToT_; /**< boundary vertex -> T vertex */
    std::vector<size_t> tTail_;                    /**< T edge -> its vertex(0) */
};

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

/**
 * The edges of `faces` (triangle indices of `tri`) lying in boundary
 * component `bcIndex`, as sorted local edge indices of that component --
 * the numbering its built triangulation uses.
 */
std::vector<size_t> boundaryEdgesOf(const regina::Triangulation<4> &tri,
                                    const std::vector<int> &faces,
                                    size_t bcIndex);

} // namespace farside

#endif
