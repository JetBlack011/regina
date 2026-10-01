//
//  preconditions.h
//
//  What a found surface must satisfy to witness anything about its row.
//

#ifndef SURFER_COBOUND_PRECONDITIONS_H
#define SURFER_COBOUND_PRECONDITIONS_H

#include <map>
#include <string>
#include <vector>

#include <triangulation/dim3.h>

#include "cobound/search/incoming.h"
#include "surfer/enumeration/surfacesearch.h"

/*! \file utils/surfer/cobound/search/preconditions.h
 *  \brief Per found surface: its boundary split into the incoming side and
 *  the others, by geometry, and its incoming curves' orientation against the
 *  row's.
 */

namespace cobordismgraph {

/* Boundary classification */

/**
 * One non-search-side ambient boundary component's identity, as classified
 * by splitBoundary().
 */
struct BoundarySide {
    std::string name;
    int components = 1;
};

/**
 * Splits a SurfaceBoundaryInfo::boundaryComponents grouping into "the
 * search side" (this row's own side) and every other ambient boundary
 * component the surface touches.
 */
struct BoundarySplit {
    size_t searchCurveCount = 0; // 0 if the search side has no boundary here
    bool searchSideRejected = false;
    /**< Set when `requiredSearchEdges` was given and the curves on
         searchSideBC are not exactly those edges -- the surface's boundary
         there is some other link, so it says nothing about this row. */
    bool unnamedSide = false;
    /**< Set when a non-search-side component carried no name at all, which
         describeBoundary_() never produces; the caller treats it as a bug. */
    std::vector<BoundarySide> otherSides;
};

/**
 * The search side is identified by geometry, never by name.
 *
 * Component `searchSideBC` is the search side. In a seeded search that is
 * all there is to it: the seed is L x {0} and no other triangle with an edge
 * in that boundary component is ever searchable, so its curves are L by
 * construction (verifyslicegenus asserts this once per row). An unseeded
 * search has no such guarantee, and passes `requiredSearchEdges` -- the
 * sorted boundary-triangulation edge indices of the row's own link -- so
 * that a surface whose search-side curves are any other edge set is
 * rejected. Identified names are deliberately not consulted: they are not
 * canonical (a census hit's "#N" varies from one identification to the
 * next), and comparing them once silently discarded whole rows.
 */
BoundarySplit
splitBoundary(const std::vector<BoundaryComponentNames> &boundaryComponents,
              size_t searchSideBC,
              const std::vector<size_t> *requiredSearchEdges = nullptr);

/** How a surface's search-side boundary compares with the row's orientation;
 * see classifyRowOrientation(). */
enum class OrientationVerdict {
    match,        /**< Some choice of orientation on each surface component
                       induces the row's own orientation. */
    mismatch,     /**< Some surface component's curves induce an orientation
                       pattern no flip of that component can fix: the surface
                       witnesses a different oriented variant of the link. */
    foreignEdge,  /**< A search-side edge is not one of the row's own. */
    incoherentCurve, /**< A single curve's edges disagree in direction, or a
                          curve's surface component is unknown. */
};

/**
 * Compares one found surface's induced boundary orientation on the
 * search side with `row`.
 *
 * `curves` come from KnottedSurface::orientedBoundaryLinks(), which orients
 * each CONNECTED COMPONENT of the surface independently and arbitrarily;
 * `surfaceComponentOf` (KnottedSurface::boundaryEdgeSurfaceComponent())
 * says which component each edge belongs to. So the curves are grouped by
 * surface component, and each group must match `row` all at once or be
 * reversed all at once; different groups are independent. That is exactly
 * the freedom the surface has: its components can be oriented separately
 * and still tubed into one oriented surface, because a split two-component
 * unlink bounds an oriented annulus for either relative orientation.
 *
 * `foreignEdge` and `incoherentCurve` cannot happen for a correctly built
 * row in a seeded search; the caller treats them as bugs.
 */
OrientationVerdict classifyRowOrientation(
    const RowOrientation &row, const std::vector<OrientedCurve> &curves,
    const std::map<const regina::Edge<3> *, size_t> &surfaceComponentOf);

} // namespace cobordismgraph

#endif // SURFER_COBOUND_PRECONDITIONS_H
