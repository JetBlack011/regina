//
//  incoming.h
//
//  The incoming link of a search, in the search's own terms.
//

#ifndef SURFER_COBOUND_INCOMING_H
#define SURFER_COBOUND_INCOMING_H

#include <unordered_map>
#include <vector>

#include <triangulation/dim3.h>

/*! \file utils/surfer/cobound/search/incoming.h
 *  \brief The incoming link of a search -- the row's own link L, on the
 *  incoming boundary T x {0} -- in the search's own terms: its edges and its
 *  PD orientation there (the row map).
 */

namespace cobordismgraph {

/**
 * A row's own PD-tagged diagram edges (knotbuilder::TriangulationWithLink's
 * `edges`/`reversed`), translated into directed pairs of *vertex indices*
 * within some triangulation combinatorially isomorphic to the row's own
 * knotbuilder triangulation -- normally the ambient search-side boundary
 * component, rebuilt via `BoundaryComponent<4>::build()` (see
 * buildRowOrientation()).
 */
struct RowOrientation {
    std::unordered_map<size_t, size_t> tailOf;
    /**< Edge index of L in the search-side triangulation -> index of that
         edge's tail vertex under the row's PD orientation. Keyed by edge,
         not by vertex pair, so two edges joining the same pair of vertices
         can never be confused. */
    std::vector<size_t> edges; /**< Sorted keys of tailOf: L's edge set. */
    std::unordered_map<size_t, size_t> rowIndexOf;
    /**< Edge index of L in the search-side triangulation -> that edge's
         position in the rowEdges given to buildRowOrientation(), i.e. which
         edge of the row's own link it is. Lets a caller tell which
         component of L a search-side curve is. */
    size_t components = 0;
    /**< How many closed curves `edges` forms, each checked to chain head to
         tail under the PD orientation (buildRowOrientation() throws
         otherwise). */
    bool divergedFromDefaultIsomorphism = false;
    /**< Whether the isomorphism isIsomorphicTo() would have returned maps L
         differently -- onto other edges, or with other directions. That is
         the map this code used to trust; see buildRowOrientation(). */
};

/**
 * Builds `rowEdges`/`rowReversed`'s RowOrientation against `searchSideTri`.
 *
 * The row's own triangulation and `searchSideTri` are isomorphic, but PD
 * triangulations have automorphisms, so "an isomorphism" does not determine
 * where L goes. When `requiredEdges` is given (a seeded search: the seed's
 * own edges in that boundary component, i.e. L x {0} exactly), the
 * isomorphism used is one taking L's edges onto exactly that set. Any such
 * choice is sound: two of them differ by an automorphism of the row's
 * triangulation preserving L, a PL homeomorphism of S^3 taking L to itself,
 * and the slice genus is invariant under homeomorphism and mirroring.
 *
 * \throws regina::InvalidArgument if `rowEdges` is empty, if the two
 * triangulations are not isomorphic, if no isomorphism takes L onto
 * `requiredEdges`, or if the image of L fails to chain into closed directed
 * curves.
 */
RowOrientation
buildRowOrientation(const std::vector<const regina::Edge<3> *> &rowEdges,
                    const std::vector<bool> &rowReversed,
                    const regina::Triangulation<3> &searchSideTri,
                    const std::vector<size_t> *requiredEdges = nullptr);
} // namespace cobordismgraph

#endif // SURFER_COBOUND_INCOMING_H
