//
//  gaussdiagram.h
//
//  Oriented link diagrams as signed Gauss data, cut into split pieces and
//  visible connected summands without losing orientation or chirality.
//

/*! \file utils/surfer/linknaming/diagrams/gaussdiagram.h
 *  \brief An oriented diagram as its signed crossings plus, per component,
 *  the crossings it passes in order; and the two decompositions a diagram can
 *  show directly.
 *
 *  This is regina::Link::fromData()'s input, and it is the whole diagram:
 *  the cyclic order of the four ends at a crossing follows from its sign and
 *  which strand is over. Every operation here keeps each component's
 *  orientation and the diagram's chirality, which is what lets a piece's name
 *  record both.
 *
 *  - splitPieces(): components joined by crossings form one piece; a diagram
 *    with several pieces is a split link, the union of the pieces.
 *  - visibleSum(): two arcs bounding the same two faces, whose removal
 *    disconnects the crossings, are the two points where a circle meets the
 *    diagram, so the link is the connected sum of the two sides, each closed
 *    along that circle. (Both arcs lie on one component, since a closed
 *    curve meets a separating circle an even number of times.)
 *
 *  Neither finds a split or sum the diagram does not show; the caller
 *  simplifies first (regina::Link::simplify() and simplifyExhaustive() never
 *  reflect or reverse, unlike rewrite()).
 */

#ifndef SURFER_LINKNAMING_GAUSSDIAGRAM_H
#define SURFER_LINKNAMING_GAUSSDIAGRAM_H

#include <optional>
#include <utility>
#include <vector>

#include <link/link.h>

namespace linknaming {

struct GaussDiagram {
    std::vector<int> signs;             /**< per crossing, +1 right-handed */
    std::vector<std::vector<long>> comps; /**< +(k+1) over crossing k, -(k+1) under */
    std::vector<size_t> origin;         /**< per component: the outgoing component it is (part of) */

    size_t crossings() const { return signs.size(); }
    size_t components() const { return comps.size(); }
    regina::Link link() const;

    /** Read a Regina link back, components in its own order. */
    static GaussDiagram of(const regina::Link &l, std::vector<size_t> origin);
    /** The same, each component its own origin (0, 1, ...). */
    static GaussDiagram of(const regina::Link &l);
};

/** The diagram's split pieces, each a group of components joined by
 *  crossings; a crossingless component is a piece of its own. */
std::vector<GaussDiagram> splitPieces(const GaussDiagram &d);

/**
 * A visible connected-sum sphere and the two summands it cuts off, or
 * nullopt if the diagram shows none. The summed component appears in both
 * summands, with the same origin.
 *
 * \pre `d` is connected (one split piece) and planar.
 */
std::optional<std::pair<GaussDiagram, GaussDiagram>> visibleSum(const GaussDiagram &d);

} // namespace linknaming

#endif
