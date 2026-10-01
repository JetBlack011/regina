// simplification.h
//
// Simplifying a diagram without losing track of its components.

/*! \file utils/surfer/linknaming/diagrams/simplification.h
 *  \brief Diagram simplification that keeps every component in its slot
 *  and every linking number: Regina's simplify(), nugatory crossings
 *  removed, and components lying over (or under) everything lifted off.
 *  The linking matrix it checks against.
 */

#pragma once

#include <optional>
#include <vector>

#include "linknaming/diagrams/gaussdiagram.h"

namespace cascade {

/// The linking matrix of a diagram, components in its own order.
std::vector<std::vector<int>> linkingMatrix(const exactnaming::GaussDiagram &d);

/**
 * Simplifies a diagram with regina::Link::simplify() (Reidemeister moves
 * only: never reflects or reverses, and keeps component indices, a zero-
 * crossing component in its slot), keeping `origin`, then removes every
 * nugatory crossing left (removeNugatoryCrossings()): a hop's row must be a
 * reduced diagram for knotbuilder's drawer to certify it. Throws std::logic_error
 * if the component count or any pairwise linking number changed: component
 * identity is what every profile is indexed by, so a violation must stop
 * everything rather than mislabel components.
 */
exactnaming::GaussDiagram simplifyKeepingComponents(const exactnaming::GaussDiagram &d);

/// A nugatory crossing of `d` (one whose removal disconnects its diagram:
/// always a self-crossing), or nullopt if `d` is reduced.
std::optional<size_t> nugatoryCrossing(const exactnaming::GaussDiagram &d);

/// `d` with every nugatory crossing removed: each by turning over the side it
/// cuts off, which deletes it, swaps over and under on that side's crossings
/// and keeps every sign. The same oriented link, components in the same
/// order, with every linking number kept.
exactnaming::GaussDiagram removeNugatoryCrossings(exactnaming::GaussDiagram d);

/// `d` with every component that is over (or under) at every crossing it
/// meets lifted off: its crossings deleted, leaving it crossingless. Such a
/// component has no self-crossing, so it is an unknot above (below)
/// everything else, and lifting it is an isotopy that splits it off. Its PD
/// code could not carry its orientation (Regina's pdAmbiguous()), so a row
/// keeping it could not be certified. Repeated until none is left.
exactnaming::GaussDiagram liftSplitComponents(exactnaming::GaussDiagram d);

} // namespace cascade
