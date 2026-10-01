//
//  edgecycles.h
//
//  Walks over a 1-subcomplex's edges: one helper per family.
//

#ifndef SURFER_LINKNAMING_EDGECYCLES_H
#define SURFER_LINKNAMING_EDGECYCLES_H

#include <cstddef>
#include <optional>
#include <vector>

#include <triangulation/dim3.h>

/*! \file utils/surfer/linknaming/complement/edgecycles.h
 *  \brief Curves made of edges: chained, oriented, counted.
 *
 *  Every walk the project makes over a set of edges is one of three kinds,
 *  and each kind has one helper here:
 *   1. **Directed chaining** (chainDirected(), countDirectedCycles()): edges
 *      already directed, joined head to tail into closed curves -- a
 *      surface's oriented boundary, a diagram's components.
 *   2. **Orienting undirected edges** (walkCurves(), walkClosedCurve()):
 *      edges joined end to end, each curve directed by its first edge -- a
 *      Link's components, one curve to be carried or measured.
 *   3. **Counting and validating** (countClosedCurves()): whether the edges
 *      are disjoint closed curves, and how many.
 *
 *  Edges are given as EdgeEnds (an id and the indices of the two vertices,
 *  tail first when directed), so one triangulation's edges, a drawing's and
 *  a boundary component's all go through the same code. Results are
 *  positions in the input. Every choice is made by input position, never by
 *  address, so the output is a function of the input sequence alone.
 */

namespace edgecycles {

/** An edge: an identity, and its two vertices (tail, head when directed). */
struct EdgeEnds {
    size_t id;
    size_t v0, v1;
};

/** Each edge's index and its vertex(0)/vertex(1) indices, in order. */
std::vector<EdgeEnds> endsOf(const std::vector<const regina::Edge<3> *> &edges);

/** What chainDirected() does when a curve cannot be continued. */
enum class OpenChain {
    stop,   /**< end that curve where it is, and go on */
    refuse, /**< give up: the result is nullopt */
};

/**
 * Chains directed edges (`v0` -> `v1`) head to tail into closed curves.
 *
 * In input order, each edge whose id is on no curve yet starts one: from its
 * tail, the walk takes the edge leaving the current vertex (the LAST such in
 * input order, should several leave it) until it is back at the start, at
 * most `edges.size() + 1` steps. When no edge leaves the current vertex,
 * `open` decides.
 *
 * \return each curve as positions in `edges`, in walk order; nullopt only
 * with OpenChain::refuse.
 */
std::optional<std::vector<std::vector<size_t>>>
chainDirected(const std::vector<EdgeEnds> &edges, OpenChain open);

/**
 * The number of closed directed curves `edges` forms, or nullopt unless
 * every vertex it touches has exactly one edge leaving it and one arriving.
 */
std::optional<size_t> countDirectedCycles(const std::vector<EdgeEnds> &edges);

/** One edge of a walked curve: its position in the input, and whether the
 *  walk runs it `v1` -> `v0`. */
struct Step {
    size_t pos;
    bool reversed;
};

/**
 * Splits undirected edges into curves, each directed by the walk: in input
 * order, each unused edge starts a curve run `v0` -> `v1`, and the walk
 * continues through the first unused edge (in input order) at the current
 * vertex, until there is none. Tolerant: an arc or a branching makes curves
 * that are not closed rather than failing.
 */
std::vector<std::vector<Step>> walkCurves(const std::vector<EdgeEnds> &edges);

/**
 * The edges as ONE simple closed curve, run from the first edge's `v0`, or
 * nullopt unless every vertex they touch meets exactly two of them (a loop
 * edge counts twice) and they are connected.
 */
std::optional<std::vector<Step>> walkClosedCurve(const std::vector<EdgeEnds> &edges);

/**
 * How many closed curves the edges form, or nullopt unless they are a
 * disjoint union of cycles: no repeated id, and every vertex they touch meets
 * exactly two of them (a loop edge counts twice).
 */
std::optional<size_t> countClosedCurves(const std::vector<EdgeEnds> &edges);

} // namespace edgecycles

#endif // SURFER_LINKNAMING_EDGECYCLES_H
