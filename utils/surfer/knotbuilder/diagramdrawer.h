//
//  diagramdrawer.h
//
//  Draws curves in knotbuilder's triangulation of S^3 as link diagrams.
//

/*! \file utils/surfer/knotbuilder/diagramdrawer.h
 *  \brief Draws edge curves of knotbuilder's triangulation of S^3 as oriented
 *  link diagrams, directly from the triangulation's own geometry.
 *
 *  knotbuilder::buildLink()'s triangulation T of S^3 is one block per
 *  crossing, each linearly the box [-1,1]^2 x [0,1] (blockgeometry.h), glued
 *  along walls into S^2 x I, plus two cones. So any edge cycle of T -- in
 *  particular every far side of a cobordism built on T, whose outgoing
 *  boundary is a copy of T -- can be drawn by projecting each edge inside
 *  its own block and reading crossings off block by block:
 *
 *    - Edges inside a block, or in its top or bottom face, project
 *      injectively into that block's square.
 *    - Every degenerate projection comes from an edge lying in a vertical
 *      wall (two blocks' shared side, or a vertical corner line). Such an
 *      edge is replaced by a small tent pushed into one adjacent block, the
 *      push growing with height so that wall points differing only in
 *      height separate. That is an isotopy of the curve: near the wall the
 *      curve has nothing else, since at any vertex it has only its own two
 *      edges and every other edge leaves the wall at once.
 *    - A region's bottom and top points both project to one corner of the
 *      squares around it. If the curve passes through both (not along the
 *      vertical edge joining them), the cyclic order of its four ends
 *      around that corner decides whether that is a crossing (top over
 *      bottom) or two passes that merely touch.
 *    - A pass through a cone apex becomes an arc over (top apex) or under
 *      (bottom apex) everything, routed through the squares, with height
 *      monotone along it so that it stays unknotted.
 *
 *  A PD code needs only the cyclic order of the four ends at each crossing,
 *  and every crossing is found inside one block (or at one corner), so no
 *  global layout of the diagram is ever computed. All arithmetic is exact
 *  (integer coordinates, 128-bit products); a degenerate projection throws
 *  Degenerate, and draw() retries with other perturbations.
 *
 *  The result is oriented: component i is traversed in the order and
 *  direction of the i-th EdgeCycle given.
 */

#ifndef SURFER_KNOTBUILDER_DIAGRAMDRAWER_H
#define SURFER_KNOTBUILDER_DIAGRAMDRAWER_H

#include <array>
#include <memory>
#include <stdexcept>
#include <vector>

#include <link/link.h>
#include <triangulation/dim3.h>

#include "blockgeometry.h"

namespace knotbuilder {

/** One edge of a curve: an edge index of T, and whether it is traversed
 *  from its vertex(1) to its vertex(0). */
struct DirectedEdge {
    size_t edge;
    bool reversed;
};

/** A closed curve: its edges in order, each one's head the next one's tail. */
using EdgeCycle = std::vector<DirectedEdge>;

/** A crossing of a drawn diagram: which components pass over and under,
 *  and its sign (+1 right-handed). */
struct CrossingInfo {
    size_t over;
    size_t under;
    int sign;
};

/**
 * An oriented link diagram. `pd` follows the KnotTheory convention (at
 * each crossing: the incoming under-strand, then counterclockwise), strand
 * labels 1..2n; component i is the i-th curve given to DiagramDrawer::draw().
 * Components with no crossing at all are not in `pd`; `crossingless` lists
 * them (each is an unknot split from the rest).
 */
struct Diagram {
    std::vector<std::array<long, 4>> pd;
    std::vector<CrossingInfo> crossings; /**< Parallel to pd. */
    size_t components = 0;
    std::vector<size_t> crossingless;

    /** The linking number of components i and j (i != j). */
    long linkingNumber(size_t i, size_t j) const;

    /** As a regina::Link, crossingless components included. */
    regina::Link link() const;
};

/** Thrown for a degenerate projection. draw() retries past it. */
class Degenerate : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

/**
 * Draws edge cycles of one knotbuilder triangulation as oriented diagrams.
 *
 * Construction precomputes the block structure (checking that `tri` really
 * is knotbuilder's: 14 tetrahedra per crossing first, every block oriented
 * alike, every vertex at one point per block); draw() is then cheap, and
 * const, so one DiagramDrawer may serve many threads.
 */
class DiagramDrawer {
  public:
    /**
     * \param tri knotbuilder::buildLink()'s triangulation, unmodified (its
     *        tetrahedron and vertex labels are what blockCoordinates()
     *        describes).
     * \param crossings the number of crossings (blocks) it was built from.
     *
     * \exception regina::InvalidArgument `tri` is not such a triangulation,
     * or its diagram is not reduced (a region meeting one crossing twice,
     * or a strand from a crossing to itself).
     */
    DiagramDrawer(const regina::Triangulation<3> &tri, size_t crossings);
    ~DiagramDrawer();
    DiagramDrawer(DiagramDrawer &&) noexcept;
    DiagramDrawer &operator=(DiagramDrawer &&) noexcept;

    /**
     * The diagram of these closed curves (disjoint, in `tri`'s 1-skeleton).
     *
     * \exception Degenerate every one of `attempts` perturbations was
     * degenerate (not seen in practice).
     * \exception regina::InvalidArgument a curve is not closed, or uses an
     * edge no block or cone contains.
     */
    Diagram draw(const std::vector<EdgeCycle> &curves, int attempts = 8) const;

    /** One attempt with a given perturbation seed; throws Degenerate. */
    Diagram drawWithSeed(const std::vector<EdgeCycle> &curves,
                         unsigned seed) const;

    /** The triangulation's own link, as buildLink() returned it (edges plus
     *  reversed flags), as EdgeCycles in component order. */
    static std::vector<EdgeCycle>
    cyclesOf(const std::vector<const regina::Edge<3> *> &edges,
             const std::vector<bool> &reversed);

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace knotbuilder

#endif
