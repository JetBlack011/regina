//
//  knotbuilder.h
//
//  Created by John Teague on 05/10/2024.
//
//  This is adapted from work of Srinivas Vadhiraj, Samantha Ward, Angela
//  Yuan, and Jingyuan Zhang performed for the Texas Experimental Geometry Lab
//  at UT Austin.

#ifndef KNOTBUILDER_H

#define KNOTBUILDER_H

#include <triangulation/dim3.h>

/*! \file utils/surfer/knotbuilder.h
 *  \brief Builds a triangulation of S^3 containing a knot/link diagram
 *  from PD code.
 */

namespace knotbuilder {
/** A planar diagram code: one 4-tuple of strand labels per crossing. */
using PDCode = std::vector<std::array<int, 4>>;

/** Parses a PD (planar diagram) code string into a PDCode. */
PDCode parsePDCode(std::string pdcode_str);

/**
 * One crossing's tetrahedral gadget, glued into the ambient triangulation
 * by buildLink() and joined to its neighbors' Blocks via glue().
 */
class Block {
  private:
    std::vector<regina::Tetrahedron<3> *> core_; /**< The crossing's central tetrahedra. */
    std::array<std::array<regina::Tetrahedron<3> *, 3>, 4> walls_;
        /**< walls_[w]: the 3 tetrahedra forming wall `w` (one of the
             crossing's 4 strand slots), to be glued to a neighboring
             Block's wall via glue(). */

  public:
    /** Builds a new crossing gadget in `tri`. */
    Block(regina::Triangulation<3> &tri);

    /** Glues this Block's wall `myWall` to `other`'s wall `otherWall`. */
    void glue(size_t myWall, Block &other, size_t otherWall);

    /** Returns the 3 edges of this crossing that trace the knot/link diagram through it. */
    const std::vector<regina::Edge<3> *> getLinkEdges() const;

    /**
     * Returns the PD-intended "reversed" flags (see TriangulationWithLink::
     * reversed) for getLinkEdges()'s 3 edges, in the same order.
     *
     * The understrand edge (getLinkEdges()[1], `core_[4]->edge(0,2)`) is
     * unconditionally not reversed: PD position 0 is always the incoming
     * understrand (true for every crossing, independent of the specific
     * PD code -- see buildLink()'s own comment), and that edge's
     * `vertex(0)` sits on the `walls_[0]` side, `vertex(1)` on the
     * `walls_[2]` side, so the PD-intended direction is always
     * `vertex(0)` -> `vertex(1)`.
     *
     * The overstrand's two edges (getLinkEdges()[0] and [2]) depend on
     * `walls1IsEntry`: whether `walls_[1]` (rather than `walls_[3]`) is
     * the PD-intended arrival end of the overstrand pass at this specific
     * crossing -- not locally fixed like the understrand, so the caller
     * must determine it via a walk over the PD code (see buildLink()).
     */
    std::vector<bool> getLinkEdgeDirections(bool walls1IsEntry) const;

    /** Returns whether `lhs` and `rhs` are the same underlying gadget (same tetrahedra). */
    friend bool operator==(const Block &lhs, const Block &rhs) {
        return lhs.core_ == rhs.core_ && lhs.walls_ == rhs.walls_;
    }
};

/** A triangulation of S^3 together with the edges tracing a knot/link diagram within it. */
struct TriangulationWithLink {
    regina::Triangulation<3> tri; /**< The triangulation. */
    std::vector<const regina::Edge<3> *> edges; /**< The edges tracing the knot/link diagram. */

    /**
     * A PD-derived direction for each of `edges`: `edges[i]` runs
     * `vertex(1)` -> `vertex(0)` (the PD-intended direction) iff
     * `reversed[i]`. Always populated, one entry per `edges` entry, same
     * order.
     *
     * This is a fixed direction on `edges`, derived once here from
     * `pdcode`'s own orientation convention (the same convention
     * `regina::Link::fromPD` uses -- see `engine/link/pd-impl.h`, only
     * referenced for its convention, never built or run against at
     * runtime). It exists so a downstream comparison can check whether a
     * *found* surface's own induced boundary direction agrees with this
     * row's own diagram, everywhere or nowhere, without recomputing
     * anything knot-theoretic at comparison time.
     */
    std::vector<bool> reversed;
};

/**
 * Builds a triangulation of S^3 from `pdcode`, along with the edges of
 * that triangulation that trace the knot/link diagram (3 per crossing).
 *
 * \throws regina::InvalidArgument if `pdcode` is not well-formed (some
 * strand label is not used in exactly two crossing slots).
 */
TriangulationWithLink buildLink(PDCode pdcode);

/**
 * Builds a copy of `tri` with every internal edge pinched away except the
 * edges in `edges` (mapped into the copy) and any loop edge (both
 * endpoint vertices coincide) -- pinching a loop drills it out instead of
 * collapsing it, leaving an ideal torus/Klein-bottle boundary component,
 * so those are always skipped.
 *
 * `tri` is left untouched. The returned edge list corresponds 1-1 (same
 * order and length, including repeats) to `edges`, mapped into the
 * returned triangulation.
 *
 * \throws regina::InvalidArgument if any edge in `edges` is a boundary
 * edge of `tri` (pinchEdge() cannot be skipped or performed on a boundary
 * edge, and buildLink()'s output is always closed, so this should not
 * arise for its output).
 */
TriangulationWithLink
reduceVertices(const regina::Triangulation<3> &tri,
               const std::vector<const regina::Edge<3> *> &edges);
} // namespace knotbuilder

#endif // KNOTBUILDER_H
