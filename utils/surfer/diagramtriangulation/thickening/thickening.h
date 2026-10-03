//
//  cobordismbuilder.h
//
//  Created by John Teague on 05/10/2024.
//
//  This is adapted from work of Srinivas Vadhiraj, Samantha Ward, Angela
//  Yuan, and Jingyuan Zhang performed for the Texas Experimental Geometry Lab
//  at UT Austin.

#ifndef COBORDISM_BUILDER_H

#define COBORDISM_BUILDER_H

#include <cassert>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include "diagramtriangulation/thickening/prism.h"
#include "diagramtriangulation/fromdiagram.h"
#include "diagramtriangulation/todiagram.h"

/*! \file utils/surfer/diagramtriangulation/thickening/thickening.h
 *  \brief Incrementally builds a cobordism triangulation one dimension up
 *  from a base triangulation.
 */

/**
 * Incrementally builds a (dim+1)-dimensional triangulation "cobordism"
 * from a dim-dimensional base triangulation, by stacking prism layers
 * (thicken()).
 *
 * \tparam dim the dimension of the base triangulation; the cobordism
 * itself lives in dimension dim+1.
 */
template <int dim>
class CobordismBuilder {
  private:
    using PrismMap = std::unordered_map<const regina::Simplex<dim> *,
                                        SimplicialPrism<dim + 1>>;

    PrismMap topPrisms_;
        /**< The most recently built thickening layer's prisms, keyed by
             base simplex. */
    bool hasPreviousLayer_ = false;
        /**< Whether at least one thicken() call has been made. */

    regina::Triangulation<dim> tri_; /**< See baseTriangulation(). */
    regina::Triangulation<dim + 1> cob_;
        /**< The cobordism built so far; see getCobordism(). */

    regina::Simplex<dim + 1> *baseBoundaryFacetSimplex_ = nullptr;
        /**< Captured on the first thicken() call only; see
             baseBoundaryComponent(). A raw Simplex<dim+1>* stays valid
             across every later thicken() call (unlike
             Face<dim,subdim>* -- see CollarBuilder's identical warning),
             so it is safe to resolve into an actual boundary component
             lazily, on demand, after construction is otherwise complete. */

  public:
    /**
     * Takes (and owns) a copy of `tri`, since thicken() requires an
     * ordered triangulation (see isOrdered()) and may need to relabel
     * vertices to achieve this.
     */
    CobordismBuilder(const regina::Triangulation<dim> &tri);

    /**
     * Returns the triangulation thicken() actually builds against.
     *
     * \warning This is *not* the same object passed to the constructor:
     * the constructor takes a copy (and may reorder it), so a
     * Simplex<dim>* or Edge<dim>* etc. obtained from the caller's own
     * triangulation is a dangling reference here, even though indices are
     * preserved across the copy and any subsequent order(). To refer to a
     * specific piece of the original triangulation after construction,
     * look it up here by index (e.g.
     * `cob.baseTriangulation().edge(origEdge->index())`).
     */
    const regina::Triangulation<dim> &baseTriangulation() const { return tri_; }

    /**
     * Returns the k-th sub-simplex of the most recently built thickening
     * layer's prism over `baseSimplex` (a simplex of baseTriangulation(),
     * not of whatever triangulation was originally passed to the
     * constructor).
     *
     * \pre At least one thicken() call has been made.
     *
     * \note Only reflects the *most recent* layer -- this data is
     * replaced wholesale on the next thicken() call, so callers needing
     * per-layer data (e.g. tracing an edge's sweep through every layer)
     * must extract it after each thicken() call, before moving to the
     * next.
     */
    regina::Simplex<dim + 1> *
    currentTopSimplex(const regina::Simplex<dim> *baseSimplex, int k) const {
        return topPrisms_.at(baseSimplex).simplex(k);
    }

    /** Returns whether `tri` is ordered -- a precondition for thicken(). */
    static bool isOrdered(const regina::Triangulation<dim> &tri);

    /**
     * Glues two boundary components of `tri` together via `iso`.
     *
     * \note Templated independently on `d` (not tied to the class's own
     * `dim`): explicit instantiation of CobordismBuilder<2> would
     * otherwise eagerly require regina::Isomorphism<1>, which doesn't
     * exist. Nothing in the codebase currently calls either this or
     * glueTriangulations() below, so no explicit instantiation of either
     * exists yet -- add one in cobordismbuilder.cpp for whatever `d` a
     * future caller needs.
     */
    template <int d>
    static regina::Triangulation<d> &
    glueBoundaries(regina::Triangulation<d> &tri, int bdryIndex1,
                   int bdryIndex2, const regina::Isomorphism<d - 1> &iso);

    /**
     * Glues a boundary component of `tri1` to one of `tri2` via `iso`,
     * producing a new triangulation.
     *
     * \note See glueBoundaries() above for why this is independently
     * templated on `d`.
     */
    template <int d>
    static regina::Triangulation<d>
    glueTriangulations(const regina::Triangulation<d> &tri1, int bdryIndex1,
                       const regina::Triangulation<d> &tri2, int bdryIndex2,
                       const regina::Isomorphism<d - 1> &iso);

    /** Adds one thickening layer (a prism per base simplex) to the cobordism. */
    inline regina::Triangulation<dim + 1> &thicken() { return thicken_(); }

    /** Adds `layers` successive thickening layers to the cobordism. */
    inline regina::Triangulation<dim + 1> &thicken(int layers) {
        for (int i = 0; i < layers; ++i) {
            thicken_();
        }

        return cob_;
    }

    /** Returns the cobordism built so far. */
    const regina::Triangulation<dim + 1> &getCobordism() const { return cob_; }

    /**
     * Returns the boundary component of getCobordism() corresponding to
     * the original base triangulation -- thicken()'s "bottom", untouched
     * by any thicken() call -- distinguishing it from the cobordism's
     * other boundary component, its "top".
     *
     * Since thicken() literally builds (base triangulation) x [0,1],
     * every boundary component this cobordism could ever have is
     * topologically -- and combinatorially -- identical to the base
     * triangulation, so isomorphism-signature matching alone can never
     * distinguish "bottom" from "top". This instead resolves an actual
     * ambient facet known, by construction, to lie on the untouched
     * original copy specifically: SimplicialPrism's own encoding
     * guarantees that simplex(0) of the very first thicken() layer holds
     * every one of the base simplex's own "bottom" vertices (at local
     * indices 0..dim) plus the "top" copy of base vertex 0 (at local
     * index dim+1) -- so the facet omitting that last local vertex is
     * exactly the base simplex's own untouched bottom facet, which stays
     * permanently unglued (thicken_()'s own wall-gluing only ever touches
     * *other* facets, and stitching/coning only ever touch the *top* of
     * whichever layer is most recent).
     *
     * \throws regina::InvalidArgument if no thicken() call has been made.
     */
    regina::BoundaryComponent<dim + 1> *baseBoundaryComponent() const;

  private:
    /** Adds a single thickening layer; see thicken(). */
    regina::Triangulation<dim + 1> &thicken_();
};

extern template class CobordismBuilder<2>;
extern template class CobordismBuilder<3>;

/**
 * One edge of a curve on a thickening's outgoing boundary -- an edge of that
 * boundary component's built triangulation -- and whether the curve runs
 * along it from its vertex(1) to its vertex(0).
 */
struct OutgoingEdge {
    const regina::Edge<3> *edge;
    bool reversed;
};

/** A curve on the outgoing boundary: its edges in order, head to tail. */
using OutgoingCurve = std::vector<OutgoingEdge>;

/**
 * The outgoing boundary of a thickening, edge by edge, as knotbuilder's T.
 *
 * A search runs in S^3 x [0,2] = T x [0,2], thickened from knotbuilder's
 * triangulation T by CobordismBuilder. Its outgoing boundary is literally
 * T x {2}: over each base tetrahedron sigma, the top prism piece
 * P_3(sigma) has sigma x {2} as a facet, with base vertex v's top copy at
 * SimplicialPrism::localVertex(v, true). OutgoingMap uses exactly that to
 * carry each edge of the outgoing boundary onto an edge of T. There is no
 * isomorphism search anywhere on the outgoing side, so none of T's
 * automorphisms -- some of which reverse orientation -- can relabel or
 * mirror the outgoing link.
 */
class OutgoingMap {
  public:
    /**
     * \param knotT knotbuilder::buildLink()'s triangulation, unmodified.
     * \param cob built from `knotT` (CobordismBuilder takes an ordered copy,
     *        relabelling vertices within tetrahedra), after its last
     *        thicken().
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
    knotbuilder::EdgeCycle carry(const OutgoingCurve &curve) const;

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
 * A link's diagram triangulated -- knotbuilder's T, with L's edges --
 * thickened into S^3 x I, and seeded with the collar L x [0, collarLayers]:
 * what a search runs in. Filled in place by buildAmbient() and never moved:
 * a namer and a search built from it hold pointers into `link.tri` and
 * `cob`.
 */
struct ThickenedLink {
    ThickenedLink() = default;
    ThickenedLink(const ThickenedLink &) = delete;
    ThickenedLink &operator=(const ThickenedLink &) = delete;

    knotbuilder::PDCode pdcode;
    knotbuilder::TriangulationWithLink link; /**< T, and L's edges in it. */
    std::optional<CobordismBuilder<3>> cob;
    regina::Triangulation<4> tri;  /**< The search's ambient. */
    std::vector<int> seedFaces;    /**< The collar, in index order; empty without one. */
    size_t searchSideBC = 0;       /**< The incoming boundary, T x {0}. */
    int componentCount = 1;
    /**< The closed curves L's edges form, counted by walking them. */
};

/**
 * Builds `out` for PD code `pdNotation`: T, and `thickenLayers` thickenings
 * with a collar through the first `collarLayers`.
 *
 * \throws regina::InvalidArgument for an unparseable or unbuildable PD code.
 */
void buildAmbient(const std::string &pdNotation, int thickenLayers,
                  int collarLayers, ThickenedLink &out);

#endif // COBORDISM_BUILDER_H
