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

#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include "simplicialprism.h"

/*! \file utils/surfer/cobordismbuilder.h
 *  \brief Incrementally builds a cobordism triangulation one dimension up
 *  from a base triangulation.
 */

/**
 * Incrementally builds a (dim+1)-dimensional triangulation "cobordism"
 * from a dim-dimensional base triangulation, by stacking prism layers
 * (thicken()) and/or capping with a cone (cone()).
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

  public:
    /**
     * Takes (and owns) a copy of `tri`, since thicken() requires an
     * ordered triangulation (see isOrdered()) and may need to relabel
     * vertices to achieve this.
     */
    CobordismBuilder(const regina::Triangulation<dim> &tri);

    /**
     * Returns the triangulation thicken()/cone() actually build against.
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

    /**
     * Caps off the current top of the cobordism with a cone on the base
     * triangulation: one new simplex per base simplex, each with a single
     * new apex vertex, glued together mirroring the base triangulation's
     * own gluings.
     *
     * If thicken() has not yet been called, the cone alone becomes the
     * whole cobordism (the cone over the base triangulation, e.g. a ball
     * when the base triangulation is a sphere). Otherwise the cone is
     * glued directly onto the most recent layer's top via
     * SimplicialPrism::capTop().
     */
    regina::Triangulation<dim + 1> &cone();

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

  private:
    /** Adds a single thickening layer; see thicken(). */
    regina::Triangulation<dim + 1> &thicken_();
};

extern template class CobordismBuilder<2>;
extern template class CobordismBuilder<3>;

#endif // COBORDISM_BUILDER_H
