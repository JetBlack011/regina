//
//  knotbuilder.h
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

template <int dim>
class CobordismBuilder {
  private:
    using PrismMap = std::unordered_map<const regina::Simplex<dim> *,
                                        SimplicialPrism<dim + 1>>;

    PrismMap topPrisms_;
    bool hasPreviousLayer_ = false;

    regina::Triangulation<dim> tri_;
    regina::Triangulation<dim + 1> cob_;

  public:
    // Takes (and owns) a copy of `tri`, since thicken() requires an ordered
    // triangulation (see isOrdered() below) and may need to relabel vertices
    // to achieve this.
    CobordismBuilder(const regina::Triangulation<dim> &tri);

    // The triangulation thicken()/cone() actually build against -- NOT the
    // same object as whatever was passed to the constructor. The
    // constructor takes a *copy* of its argument (and may reorder it), so a
    // Simplex<dim>*/Edge<dim>* etc. obtained from the caller's own
    // triangulation is a dangling reference to a different object as far as
    // this class is concerned, even though indices (tetrahedron, edge,
    // vertex) are preserved across the copy and any subsequent order().
    // Callers that need to refer to a specific piece of the original
    // triangulation (e.g. a knot edge) after construction must look it up
    // here by index (e.g. `cob.baseTriangulation().edge(origEdge->index())`)
    // rather than reusing the original pointer.
    const regina::Triangulation<dim> &baseTriangulation() const { return tri_; }

    // The k-th sub-simplex of the most-recently-built thickening layer's
    // prism over `baseSimplex` (a simplex of baseTriangulation(), NOT of
    // whatever triangulation was originally passed to the constructor --
    // see baseTriangulation()). Only valid after at least one thicken()
    // call, and only reflects the *most recent* layer -- topPrisms_ is
    // replaced wholesale on the next thicken() call, so callers that need
    // per-layer data (e.g. tracing an edge's sweep through every layer)
    // must extract it after each individual thicken() call, before moving
    // on to the next one.
    regina::Simplex<dim + 1> *
    currentTopSimplex(const regina::Simplex<dim> *baseSimplex, int k) const {
        return topPrisms_.at(baseSimplex).simplex(k);
    }

    static bool isOrdered(const regina::Triangulation<dim> &tri);

    // Kept as independently-templated members (parameter `d`, not tied to
    // the class's own `dim`) rather than plain `dim`-using members like
    // isOrdered() above: explicit instantiation of CobordismBuilder<2>
    // would otherwise eagerly require regina::Isomorphism<dim - 1>, i.e.
    // Isomorphism<1>, which doesn't exist. Lazy (on-demand) instantiation
    // avoids that; nothing in the codebase currently calls either of
    // these, so no explicit instantiation of them exists yet either --
    // add one (in cobordismbuilder.cpp) for whatever `d` a future caller
    // needs.
    template <int d>
    static regina::Triangulation<d> &
    glueBoundaries(regina::Triangulation<d> &tri, int bdryIndex1,
                   int bdryIndex2, const regina::Isomorphism<d - 1> &iso);

    template <int d>
    static regina::Triangulation<d>
    glueTriangulations(const regina::Triangulation<d> &tri1, int bdryIndex1,
                       const regina::Triangulation<d> &tri2, int bdryIndex2,
                       const regina::Isomorphism<d - 1> &iso);

    // Caps off the current top of the cobordism with a cone on tri_: one
    // new simplex per base simplex of tri_, each with a single new apex
    // vertex, glued together mirroring tri_'s own gluings. If thicken() has
    // not yet been called, the cone alone becomes the whole cobordism (a
    // triangulation of Cone(tri_), e.g. a ball when tri_ is a sphere).
    // Otherwise the cone is glued directly onto the most recent layer's top
    // via SimplicialPrism::capTop(), simplex by simplex, using the same
    // base simplex on both sides (so, as with stitchTop(), no relabelling
    // of base vertices is needed).
    regina::Triangulation<dim + 1> &cone();

    inline regina::Triangulation<dim + 1> &thicken() { return thicken_(); }

    inline regina::Triangulation<dim + 1> &thicken(int layers) {
        for (int i = 0; i < layers; ++i) {
            thicken_();
        }

        return cob_;
    }

    const regina::Triangulation<dim + 1> &getCobordism() const { return cob_; }

  private:
    regina::Triangulation<dim + 1> &thicken_();
};

extern template class CobordismBuilder<2>;
extern template class CobordismBuilder<3>;

#endif // COBORDISM_BUILDER_H
