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

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include "gluing.h"
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
    CobordismBuilder(const regina::Triangulation<dim> &tri) : tri_(tri) {
        if (isOrdered(tri_))
            return;

        if constexpr (dim == 3) {
            if (!tri_.order() || !isOrdered(tri_))
                throw regina::InvalidArgument(
                    "CobordismBuilder::CobordismBuilder(): triangulation "
                    "could not be ordered.");
        } else {
            throw regina::InvalidArgument(
                "CobordismBuilder::CobordismBuilder(): triangulation is not "
                "ordered, and automatic ordering is only implemented for "
                "dim == 3.");
        }
    }

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
    const regina::Triangulation<dim> &baseTriangulation() const {
        return tri_;
    }

    // The k-th sub-simplex of the most-recently-built thickening layer's
    // prism over `baseSimplex` (a simplex of baseTriangulation(), NOT of
    // whatever triangulation was originally passed to the constructor --
    // see baseTriangulation()). Only valid after at least one thicken()
    // call, and only reflects the *most recent* layer -- topPrisms_ is
    // replaced wholesale on the next thicken() call, so callers that need
    // per-layer data (e.g. tracing an edge's sweep through every layer)
    // must extract it after each individual thicken() call, before moving
    // on to the next one.
    regina::Simplex<dim + 1> *currentTopSimplex(
        const regina::Simplex<dim> *baseSimplex, int k) const {
        return topPrisms_.at(baseSimplex).simplex(k);
    }

    template <int d>
    static bool isOrdered(const regina::Triangulation<d> &tri) {
        for (const auto &s : tri.simplices()) {
            for (int f = 0; f <= d; ++f) {
                if (s->adjacentSimplex(f) == nullptr)
                    continue;
                regina::Perm<d + 1> g = s->adjacentGluing(f);
                std::vector<int> a;
                for (int i = 0; i < d + 1; ++i) {
                    if (i == f)
                        continue;
                    a.push_back(g[i]);
                }

                if (!std::ranges::is_sorted(a)) {
                    return false;
                }
            }
        }

        return true;
    }

    template <int d>
    static regina::Triangulation<d> &
    glueBoundaries(regina::Triangulation<d> &tri, int bdryIndex1,
                   int bdryIndex2, const regina::Isomorphism<d - 1> &iso) {
        // Defer joins to avoid modifying the triangulation while iterating
        // boundary components.
        std::vector<Gluing<d, d>> gluings;
        const regina::BoundaryComponent<d> *bdry1 =
            tri.boundaryComponent(bdryIndex1);
        const regina::BoundaryComponent<d> *bdry2 =
            tri.boundaryComponent(bdryIndex2);

        for (int i = 0; i < bdry1->size(); ++i) {
            // Since these are boundary faces, they belong to exactly 1 simplex
            const auto &emb1 = bdry1->facet(i)->front();
            const auto &emb2 = bdry2->facet(iso.simpImage(i))->front();
            regina::Simplex<d> *s1 = emb1.simplex();
            regina::Simplex<d> *s2 = emb2.simplex();

            regina::Perm<d + 1> i1 = emb1.vertices();
            regina::Perm<d + 1> i2 = emb2.vertices();
            regina::Perm<d> p = iso.facetPerm(i);
            std::array<int, d + 1> isoPerm;

            for (int j = 0; j < d; ++j) {
                isoPerm[j] = p[j];
            }
            isoPerm[d] = d;

            gluings.push_back(
                {s1, emb1.face(), s2, i2 * isoPerm * i1.inverse()});
        }

        for (const auto &gluing : gluings) {
            gluing.src->join(gluing.srcFacet, gluing.dst, gluing.gluing);
        }

        return tri;
    }

    template <int d>
    static regina::Triangulation<d>
    glueTriangulations(const regina::Triangulation<d> &tri1, int bdryIndex1,
                       const regina::Triangulation<d> &tri2, int bdryIndex2,
                       const regina::Isomorphism<d - 1> &iso) {
        regina::Triangulation<d> tri;
        tri.insertTriangulation(tri1);
        tri.insertTriangulation(tri2);

        return glueBoundaries(tri, bdryIndex1,
                              tri1.countBoundaryComponents() + bdryIndex2, iso);
    }

    // Caps off the current top of the cobordism with a cone on tri_: one
    // new simplex per base simplex of tri_, each with a single new apex
    // vertex, glued together mirroring tri_'s own gluings. If thicken() has
    // not yet been called, the cone alone becomes the whole cobordism (a
    // triangulation of Cone(tri_), e.g. a ball when tri_ is a sphere).
    // Otherwise the cone is glued directly onto the most recent layer's top
    // via SimplicialPrism::capTop(), simplex by simplex, using the same
    // base simplex on both sides (so, as with stitchTop(), no relabelling
    // of base vertices is needed).
    regina::Triangulation<dim + 1> &cone() {
        std::unordered_map<const regina::Simplex<dim> *,
                           regina::Simplex<dim + 1> *>
            coneSimplices;
        coneSimplices.reserve(tri_.size());
        for (const auto *s : tri_.simplices()) {
            coneSimplices.emplace(s, cob_.newSimplex());
        }

        for (const auto *s : tri_.simplices()) {
            regina::Simplex<dim + 1> *coneSimplex = coneSimplices.at(s);

            for (int f = 0; f <= dim; ++f) {
                const regina::Simplex<dim> *adj = s->adjacentSimplex(f);
                if (adj == nullptr || coneSimplex->adjacentSimplex(f) != nullptr)
                    continue;

                regina::Perm<dim + 1> bdryGluing = s->adjacentGluing(f);
                regina::Simplex<dim + 1> *adjConeSimplex = coneSimplices.at(adj);
                std::array<int, dim + 2> coneGluing;
                for (int i = 0; i <= dim; ++i) {
                    coneGluing[i] = bdryGluing[i];
                }
                coneGluing[dim + 1] = dim + 1;

                coneSimplex->join(f, adjConeSimplex, coneGluing);
            }
        }

        if (hasPreviousLayer_) {
            for (const auto *s : tri_.simplices()) {
                topPrisms_.at(s).capTop(coneSimplices.at(s));
            }
        }

        if (!cob_.isValid()) {
            throw regina::InvalidArgument(
                "CobordismBuilder::cone(): Resulting triangulation is not "
                "valid.");
        }

        return cob_;
    }

    inline regina::Triangulation<dim + 1> &thicken() { return thicken_(); }

    inline regina::Triangulation<dim + 1> &thicken(int layers) {
        for (int i = 0; i < layers; ++i) {
            thicken_();
        }

        return cob_;
    }

  private:
    regina::Triangulation<dim + 1> &thicken_() {
        // Make a new prism for each simplex in the triangulation. This is
        // its own layer: it gets fully glued together internally below,
        // independently of any previous layer.
        PrismMap newPrisms;
        newPrisms.reserve(tri_.size());
        for (const auto *s : tri_.simplices()) {
            newPrisms.emplace(s, cob_);
        }

        // Now glue the prisms together along their walls according to the
        // gluing of the original triangulation
        std::set<std::pair<const regina::Simplex<dim> *, int>> visited;

        for (size_t i = 0; i < tri_.size(); ++i) {
            for (int facet = 0; facet < dim + 1; ++facet) {
                const regina::Simplex<dim> *s = tri_.simplex(i);
                const regina::Simplex<dim> *adj = s->adjacentSimplex(facet);

                if (adj == nullptr || visited.contains({s, facet}))
                    continue;

                int adjFacet = s->adjacentFacet(facet);

                newPrisms.at(s).glue(facet, newPrisms.at(adj), adjFacet);

                visited.insert({s, facet});
                visited.insert({adj, adjFacet});
            }
        }

        // If there is a previous layer, stitch this layer's bottom onto its
        // top, simplex by simplex: since both layers thicken the same base
        // triangulation, the seam uses the same base simplex on both sides
        // and needs no relabelling of base vertices.
        if (hasPreviousLayer_) {
            for (const auto *s : tri_.simplices()) {
                topPrisms_.at(s).stitchTop(newPrisms.at(s));
            }
        }

        topPrisms_ = std::move(newPrisms);
        hasPreviousLayer_ = true;

        if (!cob_.isValid()) {
            throw regina::InvalidArgument(
                "CobordismBuilder::thicken(): Resulting triangulation is not "
                "valid.");
        }

        return cob_;
    }
};

#endif // COBORDISM_BUILDER_H
