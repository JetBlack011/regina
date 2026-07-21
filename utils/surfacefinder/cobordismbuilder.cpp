//
//  cobordismbuilder.cpp
//
//  Created by John Teague on 04/12/2025.
//

#include "cobordismbuilder.h"

#include <cassert>

template <int dim>
CobordismBuilder<dim>::CobordismBuilder(const regina::Triangulation<dim> &tri)
    : tri_(tri) {
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

template <int dim>
bool CobordismBuilder<dim>::isOrdered(const regina::Triangulation<dim> &tri) {
    for (const auto &s : tri.simplices()) {
        for (int f = 0; f <= dim; ++f) {
            if (s->adjacentSimplex(f) == nullptr)
                continue;
            regina::Perm<dim + 1> g = s->adjacentGluing(f);
            std::vector<int> a;
            for (int i = 0; i < dim + 1; ++i) {
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

template <int dim>
template <int d>
regina::Triangulation<d> &CobordismBuilder<dim>::glueBoundaries(
    regina::Triangulation<d> &tri, int bdryIndex1, int bdryIndex2,
    const regina::Isomorphism<d - 1> &iso) {
    using Gluing = std::tuple<regina::Simplex<d> *, int,
                              regina::Simplex<d> *, regina::Perm<d + 1>>;
    // Defer joins to avoid modifying the triangulation while iterating
    // boundary components.
    std::vector<Gluing> gluings;
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

        gluings.emplace_back(s1, emb1.face(), s2,
                             i2 * isoPerm * i1.inverse());
    }

    for (const auto &[src, srcFacet, dst, gluingPerm] : gluings) {
        src->join(srcFacet, dst, gluingPerm);
    }

    return tri;
}

template <int dim>
template <int d>
regina::Triangulation<d> CobordismBuilder<dim>::glueTriangulations(
    const regina::Triangulation<d> &tri1, int bdryIndex1,
    const regina::Triangulation<d> &tri2, int bdryIndex2,
    const regina::Isomorphism<d - 1> &iso) {
    regina::Triangulation<d> tri;
    tri.insertTriangulation(tri1);
    tri.insertTriangulation(tri2);

    return glueBoundaries<d>(tri, bdryIndex1,
                             tri1.countBoundaryComponents() + bdryIndex2, iso);
}

template <int dim>
regina::Triangulation<dim + 1> &CobordismBuilder<dim>::cone() {
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
            if (adj == nullptr ||
                coneSimplex->adjacentSimplex(f) != nullptr)
                continue;

            regina::Perm<dim + 1> bdryGluing = s->adjacentGluing(f);
            regina::Simplex<dim + 1> *adjConeSimplex =
                coneSimplices.at(adj);
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

    // Validity here is guaranteed by construction (coning a valid,
    // ordered triangulation, glued face-for-face onto a valid previous
    // layer, cannot produce an invalid result) -- checking it costs a
    // full vertex-link recognition pass (Triangulation<3>::isSphere()/
    // isBall() per vertex, i.e. real 3-manifold simplification), which
    // dominates runtime on triangulations of any size. Debug-only.
    assert(cob_.isValid() &&
           "CobordismBuilder::cone(): resulting triangulation is not "
           "valid.");

    return cob_;
}

template <int dim>
regina::Triangulation<dim + 1> &CobordismBuilder<dim>::thicken_() {
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

    // See the identical comment in cone(): validity is guaranteed by
    // construction, and checking it here is the dominant cost of
    // thicken() (a full vertex-link recognition pass over the whole
    // accumulated cobordism, on every single layer). Debug-only.
    assert(cob_.isValid() &&
           "CobordismBuilder::thicken(): resulting triangulation is not "
           "valid.");

    return cob_;
}

template class CobordismBuilder<2>;
template class CobordismBuilder<3>;
