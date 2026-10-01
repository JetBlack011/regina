//
//  simplicialprism.cpp
//
//  Created by John Teague on 07/21/2026.
//

#include "simplicialprism.h"

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

template <int dim>
SimplicialPrism<dim>::SimplicialPrism(regina::Triangulation<dim> &tri) {
    for (int i = 0; i < dim; ++i) {
        simplices_.push_back(tri.newSimplex());
    }

    regina::Perm<dim + 1> id;

    for (int i = 1; i < dim; ++i) {
        simplices_[i - 1]->join(i - 1, simplices_[i], id);
    }
}

template <int dim>
void SimplicialPrism<dim>::glue(int facet, SimplicialPrism<dim> &other,
                                int otherFacet) {
    if (!(0 <= facet && facet < dim && 0 <= otherFacet && otherFacet < dim))
        throw regina::InvalidArgument(
            "SimplicialPrism::glue(): Invalid faces");

    for (int k = 0; k < dim; ++k) {
        if (k == facet)
            continue;

        int kOther = mapExcluding_(k, facet, otherFacet);

        int myFacet = wallFacet_(k, facet);
        int otherLocalFacet = wallFacet_(kOther, otherFacet);

        std::array<int, dim + 1> image;
        for (int m = 0; m <= dim; ++m) {
            if (m == myFacet) {
                image[m] = otherLocalFacet;
                continue;
            }

            // Decode local vertex m as (base vertex, is-top-copy),
            // translate the base vertex across the (ordered) base gluing,
            // then re-encode within other.simplices_[kOther]. The top/
            // bottom level is unaffected, since the gluing only
            // identifies the base simplices and acts as the identity on
            // the interval factor.
            auto [baseVertex, isTop] = decode_(k, m);
            int otherBaseVertex =
                mapExcluding_(baseVertex, facet, otherFacet);
            image[m] = encode_(otherBaseVertex, isTop);
        }

        simplices_[k]->join(myFacet, other.simplices_[kOther],
                            regina::Perm<dim + 1>(image));
    }
}

template <int dim>
void SimplicialPrism<dim>::stitchTop(SimplicialPrism<dim> &next) {
    std::array<int, dim + 1> image;
    for (int i = 0; i <= dim; ++i)
        image[i] = (i + 1) % (dim + 1);

    simplices_[dim - 1]->join(dim - 1, next.simplices_[0],
                              regina::Perm<dim + 1>(image));
}

template <int dim>
void SimplicialPrism<dim>::capTop(regina::Simplex<dim> *coneSimplex) {
    std::array<int, dim + 1> image;
    for (int m = 0; m <= dim; ++m) {
        if (m == dim - 1) {
            image[m] = dim;
            continue;
        }
        image[m] = decode_(dim - 1, m).first;
    }

    simplices_[dim - 1]->join(dim - 1, coneSimplex,
                              regina::Perm<dim + 1>(image));
}

template <int dim>
int SimplicialPrism<dim>::mapExcluding_(int x, int a, int b) {
    // Re-rank x within {0,...,dim} with a removed (closing the gap left by
    // a), then re-open a gap at b in the target set.
    int rank = (x < a) ? x : x - 1;
    return (rank < b) ? rank : rank + 1;
}

template <int dim>
std::pair<int, bool> SimplicialPrism<dim>::decode_(int k, int m) {
    // Local vertex `dim` of every simplices_[k] is reserved for the top
    // copy of base vertex 0. Otherwise, local vertex m is the bottom copy
    // of base vertex m when m >= k, and the top copy of base vertex m + 1
    // when m < k.
    if (m == dim)
        return {0, true};
    if (k <= m)
        return {m, false};
    return {m + 1, true};
}

template <int dim>
int SimplicialPrism<dim>::encode_(int v, bool isTop) {
    // Bottom copies keep their base vertex index as the local index. The
    // top copy of base vertex 0 lives at local index dim (see decode_());
    // every other top copy shifts down by one to make room for it.
    if (!isTop)
        return v;
    return v == 0 ? dim : v - 1;
}

template <int dim>
int SimplicialPrism<dim>::wallFacet_(int k, int v) {
    // See decode_(): base vertex v's bottom copy sits at local index v
    // when k <= v, and its top copy at local index v - 1 (or dim, for
    // v == 0) when k > v. The k != v precondition documented in the header
    // guarantees exactly one of these applies.
    if (k < v)
        return v;
    return v == 0 ? dim : v - 1;
}

template class SimplicialPrism<3>;
template class SimplicialPrism<4>;
