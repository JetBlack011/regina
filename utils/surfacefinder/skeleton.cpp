//
//  skeleton.cpp
//

#include "skeleton.h"

#include <array>
#include <utility>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

template <int dim, int subdim>
Skeleton<dim, subdim>::Skeleton(const regina::Triangulation<dim> &tri)
    : tri_(&tri), numVertices_(tri.countVertices()) {
    static_assert(subdim >= 1 && subdim <= dim,
                  "Skeleton requires 1 <= subdim <= dim");
    buildSkeleton_(tri);
}

template <int dim, int subdim>
void Skeleton<dim, subdim>::buildSkeleton_(
    const regina::Triangulation<dim> &tri) {
    if constexpr (subdim == dim) {
        nodes_.reserve(tri.size());
        for (const auto &simplex : tri.simplices()) {
            nodes_.emplace_back(simplex);
        }
    } else {
        nodes_.reserve(tri.template countFaces<subdim>());
        for (const auto &face : tri.template faces<subdim>()) {
            nodes_.emplace_back(face);
        }
    }

    // Make an edge for each pair of subdim-faces that share a facet
    // (including self-folds)
    std::vector<std::vector<std::pair<size_t, int>>> buckets(
        tri.template countFaces<subdim - 1>());
    for (size_t fi = 0; fi < nodes_.size(); ++fi) {
        Face *face = nodes_[fi].face;
        for (int i = 0; i <= subdim; ++i) {
            Facet *facet = face->template face<subdim - 1>(i);
            buckets[facet->index()].emplace_back(fi, i);
        }
    }

    for (const auto &bucket : buckets) {
        for (const auto &[fi, i] : bucket) {
            for (const auto &[fj, j] : bucket) {
                // Skip only comparing an occurrence to itself; a
                // genuine self-fold (fi == fj, i != j) is kept.
                if (fi == fj && i == j) {
                    continue;
                }

                Face *src = nodes_[fi].face;
                Face *dst = nodes_[fj].face;

                regina::Perm<dim + 1> m1 =
                    src->template faceMapping<subdim - 1>(i);
                regina::Perm<dim + 1> m2 =
                    dst->template faceMapping<subdim - 1>(j);

                std::array<int, subdim + 1> gluing;
                for (int k = 0; k <= subdim; ++k) {
                    gluing[m1[k]] = m2[k];
                }

                nodes_[fi].gluings.emplace_back(
                    fi, i, fj, regina::Perm<subdim + 1>(gluing));
            }
        }
    }
}

template class Skeleton<3, 2>;
template class Skeleton<4, 2>;
