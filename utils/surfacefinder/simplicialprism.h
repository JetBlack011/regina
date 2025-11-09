#ifndef SIMPLICIALPRISM_H

#define SIMPLICIALPRISM_H

#include <algorithm>
#include <gmpxx.h>
#include <link/link.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <unistd.h>
#include <unordered_set>

#include "triangulation/forward.h"

template <int dim>
class SimplicialPrism {
  private:
    std::vector<regina::Simplex<dim> *> simplices_;

  public:
    SimplicialPrism(regina::Triangulation<dim> &tri) {
        for (int i = 0; i < dim; ++i) {
            simplices_.push_back(tri.newSimplex());
        }

        regina::Perm<dim + 1> id;

        for (int i = 1; i < dim; ++i) {
            simplices_[i - 1]->join(i - 1, simplices_[i], id);
        }
    }

    template <int facedim>
    std::array<const regina::Face<dim, facedim + 1> *, facedim + 1>
    subprism(const regina::Perm<dim> &vertices) const {
        static_assert(
            facedim >= 0 && facedim < dim - 1,
            "SimplicialPrism::faceTimesI: facedim must be in [0, dim-2]");

        std::vector<int> faceVerts;
        for (int i = 0; i < dim; ++i) {
            faceVerts.push_back(vertices[i]);
        }
        std::sort(faceVerts.begin(), faceVerts.begin() + facedim + 1);

        std::array<const regina::Face<dim, facedim + 1> *, facedim + 1> faces;

        for (int i = 0; i <= facedim; ++i) {
            std::array<int, dim + 1> facetVerts;
            std::unordered_set<int> unusedVerts;
            for (int j = 0; j <= dim; ++j) {
                unusedVerts.insert(j);
            }

            int j = 0;
            for (; j <= i; ++j) {
                facetVerts[j] = faceVerts[j] == 0 ? dim : faceVerts[j] - 1;
                unusedVerts.erase(facetVerts[j]);
            }
            for (; j <= facedim + 1; ++j) {
                facetVerts[j] = faceVerts[j - 1];
                unusedVerts.erase(facetVerts[j]);
            }
            for (; j <= dim; ++j) {
                facetVerts[j] = *unusedVerts.begin();
                unusedVerts.erase(unusedVerts.begin());
            }

            int facetIndex =
                regina::Face<dim, facedim + 1>::faceNumber(facetVerts);
            faces[i] = simplices_[faceVerts[i]]->template face<facedim + 1>(
                facetIndex);
        }

        return faces;
    }

    inline std::array<const regina::Triangle<dim> *, 2>
    edgeTimesI(const regina::Perm<dim> &vertices) const {
        return subprism<1>(vertices);
    }

    void glue(int facet, SimplicialPrism<dim> &other, int otherFacet) {
        // Positively heinous code
        if (!(0 <= facet && facet <= dim && 0 <= otherFacet &&
              otherFacet <= dim))
            throw regina::InvalidArgument(
                "SimplicialPrism::glue(): Invalid faces");

        int i = 0, j = 0;
        while (i < dim && j < dim) {
            if (i == facet)
                ++i;
            if (j == otherFacet)
                ++j;

            if (i == dim || j == dim)
                break;

            std::array<int, dim + 1> g;
            int k = 0, l = 0;
            int k0 = dim, l0 = dim;
            while (k <= i) {
                if (k == facet) {
                    k0 = k == 0 ? dim : k - 1;
                    ++k;
                }
                if (l == otherFacet) {
                    l0 = l == 0 ? dim : l - 1;
                    ++l;
                }

                g[k == 0 ? dim : k - 1] = l == 0 ? dim : l - 1;
                ++k, ++l;
            }
            k--, l--;
            while (k < dim) {
                if (k == facet)
                    k0 = k++;
                if (l == otherFacet)
                    l0 = l++;

                g[k++] = l++;
            }
            g[k0] = l0;

            std::cout << "Triangulation = "
                      << simplices_[i]->triangulation().isoSig() << "\n";
            std::cout << "facet = " << facet << ", otherFacet = " << otherFacet
                      << ", i = " << i << ", j = " << j << ", k0 = " << k0
                      << ", l0 = " << l0 << ", g = {";
            for (int m = 0; m < dim + 1; ++m) {
                std::cout << g[m];
                if (m < dim)
                    std::cout << ", ";
            }
            std::cout << "}\n";
            simplices_[i]->join(k0, other.simplices_[j], g);

            ++i, ++j;
        }
    }
};

#endif // SIMPLICIALPRISM_H
