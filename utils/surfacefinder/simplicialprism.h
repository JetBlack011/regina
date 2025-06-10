#ifndef SIMPLICIALPRISM_H

#define SIMPLICIALPRISM_H

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
            simplices_[i - 1]->join(i, simplices_[i], id);
        }
    }

    template <int facedim>
    std::array<regina::Face<dim, facedim + 1> *, facedim + 1>
    faceTimesI(const regina::Perm<dim> &vertices) {
        static_assert(
            facedim >= 0 && facedim < dim - 1,
            "SimplicialPrism::faceTimesI: facedim must be in [0, dim-2]");

        std::array<regina::Face<dim, facedim + 1> *, facedim + 1> faces;

        for (int i = 0; i < facedim + 1; ++i) {
            // Vertex nonsense
            std::unordered_set<int> unusedVertices;
            for (int j = 0; j < dim + 1; ++j) {
                unusedVertices.insert(j);
            }
            std::array<int, dim + 1> faceVertices;

            int j = 0;
            for (; j < facedim - i; ++j) {
                faceVertices[j] = vertices[j];
                unusedVertices.erase(faceVertices[j]);
            }
            for (; j < facedim + 1; ++j) {
                int vj = vertices[j - (facedim - i)];
                faceVertices[j] = vj == 0 ? dim : vj - 1;
                unusedVertices.erase(faceVertices[j]);
            }
            for (int v : unusedVertices) {
                faceVertices[j++] = v;
            }

            int faceIndex =
                regina::Face<dim, facedim + 1>::faceNumber(faceVertices);
            faces[i] = simplices_[vertices[i]]->face(faceIndex);
        }

        return faces;
    }

    void glue(int face, SimplicialPrism<dim> &other, int otherFace) {
        // Not a fan.... Must be a better way. Maybe prisms just suck.
        if (!(0 <= face && face <= dim && 0 <= otherFace && otherFace <= dim))
            throw regina::InvalidArgument(
                "SimplicialPrism::glue(): Invalid faces");

        int i = 0;
        int j = 0;

        while (i < dim && j < dim) {
            if (i == face)
                ++i;
            if (j == otherFace)
                ++j;

            if (i >= dim || j >= dim)
                break;

            int k = 0;
            int l = 0;
            int k0 = dim;
            int l0 = dim;
            std::array<int, dim + 1> gluing;

            while (k < dim + 1 && l < dim + 1) {
                if ((k <= i && k == face) || (k > i && k == face + 1))
                    k0 = k++;
                if ((l <= j && l == otherFace) || (l > j && l == otherFace + 1))
                    l0 = l++;

                if (k >= dim + 1 || l >= dim + 1)
                    break;

                gluing[k] = l;
                ++k;
                ++l;
            }

            gluing[k0] = l0;

            simplices_[i]->join(k0, other.simplices_[j], gluing);

            ++i;
            ++j;
        }
    }
};

#endif // SIMPLICIALPRISM_H
