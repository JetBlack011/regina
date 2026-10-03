// coneover.h -- the cone over a triangulation, for the local-flatness tests.
//
// CobordismBuilder<dim>::cone() built this until coning was retired from the
// search (phase 5 of the refactor); the tests that cone a knot's S³ to a
// point, to make an apex whose link is the knot, keep this small copy. One
// (dim+1)-simplex per base simplex, in the base's index order, its local
// vertices 0..dim the base simplex's own and local vertex dim+1 the shared
// apex, glued mirroring the base's gluings. The base is first ordered as
// CobordismBuilder orders it (taking an ordered copy), so the result is the
// triangulation cone() gave, simplex for simplex.

#ifndef SURFER_TESTS_CONEOVER_H
#define SURFER_TESTS_CONEOVER_H

#include <array>
#include <maths/perm.h>

#include "diagramtriangulation/thickening/thickening.h"

template <int dim>
regina::Triangulation<dim + 1> coneOver(const regina::Triangulation<dim> &tri) {
    CobordismBuilder<dim> ordered(tri);
    const regina::Triangulation<dim> &base = ordered.baseTriangulation();

    regina::Triangulation<dim + 1> cone;
    for (size_t i = 0; i < base.size(); ++i)
        cone.newSimplex();

    for (const auto *s : base.simplices()) {
        regina::Simplex<dim + 1> *apexed = cone.simplex(s->index());
        for (int f = 0; f <= dim; ++f) {
            const regina::Simplex<dim> *adj = s->adjacentSimplex(f);
            if (adj == nullptr || apexed->adjacentSimplex(f) != nullptr)
                continue;
            regina::Perm<dim + 1> gluing = s->adjacentGluing(f);
            std::array<int, dim + 2> image;
            for (int i = 0; i <= dim; ++i)
                image[i] = gluing[i];
            image[dim + 1] = dim + 1;
            apexed->join(f, cone.simplex(adj->index()),
                         regina::Perm<dim + 2>(image));
        }
    }
    return cone;
}

#endif // SURFER_TESTS_CONEOVER_H
