//
//  surfacefinder.cpp
//
//  Created by John Teague on 06/19/2024.

#include <gmpxx.h>
#include <link/link.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <unistd.h>

#include <cassert>
#include <iostream>
#include <ostream>
#include <string>

#include "triangulation/example2.h"
#include "triangulation/example3.h"
#include "triangulation/forward.h"
#include "triangulation/generic/isomorphism.h"
#include "triangulation/generic/triangulation.h"
#include "utilities/exception.h"

#include "gluing.h"

namespace {
void usage(const char *progName, const std::string &error = std::string()) {
    if (!error.empty())
        std::cerr << error << "\n\n";

    std::cerr << "Usage:\n";
    std::cerr << "    " << progName << " <isosig>\n    " << progName
              << " [ -v, --version | -h, --help ]\n\n";
    std::cerr << "    -v, --version  : Show which version of Regina is "
                 "being used\n";
    std::cerr << "    -h, --help     : Display this help\n";
    exit(1);
}

template <int dim>
regina::Triangulation<dim + 1> cone(regina::Triangulation<dim> &bdryTri) {
    regina::Triangulation<dim + 1> coneTri;

    for (auto s : bdryTri.simplices()) {
        coneTri.newSimplex();
    }

    for (int i = 0; i < bdryTri.simplices().size(); ++i) {
        auto bdrySimplex = bdryTri.simplex(i);
        auto coneSimplex = coneTri.simplex(i);

        for (int f = 0; f <= dim; ++f) {
            if (bdrySimplex->adjacentSimplex(f) == nullptr ||
                coneSimplex->adjacentSimplex(f) != nullptr)
                continue;

            regina::Perm<dim + 1> bdryGluing = bdrySimplex->adjacentGluing(f);
            regina::Simplex<dim + 1> *adjConeSimplex =
                coneTri.simplex(bdrySimplex->adjacentSimplex(f)->index());
            std::array<int, dim + 2> coneGluing;
            for (int i = 0; i < dim + 1; ++i) {
                coneGluing[i] = bdryGluing[i];
            }
            coneGluing[dim + 1] = dim + 1;

            coneSimplex->join(f, adjConeSimplex, coneGluing);
        }
    }

    return coneTri;
}

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

    void glue(int face, SimplicialPrism<dim> &other, int otherFace) {
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

template <int dim>
regina::Triangulation<dim + 1> thicken(const regina::Triangulation<dim> &tri) {
    regina::Triangulation<dim + 1> thickenedTri;

    std::vector<SimplicialPrism<dim + 1>> prisms;

    for (auto s : tri.simplices()) {
        prisms.emplace_back(SimplicialPrism(thickenedTri));
    }

    std::set<std::pair<const regina::Simplex<dim> *, int>> visited;

    for (int i = 0; i < prisms.size(); ++i) {
        for (int f = 0; f < dim + 1; ++f) {
            const regina::Simplex<dim> *s = tri.simplex(i);
            const regina::Simplex<dim> *adj = s->adjacentSimplex(f);

            if (adj == nullptr || visited.contains({s, f}))
                continue;

            int adjFacet = s->adjacentFacet(f);

            prisms[i].glue(f, prisms[adj->index()], adjFacet);

            visited.insert({s, f});
            visited.insert({adj, adjFacet});
        }
    }

    return thickenedTri;
}

template <int dim>
regina::Triangulation<dim + 1> glue(const regina::Triangulation<dim> &tri1,
                                    const regina::Triangulation<dim> &tri2,
                                    const regina::Isomorphism<dim> &iso) {

}
} // namespace

int main(int argc, char *argv[]) {
    // Check for standard arguments:
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help")
            usage(argv[0]);
        if (arg == "-v" || arg == "--version") {
            if (argc != 2)
                usage(argv[0], "Option --version cannot be used with "
                               "any other arguments.");
            std::cout << PACKAGE_BUILD_STRING << "\n";
            exit(0);
        }
    }

    if (argc != 2) {
        usage(argv[0], "Please provide an isomorphism signature.");
    }

    std::string isoSig = argv[1];
    regina::Triangulation<3> tri1(isoSig);
    regina::Triangulation<3> tri2(isoSig);

    std::vector<regina::Isomorphism<3>> isos;

    tri1.findAllIsomorphisms(tri2, [&isos](const regina::Isomorphism<3> &iso) {
        isos.push_back(iso);
        return false;
    });

    for (const auto &iso : isos) {
        std::cout << iso << "\n\n";
    }

    ////regina::Triangulation<2> tri = regina::Example<2>::torus();

    // auto coneTri = cone(tri);

    // std::cout << "[+] Number of simplices = " << tri.simplices().size() <<
    // "\n"; std::cout << "[+] Number of cone simplices = " <<
    // coneTri.simplices().size()
    //           << "\n\n\n";

    // auto thickenedTri = thicken(tri);

    // std::cout << "[+] Number of thickened simplices = "
    //           << thickenedTri.simplices().size() << "\n\n";
    // std::cout << "Cone isosig = " << coneTri.isoSig() << "\n";
    // std::cout << "Thickened isosig = " << thickenedTri.isoSig() << "\n\n";

    // std::cout << "Thickened boundaries: ";

    // for (auto comp : thickenedTri.boundaryComponents()) {
    //     std::cout << comp->build().isoSig() << "\n\n";
    // }

    // regina::Triangulation<4> gluedTri;

    // regina::Perm<5> id;

    // coneTri.moveContentsTo(gluedTri);
    // thickenedTri.moveContentsTo(gluedTri);

    // int totalSimplices = tri.simplices().size();
    // for (int i = 0; i < totalSimplices; ++i) {
    //     gluedTri.simplex(i)->join(4, gluedTri.simplex(totalSimplices + 3 + i
    //     * 4), id);
    // }

    // std::cout << "ISO SIG = " << gluedTri.isoSig() << "\n";

    return 0;
}
