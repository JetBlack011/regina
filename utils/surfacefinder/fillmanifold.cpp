//
//  surfacefinder.cpp
//
//  Created by John Teague on 06/19/2024.

#include <gmpxx.h>
#include <link/link.h>
#include <optional>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <unistd.h>

#include <cassert>
#include <iostream>
#include <ostream>
#include <string>

#include "gluing.h"
#include "triangulation/example2.h"
#include "triangulation/example3.h"
#include "triangulation/forward.h"
#include "triangulation/generic/boundarycomponent.h"
#include "triangulation/generic/isomorphism.h"
#include "triangulation/generic/triangulation.h"
#include "utilities/exception.h"

namespace {
void usage(const char *progName, const std::string &error = std::string()) {
    if (!error.empty())
        std::cerr << error << "\n\n";

    std::cerr << "Usage:\n";
    std::cerr << "    " << progName << " {dim} <isosig>\n"
              << "    dim  : Specify the dimension of the given isomorphism "
                 "signature (optional, default is 3)\n"
              << "    " << progName << " [ -v, --version | -h, --help ]\n\n";
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

template <int dim>
std::optional<regina::Triangulation<dim + 1>>
thickenHelper(regina::Triangulation<dim> &tri) {
    regina::Triangulation<dim + 1> thickenedTri;
    std::vector<SimplicialPrism<dim + 1>> prisms;
    prisms.reserve(tri.size());
    for (int i = 0; i < tri.size(); ++i) {
        prisms.emplace_back(thickenedTri);
    }
    std::set<std::pair<const regina::Simplex<dim> *, int>> visited;

    int counter = 0;
    for (int i = 0; i < prisms.size(); ++i) {
        for (int f = 0; f < dim + 1; ++f) {
            const regina::Simplex<dim> *s = tri.simplex(i);
            const regina::Simplex<dim> *adj = s->adjacentSimplex(f);

            if (adj == nullptr || visited.contains({s, f}))
                continue;

            int adjFacet = s->adjacentFacet(f);

            prisms[i].glue(f, prisms[adj->index()], adjFacet);
            ++counter;

            visited.insert({s, f});
            visited.insert({adj, adjFacet});
        }
    }

    return thickenedTri;
}
template <int dim>
std::optional<regina::Triangulation<dim + 1>>
thicken(regina::Triangulation<dim> &tri) {
    thickenHelper<dim>(tri);
}

template <>
std::optional<regina::Triangulation<3>> thicken(regina::Triangulation<2> &tri) {
    //std::cout << "Orienting...\n";
    //tri.orient();
    for (int i = 0; i < tri.simplices().size(); ++i) {
        for (int f = 0; f < 3; ++f) {
            if (tri.simplex(i)->adjacentSimplex(f) == nullptr)
                continue;
            std::cout << "Simplex " << i << " glued along face " << f << " by " << tri.simplex(i)->adjacentGluing(f) << " to " << tri.simplex(i)->adjacentSimplex(f)->index() << "\n";
        }
    }
    return thickenHelper(tri);
}

template <>
std::optional<regina::Triangulation<4>> thicken(regina::Triangulation<3> &tri) {
    if (!tri.order())
        return std::nullopt;
    return thickenHelper(tri);
}

template <int dim>
regina::Triangulation<dim>
glue(const regina::Triangulation<dim> &tri1, int bdryIndex1,
     const regina::Triangulation<dim> &tri2, int bdryIndex2,
     const regina::Isomorphism<dim - 1> &iso) {
    regina::Triangulation<dim> tri;
    tri.insertTriangulation(tri1);
    tri.insertTriangulation(tri2);
    std::vector<Gluing<dim, dim>> gluings;
    const regina::BoundaryComponent<dim> *bdry1 =
        tri.boundaryComponent(bdryIndex1);
    const regina::BoundaryComponent<dim> *bdry2 =
        tri.boundaryComponent(tri1.countBoundaryComponents() + bdryIndex2);

    for (int i = 0; i < bdry1->size(); ++i) {
        // Since these are boundary faces, they belong to exactly 1 simplex
        const auto &emb1 = bdry1->facet(i)->front();
        const auto &emb2 = bdry2->facet(iso.simpImage(i))->front();
        regina::Simplex<dim> *s1 = emb1.simplex();
        regina::Simplex<dim> *s2 = emb2.simplex();

        regina::Perm<dim + 1> i1 = emb1.vertices();
        regina::Perm<dim + 1> i2 = emb2.vertices();
        regina::Perm<dim> p = iso.facetPerm(i);
        std::array<int, dim + 1> isoPerm;

        for (int j = 0; j < dim; ++j) {
            isoPerm[j] = p[j];
        }
        isoPerm[dim] = dim;

        gluings.push_back({s1, emb1.face(), s2, i2 * isoPerm * i1.inverse()});

        // std::cout << "bdryIndex = " << i << " to " << iso.simpImage(i) << ",
        // IsoPerm = " << p
        //           << ", i1 = " << i1 << ", i1^-1 = " << i1.inverse()
        //           << ", i2 = " << i2 << "\n";
    }

    for (const auto &gluing : gluings) {
        // std::cout << "Gluing " << gluing.src->index() << " and "
        //           << gluing.dst->index() << " along " << gluing.srcFacet <<
        //           ", "
        //           << gluing.gluing << "\n";
        gluing.src->join(gluing.srcFacet, gluing.dst, gluing.gluing);
    }

    return tri;
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

    std::string isoSig;
    int dim = 3;

    if (argc == 2) {
        isoSig = argv[1];
    } else if (argc == 3) {
        isoSig = argv[2];
        dim = std::stoi(argv[1]);
    } else {
        usage(argv[0],
              "Please provide an isomorphism signature and/or dimension.");
    }

    regina::Triangulation<3> threeMfld = regina::Example<3>::sphere600();
    std::cout << "Original = " << threeMfld.isoSig() << "\n";
    regina::Triangulation<4> tri1 = thicken(threeMfld).value();
    regina::Triangulation<4> tri2 = cone(threeMfld);
    // regina::Triangulation tri = regina::Example<2>::sphere();
    // tri.subdivide();
    //regina::Triangulation<2> sphere = regina::Example<2>::sphere();
    //regina::Triangulation<3> tri1 = thicken(sphere).value();
    //tri1.subdivide();
    //regina::Triangulation<3> tri2(tri1);

    //std::cout << "Tri1 = " << tri1.isoSig() << "\n";

    // Find all isomorphisms between the two triangulations
    std::vector<regina::Isomorphism<3>> isos;
    tri1.boundaryComponent(0)->build().findAllIsomorphisms(
        tri2.boundaryComponent(0)->build(),
        [&isos](const regina::Isomorphism<3> &iso) {
            isos.push_back(iso);
            return true;
        });
    std::cout << "[+] Found " << isos.size() << " isomorphisms\n";
    for (const auto &iso : isos) {
        std::cout << iso << "\n";
    }
    std::cout << "\n";
    regina::Triangulation<4> tri = glue(tri1, 0, tri2, 0, isos[0]);

    //std::cout << "Tri = " << tri.isoSig() << "\n";
    //if (tri.isValid()) {
    //    std::cout << "[+] Triangulation is valid\n";
    //} else {
    //    std::cout << "[!] Triangulation is NOT valid\n";
    //}

    //isos.clear();
    //tri.boundaryComponent(0)->build().findAllIsomorphisms(
    //    tri2.boundaryComponent(0)->build(),
    //    [&isos](const regina::Isomorphism<3> &iso) {
    //        isos.push_back(iso);
    //        return true;
    //    });
    //std::cout << "[+] Found " << isos.size() << " isomorphisms\n";
    //for (const auto &iso : isos) {
    //    std::cout << iso << "\n";
    //}
    //std::cout << "\n";
    //regina::Triangulation<4> thickerTri = glue(tri, 0, tri2, 0, isos[0]);

    //std::cout << "Tri = " << thickerTri.isoSig() << "\n";
    //if (tri.isValid()) {
    //    std::cout << "[+] Triangulation is valid\n";
    //} else {
    //    std::cout << "[!] Triangulation is NOT valid\n";
    //}

    // auto conedTri = cone(tri);

    // std::cout << "[+] Number of simplices = " << tri.simplices().size() <<
    // "\n"; std::cout << "[+] Number of cone simplices = "
    //           << conedTri.simplices().size() << "\n";

    // std::vector<regina::Isomorphism<2>> isos;

    // auto res = thicken(tri);
    // if (!res.has_value()) {
    //     std::cout << "[!] Given triangulation could not be thickened! Make "
    //                  "sure it's possible to order it (e.g. make sure the link
    //                  " "diagram you gave to knotbuilder is alternating).\n";
    //     return -1;
    // }
    // const auto &thickTri = res.value();
    // regina::Triangulation evenThickerTri(thickTri);

    // evenThickerTri.boundaryComponent(0)->build().findAllIsomorphisms(
    //     thickTri.boundaryComponent(1)->build(),
    //     [&isos](const regina::Isomorphism<2> &iso) {
    //         isos.push_back(iso);
    //         return true;
    //     });

    // for (int i = 0; i < 2; ++i) {
    //     std::cout << isos[0] << "\n";
    //     evenThickerTri = glue(evenThickerTri, 0, thickTri, 1, isos[0]);
    // }

    // isos.clear();
    // evenThickerTri.boundaryComponent(0)->build().findAllIsomorphisms(
    //     conedTri.boundaryComponent(0)->build(),
    //     [&isos](const regina::Isomorphism<2> &iso) {
    //         isos.push_back(iso);
    //         return false;
    //     });

    //// regina::Triangulation<4> thickenedTri =
    //// regina::Example<4>::iBundle(tri);

    // std::cout << "[+] Number of thickened simplices = "
    //           << thickTri.simplices().size() << "\n";
    // std::cout << "[+] Number of even more thickened simplices = "
    //           << evenThickerTri.simplices().size() << "\n";
    // std::cout << "[+] Total simplices after gluing = "
    //           << conedTri.simplices().size() + thickTri.simplices().size()
    //           << "\n\n";
    // if (thickTri.isValid()) {
    //     std::cout << "[+] Thickened triangulation is valid\n";
    // } else {
    //     std::cout << "[!] Thickened triangulation is NOT valid\n";
    // }
    // if (evenThickerTri.isValid()) {
    //     std::cout << "[+] Even thicker triangulation is valid\n";
    // } else {
    //     std::cout << "[!] Even thicker triangulation is NOT valid\n";
    // }

    // std::cout << "[+] Cone isosig = " << conedTri.isoSig() << "\n\n";
    // std::cout << "[+] Thickened isosig = " << thickTri.isoSig() << "\n\n";
    // std::cout << "[+] Even thicker isosig = " << evenThickerTri.isoSig()
    //           << "\n\n";

    // for (const auto &iso : isos) {
    //     auto gluedTri = glue(evenThickerTri, 0, conedTri, 0, iso);

    //    if (gluedTri.isValid()) {
    //        std::cout << "[+] Thickened + coned triangulation is valid\n";
    //    } else {
    //        std::cout << "[!] Thickened + coned triangulation is NOT valid\n";
    //    }
    //    std::cout << "ISO SIG = " << gluedTri.isoSig() << "\n";
    //}
    return 0;
}
