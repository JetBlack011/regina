//
//  surfacefinder.cpp
//
//  Created by John Teague on 06/19/2024.
//

#include <gmpxx.h>
#include <link/link.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <unistd.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "triangulation/forward.h"
#include "triangulation/generic/boundarycomponent.h"
#include "triangulation/generic/triangulation.h"

#include "gluinggraph.h"
#include "knottedsurfaces.h"

namespace {
/* Constants */
static const regina::Triangulation<2> GENUS_2_SURFACE =
    regina::Triangulation<2>::fromGluings(6, {{0, 2, 5, {2, 1, 0}},
                                              {0, 1, 1, {0, 2, 1}},
                                              {0, 0, 5, {1, 0, 2}},
                                              {1, 1, 2, {0, 2, 1}},
                                              {1, 0, 3, {0, 2, 1}},
                                              {2, 1, 3, {0, 2, 1}},
                                              {2, 0, 4, {0, 2, 1}},
                                              {3, 1, 4, {0, 2, 1}},
                                              {4, 1, 5, {0, 2, 1}}});

void usage(const char *progName, const std::string &error = std::string()) {
    if (!error.empty())
        std::cerr << error << "\n\n";

    std::cerr << "Usage:\n";
    std::cerr << "    " << progName
              << " { -a, --all | -b, --boundary | -c, --closed | -l, --links } "
                 "<isosig>\n"
                 "    "
              << progName << " [ -v, --version | -h, --help ]\n\n";
    std::cerr << "    -a, --all      : Find all surfaces, regardless of "
                 "boundary conditions\n";
    std::cerr << "    -b, --boundary : Find surfaces such that their boundary is "
                 "contained entirely in\n"
                 "                     the boundary of the given 4-manifold (Note, "
                 "if the 4-manifold\n                     is closed, this is "
                 "equivalent to --closed)\n";
    std::cerr << "    -c, --closed   : Find only closed surfaces in the given "
                 "4-manifold\n";
    std::cerr << "    -l, --links    : Same as --boundary, but also gives "
                 "a list of link types which\n"
                 "                     are the boundaries of the surfaces "
                 "we find\n\n";
    std::cerr << "    -v, --version  : Show which version of Regina is "
                 "being used\n";
    std::cerr << "    -h, --help     : Display this help\n";
    exit(1);
}

template <int dim>
void surfacesDetail(std::set<KnottedSurface<dim>> &surfaces, SurfaceCondition cond) {
    std::cout << "--- "
              << (cond == SurfaceCondition::all
                      ? ""
                      : (cond == SurfaceCondition::boundary
                             ? "Proper "
                             : (cond == SurfaceCondition::closed ? "Closed " : "")))
              << "Surfaces ---\n";
    int surfaceCount = 0;
    int closedCount = 0;

    std::map<std::string, int> countMap;

    for (const auto &surface : surfaces) {
        ++surfaceCount;
        if (surface.surface().isClosed()) {
            ++closedCount;
        }

        // For checking if properness works
        bool isProper = true;
        for (const regina::BoundaryComponent<2> *comp : surface.surface().boundaryComponents()) {
            for (const regina::Edge<2> *edge : comp->edges()) {
                if (edge->isBoundary() && !surface.image(edge)->isBoundary()) {
                    std::cout << "NOTE: " << surface.detail() << " is not proper!!\n";
                    isProper = false;
                    return;
                }
            }
        }

        // if (isProper) {
        //     std::cout << surface.detail() << " is PROPER!"
        //               << "\n";
        // }

        auto search = countMap.find(surface.detail());
        if (search == countMap.end()) {
            countMap.insert({surface.detail(), 1});
        } else {
            ++search->second;
        }
    }

    std::cout << "\n";

    std::vector<std::string> descriptions;
    descriptions.reserve(countMap.size());
    for (const auto &[description, count] : countMap) {
        descriptions.push_back(description);
    }

    std::sort(descriptions.begin(), descriptions.end(),
              [](const std::string &first, const std::string &second) {
                  return first.size() < second.size() ||
                         (first.size() == second.size() && first < second);
              });

    for (std::string &description : descriptions) {
        int count = countMap.find(description)->second;
        std::cout << "- " << count << (count == 1 ? " copy of " : " copies of ") << description
                  << "\n";
    }

    std::cout << "\n";
    if (cond != SurfaceCondition::closed) {
        std::cout << "Total admissible surfaces = " << surfaceCount << "\n";
    }
    std::cout << "Total closed surfaces = " << closedCount << "\n\n";
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

    if (argc != 3) {
        usage(argv[0], "Please specify a surface condition and provide an "
                       "isomorphism signature.");
    }

    std::string arg = argv[1];
    std::string isoSig = argv[2];
    SurfaceCondition cond;

    if (arg == "-a" || arg == "--all") {
        cond = SurfaceCondition::all;
    } else if (arg == "-b" || arg == "--boundary" || arg == "-l" || arg == "--links") {
        cond = SurfaceCondition::boundary;
    } else if (arg == "-c" || arg == "--closed") {
        cond = SurfaceCondition::closed;
    } else {
        usage(argv[0], "Please specify a valid surface condition.");
    }

    //    regina::Triangulation<4> tri;
    //    tri.newSimplex();
    //    tri.pachner(tri.pentachoron(0));
    //
    // regina::Triangulation<3> tri("caba");
    // for (const regina::Tetrahedron<3> *tet : tri.tetrahedra()) {
    //     std::cout << "Tetrahedron " << tet->index() << "\n";
    //     for (int i = 0; i < 4; ++i) {
    //         std::cout << "Vertex " << i << " = " << tet->vertex(i)->index()
    //                   << "\n";
    //     }
    //     std::cout << "\n";
    // }
    // std::vector<int> edgeIndices = {2, 0, 7, 8, 5};
    // Curve<3> c;
    // for (int i : edgeIndices) {
    //     c.push_back(tri.edge(i));
    // }
    // Knot k(tri, c);
    // std::cout << k << "\n";
    // k.simplify();
    // std::cout << k;

    regina::Triangulation<4> tri(isoSig);

    // std::cout << "Original = " << tri.isoSig() << "\n";
    // Knot k = {tri, {}};
    // std::vector<regina::Tetrahedron<3> *> tets;
    // for (regina::Tetrahedron<3> *tet : k.tri_.tetrahedra()) {
    //     tets.push_back(tet);
    // }

    // k.subdivideSharedVertexSequence_(tets);
    // std::cout << "Subdivided = " << k.tri_.isoSig() << "\n";
    //  tri.newSimplex();

    // regina::Triangulation<2> tri = regina::Example<2>::orientable(6, 0);
    // tri.newSimplex();
    // tri.subdivide();
    // tri.subdivide();
    //    tri.subdivide();

    //auto tri = regina::Triangulation<3>("caba");
    //std::cout << tri.detail() << "\n";

    //auto k = Knot(tri, {tri.edge(2), tri.edge(6), tri.edge(7), tri.edge(3), tri.edge(5)});
    //std::cout << k << "\n";
    //k.shrink_();
    //std::cout << k << "\n";

    GluingGraph graph(tri, cond);
    auto &surfaces = graph.findSurfaces();

    surfacesDetail(surfaces, cond);

    for (const auto *comp : tri.boundaryComponents()) {
        for (const auto *v : comp->vertices()) {
            std::cout << v->index() << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    int numNonUnlinks = 0;
    for (auto surface : surfaces) {
        if (surface.surface().isClosed())
            continue;
        Link l = surface.boundary();
        // std::cout << l << "\n";
        // l.simplify();
        if (!l.isUnlink()) {
            ++numNonUnlinks;
            std::cout << "DETAIL: " << surface.detail() << "\n";
            std::cout << l << "\n";
        }
    }

    std::cout << "NUMBER OF NONTRIVIAL LINKS = " << numNonUnlinks;

    return 0;
}
