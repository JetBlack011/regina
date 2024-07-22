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

#include <cstring>
#include <iostream>
#include <ostream>
#include <string>
#include <tuple>
#include <utility>

#include "surfaceknots.h"
#include "gluinggraph.h"
#include "triangulation/example3.h"
#include "triangulation/example4.h"
#include "triangulation/forward.h"
#include "triangulation/generic/boundarycomponent.h"
#include "triangulation/generic/triangulation.h"

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




// template <int n>
// regina::Triangulation<3> knotComplement(const regina::Triangulation<3> t,
//                                         const Knot &k) {
//     return linkComplement(t, {k});
// }

template <typename T, typename D>
std::ostream &operator<<(std::ostream &os, const std::pair<T, D> &p) {
    os << "(" << p.first << ", " << p.second << ")";
    return os;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const std::tuple<T, T, T> &p) {
    os << "(" << std::get<0>(p) << ", " << std::get<1>(p) << ", "
       << std::get<2>(p) << ")";
    return os;
}

void usage(const char *progName, const std::string &error = std::string()) {
    if (!error.empty()) std::cerr << error << "\n\n";

    std::cerr << "Usage:\n";
    std::cerr << "    " << progName
              << " [ -c, --closed | -b, --boundary | -l, --links ] <isosig>\n"
                 "    "
              << progName << " [ -v, --version | -h, --help ]\n\n";
    std::cerr << "    -c, --closed   : Find closed surfaces in the given "
                 "4-manifold\n";
    std::cerr << "    -b, --boundary : Find surfaces such that its boundary is "
                 "contained entirely in\n"
                 "                     the boundary of the given 4-manifold\n";
    std::cerr
        << "    -l, --links    : Same as --boundary, but also gives a list of "
           "link types which\n"
           "                     are the boundaries of the surfaces we find\n";
    std::cerr
        << "    -v, --version  : Show which version of Regina is being used\n";
    std::cerr << "    -?, --help     : Display this help\n";
    exit(1);
}

int main(int argc, char *argv[]) {
    // Check for standard arguments:
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
            usage(argv[0]);
        if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--version") == 0) {
            if (argc != 2)
                usage(argv[0],
                      "Option --version cannot be used with "
                      "any other arguments.");
            std::cout << PACKAGE_BUILD_STRING << "\n";
            exit(0);
        }
    }

    if (argc != 3) {
        usage(argv[0],
              "Please specify a surface condition and provide an "
              "isomorphism signature.");
    }

    std::string isoSig = argv[2];
    SurfaceCondition cond;

    if (strcmp(argv[1], "-c") == 0 || strcmp(argv[1], "--closed") == 0) {
        cond = SurfaceCondition::closed;
    } else if (strcmp(argv[1], "-b") == 0 ||
               strcmp(argv[1], "--boundary") == 0) {
        cond = SurfaceCondition::boundary;
    } else if (strcmp(argv[1], "-l") == 0 || strcmp(argv[1], "--links") == 0) {
        cond = SurfaceCondition::links;
    } else {
        usage(argv[0], "Please specify a valid surface condition.");
    }

    // regina::Triangulation<4> tri(isoSig);

    // regina::Triangulation<3> bdry = tri.boundaryComponent(0)->build();

    // auto tri = regina::Example<4>::fourSphere();
    regina::Triangulation<4> tri(isoSig);
    //regina::Triangulation<3> tri("bkaagj");
    std::cout << tri.detail() << "\n\n";
    // Knot<3> k;

    // knotComplement(tri, k);

    std::cout << "Boundary = ";
    for (const auto comp : tri.boundaryComponents()) {
        for (const auto edge : comp->edges()) {
            std::cout << edge->index() << " ";
        }
    }
    std::cout << "\n\n";

    std::cout << "--- Closed Surfaces ---\n";
    int surfaceCount = 0;
    int closedCount = 0;
    GluingGraph graph(&tri, cond);

    auto surfaceKnots = graph.findSurfaceKnots();
    std::map<std::string, int> countMap;

    for (const auto &surfaceKnot : surfaceKnots) {
        ++surfaceCount;
        if (surfaceKnot.surface().isClosed()) {
            ++closedCount;
        } else {
            // std::cout << "Surface boundary = ";
            for (const auto comp :
            surfaceKnot.surface().boundaryComponents())
            {
                for (const auto edge : comp->edges()) {
                    // std::cout << surfaceKnot.findEdgeInTri(edge)->index() << "
                    // ";
                }
            }
            // regina::snappea::DualOneSkeletonCurve curve;
            // std::cout << "\n";
        }
        // std::cout << " - " << surfaceKnot << "\n";
        auto search = countMap.find(surfaceKnot.detail());
        if (search == countMap.end()) {
            countMap.insert({surfaceKnot.detail(), 1});
        } else {
            ++search->second;
        }
    }

    graph.reset();

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
        std::cout << "- " << count << (count == 1 ? " copy of " : " copies of
        ")
                  << description << "\n";
    }

    std::cout << "\n";
    std::cout << "Total admissible surfaces = " << surfaceCount
              << "\nTotal closed surfaces = " << closedCount << "\n\n";

    return 0;
}
