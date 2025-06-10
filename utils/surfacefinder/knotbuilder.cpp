//
//  knotbuilder.cpp
//
//  Created by John Teague on 05/10/2025.
//
//  Algorithm adapted from work of Srinivas Vadhiraj, Samantha Ward, Angela
//  Yuan, and Jingyuan Zhang in the Texas Experimental Geometry Lab.

#include <gmpxx.h>
#include <link/link.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <unistd.h>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <string>

#include "triangulation/forward.h"
#include "triangulation/generic/triangulation.h"
#include "utilities/exception.h"

#include "knotbuilder.h"

knotbuilder::PDCode knotbuilder::parsePDCode(std::string pdcode_str) {
    std::vector<std::array<int, 4>> pdcode;

    for (char &c : pdcode_str) {
        if (!std::isdigit(c)) {
            c = ' ';
        }
    }

    std::stringstream ss(pdcode_str);
    std::vector<int> pdlist;
    int token;
    while (ss >> token) {
        pdlist.push_back(token);
    }

    bool isZeroIndexed = false;
    if (std::ranges::find(pdlist, 0) != pdlist.end()) {
        isZeroIndexed = true;
    }

    if (!isZeroIndexed) {
        for (int &i : pdlist) {
            --i;
        }
    }

    for (int i = 0; i < pdlist.size(); i += 4) {
        std::array<int, 4> crossing;
        for (int j = 0; j < 4; j++) {
            crossing[j] = pdlist[i + j];
        }
        pdcode.push_back(crossing);
    }

    return pdcode;
}

knotbuilder::Block::Block(regina::Triangulation<3> &tri) {
    const regina::Perm<4> id;
    core_.reserve(6);
    for (int i = 0; i < 6; ++i) {
        core_.push_back(tri.newTetrahedron());
    }

    core_[0]->join(0, core_[4], id);
    core_[1]->join(3, core_[5], {1, 2, 3, 0});
    core_[2]->join(0, core_[5], {0, 1});
    core_[3]->join(3, core_[4], {0, 2, 3, 1});
    core_[4]->join(3, core_[5], id);

    std::vector<regina::Tetrahedron<3> *> wallTets;
    wallTets.reserve(8);
    for (int i = 0; i < 8; ++i) {
        wallTets.push_back(tri.newTetrahedron());
    }

    wallTets[0]->join(1, core_[0], {1, 2});
    wallTets[1]->join(3, core_[0], id);
    wallTets[1]->join(0, core_[1], {2, 0, 1, 3});
    wallTets[2]->join(0, core_[1], {0, 1});
    wallTets[3]->join(0, core_[1], id);
    wallTets[3]->join(3, core_[2], {0, 2, 3, 1});
    wallTets[4]->join(1, core_[2], {1, 2});
    wallTets[5]->join(3, core_[2], id);
    wallTets[5]->join(0, core_[3], {2, 0, 1, 3});
    wallTets[6]->join(0, core_[3], {0, 1});
    wallTets[7]->join(0, core_[3], id);
    wallTets[7]->join(3, core_[0], {0, 2, 3, 1});

    for (size_t i = 0; i < 4; ++i) {
        walls_[i] = {wallTets[2 * i], wallTets[2 * i + 1],
                     wallTets[(2 * i + 2) % 8]};
    }

    //std::cout << "New block = " << tri.isoSig() << "\n";
}

void knotbuilder::Block::glue(size_t myWall, Block &other, size_t otherWall) {
    if (*this == other &&
        std::max(myWall, otherWall) - std::min(myWall, otherWall) == 2)
        throw regina::InvalidArgument("Invalid block gluing! A block cannot be "
                                      "glued to itself along opposite faces.");

    int f0, f1, f2;
    regina::Perm<4> g0, g1, g2;

    if (myWall % 2 == 0 && otherWall % 2 == 0) {
        f0 = 3;
        f1 = 2;
        f2 = 2;
    } else if (myWall % 2 == 0 && otherWall % 2 == 1) {
        f0 = 3;
        g0 = {2, 3};
        f1 = 2;
        g1 = {1, 2};
        f2 = 2;
        g2 = {1, 2};
    } else if (myWall % 2 == 1 && otherWall % 2 == 0) {
        f0 = 1;
        g0 = {1, 2};
        f1 = 1;
        g1 = {1, 2};
        f2 = 2;
        g2 = {2, 3};
    } else if (myWall % 2 == 1 && otherWall % 2 == 1) {
        f0 = 1;
        f1 = 1;
        f2 = 2;
    }

    if (myWall % 2 == otherWall % 2) {
        walls_[myWall][0]->join(f0, other.walls_[otherWall][0], g0);
        walls_[myWall][1]->join(f1, other.walls_[otherWall][1], g1);
        walls_[myWall][2]->join(f2, other.walls_[otherWall][2], g2);
    } else {
        walls_[myWall][0]->join(f0, other.walls_[otherWall][2], g0);
        walls_[myWall][1]->join(f1, other.walls_[otherWall][1], g1);
        walls_[myWall][2]->join(f2, other.walls_[otherWall][0], g2);
    }
}

template <int dim>
bool isOrdered(const regina::Triangulation<dim> &tri) {
    for (const auto &s : tri.simplices()) {
        for (int f = 0; f < dim; ++f) {
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

bool knotbuilder::buildLink(regina::Triangulation<3> &tri,
                            std::string &pdcode_str,
                            std::vector<const regina::Edge<3> *> &edges) {
    PDCode pdcode = parsePDCode(pdcode_str);
    size_t numCrossings = pdcode.size();

    std::cout << "Building link with " << numCrossings << " crossings.\n";

    std::vector<Block> blocks;

    blocks.reserve(numCrossings);
    for (size_t i = 0; i < numCrossings; ++i) {
        blocks.emplace_back(tri);
    }

    // std::cout << "Building strands.\n";
    std::vector<std::vector<std::pair<int, int>>> strands(2 * numCrossings);

    for (int n = 0; n < numCrossings; ++n) {
        strands[pdcode[n][0]].emplace_back(n, 0);
        strands[pdcode[n][1]].emplace_back(n, 1);
        strands[pdcode[n][2]].emplace_back(n, 2);
        strands[pdcode[n][3]].emplace_back(n, 3);
    }

    // std::cout << "Strands:\n";
    // for (int i = 0; i < strands.size(); ++i) {
    //     std::cout << "Strand " << i + 1 << ": ";
    //     for (const auto &pair : strands[i]) {
    //         std::cout << "(" << pair.first << ", " << pair.second << ") ";
    //     }
    //     std::cout << "\n";
    // }

    // std::cout << tri.isoSig() << "\n";
    // std::cout << "Gluing blocks.\n";
    for (const auto &strand : strands) {
        if (strand.size() != 2)
            return false;

        Block &block1 = blocks[strand[0].first];
        Block &block2 = blocks[strand[1].first];
        int wall1 = strand[0].second;
        int wall2 = strand[1].second;

        // std::cout << "Gluing blocks " << strand[0].first << " and "
        //           << strand[1].first << " on walls " << wall1 << " and "
        //           << wall2 << "\n";

        block1.glue(wall1, block2, wall2);
    }

    tri.finiteToIdeal();

    for (const auto &block : blocks) {
        const auto blockEdges = block.getEdges();
        edges.insert(edges.begin(), blockEdges.begin(), blockEdges.end());
    }

    return true;
}

const std::vector<regina::Edge<3> *> knotbuilder::Block::getEdges() const {
    return {core_[4]->edge(5), core_[5]->edge(5), core_[4]->edge(0)};
}

namespace {
void usage(const char *progName, const std::string &error = std::string()) {
    if (!error.empty())
        std::cerr << error << "\n\n";

    std::cerr << "Usage:\n";
    std::cerr << "    " << progName << " <PD Code>\n    " << progName
              << " [ -v, --version | -h, --help ]\n\n";
    std::cerr << "    -v, --version  : Show which version of Regina is "
                 "being used\n";
    std::cerr << "    -h, --help     : Display this help\n";
    exit(1);
}
} // namespace

// int main(int argc, char *argv[]) {
//     // Check for standard arguments:
//     for (int i = 1; i < argc; ++i) {
//         std::string arg = argv[i];
//         if (arg == "-h" || arg == "--help")
//             usage(argv[0]);
//         if (arg == "-v" || arg == "--version") {
//             if (argc != 2)
//                 usage(argv[0], "Option --version cannot be used with "
//                                "any other arguments.");
//             std::cout << PACKAGE_BUILD_STRING << "\n";
//             exit(0);
//         }
//     }
//
//     if (argc != 2) {
//         usage(argv[0], "Please provide a valid PD Code.");
//     }
//
//     regina::Triangulation<3> tri;
//     std::vector<const regina::Edge<3> *> edges;
//
//     std::string pdcode_str = argv[1];
//     if (pdcode_str.length() <= 4 || !knotbuilder::buildLink(tri, pdcode_str,
//     edges)) {
//         usage(argv[0], "Please provide a valid PD Code.");
//     }
//
//     std::cout << "\nTriangulation of S^3 = " << tri.isoSig() << "\n\n";
//
//     Link l(tri, edges);
//     std::cout << "Triangulation of the complement = "
//               << l.buildComplement().isoSig() << "\n";
//
//     return 0;
// }
