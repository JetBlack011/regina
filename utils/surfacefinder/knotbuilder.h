//
//  knotbuilder.h 
//
//  Created by John Teague on 05/10/2024.
//
//  This is adapted from work of Srinivas Vadhiraj, Samantha Ward, Angela
//  Yuan, and Jingyuan Zhang performed for the Texas Experimental Geometry Lab
//  at UT Austin.

#include <algorithm>
#include <gmpxx.h>
#include <link/link.h>
#include <sstream>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <unistd.h>

#include <cassert>
#include <iostream>
#include <ostream>
#include <string>

#include "knottedsurfaces.h"
#include "triangulation/example2.h"
#include "triangulation/example3.h"
#include "triangulation/forward.h"
#include "triangulation/generic/triangulation.h"
#include "utilities/exception.h"

#include "gluing.h"

std::vector<std::vector<int>> parsePDCode(std::string pdcode_str) {
    std::vector<std::vector<int>> pdcode;

    // remove the brackets
    pdcode_str = pdcode_str.substr(2, pdcode_str.length() - 4);

    // Replace delineating characters with spaces
    std::string target = "),(";
    std::string replacement = " ";

    size_t pos = 0;
    while ((pos = pdcode_str.find(target, pos)) != std::string::npos) {
        pdcode_str.replace(pos, target.length(), replacement);
        pos += replacement.length(); // move past the replacement
    }

    // Split into crossings
    std::istringstream pdcode_iss(pdcode_str);
    std::string token;
    std::vector<std::istringstream> crossings;

    while (pdcode_iss >> token) {
        crossings.emplace_back(token);
    }

    // Convert to final pdcode
    for (std::istringstream &crossing_iss : crossings) {
        std::vector<int> crossing;
        std::string token;
        while (std::getline(crossing_iss, token, ',')) {
            crossing.push_back(std::stoi(token));
        }
        pdcode.push_back(crossing);
    }

    return pdcode;
}

class Block {
  private:
    std::vector<regina::Tetrahedron<3> *> core_;
    std::array<std::array<regina::Tetrahedron<3> *, 3>, 4> walls_;

  public:
    Block(regina::Triangulation<3> &tri) {
        const regina::Perm<4> id;
        core_.reserve(6);
        for (int i = 0; i < 6; ++i) {
            core_.push_back(tri.newTetrahedron());
        }

        core_[0]->join(0, core_[4], id);
        core_[1]->join(0, core_[5], id);
        core_[2]->join(0, core_[5], {0, 1});
        core_[3]->join(0, core_[4], {0, 1});
        core_[4]->join(2, core_[5], id);

        std::vector<regina::Tetrahedron<3> *> wallTets;
        wallTets.reserve(8);
        for (int i = 0; i < 8; ++i) {
            wallTets.push_back(tri.newTetrahedron());
        }

        wallTets[0]->join(3, core_[0], id);
        wallTets[1]->join(2, core_[0], id);
        wallTets[1]->join(0, core_[1], {0, 2});
        wallTets[2]->join(3, core_[1], id);
        wallTets[3]->join(1, core_[1], id);
        wallTets[3]->join(0, core_[2], {0, 1});
        wallTets[4]->join(3, core_[2], id);
        wallTets[5]->join(2, core_[2], id);
        wallTets[5]->join(0, core_[3], {0, 2});
        wallTets[6]->join(3, core_[3], id);
        wallTets[7]->join(1, core_[3], id);
        wallTets[7]->join(0, core_[0], {0, 1});

        for (size_t i = 0; i < 4; ++i) {
            walls_[i] = {wallTets[2 * i], wallTets[2 * i + 1],
                         wallTets[(2 * i + 2) % 8]};
        }
    }

    void glue(size_t myWall, Block &other, size_t otherWall) {
        int minWall, maxWall;
        Block *minBlock, *maxBlock;
        if (myWall <= otherWall) {
            minWall = myWall;
            minBlock = this;
            maxWall = otherWall;
            maxBlock = &other;
        } else {
            minWall = otherWall;
            minBlock = &other;
            maxWall = myWall;
            maxBlock = this;
        }

        int minFace = (minWall % 2 == 0) ? 2 : 1;
        regina::Perm<4> g0;
        regina::Perm<4> g1;
        regina::Perm<4> g2;

        if ((maxWall - minWall) % 2 == 0) {
            g1 = {0, minFace};
        } else if (minFace == 2) {
            g0 = {1, 2};
            g1 = {1, 2, 0, 3};
            g2 = {1, 2};
        } else if (minFace == 1) {
            g0 = {1, 2};
            g1 = {2, 0, 1, 3};
            g2 = {1, 2};
        }

        minBlock->walls_[minWall][0]->join(minFace,
                                           maxBlock->walls_[maxWall][2], g0);
        minBlock->walls_[minWall][1]->join(3, maxBlock->walls_[maxWall][1], g1);
        minBlock->walls_[minWall][2]->join(minFace,
                                           maxBlock->walls_[maxWall][0], g2);
    }

    const std::vector<regina::Edge<3> *> getEdges() const {
        return {core_[4]->edge(5), core_[5]->edge(5), core_[4]->edge(0)};
    }
};

bool buildKnot(regina::Triangulation<3> &tri, std::string &pdcode_str,
               std::vector<regina::Edge<3> *> &knotEdges) {
    std::vector<std::vector<int>> pdcode = parsePDCode(pdcode_str);
    size_t numCrossings = pdcode.size();

    std::cout << "Building knot with " << numCrossings << " crossings.\n";
    std::vector<Block> blocks;

    blocks.reserve(numCrossings);
    for (size_t i = 0; i < numCrossings; ++i) {
        blocks.emplace_back(tri);
    }

    // std::cout << "Building strands.\n";
    std::vector<std::vector<std::pair<int, int>>> strands(2 * numCrossings);

    for (int n = 0; n < numCrossings; ++n) {
        strands[pdcode[n][0] - 1].emplace_back(n, 0);
        strands[pdcode[n][1] - 1].emplace_back(n, 1);
        strands[pdcode[n][2] - 1].emplace_back(n, 2);
        strands[pdcode[n][3] - 1].emplace_back(n, 3);
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
        knotEdges.insert(knotEdges.end(), blockEdges.begin(), blockEdges.end());
    }

    return true;
}
