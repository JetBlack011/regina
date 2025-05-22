//
//  knotbuilder.h
//
//  Created by John Teague on 05/10/2024.
//
//  This is adapted from work of Srinivas Vadhiraj, Samantha Ward, Angela
//  Yuan, and Jingyuan Zhang performed for the Texas Experimental Geometry Lab
//  at UT Austin.

#ifndef KNOTBUILDER_H

#define KNOTBUILDER_H

#include <gmpxx.h>
#include <link/link.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <unistd.h>

#include <cassert>
#include <string>

#include "triangulation/forward.h"

using PDCode = std::vector<std::array<int, 4>>;

PDCode parsePDCode(std::string pdcode_str);

class Block {
  private:
    std::vector<regina::Tetrahedron<3> *> core_;
    std::array<std::array<regina::Tetrahedron<3> *, 3>, 4> walls_;

  public:
    Block(regina::Triangulation<3> &tri);

    void glue(size_t myWall, Block &other, size_t otherWall);

    const std::vector<regina::Edge<3> *> getEdges() const;

    friend bool operator==(const Block &lhs, const Block &rhs) {
        return lhs.core_ == rhs.core_ && lhs.walls_ == rhs.walls_;
    }
};

bool buildLink(regina::Triangulation<3> &tri, std::string &pdcode_str,
               std::vector<const regina::Edge<3> *> &edges);

#endif
