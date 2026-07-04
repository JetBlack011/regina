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

#include <triangulation/dim3.h>

namespace knotbuilder {
using PDCode = std::vector<std::array<int, 4>>;

PDCode parsePDCode(std::string pdcode_str);

class Block {
  private:
    std::vector<regina::Tetrahedron<3> *> core_;
    std::array<std::array<regina::Tetrahedron<3> *, 3>, 4> walls_;

  public:
    Block(regina::Triangulation<3> &tri);

    void glue(size_t myWall, Block &other, size_t otherWall);

    const std::vector<regina::Edge<3> *> getLinkEdges() const;

    friend bool operator==(const Block &lhs, const Block &rhs) {
        return lhs.core_ == rhs.core_ && lhs.walls_ == rhs.walls_;
    }
};

bool buildLink(regina::Triangulation<3> &tri, PDCode pdcode,
               std::vector<const regina::Edge<3> *> &edges);
} // namespace knotbuilder

#endif // KNOTBUILDER_H
