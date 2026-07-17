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

struct TriangulationWithEdges {
    regina::Triangulation<3> tri;
    std::vector<const regina::Edge<3> *> edges;
};

// Builds a triangulation of S^3 from `pdcode`, along with the edges of that
// triangulation that trace the knot/link diagram (3 per crossing).
//
// Throws regina::InvalidArgument if `pdcode` is not well-formed (some strand
// label is not used in exactly two crossing slots).
TriangulationWithEdges buildLink(PDCode pdcode);

// Builds a copy of `tri` with every internal edge pinched away except:
//   - the edges in `edges` (mapped into the copy), and
//   - any edge that is a loop (both endpoint vertices coincide) -- pinching
//     a loop edge drills it out instead of collapsing it, leaving an ideal
//     torus/Klein-bottle boundary component, so those are always skipped.
//
// `tri` is left untouched. The returned edge list corresponds 1-1 (same
// order, same length, including repeats) to `edges`, mapped into the
// returned triangulation.
//
// Throws regina::InvalidArgument if any edge in `edges` is a boundary edge
// of `tri` (pinchEdge() cannot be skipped or performed on a boundary edge,
// and buildLink()'s output is always closed, so this should not arise for
// its output).
TriangulationWithEdges
reduceVertices(const regina::Triangulation<3> &tri,
               const std::vector<const regina::Edge<3> *> &edges);
} // namespace knotbuilder

#endif // KNOTBUILDER_H
