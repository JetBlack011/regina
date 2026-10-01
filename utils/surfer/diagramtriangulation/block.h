//
//  blockgeometry.h
//
//  The linear geometry of knotbuilder's crossing block.
//

/*! \file utils/surfer/knotbuilder/blockgeometry.h
 *  \brief knotbuilder's 14-tetrahedron crossing block, as the box
 *  [-1,1]^2 x [0,1].
 *
 *  knotbuilder::buildLink() builds its triangulation T of S^3 from a
 *  diagram: one Block per crossing (tetrahedra 14k..14k+13, in Block's own
 *  order: cores 0-5, then walls 0-7), glued wall to wall into S^2 x I, then
 *  both boundary spheres coned off (finiteToIdeal()). Each block is
 *  linearly the box [-1,1]^2 x [0,1] (in units of BLOCK_SCALE), with all 13
 *  of its vertices on the box's surface:
 *
 *    - the 8 box corners;
 *    - the 4 strand points, at the midpoints of the bottom edges (one per
 *      wall: where a strand of the diagram leaves the crossing);
 *    - the centre of the top face, where the over-strand peaks (nudged off
 *      centre, so that the peak does not project onto the under-strand).
 *
 *  The under-strand runs straight across the bottom face; the over-strand
 *  arches from one strand point up through the top centre and down to the
 *  opposite one. Each wall is a vertical rectangle fanned from its strand
 *  point, so neighbouring blocks' coordinates agree on the walls they share.
 *
 *  blockCoordinates() gives every (tetrahedron, vertex) of the block its
 *  point; verifyBlockModel() proves the embedding exactly. See
 *  diagramdrawer.h for what this is used for.
 */

#ifndef SURFER_KNOTBUILDER_BLOCKGEOMETRY_H
#define SURFER_KNOTBUILDER_BLOCKGEOMETRY_H

#include <array>
#include <cstdint>
#include <string>

namespace knotbuilder {

/** An integer coordinate. The block is [-BLOCK_SCALE, BLOCK_SCALE]^2 x [0, BLOCK_SCALE]. */
using BlockCoord = std::int64_t;
constexpr BlockCoord BLOCK_SCALE = BlockCoord(1) << 20;

/** A point of one block's box, in that block's own frame. */
struct BlockPoint {
    BlockCoord x, y, z;
    bool operator==(const BlockPoint &) const = default;
};

/**
 * The box coordinates of knotbuilder's crossing block: entry [j][i] is the
 * point of vertex i of the block's j-th tetrahedron (Block's order: cores
 * 0-5, then walls 0-7). Derived from a freshly built Block by the vertices'
 * combinatorial roles (strand points from getLinkEdges(), top corners as
 * the top centre's neighbours, bottom corners as the rim vertices between
 * two strand points), never from Regina's vertex numbering.
 */
const std::array<std::array<BlockPoint, 4>, 14> &blockCoordinates();

/**
 * Checks blockCoordinates() exactly against a freshly built Block: every
 * tetrahedron nondegenerate and oriented consistently with the block, their
 * volumes summing to the box's, and every boundary triangle lying in a face
 * of the box with each face covered exactly. With the boundary carried onto
 * the box's boundary with degree one, positive volumes summing to the box's
 * mean every point of the box lies in exactly one tetrahedron: the block
 * tiles the box. Returns an empty string if so, else why not.
 */
std::string verifyBlockModel();

} // namespace knotbuilder

#endif
