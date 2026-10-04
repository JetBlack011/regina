//
//  block.cpp
//

#include "diagramtriangulation/block.h"

#include <map>
#include <optional>
#include <set>
#include <vector>

#include <triangulation/dim3.h>

#include "diagramtriangulation/fromdiagram.h"

namespace knotbuilder {

namespace {
using Wide = __int128;
constexpr BlockCoord S = BLOCK_SCALE;
struct UV {
    BlockCoord x, y;
};
Wide cross2(BlockCoord ax, BlockCoord ay, BlockCoord bx, BlockCoord by) {
    return Wide(ax) * by - Wide(ay) * bx;
}
Wide det(const BlockPoint &a, const BlockPoint &b, const BlockPoint &c,
         const BlockPoint &d) {
    Wide ux = b.x - a.x, uy = b.y - a.y, uz = b.z - a.z;
    Wide vx = c.x - a.x, vy = c.y - a.y, vz = c.z - a.z;
    Wide wx = d.x - a.x, wy = d.y - a.y, wz = d.z - a.z;
    return ux * (vy * wz - vz * wy) - uy * (vx * wz - vz * wx) +
           uz * (vx * wy - vy * wx);
}
} // namespace

const std::array<std::array<BlockPoint, 4>, 14> &blockCoordinates() {
    static const std::array<std::array<BlockPoint, 4>, 14> coords = [] {
        regina::Triangulation<3> b;
        Block block(b);
        auto *core4 = b.tetrahedron(4);
        auto *core5 = b.tetrahedron(5);
        const size_t top = core4->vertex(3)->index();      // over-strand apex
        const size_t under0 = core4->vertex(0)->index();   // under-strand, PD in
        const size_t under1 = core4->vertex(2)->index();   // under-strand, PD out
        const size_t over0 = core4->vertex(1)->index();    // over-strand, walls_[3] side
        const size_t over1 = core5->vertex(1)->index();    // over-strand, walls_[1] side
        if (core5->vertex(3)->index() != top)
            throw regina::InvalidArgument(
                "blockCoordinates(): the over-strand's two edges do not meet");

        std::map<size_t, BlockPoint> at;
        at[under0] = {0, -S, 0};
        at[under1] = {0, S, 0};
        at[over0] = {-S, 0, 0};
        at[over1] = {S, 0, 0};
        // Nudged off the square's centre, so that the over-strand's peak
        // does not project onto the under-strand.
        at[top] = {S / 7, S / 11, S};

        // Boundary adjacency: bottom corners are the rim vertices next to
        // two strand points; top corners are the other neighbours of the
        // top centre, each above the bottom corner it shares an edge with.
        std::map<size_t, std::set<size_t>> adj;
        for (auto *t : b.triangles()) {
            if (!t->isBoundary()) continue;
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    if (i != j)
                        adj[t->vertex(i)->index()].insert(t->vertex(j)->index());
        }
        std::set<size_t> strandPts{under0, under1, over0, over1};
        for (size_t v = 0; v < b.countVertices(); ++v) {
            if (at.contains(v) || adj[v].contains(top)) continue;
            std::vector<size_t> nearStrands;
            for (size_t w : adj[v])
                if (strandPts.contains(w)) nearStrands.push_back(w);
            if (nearStrands.size() != 2) continue;
            BlockCoord x = 0, y = 0;
            for (size_t w : nearStrands) {
                x += at[w].x;
                y += at[w].y;
            }
            at[v] = {x, y, 0};
        }
        for (size_t v : adj[top]) {
            if (at.contains(v)) continue;
            std::optional<BlockPoint> below;
            for (size_t w : adj[v])
                if (at.contains(w) && at[w].z == 0 && std::abs(at[w].x) == S &&
                    std::abs(at[w].y) == S)
                    below = at[w];
            if (!below)
                throw regina::InvalidArgument(
                    "blockCoordinates(): a top corner has no bottom corner below it");
            at[v] = {below->x, below->y, S};
        }
        if (at.size() != b.countVertices())
            throw regina::InvalidArgument(
                "blockCoordinates(): the block's vertices are not where expected");

        std::array<std::array<BlockPoint, 4>, 14> out;
        for (int j = 0; j < 14; ++j)
            for (int i = 0; i < 4; ++i)
                out[j][i] = at.at(b.tetrahedron(j)->vertex(i)->index());
        return out;
    }();
    return coords;
}

std::string verifyBlockModel() {
    const auto &c = blockCoordinates();
    regina::Triangulation<3> b;
    Block block(b);
    if (!b.isOrientable()) return "block not orientable";
    Wide total = 0;
    std::optional<bool> positive;
    for (int j = 0; j < 14; ++j) {
        Wide d = det(c[j][0], c[j][1], c[j][2], c[j][3]);
        if (d == 0) return "degenerate tetrahedron " + std::to_string(j);
        // Signs must agree once each tetrahedron's own orientation is
        // factored in, i.e. the map is orientation-consistent.
        bool pos = (d > 0) == (b.tetrahedron(j)->orientation() > 0);
        if (positive && *positive != pos)
            return "tetrahedra oriented inconsistently";
        positive = pos;
        total += d < 0 ? -d : d;
    }
    // |det| is 6 volume; the box has volume (2S)(2S)(S).
    if (total != Wide(24) * S * S * S) return "volumes do not sum to the box";

    // Boundary: each triangle in a face of the box, faces covered exactly.
    std::map<std::pair<int, BlockCoord>, Wide> area2; // (axis, value) -> 2 x area
    for (int j = 0; j < 14; ++j)
        for (int f = 0; f < 4; ++f) {
            if (b.tetrahedron(j)->adjacentTetrahedron(f)) continue;
            std::vector<BlockPoint> p;
            for (int i = 0; i < 4; ++i)
                if (i != f) p.push_back(c[j][i]);
            bool placed = false;
            for (int axis = 0; axis < 3 && !placed; ++axis) {
                auto get = [axis](const BlockPoint &q) {
                    return axis == 0 ? q.x : axis == 1 ? q.y : q.z;
                };
                BlockCoord v = get(p[0]);
                bool extreme = axis == 2 ? (v == 0 || v == S) : (v == S || v == -S);
                if (!extreme || get(p[1]) != v || get(p[2]) != v) continue;
                auto uv = [axis](const BlockPoint &q) -> UV {
                    return axis == 0 ? UV{q.y, q.z}
                                     : axis == 1 ? UV{q.x, q.z} : UV{q.x, q.y};
                };
                UV a = uv(p[0]), bb = uv(p[1]), cc = uv(p[2]);
                Wide a2 = cross2(bb.x - a.x, bb.y - a.y, cc.x - a.x, cc.y - a.y);
                area2[{axis, v}] += a2 < 0 ? -a2 : a2;
                placed = true;
            }
            if (!placed) return "a boundary triangle is off the box";
        }
    for (const auto &[face, a2] : area2) {
        Wide want = face.first == 2 ? Wide(8) * S * S : Wide(4) * S * S;
        if (a2 != want) return "a face of the box is not covered exactly";
    }
    if (area2.size() != 6) return "the boundary does not cover all six faces";
    return {};
}


} // namespace knotbuilder
