//
//  knotbuilder.cpp
//
//  Created by John Teague on 05/10/2025.
//
//  Algorithm adapted from work of Srinivas Vadhiraj, Samantha Ward, Angela
//  Yuan, and Jingyuan Zhang in the Texas Experimental Geometry Lab.

#include "knotbuilder.h"

#include <unordered_map>
#include <unordered_set>

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

    core_[0]->join(2, core_[4], {0, 1, 2, 3});
    core_[1]->join(2, core_[5], {0, 1, 2, 3});
    core_[2]->join(2, core_[5], {1, 2, 0, 3});
    core_[3]->join(2, core_[4], {1, 2, 0, 3});
    core_[4]->join(1, core_[5], {0, 1, 2, 3});

    std::vector<regina::Tetrahedron<3> *> wallTets;
    wallTets.reserve(8);
    for (int i = 0; i < 8; ++i) {
        wallTets.push_back(tri.newTetrahedron());
    }

    wallTets[0]->join(0, core_[0], {3, 0, 1, 2});
    wallTets[1]->join(2, core_[0], {0, 2, 1, 3});
    wallTets[1]->join(1, core_[1], {0, 1, 2, 3});
    wallTets[2]->join(0, core_[1], {3, 0, 1, 2});
    wallTets[3]->join(1, core_[1], {1, 0, 2, 3});
    wallTets[3]->join(2, core_[2], {0, 2, 1, 3});
    wallTets[4]->join(0, core_[2], {3, 0, 1, 2});
    wallTets[5]->join(2, core_[2], {1, 2, 0, 3});
    wallTets[5]->join(1, core_[3], {1, 0, 2, 3});
    wallTets[6]->join(0, core_[3], {3, 0, 1, 2});
    wallTets[7]->join(1, core_[3], {0, 1, 2, 3});
    wallTets[7]->join(2, core_[0], {1, 2, 0, 3});

    for (size_t i = 0; i < 4; ++i) {
        walls_[i] = {wallTets[2 * i], wallTets[2 * i + 1],
                     wallTets[(2 * i + 2) % 8]};
    }

    // std::cout << "New block = " << tri.isoSig() << "\n";
}

namespace {
// Local-vertex roles for a wall's corner tetrahedron in its "before" role
// (walls_[i][0], touching the core "before" this wall) or "after" role
// (walls_[i][2], touching the core "after" this wall): which local vertex
// is excluded (the facet to join on), which is the corner tet's own private
// apex, which is this wall's shared anchor vertex (the point common to all
// three of the wall's tetrahedra), and which is the vertex shared with the
// adjacent core.
struct CornerRole {
    int excl, own, anchor, priv;
};

// Local-vertex roles for a wall's middle tetrahedron (walls_[i][1]): the
// excluded vertex, the wall's anchor vertex, and the two vertices shared
// with the "before" and "after" cores respectively.
struct MiddleRole {
    int excl, anchor, privBefore, privAfter;
};

constexpr CornerRole beforeRole[4] = {
    {2, 0, 1, 3},
    {1, 0, 2, 3},
    {1, 0, 2, 3},
    {2, 0, 1, 3},
};

constexpr CornerRole afterRole[4] = {
    {2, 0, 1, 3},
    {2, 0, 1, 3},
    {1, 0, 2, 3},
    {1, 0, 2, 3},
};

constexpr MiddleRole midRole[4] = {
    {3, 0, 1, 2},
    {3, 0, 2, 1},
    {3, 0, 1, 2},
    {3, 0, 2, 1},
};

regina::Perm<4> matchCorner(const CornerRole &mine, const CornerRole &theirs) {
    std::array<int, 4> image;
    image[mine.excl] = theirs.excl;
    image[mine.own] = theirs.own;
    image[mine.anchor] = theirs.anchor;
    image[mine.priv] = theirs.priv;
    return regina::Perm<4>(image);
}

regina::Perm<4> matchMiddle(const MiddleRole &mine, const MiddleRole &theirs) {
    std::array<int, 4> image;
    image[mine.excl] = theirs.excl;
    image[mine.anchor] = theirs.anchor;
    image[mine.privBefore] = theirs.privAfter;
    image[mine.privAfter] = theirs.privBefore;
    return regina::Perm<4>(image);
}
} // namespace

void knotbuilder::Block::glue(size_t myWall, Block &other, size_t otherWall) {
    if (*this == other &&
        std::max(myWall, otherWall) - std::min(myWall, otherWall) == 2)
        throw regina::InvalidArgument("Invalid block gluing! A block cannot be "
                                      "glued to itself along opposite faces.");

    // This wall's "before" corner meets the other wall's "after" corner
    // (they are the two ends of the same diagram edge, and the region
    // trailing this wall matches the region leading into the other), and
    // symmetrically for "after" meeting "before".
    const CornerRole &myBefore = beforeRole[myWall];
    const CornerRole &theirAfter = afterRole[otherWall];
    walls_[myWall][0]->join(myBefore.excl, other.walls_[otherWall][2],
                            matchCorner(myBefore, theirAfter));

    const CornerRole &myAfter = afterRole[myWall];
    const CornerRole &theirBefore = beforeRole[otherWall];
    walls_[myWall][2]->join(myAfter.excl, other.walls_[otherWall][0],
                            matchCorner(myAfter, theirBefore));

    const MiddleRole &myMid = midRole[myWall];
    const MiddleRole &theirMid = midRole[otherWall];
    walls_[myWall][1]->join(myMid.excl, other.walls_[otherWall][1],
                            matchMiddle(myMid, theirMid));
}

knotbuilder::TriangulationWithEdges knotbuilder::buildLink(PDCode pdcode) {
    size_t numCrossings = pdcode.size();

    std::vector<std::vector<std::pair<int, int>>> strands(2 * numCrossings);

    for (int n = 0; n < numCrossings; ++n) {
        strands[pdcode[n][0]].emplace_back(n, 0);
        strands[pdcode[n][1]].emplace_back(n, 1);
        strands[pdcode[n][2]].emplace_back(n, 2);
        strands[pdcode[n][3]].emplace_back(n, 3);
    }

    for (const auto &strand : strands) {
        if (strand.size() != 2)
            throw regina::InvalidArgument(
                "buildLink(): every PD code strand label must appear in "
                "exactly two crossing slots");
    }

    regina::Triangulation<3> tri;

    std::vector<Block> blocks;
    blocks.reserve(numCrossings);
    for (size_t i = 0; i < numCrossings; ++i) {
        blocks.emplace_back(tri);
    }

    for (const auto &strand : strands) {
        Block &block1 = blocks[strand[0].first];
        Block &block2 = blocks[strand[1].first];
        int wall1 = strand[0].second;
        int wall2 = strand[1].second;

        block1.glue(wall1, block2, wall2);
    }

    tri.finiteToIdeal();

    std::vector<const regina::Edge<3> *> edges;
    for (const auto &block : blocks) {
        const auto blockEdges = block.getLinkEdges();
        edges.insert(edges.begin(), blockEdges.begin(), blockEdges.end());
    }

    return {std::move(tri), std::move(edges)};
}

const std::vector<regina::Edge<3> *> knotbuilder::Block::getLinkEdges() const {
    return {core_[4]->edge(1, 3), core_[4]->edge(0, 2), core_[5]->edge(1, 3)};
}

namespace {
// A tetrahedron + local edge/vertex number stays valid across repeated
// pinchEdge() calls (which only append new tetrahedra), unlike Edge<3>*/
// Vertex<3>* themselves (pinchEdge() fully rebuilds the skeleton, so every
// Edge<3>*/Vertex<3>* is invalidated on every call -- see
// engine/triangulation/dim3/moves.cpp).
struct EdgeDescriptor {
    regina::Tetrahedron<3> *tet;
    int localEdge;

    regina::Edge<3> *resolve() const { return tet->edge(localEdge); }
};

struct VertexDescriptor {
    regina::Tetrahedron<3> *tet;
    int localVertex;

    regina::Vertex<3> *resolve() const { return tet->vertex(localVertex); }
};
} // namespace

knotbuilder::TriangulationWithEdges knotbuilder::reduceVertices(
    const regina::Triangulation<3> &tri,
    const std::vector<const regina::Edge<3> *> &edges) {
    regina::Triangulation<3> newTri(tri);

    // Preserved edges themselves are never pinch candidates, but that alone
    // isn't enough: pinchEdge() on some other, unrelated edge merges its two
    // endpoint vertices, and if either of those happens to also be an
    // endpoint of a preserved edge, the preserved edge's identity (and
    // possibly its very shape -- it could collapse into a loop) changes as
    // a side effect. So every vertex touched by a preserved edge must also
    // be off-limits as a pinch-candidate endpoint.
    std::unordered_map<const regina::Edge<3> *, EdgeDescriptor> descByOrig;
    std::vector<VertexDescriptor> protectedVertices;
    std::unordered_set<regina::Vertex<3> *> seenVertices;

    for (const regina::Edge<3> *e : edges) {
        if (descByOrig.contains(e))
            continue;

        regina::Edge<3> *mapped = newTri.edge(e->index());
        if (mapped->isBoundary())
            throw regina::InvalidArgument(
                "reduceVertices(): a preserved edge is a boundary edge");

        auto emb = mapped->front();
        descByOrig[e] = {emb.tetrahedron(), emb.edge()};

        for (int i = 0; i < 2; ++i) {
            regina::Vertex<3> *v = mapped->vertex(i);
            if (!seenVertices.insert(v).second)
                continue;
            auto vemb = v->front();
            protectedVertices.push_back({vemb.tetrahedron(), vemb.vertex()});
        }
    }

    // Repeatedly pinch one eligible edge at a time, rescanning from scratch
    // every iteration since pinchEdge() invalidates every Edge<3>*/
    // Vertex<3>* in newTri (see the descriptor structs above for how the
    // preserved/protected sets survive this).
    while (true) {
        std::unordered_set<const regina::Edge<3> *> preservedNow;
        for (const auto &[orig, desc] : descByOrig)
            preservedNow.insert(desc.resolve());

        std::unordered_set<const regina::Vertex<3> *> protectedNow;
        for (const auto &desc : protectedVertices)
            protectedNow.insert(desc.resolve());

        regina::Edge<3> *candidate = nullptr;
        for (regina::Edge<3> *e : newTri.edges()) {
            if (preservedNow.contains(e))
                continue;
            if (e->vertex(0) == e->vertex(1))
                continue; // pinching a loop would leave ideal boundary
            if (protectedNow.contains(e->vertex(0)) ||
                protectedNow.contains(e->vertex(1)))
                continue; // would disturb a preserved edge's endpoint
            candidate = e;
            break;
        }

        if (!candidate)
            break;

        newTri.pinchEdge(candidate);
    }

    std::vector<const regina::Edge<3> *> newEdges;
    newEdges.reserve(edges.size());
    for (const regina::Edge<3> *e : edges)
        newEdges.push_back(descByOrig.at(e).resolve());

    return {std::move(newTri), std::move(newEdges)};
}
