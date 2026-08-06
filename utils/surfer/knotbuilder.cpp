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

    // PD codes conventionally 1-index strand labels; detect 0-indexed
    // input (a literal 0 can only appear in that case) and normalize to
    // 0-indexed internally either way.
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
    // Fixed 14-tetrahedron gluing pattern for one crossing: 6 "core"
    // tetrahedra encoding the over/under strand crossing itself, plus 8
    // "wall" tetrahedra (2 per side) forming the 4 walls glue() joins to
    // neighboring Blocks.
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

// Builds the permutation identifying `mine`'s local vertices with
// `theirs`', role-for-role.
regina::Perm<4> matchCorner(const CornerRole &mine, const CornerRole &theirs) {
    std::array<int, 4> image;
    image[mine.excl] = theirs.excl;
    image[mine.own] = theirs.own;
    image[mine.anchor] = theirs.anchor;
    image[mine.priv] = theirs.priv;
    return regina::Perm<4>(image);
}

// As matchCorner(), for a wall's middle tetrahedron.
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

namespace {
using Strands = std::vector<std::vector<std::pair<int, int>>>;
using PDPos = std::pair<int, int>;

// Mirrors regina::Link::fromPD's own dir[]/occ[] walk
// (engine/link/pd-impl.h), adapted to run directly against knotbuilder's
// own already-computed `strands` occurrences rather than building a
// regina::Link. strands[label] plays the role of pd-impl.h's occ[label]
// (occ[label].first/.second are strands[label][0]/[1]).
//
// dir[label] is 1 if the strand flows strands[label][0] -> strands[label][1]
// (i.e. its first occurrence is where it *leaves* that strand, its second
// where it *arrives*), or -1 for the reverse. Always populated (every
// strand label appears exactly twice, checked by buildLink() before this
// runs), since every component either passes under some crossing (handled
// by the first loop below, seeded from that crossing's own predetermined
// arrival at position 0) or -- for the rare case of a component that is
// always the overcrossing, so the PD code doesn't define its own
// orientation -- gets an arbitrary but still self-consistent direction
// from the second loop, exactly as pd-impl.h does for the same case.
std::vector<int> computeStrandDirections(const knotbuilder::PDCode &pdcode,
                                         const Strands &strands) {
    std::vector<int> dir(strands.size(), 0);

    auto walk = [&](int start) {
        int s = start;
        PDPos pos;
        while (true) {
            pos = (dir[s] > 0) ? strands[s][1] : strands[s][0];
            pos.second ^= 2; // same crossing's other half of this pass
                              // (understrand: 0<->2, overstrand: 1<->3)
            s = pdcode[pos.first][pos.second];
            if (s == start)
                break;
            dir[s] = (strands[s][0] == pos) ? 1 : -1;
        }
    };

    // Position 0 is always the incoming understrand -- the only locally
    // predetermined direction (see the comment on buildLink()'s own use of
    // this fact). Seed every component reachable this way.
    for (int n = 0; n < static_cast<int>(pdcode.size()); ++n) {
        int start = pdcode[n][0];
        if (dir[start])
            continue;
        PDPos here{n, 0};
        dir[start] = (strands[start][0] == here) ? -1 : 1;
        walk(start);
    }

    // Any label still unresolved belongs to a component that never
    // provides a position-0 anchor (every crossing it's involved in, it's
    // always the overstrand) -- arbitrary but consistent, as pd-impl.h
    // does for the same case.
    for (int n = 0; n < static_cast<int>(pdcode.size()); ++n) {
        int start = pdcode[n][1];
        if (dir[start])
            continue;
        dir[start] = 1;
        walk(start);
    }

    return dir;
}

// Whether position `pos` of crossing `n` is the PD-intended arrival point
// of whichever strand label occupies it -- i.e. whether that strand's
// directed span (per `dir`, computed above) ends rather than begins here.
bool isEntry(const knotbuilder::PDCode &pdcode, const Strands &strands,
            const std::vector<int> &dir, int n, int pos) {
    int label = pdcode[n][pos];
    PDPos here{n, pos};
    int occIdx = (strands[label][0] == here) ? 0 : 1;
    int arrivalIdx = (dir[label] > 0) ? 1 : 0;
    return occIdx == arrivalIdx;
}
} // namespace

knotbuilder::TriangulationWithLink knotbuilder::buildLink(PDCode pdcode) {
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

    // PD-intended direction per crossing's overstrand pass: not locally
    // fixed like the understrand (whose position-0-is-arrival rule holds
    // unconditionally for every crossing -- see computeStrandDirections()'s
    // own comment), so this walk is required to resolve it per crossing.
    std::vector<int> dir = computeStrandDirections(pdcode, strands);

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
    std::vector<bool> reversed;
    for (size_t n = 0; n < blocks.size(); ++n) {
        const auto blockEdges = blocks[n].getLinkEdges();
        bool walls1IsEntry =
            isEntry(pdcode, strands, dir, static_cast<int>(n), 1);
        const auto blockReversed =
            blocks[n].getLinkEdgeDirections(walls1IsEntry);
        edges.insert(edges.begin(), blockEdges.begin(), blockEdges.end());
        reversed.insert(reversed.begin(), blockReversed.begin(),
                        blockReversed.end());
    }

    return {std::move(tri), std::move(edges), std::move(reversed)};
}

const std::vector<regina::Edge<3> *> knotbuilder::Block::getLinkEdges() const {
    return {core_[4]->edge(1, 3), core_[4]->edge(0, 2), core_[5]->edge(1, 3)};
}

std::vector<bool>
knotbuilder::Block::getLinkEdgeDirections(bool walls1IsEntry) const {
    // Order matches getLinkEdges(): {core_[4]->edge(1,3), core_[4]->edge(0,2),
    // core_[5]->edge(1,3)}.
    //
    // core_[4]->edge(0,2) (the understrand, == core_[5]->edge(0,2)):
    // unconditionally not reversed -- see this method's own doc comment.
    //
    // Overstrand: the PD-intended traversal is walls_[1] -> apex ->
    // walls_[3] when walls1IsEntry, or walls_[3] -> apex -> walls_[1]
    // otherwise. core_[5]->edge(1,3) runs walls_[1]-side (vertex(0)) ->
    // apex (vertex(1)); core_[4]->edge(1,3) runs walls_[3]-side
    // (vertex(0)) -> apex (vertex(1)). So when walls1IsEntry,
    // core_[5]->edge(1,3) already runs vertex(0)->vertex(1) (not
    // reversed) while core_[4]->edge(1,3) runs against its own
    // vertex(0)->vertex(1) (reversed) -- and vice versa otherwise.
    return {/* core_[4]->edge(1,3) */ walls1IsEntry,
            /* core_[4]->edge(0,2) */ false,
            /* core_[5]->edge(1,3) */ !walls1IsEntry};
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

knotbuilder::TriangulationWithLink
knotbuilder::reduceVertices(const regina::Triangulation<3> &tri,
                            const std::vector<const regina::Edge<3> *> &edges) {
    regina::Triangulation<3> newTri(tri);

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
        descByOrig[e] = {.tet = emb.tetrahedron(), .localEdge = emb.edge()};

        for (int i = 0; i < 2; ++i) {
            regina::Vertex<3> *v = mapped->vertex(i);
            if (!seenVertices.insert(v).second)
                continue;
            auto vemb = v->front();
            protectedVertices.push_back({vemb.tetrahedron(), vemb.vertex()});
        }
    }

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

    return {.tri = std::move(newTri), .edges = std::move(newEdges)};
}
