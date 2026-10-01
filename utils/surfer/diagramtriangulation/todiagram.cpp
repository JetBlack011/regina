//
//  diagramdrawer.cpp
//
//  See diagramdrawer.h for the construction and why it is an isotopy.
//

#include "diagramtriangulation/todiagram.h"

#include <algorithm>
#include <map>
#include <numeric>
#include <optional>
#include <queue>
#include <set>
#include <unordered_map>

namespace knotbuilder {

namespace {

using Wide = __int128;

constexpr BlockCoord S = BLOCK_SCALE;
// Route heights: far above (below) every block, rising (falling) along the
// route so that it stays an unknotted arc.
constexpr BlockCoord ROUTE_HEIGHT = BlockCoord(1) << 40;
constexpr BlockCoord ROUTE_STEP = BlockCoord(1) << 8;

struct XY {
    BlockCoord x, y;
    bool operator==(const XY &) const = default;
    bool operator<(const XY &o) const { return x != o.x ? x < o.x : y < o.y; }
};

Wide cross(BlockCoord ax, BlockCoord ay, BlockCoord bx, BlockCoord by) {
    return Wide(ax) * by - Wide(ay) * bx;
}

// The side a wall lies on, in a block's own frame: side 0 is y = -S, 1 is
// x = +S, 2 is y = +S, 3 is x = -S. At a corner, the quarter of directions
// into the square runs counterclockwise from its "start" side to its "end"
// side; crossing the end side leads to the next block counterclockwise
// around that corner.
struct CornerSides {
    int start, end;
};
CornerSides cornerSides(XY c) {
    if (c.x < 0 && c.y < 0) return {0, 3};
    if (c.x > 0 && c.y < 0) return {1, 0};
    if (c.x > 0 && c.y > 0) return {2, 1};
    return {3, 2};
}
bool onSide(const BlockPoint &p, int side) {
    switch (side) {
        case 0: return p.y == -S;
        case 1: return p.x == S;
        case 2: return p.y == S;
        default: return p.x == -S;
    }
}
// Unit direction (in units of S) along `side`, pointing away from corner c.
XY sideDirectionFrom(XY c, int side) {
    if (side == 0 || side == 2) return {c.x < 0 ? S : -S, 0};
    return {0, c.y < 0 ? S : -S};
}

} // namespace

namespace {
Wide det(const BlockPoint &a, const BlockPoint &b, const BlockPoint &c,
         const BlockPoint &d) {
    Wide ux = b.x - a.x, uy = b.y - a.y, uz = b.z - a.z;
    Wide vx = c.x - a.x, vy = c.y - a.y, vz = c.z - a.z;
    Wide wx = d.x - a.x, wy = d.y - a.y, wz = d.z - a.z;
    return ux * (vy * wz - vz * wy) - uy * (vx * wz - vz * wx) +
           uz * (vx * wy - vy * wx);
}
} // namespace

// ---------------------------------------------------------------- diagrams

long Diagram::linkingNumber(size_t i, size_t j) const {
    long twice = 0;
    for (const CrossingInfo &c : crossings)
        if ((c.over == i && c.under == j) || (c.over == j && c.under == i))
            twice += c.sign;
    return twice / 2;
}

regina::Link Diagram::link() const {
    std::vector<int> signs;
    signs.reserve(crossings.size());
    for (const CrossingInfo &c : crossings) signs.push_back(c.sign);
    return regina::Link::fromData(signs.begin(), signs.end(), gauss.begin(), gauss.end());
}

// ---------------------------------------------------------------- the drawer

struct DiagramDrawer::Impl {
    regina::Triangulation<3> tri;
    size_t blocks;
    std::vector<std::unordered_map<size_t, BlockPoint>> pos;   // block -> vertex -> point
    std::vector<std::map<XY, std::pair<long, long>>> cornerVerts; // block -> corner -> (bottom, top)
    std::unordered_map<size_t, std::vector<size_t>> blocksOf;
    std::unordered_map<size_t, int> apex;                 // vertex -> +1 top, -1 bottom
    std::vector<std::array<size_t, 4>> nbr;               // block, side -> block
    std::vector<std::vector<size_t>> edgeBlocks;          // edge -> blocks containing it

    struct Piece {
        size_t block;
        BlockPoint p0, p1;
        long v0 = -1, v1 = -1; // curve vertex at an end, or -1
        size_t comp = 0, order = 0;
    };

    // A crossing event's position along a piece: num / den in [0, 1], or
    // (2, 1) for "at the piece's far end, after everything else".
    struct Param {
        Wide num, den;
        bool operator<(const Param &o) const { return num * o.den < o.num * den; }
    };

    Impl(const regina::Triangulation<3> &t, size_t crossings);
    Diagram draw(const std::vector<EdgeCycle> &curves, unsigned seed) const;
    void route(std::vector<Piece> &out, size_t p, size_t r, int sign,
               unsigned seed) const;
    std::pair<XY, XY> gate(size_t b, int side, size_t nb, int which) const;
};

DiagramDrawer::Impl::Impl(const regina::Triangulation<3> &t, size_t crossings)
    : tri(t), blocks(crossings) {
    const auto &c = blockCoordinates();
    if (tri.size() < 14 * blocks || !tri.isOrientable())
        throw regina::InvalidArgument(
            "DiagramDrawer: not a knotbuilder triangulation of this many crossings");
    pos.resize(blocks);
    cornerVerts.resize(blocks);
    std::optional<bool> positive;
    for (size_t k = 0; k < blocks; ++k)
        for (int j = 0; j < 14; ++j) {
            const regina::Tetrahedron<3> *tet = tri.tetrahedron(14 * k + j);
            Wide d = det(c[j][0], c[j][1], c[j][2], c[j][3]);
            bool pos1 = (d > 0) == (tet->orientation() > 0);
            if (positive && *positive != pos1)
                throw regina::InvalidArgument(
                    "DiagramDrawer: blocks are not all oriented alike");
            positive = pos1;
            for (int i = 0; i < 4; ++i) {
                size_t v = tet->vertex(i)->index();
                auto [it, fresh] = pos[k].emplace(v, c[j][i]);
                if (!fresh && !(it->second == c[j][i]))
                    throw regina::InvalidArgument(
                        "DiagramDrawer: a vertex sits at two points of one block "
                        "(the diagram is not reduced)");
            }
        }
    for (size_t k = 0; k < blocks; ++k)
        for (const auto &[v, p] : pos[k]) {
            blocksOf[v].push_back(k);
            if (std::abs(p.x) == S && std::abs(p.y) == S) {
                auto &slot = cornerVerts[k].try_emplace(XY{p.x, p.y}, -1, -1).first->second;
                (p.z == 0 ? slot.first : slot.second) = static_cast<long>(v);
            }
        }
    // Cone apexes: the vertex of a cone tetrahedron opposite its face on a
    // block, top or bottom by that face's height.
    for (size_t i = 14 * blocks; i < tri.size(); ++i) {
        const regina::Tetrahedron<3> *tet = tri.tetrahedron(i);
        for (int f = 0; f < 4; ++f) {
            const regina::Tetrahedron<3> *adj = tet->adjacentTetrahedron(f);
            if (!adj || adj->index() >= 14 * blocks) continue;
            size_t k = adj->index() / 14;
            size_t someBase = tet->vertex(f == 0 ? 1 : 0)->index();
            apex[tet->vertex(f)->index()] = pos[k].at(someBase).z == S ? 1 : -1;
        }
    }
    // Neighbours: across side s, the other block holding that side's strand point.
    nbr.resize(blocks);
    for (size_t k = 0; k < blocks; ++k)
        for (int s = 0; s < 4; ++s) {
            long mid = -1;
            for (const auto &[v, p] : pos[k])
                if (p.z == 0 && onSide(p, s) && (s % 2 == 0 ? p.x == 0 : p.y == 0))
                    mid = static_cast<long>(v);
            if (mid < 0)
                throw regina::InvalidArgument("DiagramDrawer: a side has no strand point");
            const auto &bs = blocksOf.at(static_cast<size_t>(mid));
            if (bs.size() != 2 || bs[0] == bs[1])
                throw regina::InvalidArgument(
                    "DiagramDrawer: a strand runs from a crossing to itself "
                    "(the diagram is not reduced)");
            nbr[k][s] = bs[0] == k ? bs[1] : bs[0];
        }
    edgeBlocks.resize(tri.countEdges());
    for (size_t e = 0; e < tri.countEdges(); ++e) {
        std::set<size_t> bs;
        for (const auto &emb : *tri.edge(e))
            if (emb.tetrahedron()->index() < 14 * blocks)
                bs.insert(emb.tetrahedron()->index() / 14);
        edgeBlocks[e].assign(bs.begin(), bs.end());
    }
}

std::pair<XY, XY> DiagramDrawer::Impl::gate(size_t b, int side, size_t nb, int which) const {
    // A point of the side shared by b and nb, the same point seen from both
    // frames: a fixed fraction of the way between the side's two bottom
    // corners, taken in vertex-index order so both frames agree.
    std::vector<size_t> corners;
    for (const auto &[v, p] : pos[b])
        if (p.z == 0 && onSide(p, side) && std::abs(p.x) == S && std::abs(p.y) == S)
            corners.push_back(v);
    if (corners.size() != 2) throw Degenerate("gate: side without two corners");
    std::sort(corners.begin(), corners.end());
    // 2/7 for top routes, 3/7 for bottom ones: never a strand point (1/2).
    const BlockCoord num = which > 0 ? 2 : 3, den = 7;
    auto at = [&](size_t k) {
        const BlockPoint &a = pos[k].at(corners[0]), &z = pos[k].at(corners[1]);
        return XY{a.x + (z.x - a.x) * num / den, a.y + (z.y - a.y) * num / den};
    };
    return {at(b), at(nb)};
}

void DiagramDrawer::Impl::route(std::vector<Piece> &out, size_t p, size_t r, int sign,
                         unsigned seed) const {
    // Breadth-first over blocks, from one holding p to one holding r.
    const auto &from = blocksOf.at(p);
    const auto &to = blocksOf.at(r);
    std::unordered_map<size_t, std::pair<long, int>> prev; // block -> (block, side)
    std::queue<size_t> q;
    for (size_t b : from) {
        prev[b] = {-1, -1};
        q.push(b);
    }
    long goal = -1;
    while (!q.empty()) {
        size_t b = q.front();
        q.pop();
        if (std::find(to.begin(), to.end(), b) != to.end()) {
            goal = static_cast<long>(b);
            break;
        }
        for (int s = 0; s < 4; ++s) {
            size_t n = nbr[b][s];
            if (!prev.contains(n)) {
                prev[n] = {static_cast<long>(b), s};
                q.push(n);
            }
        }
    }
    std::vector<size_t> path{static_cast<size_t>(goal)};
    while (prev.at(path.back()).first >= 0)
        path.push_back(static_cast<size_t>(prev.at(path.back()).first));
    std::reverse(path.begin(), path.end());

    unsigned salt = seed * 2654435761u + static_cast<unsigned>(p * 40503 + r);
    auto interior = [&]() {
        salt = salt * 1103515245u + 12345u;
        BlockCoord x = static_cast<BlockCoord>((salt >> 8) % 1201) - 600;
        salt = salt * 1103515245u + 12345u;
        BlockCoord y = static_cast<BlockCoord>((salt >> 8) % 1201) - 600;
        return XY{x * S / 971, y * S / 977};
    };
    BlockCoord z = sign * ROUTE_HEIGHT;
    auto next = [&]() { z += sign * ROUTE_STEP; return z; };
    const BlockPoint &ps = pos[path[0]].at(p);
    BlockPoint cur{ps.x, ps.y, next()};
    long curV = static_cast<long>(p);
    for (size_t i = 0; i < path.size(); ++i) {
        size_t b = path[i];
        XY m = interior();
        BlockPoint mid{m.x, m.y, next()};
        out.push_back({b, cur, mid, curV, -1});
        if (i + 1 < path.size()) {
            size_t n = path[i + 1];
            int side = -1;
            for (int s = 0; s < 4; ++s)
                if (nbr[b][s] == n) side = s;
            auto [gb, gn] = gate(b, side, n, sign);
            BlockCoord zg = next();
            out.push_back({b, mid, {gb.x, gb.y, zg}, -1, -1});
            cur = {gn.x, gn.y, zg};
            curV = -1;
        } else {
            const BlockPoint &pr = pos[b].at(r);
            out.push_back({b, mid, {pr.x, pr.y, next()}, -1, static_cast<long>(r)});
        }
    }
}

Diagram DiagramDrawer::Impl::draw(const std::vector<EdgeCycle> &curves,
                           unsigned seed) const {
    std::vector<Piece> pieces;
    std::vector<std::vector<size_t>> cyclesV;          // vertices, per component
    const BlockCoord eps = S / 200;
    const XY tgt{S / 5 + static_cast<BlockCoord>(seed % 13) * S / 311,
                 -S / 7 + static_cast<BlockCoord>(seed % 11) * S / 293};

    for (size_t c = 0; c < curves.size(); ++c) {
        const EdgeCycle &ec = curves[c];
        const size_t m = ec.size();
        if (m < 2) throw regina::InvalidArgument("draw: a curve has fewer than two edges");
        std::vector<size_t> vs(m);
        for (size_t k = 0; k < m; ++k) {
            const regina::Edge<3> *e = tri.edge(ec[k].edge);
            size_t tail = e->vertex(ec[k].reversed ? 1 : 0)->index();
            size_t head = e->vertex(ec[k].reversed ? 0 : 1)->index();
            vs[k] = tail;
            const regina::Edge<3> *nx = tri.edge(ec[(k + 1) % m].edge);
            size_t nextTail = nx->vertex(ec[(k + 1) % m].reversed ? 1 : 0)->index();
            if (head != nextTail)
                throw regina::InvalidArgument("draw: a curve is not a closed path");
        }
        size_t start = 0;
        while (start < m && apex.contains(vs[start])) ++start;
        if (start == m) throw regina::InvalidArgument("draw: a curve lies in a cone");
        std::rotate(vs.begin(), vs.begin() + start, vs.end());
        std::vector<DirectedEdge> es(ec.begin(), ec.end());
        std::rotate(es.begin(), es.begin() + start, es.end());
        cyclesV.push_back(vs);

        size_t order = 0;
        auto push = [&](Piece pc) {
            pc.comp = c;
            pc.order = order++;
            pieces.push_back(pc);
        };
        for (size_t k = 0; k < m;) {
            size_t p = vs[k], q = vs[(k + 1) % m];
            if (auto a = apex.find(q); a != apex.end()) {
                size_t r = vs[(k + 2) % m];
                std::vector<Piece> rp;
                route(rp, p, r, a->second, seed);
                for (Piece &pc : rp) push(pc);
                k += 2;
                continue;
            }
            const auto &ks = edgeBlocks[es[k].edge];
            if (ks.empty())
                throw regina::InvalidArgument("draw: an edge lies in no block");
            size_t b = ks.front();
            const BlockPoint &P0 = pos[b].at(p), &P1 = pos[b].at(q);
            if (ks.size() == 1) {
                push({b, P0, P1, static_cast<long>(p), static_cast<long>(q)});
                ++k;
                continue;
            }
            // A wall or corner edge: a tent into block b, the push growing
            // with height (see the header).
            auto lift = [&](BlockCoord x, BlockCoord y, BlockCoord zz) {
                Wide f = Wide(S) + zz; // S..2S
                BlockCoord dx = static_cast<BlockCoord>(Wide(tgt.x - x) * f * eps / (Wide(S) * S));
                BlockCoord dy = static_cast<BlockCoord>(Wide(tgt.y - y) * f * eps / (Wide(S) * S));
                BlockCoord z2 = static_cast<BlockCoord>(Wide(zz) * (S - 2 * eps) / S + eps);
                return BlockPoint{x + dx, y + dy, z2};
            };
            std::vector<BlockPoint> pts{P0};
            if (P0.x != P1.x || P0.y != P1.y) {
                pts.push_back(lift((P0.x + P1.x) / 2, (P0.y + P1.y) / 2, (P0.z + P1.z) / 2));
            } else {
                pts.push_back(lift(P0.x, P0.y, (2 * P0.z + P1.z) / 3));
                XY t2{tgt.x + S / 3, tgt.y - S / 4};
                BlockCoord zz = (P0.z + 2 * P1.z) / 3;
                pts.push_back({P0.x + (t2.x - P0.x) * eps / S,
                               P0.y + (t2.y - P0.y) * eps / S,
                               static_cast<BlockCoord>(Wide(zz) * (S - 2 * eps) / S + eps)});
            }
            pts.push_back(P1);
            for (size_t j = 0; j + 1 < pts.size(); ++j)
                push({b, pts[j], pts[j + 1], j == 0 ? static_cast<long>(p) : -1,
                      j + 2 == pts.size() ? static_cast<long>(q) : -1});
            ++k;
        }
    }

    // ---- crossings inside blocks
    std::vector<size_t> compLen(curves.size(), 0);
    for (const Piece &pc : pieces) ++compLen[pc.comp];
    auto consecutive = [&](const Piece &a, const Piece &b) {
        if (a.comp != b.comp) return false;
        size_t d = a.order > b.order ? a.order - b.order : b.order - a.order;
        return d == 1 || d + 1 == compLen[a.comp];
    };
    auto sameCornerOtherVertex = [](const Piece &a, const Piece &b) {
        const std::pair<BlockPoint, long> ea[2] = {{a.p0, a.v0}, {a.p1, a.v1}};
        const std::pair<BlockPoint, long> eb[2] = {{b.p0, b.v0}, {b.p1, b.v1}};
        for (const auto &[p, v] : ea)
            for (const auto &[q, w] : eb)
                if (v >= 0 && w >= 0 && v != w && p.x == q.x && p.y == q.y)
                    return true;
        return false;
    };

    struct Event {
        Param at;
        size_t crossing;
        bool over;
    };
    struct Cross {
        bool corner;
        XY dOver, dUnder;                           // plain
        std::vector<std::pair<long, bool>> ccw;     // corner: (vertex, incoming?) ccw
        long top = -1, bottom = -1;
        size_t overComp = 0, underComp = 0;
    };
    std::vector<std::vector<Event>> events(pieces.size());
    std::vector<Cross> xs;

    std::unordered_map<size_t, std::vector<size_t>> byBlock;
    for (size_t i = 0; i < pieces.size(); ++i) byBlock[pieces[i].block].push_back(i);
    for (const auto &[blk, idxs] : byBlock) {
        for (size_t x = 0; x < idxs.size(); ++x)
            for (size_t y = x + 1; y < idxs.size(); ++y) {
                const Piece &A = pieces[idxs[x]], &B = pieces[idxs[y]];
                if (consecutive(A, B) || sameCornerOtherVertex(A, B)) continue;
                BlockCoord dx = A.p1.x - A.p0.x, dy = A.p1.y - A.p0.y;
                BlockCoord ex = B.p1.x - B.p0.x, ey = B.p1.y - B.p0.y;
                BlockCoord rx = B.p0.x - A.p0.x, ry = B.p0.y - A.p0.y;
                Wide den = cross(dx, dy, ex, ey);
                if (den == 0) {
                    if (cross(rx, ry, dx, dy) == 0 && (dx || dy)) {
                        Wide L = Wide(dx) * dx + Wide(dy) * dy;
                        Wide a = Wide(rx) * dx + Wide(ry) * dy;
                        Wide b = Wide(B.p1.x - A.p0.x) * dx + Wide(B.p1.y - A.p0.y) * dy;
                        if (std::max(std::min(a, b), Wide(0)) <= std::min(std::max(a, b), L))
                            throw Degenerate("collinear projections overlap");
                    }
                    continue;
                }
                Wide sN = cross(rx, ry, ex, ey), tN = cross(rx, ry, dx, dy);
                if (den < 0) { den = -den; sN = -sN; tN = -tN; }
                bool inside = sN > 0 && sN < den && tN > 0 && tN < den;
                bool closed = sN >= 0 && sN <= den && tN >= 0 && tN <= den;
                if (!closed) continue;
                if (!inside) throw Degenerate("projections touch at an end");
                Wide za = Wide(A.p0.z) * den + sN * (A.p1.z - A.p0.z);
                Wide zb = Wide(B.p0.z) * den + tN * (B.p1.z - B.p0.z);
                if (za == zb) throw Degenerate("pieces meet in space");
                bool aOver = za > zb;
                const Piece &O = aOver ? A : B, &U = aOver ? B : A;
                size_t oi = aOver ? idxs[x] : idxs[y], ui = aOver ? idxs[y] : idxs[x];
                Cross cr;
                cr.corner = false;
                cr.dOver = {O.p1.x - O.p0.x, O.p1.y - O.p0.y};
                cr.dUnder = {U.p1.x - U.p0.x, U.p1.y - U.p0.y};
                cr.overComp = O.comp;
                cr.underComp = U.comp;
                size_t id = xs.size();
                xs.push_back(cr);
                events[oi].push_back({{aOver ? sN : tN, den}, id, true});
                events[ui].push_back({{aOver ? tN : sN, den}, id, false});
            }
    }

    // ---- crossings at corners
    std::unordered_map<long, std::pair<long, long>> atVertex; // vertex -> (incoming, outgoing piece)
    for (size_t i = 0; i < pieces.size(); ++i) {
        if (pieces[i].v1 >= 0) {
            auto &slot = atVertex.try_emplace(pieces[i].v1, -1, -1).first->second;
            slot.first = static_cast<long>(i);
        }
        if (pieces[i].v0 >= 0) {
            auto &slot = atVertex.try_emplace(pieces[i].v0, -1, -1).first->second;
            slot.second = static_cast<long>(i);
        }
    }
    std::set<std::pair<long, long>> checked;
    for (const auto &[v, io] : atVertex) {
        const Piece &out0 = pieces[io.second];
        const BlockPoint &P = out0.p0;
        if (std::abs(P.x) != S || std::abs(P.y) != S) continue;
        const auto &cv = cornerVerts[out0.block].at({P.x, P.y});
        long w = cv.first == v ? cv.second : cv.first;
        if (w < 0 || !atVertex.contains(w)) continue;
        std::pair<long, long> key{std::min(v, w), std::max(v, w)};
        if (checked.contains(key)) continue;
        checked.insert(key);
        // This includes a curve running up the corner's own vertical edge
        // (v and w consecutive). Its tent leaves the corner point and comes
        // back to it, so the shadow passes the corner twice -- once on the
        // way in at v, once on the way out at w -- and those two passes
        // cross whenever their ends alternate, exactly as for any other
        // two passes. (Skipping this case once lost that crossing, which
        // left a non-planar diagram.)

        // The blocks around this corner, counterclockwise, starting from out0's.
        const long regionBottom = cv.first;
        std::map<std::pair<size_t, XY>, size_t> rank;
        size_t b = out0.block;
        XY cr{P.x, P.y};
        for (size_t guard = 0; guard < 4 * blocks + 4; ++guard) {
            if (rank.contains({b, cr})) break;
            rank[{b, cr}] = rank.size();
            size_t nb = nbr[b][cornerSides(cr).end];
            const BlockPoint &np = pos[nb].at(static_cast<size_t>(regionBottom));
            b = nb;
            cr = {np.x, np.y};
        }
        struct End {
            long vertex;
            bool incoming;
            size_t block;
            XY corner, dir;
        };
        std::vector<End> ends;
        for (long vv : {v, w}) {
            auto [inc, outg] = atVertex.at(vv);
            const Piece &pi = pieces[inc], &po = pieces[outg];
            ends.push_back({vv, true, pi.block, {pi.p1.x, pi.p1.y},
                            {pi.p0.x - pi.p1.x, pi.p0.y - pi.p1.y}});
            ends.push_back({vv, false, po.block, {po.p0.x, po.p0.y},
                            {po.p1.x - po.p0.x, po.p1.y - po.p0.y}});
        }
        for (const End &e : ends) {
            if (!rank.contains({e.block, e.corner}))
                throw Degenerate("an end's block is not around its corner");
            CornerSides cs = cornerSides(e.corner);
            XY s0 = sideDirectionFrom(e.corner, cs.start);
            XY s1 = sideDirectionFrom(e.corner, cs.end);
            if (cross(s0.x, s0.y, e.dir.x, e.dir.y) <= 0 ||
                cross(e.dir.x, e.dir.y, s1.x, s1.y) <= 0)
                throw Degenerate("an end runs along a side at its corner");
        }
        std::vector<size_t> ord(4);
        std::iota(ord.begin(), ord.end(), 0);
        bool tie = false;
        std::sort(ord.begin(), ord.end(), [&](size_t a, size_t bb) {
            size_t ra = rank.at({ends[a].block, ends[a].corner});
            size_t rb = rank.at({ends[bb].block, ends[bb].corner});
            if (ra != rb) return ra < rb;
            Wide c2 = cross(ends[a].dir.x, ends[a].dir.y, ends[bb].dir.x, ends[bb].dir.y);
            if (c2 == 0 && a != bb) tie = true;
            return c2 > 0;
        });
        if (tie) throw Degenerate("two ends leave a corner in one direction");
        if (!(ends[ord[0]].vertex == ends[ord[2]].vertex &&
              ends[ord[1]].vertex == ends[ord[3]].vertex))
            continue; // the two passes touch without crossing
        long top = pos[out0.block].at(static_cast<size_t>(v)).z == S ? v : w;
        long bottom = top == v ? w : v;
        Cross c;
        c.corner = true;
        for (size_t k : ord) c.ccw.push_back({ends[k].vertex, ends[k].incoming});
        c.top = top;
        c.bottom = bottom;
        c.overComp = pieces[atVertex.at(top).first].comp;
        c.underComp = pieces[atVertex.at(bottom).first].comp;
        size_t id = xs.size();
        xs.push_back(c);
        events[atVertex.at(top).first].push_back({{2, 1}, id, true});
        events[atVertex.at(bottom).first].push_back({{2, 1}, id, false});
    }

    // ---- arcs, PD code, signs
    struct EndLabels {
        long underIn = 0, underOut = 0, overIn = 0, overOut = 0;
    };
    std::vector<EndLabels> lab(xs.size());
    Diagram out;
    out.components = curves.size();
    long label = 0;
    for (size_t c = 0; c < curves.size(); ++c) {
        std::vector<size_t> idx;
        for (size_t i = 0; i < pieces.size(); ++i)
            if (pieces[i].comp == c) idx.push_back(i);
        std::sort(idx.begin(), idx.end(),
                  [&](size_t a, size_t b) { return pieces[a].order < pieces[b].order; });
        std::vector<std::pair<size_t, bool>> seq;
        for (size_t i : idx) {
            auto evs = events[i];
            std::sort(evs.begin(), evs.end(),
                      [](const Event &a, const Event &b) { return a.at < b.at; });
            // Two crossings at one point of a piece mean three pieces through
            // one point of the plane: their order along each strand is then
            // arbitrary, and a wrong order can still be planar. Never guess.
            for (size_t k = 1; k < evs.size(); ++k)
                if (!(evs[k - 1].at < evs[k].at))
                    throw Degenerate("three pieces cross at one point");
            for (const Event &e : evs) seq.push_back({e.crossing, e.over});
        }
        std::vector<long> &g = out.gauss.emplace_back();
        for (auto [id, over] : seq)
            g.push_back(over ? static_cast<long>(id) + 1 : -static_cast<long>(id) - 1);
        if (seq.empty()) {
            out.crossingless.push_back(c);
            continue;
        }
        const long n = static_cast<long>(seq.size());
        for (long k = 0; k < n; ++k) {
            auto [id, over] = seq[k];
            long in = label + k + 1, outL = label + (k + 1) % n + 1;
            if (over) { lab[id].overIn = in; lab[id].overOut = outL; }
            else { lab[id].underIn = in; lab[id].underOut = outL; }
        }
        label += n;
    }
    for (size_t id = 0; id < xs.size(); ++id) {
        const Cross &c = xs[id];
        const EndLabels &e = lab[id];
        bool overFromRight; // ccw from under-in: over-in comes next
        if (!c.corner) {
            overFromRight = cross(c.dUnder.x, c.dUnder.y, c.dOver.x, c.dOver.y) > 0;
        } else {
            size_t k0 = 0;
            while (!(c.ccw[k0].first == c.bottom && c.ccw[k0].second)) ++k0;
            overFromRight = c.ccw[(k0 + 1) % 4].second; // next end is over, incoming
        }
        if (overFromRight)
            out.pd.push_back({e.underIn, e.overIn, e.underOut, e.overOut});
        else
            out.pd.push_back({e.underIn, e.overOut, e.underOut, e.overIn});
        // [under-in, over-in, under-out, over-out] is a left-handed crossing.
        out.crossings.push_back({c.overComp, c.underComp, overFromRight ? -1 : 1});
    }
    // Every closed curve in S^3 has a planar diagram, so a drawing that is
    // not one is a defect here, and must never reach a name: a crossing
    // missed or misordered leaves a virtual diagram, which simplify() can
    // carry anywhere.
    if (!out.pd.empty() && !out.link().isClassical())
        throw NonPlanar("the drawing is not a planar diagram");
    return out;
}

DiagramDrawer::DiagramDrawer(const regina::Triangulation<3> &tri, size_t crossings)
    : impl_(std::make_unique<Impl>(tri, crossings)) {}
DiagramDrawer::~DiagramDrawer() = default;
DiagramDrawer::DiagramDrawer(DiagramDrawer &&) noexcept = default;
DiagramDrawer &DiagramDrawer::operator=(DiagramDrawer &&) noexcept = default;

Diagram DiagramDrawer::drawWithSeed(const std::vector<EdgeCycle> &curves,
                             unsigned seed) const {
    return impl_->draw(curves, seed);
}

Diagram DiagramDrawer::draw(const std::vector<EdgeCycle> &curves, int attempts) const {
    for (int a = 0;; ++a) {
        try {
            return impl_->draw(curves, static_cast<unsigned>(a) + 1);
        } catch (const Degenerate &) {
            if (a + 1 >= attempts) throw;
        }
    }
}

std::vector<EdgeCycle>
DiagramDrawer::cyclesOf(const std::vector<const regina::Edge<3> *> &edges,
                 const std::vector<bool> &reversed) {
    std::unordered_map<size_t, size_t> fromTail; // tail vertex -> position in edges
    for (size_t i = 0; i < edges.size(); ++i)
        fromTail[edges[i]->vertex(reversed[i] ? 1 : 0)->index()] = i;
    std::vector<bool> used(edges.size(), false);
    std::vector<EdgeCycle> out;
    for (size_t s = 0; s < edges.size(); ++s) {
        if (used[s]) continue;
        EdgeCycle cyc;
        size_t i = s;
        while (!used[i]) {
            used[i] = true;
            cyc.push_back({edges[i]->index(), static_cast<bool>(reversed[i])});
            size_t head = edges[i]->vertex(reversed[i] ? 0 : 1)->index();
            i = fromTail.at(head);
        }
        out.push_back(std::move(cyc));
    }
    return out;
}

} // namespace knotbuilder
