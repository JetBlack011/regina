//
//  todiagram_test.cpp
//
//  knotbuilder::DiagramDrawer, against independent answers:
//
//    1. The block model is an exact embedding of knotbuilder's block.
//    2. Drawing knotbuilder's own link from its triangulation gives back
//       the input diagram exactly -- same signature, with neither mirror
//       nor reversal allowed -- for knots and for several orientations of
//       links, and the drawn linking numbers are the input's.
//    3. Random closed curves (many through walls, corners and cone points):
//       the drawing simplifies to an unknot exactly when the complement
//       drilled from the triangulation itself is a solid torus.
//
//  An optional argument sweeps a whole table instead (slow; not in ctest):
//
//    ./diagramdrawer_test ../../../../../cobordism-atlas/data/<table>.csv [max crossings]
//

#include <fstream>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <functional>
#include <map>
#include <optional>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include <link/link.h>
#include <triangulation/dim3.h>

#include "diagramtriangulation/todiagram.h"
#include "diagramtriangulation/fromdiagram.h"
#include "linknaming/complement/linkcomplement.h"

using namespace knotbuilder;

static int passed = 0;
static int failed_count = 0;

#define EXPECT_EQ(actual, expected, desc)                                     \
    do {                                                                      \
        auto _a = (actual);                                                   \
        auto _e = (expected);                                                 \
        if (_a == _e) {                                                       \
            ++passed;                                                         \
        } else {                                                              \
            std::cout << "  FAIL: " << (desc) << "\n        expected " << _e \
                      << ", got " << _a << "\n";                              \
            ++failed_count;                                                   \
        }                                                                     \
    } while (0)

namespace {

std::vector<std::array<long, 4>> parsePD(const std::string &s) {
    std::vector<long> n;
    long cur = -1;
    for (char ch : s) {
        if (ch >= '0' && ch <= '9') cur = (cur < 0 ? 0 : cur * 10) + (ch - '0');
        else if (cur >= 0) { n.push_back(cur); cur = -1; }
    }
    if (cur >= 0) n.push_back(cur);
    std::vector<std::array<long, 4>> out;
    for (size_t i = 0; i + 3 < n.size(); i += 4)
        out.push_back({n[i], n[i + 1], n[i + 2], n[i + 3]});
    return out;
}

// The signed linking matrix of a regina::Link, components in Regina's order.
std::vector<std::vector<long>> linkingMatrix(const regina::Link &l) {
    std::vector<std::vector<long>> m(l.countComponents(),
                                     std::vector<long>(l.countComponents(), 0));
    std::unordered_map<const regina::Crossing *, std::vector<size_t>> on;
    for (size_t c = 0; c < l.countComponents(); ++c) {
        regina::StrandRef start = l.component(c);
        if (!start) continue;
        regina::StrandRef s = start;
        do {
            on[s.crossing()].push_back(c);
            s = s.next();
        } while (s != start);
    }
    for (const auto &[x, comps] : on)
        if (comps.size() == 2 && comps[0] != comps[1]) {
            m[comps[0]][comps[1]] += x->sign();
            m[comps[1]][comps[0]] += x->sign();
        }
    for (auto &row : m)
        for (long &v : row) v /= 2;
    return m;
}

// Whether some relabelling of components carries matrix a onto b exactly.
bool sameUpToRelabelling(const std::vector<std::vector<long>> &a,
                         const std::vector<std::vector<long>> &b) {
    if (a.size() != b.size()) return false;
    std::vector<size_t> perm(a.size());
    std::iota(perm.begin(), perm.end(), 0);
    do {
        bool ok = true;
        for (size_t i = 0; i < a.size() && ok; ++i)
            for (size_t j = 0; j < a.size() && ok; ++j)
                ok = a[i][j] == b[perm[i]][perm[j]];
        if (ok) return true;
    } while (std::next_permutation(perm.begin(), perm.end()));
    return false;
}

// Returns whether the drawing reproduced the input exactly: the same
// oriented diagram (signature with neither mirror nor reversal allowed),
// and -- independently of that signature -- the same signed linking
// matrix up to relabelling components, so that no component is drawn
// backwards.
bool reproduce(const std::string &name, const std::string &pd, bool verbose) {
    auto input = parsePD(pd);
    auto [tri, edges, reversed] = knotbuilder::buildLink(knotbuilder::parsePDCode(pd));
    DiagramDrawer drawer(tri, input.size());
    Diagram d = drawer.draw(DiagramDrawer::cyclesOf(edges, reversed));
    regina::Link in = regina::Link::fromPD(input.begin(), input.end());
    regina::Link out = d.link();
    if (in.sig<2>(false, false, true) != out.sig<2>(false, false, true)) {
        if (verbose)
            std::cout << "  " << name << ": drawn diagram differs from the input\n";
        return false;
    }
    std::vector<std::vector<long>> drawnLk(d.components, std::vector<long>(d.components, 0));
    for (size_t i = 0; i < d.components; ++i)
        for (size_t j = 0; j < d.components; ++j)
            if (i != j) drawnLk[i][j] = d.linkingNumber(i, j);
    if (!sameUpToRelabelling(linkingMatrix(in), drawnLk) ||
        !sameUpToRelabelling(linkingMatrix(in), linkingMatrix(out))) {
        if (verbose)
            std::cout << "  " << name << ": linking matrices differ\n";
        return false;
    }
    return true;
}

void test_block_model() {
    EXPECT_EQ(verifyBlockModel(), std::string(),
              "the block is exactly the box: tetrahedra tile it, boundary "
              "triangles tile its faces");
}

void test_reproduce_battery() {
    const std::vector<std::pair<std::string, std::string>> rows = {
        {"3_1", "[[1;5;2;4];[3;1;4;6];[5;3;6;2]]"},
        {"4_1", "[[4;2;5;1];[8;6;1;5];[6;3;7;4];[2;7;3;8]]"},
        {"5_1", "[[2;8;3;7];[4;10;5;9];[6;2;7;1];[8;4;9;3];[10;6;1;5]]"},
        {"8_20", "[[1;7;2;6];[4;13;5;14];[5;9;6;8];[7;3;8;2];[10;15;11;16];"
                 "[12;9;13;10];[14;3;15;4];[16;11;1;12]]"},
        {"L2a1{0}", "PD[X[4; 1; 3; 2]; X[2; 3; 1; 4]]"},
        {"L2a1{1}", "PD[X[4; 2; 3; 1]; X[2; 4; 1; 3]]"},
        {"L4a1{0}", "PD[X[6; 1; 7; 2]; X[8; 3; 5; 4]; X[2; 5; 3; 6]; X[4; 7; 1; 8]]"},
        {"L4a1{1}", "PD[X[6; 2; 7; 1]; X[8; 4; 5; 3]; X[2; 8; 3; 7]; X[4; 6; 1; 5]]"},
        {"L6a3{0}", "PD[X[8; 1; 9; 2]; X[2; 9; 3; 10]; X[10; 3; 11; 4]; "
                    "X[12; 5; 7; 6]; X[6; 7; 1; 8]; X[4; 11; 5; 12]]"},
        {"L6a3{1}", "PD[X[10; 2; 11; 1]; X[2; 10; 3; 9]; X[8; 4; 9; 3]; "
                    "X[12; 6; 7; 5]; X[6; 12; 1; 11]; X[4; 8; 5; 7]]"},
        {"L6a4{0;0}", "PD[X[6; 1; 7; 2]; X[12; 8; 9; 7]; X[4; 12; 1; 11]; "
                      "X[10; 5; 11; 6]; X[8; 4; 5; 3]; X[2; 9; 3; 10]]"},
        {"L6a4{1;1}", "PD[X[6; 2; 7; 1]; X[12; 6; 9; 5]; X[4; 9; 1; 10]; "
                      "X[10; 7; 11; 8]; X[8; 3; 5; 4]; X[2; 12; 3; 11]]"},
        {"L6n1{1;0}", "PD[X[6; 2; 7; 1]; X[12; 5; 9; 6]; X[4; 12; 1; 11]; "
                      "X[7; 10; 8; 11]; X[3; 5; 4; 8]; X[9; 3; 10; 2]]"},
    };
    for (const auto &[name, pd] : rows)
        EXPECT_EQ(reproduce(name, pd, true), true,
                  name + ": drawing knotbuilder's own link reproduces the "
                         "input diagram (no mirror, no reversal)");
}

// The orientation checks can fail: reversing one component of the Hopf link
// must change both the oriented diagram and its linking number.
void test_reversed_component_is_caught() {
    const std::string pd = "PD[X[4; 1; 3; 2]; X[2; 3; 1; 4]]"; // L2a1{0}
    auto input = parsePD(pd);
    auto [tri, edges, reversed] = knotbuilder::buildLink(knotbuilder::parsePDCode(pd));
    DiagramDrawer drawer(tri, input.size());
    auto cycles = DiagramDrawer::cyclesOf(edges, reversed);
    Diagram asIs = drawer.draw(cycles);
    std::reverse(cycles[1].begin(), cycles[1].end());
    for (auto &de : cycles[1]) de.reversed = !de.reversed;
    Diagram flipped = drawer.draw(cycles);
    regina::Link in = regina::Link::fromPD(input.begin(), input.end());
    EXPECT_EQ(in.sig<2>(false, false, true) == flipped.link().sig<2>(false, false, true), false,
              "reversing one component changes the oriented diagram's signature");
    EXPECT_EQ(flipped.linkingNumber(0, 1), -asIs.linkingNumber(0, 1),
              "and negates the linking number");
    EXPECT_EQ(asIs.linkingNumber(0, 1) != 0, true, "which is nonzero for the Hopf link");
}

// Random simple closed edge paths, by randomized depth-first extension.
std::vector<EdgeCycle> randomCycle(const regina::Triangulation<3> &tri,
                                   std::mt19937 &rng, size_t maxLen) {
    std::unordered_map<size_t, std::vector<std::pair<size_t, size_t>>> adj;
    for (const regina::Edge<3> *e : tri.edges()) {
        size_t a = e->vertex(0)->index(), b = e->vertex(1)->index();
        if (a == b) continue;
        adj[a].push_back({b, e->index()});
        adj[b].push_back({a, e->index()});
    }
    for (int attempt = 0; attempt < 200; ++attempt) {
        size_t s = rng() % tri.countVertices();
        std::vector<size_t> path{s};
        std::vector<size_t> epath;
        std::set<size_t> seen{s};
        while (path.size() <= maxLen) {
            size_t cur = path.back();
            std::vector<size_t> closers;
            for (auto [w, e] : adj[cur])
                if (w == s && path.size() >= 3 && (epath.empty() || e != epath.back()))
                    closers.push_back(e);
            if (!closers.empty() && rng() % 3 == 0) {
                epath.push_back(closers[rng() % closers.size()]);
                // direct each edge along the walk
                EdgeCycle cyc;
                for (size_t k = 0; k < epath.size(); ++k) {
                    size_t from = path[k];
                    const regina::Edge<3> *e = tri.edge(epath[k]);
                    cyc.push_back({epath[k], e->vertex(0)->index() != from});
                }
                return {cyc};
            }
            std::vector<std::pair<size_t, size_t>> opts;
            for (auto [w, e] : adj[cur])
                if (!seen.contains(w)) opts.push_back({w, e});
            if (opts.empty()) break;
            auto [w, e] = opts[rng() % opts.size()];
            path.push_back(w);
            epath.push_back(e);
            seen.insert(w);
        }
    }
    return {};
}

void test_random_cycles_against_drilling() {
    const std::vector<std::pair<std::string, std::string>> rows = {
        {"3_1", "[[1;5;2;4];[3;1;4;6];[5;3;6;2]]"},
        {"6_1", "[[1;7;2;6];[3;10;4;11];[5;3;6;2];[7;1;8;12];[9;4;10;5];[11;9;12;8]]"},
    };
    std::mt19937 rng(20260926);
    int agree = 0, total = 0, throughWalls = 0;
    for (const auto &[name, pd] : rows) {
        auto [tri, edges, reversed] = knotbuilder::buildLink(knotbuilder::parsePDCode(pd));
        DiagramDrawer drawer(tri, parsePD(pd).size());
        for (int i = 0; i < 40; ++i) {
            auto cycles = randomCycle(tri, rng, 6 + rng() % 18);
            if (cycles.empty()) continue;
            Diagram d = drawer.draw(cycles);
            regina::Link drawn = d.link();
            drawn.simplify();
            bool drawnUnknot = drawn.size() == 0;

            std::vector<const regina::Edge<3> *> es;
            for (const DirectedEdge &de : cycles[0]) es.push_back(tri.edge(de.edge));
            regina::Triangulation<3> comp = Link(tri, es).buildComplement();
            bool drilledSolidTorus = comp.isSolidTorus();
            ++total;
            if (drawnUnknot == drilledSolidTorus) ++agree;
            // Curves through a wall, a corner line or a cone: an edge in
            // two or more blocks, or in none.
            const size_t nb = parsePD(pd).size();
            for (const DirectedEdge &de : cycles[0]) {
                std::set<size_t> blocksWith;
                for (const auto &emb : *tri.edge(de.edge))
                    if (emb.tetrahedron()->index() < 14 * nb)
                        blocksWith.insert(emb.tetrahedron()->index() / 14);
                if (blocksWith.size() != 1) { ++throughWalls; break; }
            }
        }
    }
    EXPECT_EQ(agree, total,
              "random curves: drawn unknot <=> drilled complement is a solid torus");
    EXPECT_EQ(total >= 60, true, "enough random curves were drawn");
    std::cout << "  (" << total << " random curves, " << throughWalls
              << " through a wall, corner line or cone point)\n";
}

// Every simple closed edge path of 2..maxLen edges, once each (from its least
// vertex, in one direction), directed along the walk.
std::vector<EdgeCycle> allShortCycles(const regina::Triangulation<3> &tri, size_t maxLen) {
    std::vector<std::vector<std::pair<size_t, size_t>>> adj(tri.countVertices());
    for (const regina::Edge<3> *e : tri.edges()) {
        size_t a = e->vertex(0)->index(), b = e->vertex(1)->index();
        if (a == b) continue;
        adj[a].push_back({b, e->index()});
        adj[b].push_back({a, e->index()});
    }
    std::vector<EdgeCycle> out;
    std::vector<size_t> path, epath;
    std::vector<bool> on(tri.countVertices(), false);
    auto directed = [&](size_t e, size_t from) {
        return DirectedEdge{e, tri.edge(e)->vertex(0)->index() != from};
    };
    std::function<void(size_t)> extend = [&](size_t s) {
        const size_t cur = path.back();
        for (auto [w, e] : adj[cur]) {
            if (w == s && !epath.empty() && (epath.size() > 1 || e != epath[0]) &&
                epath.front() < e) {
                EdgeCycle cyc;
                for (size_t k = 0; k < epath.size(); ++k) cyc.push_back(directed(epath[k], path[k]));
                cyc.push_back(directed(e, cur));
                out.push_back(std::move(cyc));
            }
            if (w > s && !on[w] && path.size() < maxLen) {
                path.push_back(w);
                epath.push_back(e);
                on[w] = true;
                extend(s);
                on[w] = false;
                path.pop_back();
                epath.pop_back();
            }
        }
    };
    for (size_t s = 0; s < tri.countVertices(); ++s) {
        path = {s};
        epath.clear();
        on[s] = true;
        extend(s);
        on[s] = false;
    }
    return out;
}

EdgeCycle reversedCycle(EdgeCycle c) {
    std::reverse(c.begin(), c.end());
    for (DirectedEdge &de : c) de.reversed = !de.reversed;
    return c;
}

// Exhaustive where the random curves above only sample: every short closed
// curve in T, drawn with every perturbation draw() may try and in both
// directions, must come out as a PLANAR diagram. A short curve stays near
// one wall or corner, so this covers every local configuration a curve of
// that length can make -- in particular a curve running up a vertical
// corner edge, whose tent projects to a loop through the corner point. A
// drawing that is not classical is a defect of the drawer, never of the
// curve (a closed curve in S^3 always has a planar diagram).
void test_every_short_cycle_draws_planar() {
    const std::vector<std::pair<std::string, std::string>> rows = {
        {"3_1", "[[1;5;2;4];[3;1;4;6];[5;3;6;2]]"},
        {"4_1", "[[4;2;5;1];[8;6;1;5];[6;3;7;4];[2;7;3;8]]"},
        {"5_2", "[[1;5;2;4];[3;9;4;8];[5;1;6;10];[7;3;8;2];[9;7;10;6]]"},
        {"L4a1{0}", "PD[X[6; 1; 7; 2]; X[8; 3; 5; 4]; X[2; 5; 3; 6]; X[4; 7; 1; 8]]"},
    };
    // Longer curves, or drilling every curve, for a manual deep run:
    //   SHORT_CYCLES_MAX_LEN=7 SHORT_CYCLES_DRILL_EVERY=1 ./diagramdrawer_test
    auto envOr = [](const char *var, size_t dflt) {
        const char *v = std::getenv(var);
        return v ? static_cast<size_t>(std::stoul(v)) : dflt;
    };
    const size_t maxLen = envOr("SHORT_CYCLES_MAX_LEN", 5);
    const size_t drillEvery = envOr("SHORT_CYCLES_DRILL_EVERY", 25);
    const auto start = std::chrono::steady_clock::now();
    long cycles = 0, drawn = 0, degenerate = 0, invalid = 0, nonPlanar = 0;
    long drilled = 0, drilledAgree = 0, drillMs = 0;
    for (const auto &[name, pd] : rows) {
        auto [tri, edges, reversed] = knotbuilder::buildLink(knotbuilder::parsePDCode(pd));
        DiagramDrawer drawer(tri, parsePD(pd).size());
        const std::vector<EdgeCycle> all = allShortCycles(tri, maxLen);
        std::cout << "  " << name << ": " << all.size() << " closed curves\n";
        for (const EdgeCycle &c : all) {
            ++cycles;
            bool firstDrawing = cycles % static_cast<long>(drillEvery) == 0;
            for (const EdgeCycle &dir : {c, reversedCycle(c)})
                for (unsigned seed = 1; seed <= 8; ++seed) {
                    try {
                        Diagram d = drawer.drawWithSeed({dir}, seed);
                        ++drawn;
                        regina::Link l = d.link();
                        if (!l.isClassical()) throw NonPlanar("not classical");
                        // The knot type, independently, once per curve: a
                        // drawing that is planar but wrong would show here.
                        if (firstDrawing) {
                            firstDrawing = false;
                            const auto t0 = std::chrono::steady_clock::now();
                            l.simplify();
                            std::vector<const regina::Edge<3> *> es;
                            for (const DirectedEdge &de : dir) es.push_back(tri.edge(de.edge));
                            ++drilled;
                            if ((l.size() == 0) == Link(tri, es).buildComplement().isSolidTorus())
                                ++drilledAgree;
                            drillMs += std::chrono::duration_cast<std::chrono::milliseconds>(
                                           std::chrono::steady_clock::now() - t0).count();
                        }
                    } catch (const NonPlanar &) {
                        if (++nonPlanar <= 5) {
                            std::cout << "  " << name << ": non-planar drawing, seed " << seed
                                      << ", edges";
                            for (const DirectedEdge &de : dir)
                                std::cout << ' ' << (de.reversed ? "-" : "+") << de.edge;
                            std::cout << "\n";
                        }
                    } catch (const Degenerate &) {
                        ++degenerate;
                    } catch (const regina::InvalidArgument &) {
                        ++invalid; // a curve through a cone point only
                    }
                }
        }
    }
    EXPECT_EQ(nonPlanar, 0L, "every short closed curve draws as a planar diagram");
    EXPECT_EQ(drilledAgree, drilled,
              "short curves: drawn unknot <=> drilled complement is a solid torus");
    EXPECT_EQ(cycles > 1000, true, "enough short curves were enumerated");
    std::cout << "  (" << cycles << " closed curves of <= " << maxLen << " edges; " << drawn
              << " drawings, " << degenerate << " degenerate perturbations, " << invalid
              << " invalid; " << drilled << " checked against drilling in " << drillMs
              << " ms; "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     std::chrono::steady_clock::now() - start).count()
              << " ms in all)\n";
}

// The one place two DIFFERENT curves meet across blocks: one through a
// region's bottom point, the other through its top point, which project to
// the same corner. Every such pair of short curves, drawn with every
// perturbation, must be planar and must agree -- HOMFLY polynomial and
// linking number -- across perturbations, across the order of the two
// components, and after reversing both (which fixes the HOMFLY polynomial of
// an oriented link). Different perturbations move tents and cone routes, so
// disagreement would expose a defect even where planarity cannot.
void test_two_curves_at_corners() {
    const std::vector<std::pair<std::string, std::string>> rows = {
        {"3_1", "[[1;5;2;4];[3;1;4;6];[5;3;6;2]]"},
        {"4_1", "[[4;2;5;1];[8;6;1;5];[6;3;7;4];[2;7;3;8]]"},
        {"L4a1{0}", "PD[X[6; 1; 7; 2]; X[8; 3; 5; 4]; X[2; 5; 3; 6]; X[4; 7; 1; 8]]"},
    };
    const size_t maxLen = 5, perCorner = 400;
    const auto start = std::chrono::steady_clock::now();
    long pairs = 0, drawings = 0, degenerate = 0, nonPlanar = 0, disagree = 0;
    for (const auto &[name, pd] : rows) {
        const size_t nb = parsePD(pd).size();
        auto [tri, edges, reversed] = knotbuilder::buildLink(knotbuilder::parsePDCode(pd));
        DiagramDrawer drawer(tri, nb);
        // Each region's (bottom, top) vertex pair, from the block model.
        std::set<std::pair<size_t, size_t>> corners;
        const auto &coords = blockCoordinates();
        for (size_t k = 0; k < nb; ++k) {
            std::map<std::pair<long, long>, std::pair<long, long>> at; // (x, y) -> (bottom, top)
            for (int j = 0; j < 14; ++j)
                for (int i = 0; i < 4; ++i) {
                    const BlockPoint &p = coords[j][i];
                    if (std::abs(p.x) != BLOCK_SCALE || std::abs(p.y) != BLOCK_SCALE) continue;
                    long v = static_cast<long>(tri.tetrahedron(14 * k + j)->vertex(i)->index());
                    auto &slot = at.try_emplace({p.x, p.y}, -1, -1).first->second;
                    (p.z == 0 ? slot.first : slot.second) = v;
                }
            for (const auto &[xy, bt] : at)
                if (bt.first >= 0 && bt.second >= 0)
                    corners.insert({static_cast<size_t>(bt.first), static_cast<size_t>(bt.second)});
        }
        const std::vector<EdgeCycle> all = allShortCycles(tri, maxLen);
        auto vertices = [&](const EdgeCycle &c) {
            std::set<size_t> vs;
            for (const DirectedEdge &de : c)
                for (int e = 0; e < 2; ++e) vs.insert(tri.edge(de.edge)->vertex(e)->index());
            return vs;
        };
        std::vector<std::set<size_t>> vsOf;
        for (const EdgeCycle &c : all) vsOf.push_back(vertices(c));
        for (auto [bottom, top] : corners) {
            std::vector<size_t> viaBottom, viaTop;
            for (size_t i = 0; i < all.size(); ++i) {
                bool b = vsOf[i].contains(bottom), t = vsOf[i].contains(top);
                if (b && !t) viaBottom.push_back(i);
                if (t && !b) viaTop.push_back(i);
            }
            std::vector<std::pair<size_t, size_t>> cand;
            for (size_t a : viaBottom)
                for (size_t b : viaTop) {
                    bool disjoint = std::none_of(vsOf[a].begin(), vsOf[a].end(),
                                                 [&](size_t v) { return vsOf[b].contains(v); });
                    if (disjoint) cand.push_back({a, b});
                }
            // A deterministic spread of at most perCorner pairs per corner.
            const size_t stride = std::max<size_t>(1, cand.size() / perCorner);
            for (size_t p = 0; p < cand.size(); p += stride) {
                const EdgeCycle &A = all[cand[p].first], &B = all[cand[p].second];
                ++pairs;
                std::optional<std::string> homfly;
                std::optional<long> lk;
                bool bad = false;
                auto check = [&](const std::vector<EdgeCycle> &curves, unsigned seed,
                                 long lkSign) {
                    try {
                        Diagram d = drawer.drawWithSeed(curves, seed);
                        ++drawings;
                        regina::Link rl = d.link();
                        // The link must keep every drawn orientation, even of
                        // a component passing over everything it meets.
                        if (rl.linking() != d.linkingNumber(0, 1)) bad = true;
                        std::string h = rl.homfly().str();
                        long l = lkSign * d.linkingNumber(0, 1);
                        if (!homfly) { homfly = h; lk = l; }
                        else if (*homfly != h || *lk != l) bad = true;
                    } catch (const NonPlanar &) {
                        ++nonPlanar;
                        bad = true;
                    } catch (const Degenerate &) {
                        ++degenerate;
                    }
                };
                for (unsigned seed = 1; seed <= 8; ++seed) check({A, B}, seed, 1);
                check({B, A}, 1, 1);
                check({reversedCycle(A), reversedCycle(B)}, 1, 1);
                if (bad && ++disagree <= 5)
                    std::cout << "  " << name << ": curves at corner (" << bottom << ", " << top
                              << ") draw inconsistently\n";
            }
        }
    }
    EXPECT_EQ(nonPlanar, 0L, "two curves at a corner always draw as a planar diagram");
    EXPECT_EQ(disagree, 0L,
              "two curves at a corner: every perturbation, component order and "
              "global reversal gives the same HOMFLY polynomial and linking number");
    EXPECT_EQ(pairs > 500, true, "enough pairs of curves at corners were drawn");
    std::cout << "  (" << pairs << " pairs of curves at corners, " << drawings << " drawings, "
              << degenerate << " degenerate perturbations; "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     std::chrono::steady_clock::now() - start).count()
              << " ms)\n";
}

void sweepTable(const std::string &path, int maxCrossings) {
    std::ifstream in(path);
    std::string line;
    std::getline(in, line);
    int rows = 0, same = 0;
    while (std::getline(in, line)) {
        std::string name = line.substr(0, line.find(','));
        std::string rest = line.substr(line.find(',') + 1);
        std::string pd = rest.substr(0, rest.rfind(','));
        int crossings = static_cast<int>(parsePD(pd).size());
        if (maxCrossings > 0 && crossings > maxCrossings) continue;
        ++rows;
        try {
            if (reproduce(name, pd, true)) ++same;
        } catch (const std::exception &e) {
            std::cout << "  " << name << ": " << e.what() << "\n";
        }
    }
    std::cout << same << "/" << rows << " rows of " << path
              << " reproduced exactly\n";
    passed += same;
    failed_count += rows - same;
}

} // namespace

int main(int argc, char **argv) {
    if (argc >= 2) {
        sweepTable(argv[1], argc >= 3 ? std::stoi(argv[2]) : 0);
    } else {
        test_block_model();
        test_reproduce_battery();
        test_reversed_component_is_caught();
        test_random_cycles_against_drilling();
        test_every_short_cycle_draws_planar();
        test_two_curves_at_corners();
    }
    std::cout << passed << " passed, " << failed_count << " failed\n";
    return failed_count == 0 ? 0 : 1;
}
