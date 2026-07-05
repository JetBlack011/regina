// surfacefinder_test.cpp
// Sanity checks for SurfaceFinder: structure (nodes/edges) and surface counts
// on small triangulations where the answer is known by hand.

#include <iostream>
#include <sstream>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <triangulation/example3.h>
#include <triangulation/example4.h>
#include <unistd.h>
#include <unordered_set>

#include "cobordismbuilder.h"
#include "surfacefinder.h"

static int passed = 0, failed_count = 0;

namespace {
bool colorEnabled() {
    static bool enabled = isatty(fileno(stdout));
    return enabled;
}
std::ostream &green(std::ostream &os) {
    return colorEnabled() ? os << "\033[32m" : os;
}
std::ostream &red(std::ostream &os) {
    return colorEnabled() ? os << "\033[31m" : os;
}
std::ostream &bold(std::ostream &os) {
    return colorEnabled() ? os << "\033[1m" : os;
}
std::ostream &resetColor(std::ostream &os) {
    return colorEnabled() ? os << "\033[0m" : os;
}
} // namespace

#define EXPECT_EQ(actual, expected, desc)                                      \
    do {                                                                       \
        auto _a = (actual);                                                    \
        auto _e = (expected);                                                  \
        if (_a == _e) {                                                        \
            std::cout << green << "  PASS: " << resetColor << (desc) << "\n";  \
            ++passed;                                                          \
        } else {                                                               \
            std::cout << red << "  FAIL: " << (desc) << "\n"                   \
                      << "        expected " << _e << ", got " << _a           \
                      << resetColor << "\n";                                   \
            ++failed_count;                                                    \
        }                                                                      \
    } while (0)

#define EXPECT_GE(actual, expected, desc)                                      \
    do {                                                                       \
        auto _a = (actual);                                                    \
        auto _e = (expected);                                                  \
        if (_a >= _e) {                                                        \
            std::cout << green << "  PASS: " << resetColor << (desc) << "\n";  \
            ++passed;                                                          \
        } else {                                                               \
            std::cout << red << "  FAIL: " << (desc) << "\n"                   \
                      << "        expected >= " << _e << ", got " << _a        \
                      << resetColor << "\n";                                   \
            ++failed_count;                                                    \
        }                                                                      \
    } while (0)

// ────────────────────────────────────────────────────────────────────
// Single tetrahedron = B^3.
//
// The four boundary triangles form a gluing graph K_4 (every pair shares
// exactly one edge).  Every non-empty subset of the four triangles is a
// valid embedded proper surface:
//
//   1-triangle subsets (4): single triangle = disc
//   2-triangle subsets (6): two triangles glued along one edge = disc
//   3-triangle subsets (4): "cap" of the tetrahedron = disc
//   4-triangle subset  (1): full boundary S^2 = sphere
//
// Expected totals:
//   --all      : 15
//   --boundary : 15  (all are proper: every edge lies on ∂B^3 = S^2)
//   --closed   :  1  (only the full S^2)
// ────────────────────────────────────────────────────────────────────
void test_single_tetrahedron() {
    std::cout << "\n--- single tetrahedron (B^3) ---\n";

    regina::Triangulation<3> tri;
    tri.newSimplex();

    EXPECT_EQ(tri.countTriangles(), 4, "triangulation has 4 triangles");
    EXPECT_EQ(tri.isClosed(), false, "triangulation has boundary");

    // Graph structure
    {
        SurfaceFinder<3> g(tri, SurfaceCondition::all);
        EXPECT_EQ(g.countNodes(), 4, "K4: 4 nodes");
        // K4: each of 4 nodes has 3 neighbours → 12 directed edges total
        EXPECT_EQ(g.countEdges(), 12, "K4: 12 directed adjacency entries");
    }

    // Surface counts
    {
        SurfaceFinder<3> g(tri, SurfaceCondition::all);
        auto &s = g.findSurfaces();
        EXPECT_EQ((int)s.size(), 15, "--all: 15 surfaces");
    }
    {
        SurfaceFinder<3> g(tri, SurfaceCondition::boundary);
        auto &s = g.findSurfaces();
        EXPECT_EQ((int)s.size(), 15, "--boundary: 15 proper surfaces");
    }
    {
        SurfaceFinder<3> g(tri, SurfaceCondition::closed);
        auto &s = g.findSurfaces();
        EXPECT_EQ((int)s.size(), 1, "--closed: 1 closed surface (S^2)");
    }
}

// ────────────────────────────────────────────────────────────────────
// findSurfaces(startingTriangles): search with a fixed base of triangles
// that must already be part of the surface (e.g. an annulus whose boundary
// is a known knot/link -- the intended use case). On the tetrahedron's K4
// graph, every subset is connected, so "surfaces containing a fixed base
// set B" is exactly "B plus any subset of the remaining triangles":
//   |B| = 1 (say {n0}): 2^3 = 8 supersets  (1 + 3 + 3 + 1, by size)
//   |B| = 2 (say {n0,n1}): 2^2 = 4 supersets ({n0,n1}, plus each of the
//     other two singly, plus both)
// All of these are proper (as in test_single_tetrahedron), so --boundary
// matches --all; only the full 4-triangle surface is closed.
//
// This also regression-tests the bug this overload used to have: it only
// ever added each starting triangle as an independent single-triangle seed
// and blocked the others, so it could never find a surface actually using
// two or more starting triangles together -- exactly the annulus use case.
// ────────────────────────────────────────────────────────────────────
void test_starting_triangles_overload() {
    std::cout << "\n--- findSurfaces(startingTriangles) on B^3 (K4) ---\n";

    regina::Triangulation<3> tri;
    tri.newSimplex();

    regina::Triangle<3> *t0 = tri.triangle(0);
    regina::Triangle<3> *t1 = tri.triangle(1);

    // Throws on a triangle that isn't part of this triangulation at all.
    {
        regina::Triangulation<3> other;
        other.newSimplex();
        std::unordered_set<regina::Triangle<3> *> foreign{other.triangle(0)};

        SurfaceFinder<3> g(tri, SurfaceCondition::all);
        bool threw = false;
        try {
            g.findSurfaces(foreign);
        } catch (const regina::InvalidArgument &) {
            threw = true;
        }
        EXPECT_EQ(threw, true,
                  "throws on a starting triangle outside the graph");
    }

    // Single-triangle base: every superset of {t0} within K4.
    {
        std::unordered_set<regina::Triangle<3> *> base{t0};

        SurfaceFinder<3> g(tri, SurfaceCondition::all);
        EXPECT_EQ((int)g.findSurfaces(base).size(), 8,
                  "base {t0}, --all: 8 surfaces contain t0");
    }
    {
        std::unordered_set<regina::Triangle<3> *> base{t0};
        SurfaceFinder<3> g(tri, SurfaceCondition::boundary);
        EXPECT_EQ((int)g.findSurfaces(base).size(), 8,
                  "base {t0}, --boundary: 8 proper surfaces contain t0");
    }
    {
        std::unordered_set<regina::Triangle<3> *> base{t0};
        SurfaceFinder<3> g(tri, SurfaceCondition::closed);
        EXPECT_EQ((int)g.findSurfaces(base).size(), 1,
                  "base {t0}, --closed: only the full sphere contains t0");
    }

    // Two-triangle base: every superset of {t0, t1} -- the case the old
    // implementation could never find, since it required both to be
    // present simultaneously rather than as alternative independent seeds.
    {
        std::unordered_set<regina::Triangle<3> *> base{t0, t1};
        SurfaceFinder<3> g(tri, SurfaceCondition::all);
        EXPECT_EQ((int)g.findSurfaces(base).size(), 4,
                  "base {t0,t1}, --all: 4 surfaces contain both");
    }
    {
        std::unordered_set<regina::Triangle<3> *> base{t0, t1};
        SurfaceFinder<3> g(tri, SurfaceCondition::closed);
        EXPECT_EQ((int)g.findSurfaces(base).size(), 1,
                  "base {t0,t1}, --closed: only the full sphere contains "
                  "both");
    }
}

// ────────────────────────────────────────────────────────────────────
// findSurfaces(startingTriangles) must reject a base that isn't already a
// single connected partial surface -- using two tetrahedra sharing a face,
// where a_i is adjacent to b_i (same index, via the shared face's edges)
// but *not* to b_j for j != i, so {a0, b1} shares no edge at all.
// ────────────────────────────────────────────────────────────────────
void test_starting_triangles_disconnected_base() {
    std::cout << "\n--- findSurfaces(startingTriangles) rejects a "
                 "disconnected base ---\n";

    regina::Triangulation<3> tri;
    auto *t0 = tri.newSimplex();
    auto *t1 = tri.newSimplex();
    t0->join(3, t1, regina::Perm<4>());

    // a0 = t0's face opposite vertex 0; b1 = t1's face opposite vertex 1.
    regina::Triangle<3> *a0 = t0->triangle(0);
    regina::Triangle<3> *b1 = t1->triangle(1);
    std::unordered_set<regina::Triangle<3> *> base{a0, b1};

    SurfaceFinder<3> g(tri, SurfaceCondition::all);
    bool threw = false;
    try {
        g.findSurfaces(base);
    } catch (const regina::InvalidArgument &) {
        threw = true;
    }
    EXPECT_EQ(threw, true, "throws on a disconnected starting set");
}

// ────────────────────────────────────────────────────────────────────
// surfaces() accessor and operator<<: exercised together on the
// tetrahedron, since both are thin wrappers around already-tested state.
// ────────────────────────────────────────────────────────────────────
void test_accessor_and_stream_output() {
    std::cout << "\n--- surfaces() accessor and operator<< ---\n";

    regina::Triangulation<3> tri;
    tri.newSimplex();

    SurfaceFinder<3> g(tri, SurfaceCondition::all);
    auto &found = g.findSurfaces();
    EXPECT_EQ(&g.surfaces(), &found,
              "surfaces() returns the same set findSurfaces() returned");
    EXPECT_EQ((int)g.surfaces().size(), 15,
              "surfaces() reflects the most recent findSurfaces() call");

    std::ostringstream out;
    out << g;
    std::string s = out.str();
    EXPECT_EQ(s.find("SurfaceFinder<3>") != std::string::npos, true,
              "operator<< names the class and dimension");
    EXPECT_EQ(s.find("nodes = 4") != std::string::npos, true,
              "operator<< reports the node count");
    EXPECT_EQ(s.find("edges = 12") != std::string::npos, true,
              "operator<< reports the edge count");
}

// ────────────────────────────────────────────────────────────────────
// Two tetrahedra sharing one face → solid that is still B^3.
//
// Triangles: 4 + 4 − 1 = 7 total (6 boundary, 1 interior, call it M).
// The interior triangle is improper: its three edges are interior manifold
// edges, so any surface whose boundary includes one of them fails isProper().
//
// The gluing graph here is *not* simply "two K4 cliques sharing the vertex
// M": each of M's three edges is an ambient edge of degree 3 (touched by M
// and by one triangle from *each* tetrahedron), so the corresponding
// triangle from tetrahedron A is also directly adjacent to the
// corresponding triangle from tetrahedron B, not just to M. This makes exact
// hand-enumeration error-prone, so the expected counts below are pinned to
// what the (independently-validated-on-K4) exhaustive search reports,
// rather than derived by hand:
//   --all      : 72
//   --boundary : 72  (every embeddable subset here turns out proper)
//   --closed   : 3   (the boundary S^2 of each tetrahedron alone, plus one
//                      further closed surface mixing triangles from both)
// This test previously only checked loose lower bounds (>= 7 / >= 6 / >= 1)
// because the old triangle-chain DFS could not find surfaces requiring a
// triangle to have two or more simultaneously-attached branches -- see
// SurfaceFinder::extend_'s comment for why.
// ────────────────────────────────────────────────────────────────────
void test_two_tetrahedra() {
    std::cout << "\n--- two tetrahedra sharing one face ---\n";

    regina::Triangulation<3> tri;
    auto *t0 = tri.newSimplex();
    auto *t1 = tri.newSimplex();
    t0->join(3, t1, regina::Perm<4>());

    EXPECT_EQ(tri.countTriangles(), 7, "7 distinct triangles");

    {
        SurfaceFinder<3> g(tri, SurfaceCondition::all);
        EXPECT_EQ(g.countNodes(), 7, "7 nodes in gluing graph");
        auto &s = g.findSurfaces();
        EXPECT_EQ((int)s.size(), 72, "--all: 72 surfaces");
    }
    {
        SurfaceFinder<3> g(tri, SurfaceCondition::boundary);
        auto &s = g.findSurfaces();
        EXPECT_EQ((int)s.size(), 72, "--boundary: 72 proper surfaces");
    }
    {
        SurfaceFinder<3> g(tri, SurfaceCondition::closed);
        auto &s = g.findSurfaces();
        EXPECT_EQ((int)s.size(), 3, "--closed: 3 closed surfaces");
    }
}

// ────────────────────────────────────────────────────────────────────
// S^3 (closed 3-manifold).
//
// In a closed manifold isProper() ≡ isClosed(), so --boundary and
// --closed must agree.
//
// We test two triangulations:
//   (a) Regina's minimal 1-tetrahedron S^3 – very coarse; the 2
//       distinct triangles are likely too few / self-intersecting for
//       any embedded closed surface to appear.
//   (b) The boundary of a 4-simplex: 5 tetrahedra, 10 triangles,
//       guaranteed to contain at least one S^2.
// ────────────────────────────────────────────────────────────────────
void test_three_sphere() {
    std::cout << "\n--- S^3 (closed 3-manifold, minimal) ---\n";

    {
        auto tri = regina::Example<3>::threeSphere();
        EXPECT_EQ(tri.isClosed(), true, "minimal S^3 is closed");
        std::cout << "    " << tri.size() << " tetrahedra, "
                  << tri.countTriangles() << " triangles\n";

        size_t closed_count, boundary_count;
        {
            SurfaceFinder<3> g(tri, SurfaceCondition::closed);
            closed_count = g.findSurfaces().size();
            std::cout << "    --closed:   " << closed_count << "\n";
        }
        {
            SurfaceFinder<3> g(tri, SurfaceCondition::boundary);
            boundary_count = g.findSurfaces().size();
            std::cout << "    --boundary: " << boundary_count << "\n";
        }
        EXPECT_EQ(closed_count, boundary_count,
                  "--closed == --boundary in a closed manifold");
    }

    std::cout << "\n--- S^3 (boundary of 4-simplex, 5 tetrahedra) ---\n";
    {
        // Build ∂Δ^4: one 4-simplex whose boundary gives a 5-tetrahedra
        // triangulation of S^3.  This is fine enough to contain an S^2.
        regina::Triangulation<4> fourBall;
        fourBall.newSimplex();
        auto tri = fourBall.boundaryComponent(0)->build();

        EXPECT_EQ(tri.isClosed(), true, "∂Δ^4 is closed");
        std::cout << "    " << tri.size() << " tetrahedra, "
                  << tri.countTriangles() << " triangles\n";

        size_t closed_count, boundary_count;
        {
            SurfaceFinder<3> g(tri, SurfaceCondition::closed);
            closed_count = g.findSurfaces().size();
            std::cout << "    --closed:   " << closed_count << "\n";
        }
        {
            SurfaceFinder<3> g(tri, SurfaceCondition::boundary);
            boundary_count = g.findSurfaces().size();
            std::cout << "    --boundary: " << boundary_count << "\n";
        }
        EXPECT_EQ(closed_count, boundary_count,
                  "--closed == --boundary in a closed manifold");
        // 15 = 5 tetrahedron-boundary spheres (one per facet of Δ^4) + 10
        // further closed surfaces mixing triangles across facets, confirmed
        // by the exhaustive search rather than hand-enumerated.
        EXPECT_EQ((int)closed_count, 15,
                  "∂Δ^4 triangulation has 15 embedded closed surfaces");
    }
}

// ────────────────────────────────────────────────────────────────────
// Bounded pipeline smoke test: CobordismBuilder<3>::cone() -> SurfaceFinder<4>
// -> KnottedSurface<4>::boundary(), on a small hand-built (not knotbuilder)
// S^3 with a known unknotted loop, checked end to end.
//
// A full *unconstrained* findSurfaces() on real knotbuilder output is not
// tractable: even the smallest realistic case (Hopf link PD code, cone()
// with zero thickening layers) still had 10,000+ surfaces after 2 million
// search calls with no end in sight. The fix for that is exactly the
// seeded search exercised here: give findSurfaces(startingTriangles) the
// triangles that already trace out the known loop, so the search only has
// to find the cap, not explore the whole graph. Wiring that up against
// knotbuilder's *actual* Block-gadget output -- i.e. identifying which
// triangles trace a specific knot/link edge's sweep through thicken()/
// cone() -- is real infrastructure that doesn't exist yet; this test
// covers the same data flow (Link extraction, CobordismBuilder capping,
// seeded SurfaceFinder search, boundary Link comparison) on a manifold
// small enough to make the point without it.
// ────────────────────────────────────────────────────────────────────
void test_cone_pipeline_finds_disc_bounding_known_unknot() {
    std::cout << "\n--- cone() pipeline: seeded search finds discs "
                 "bounding a known unknot ---\n";

    // ∂Δ^4: S^3, 5 tetrahedra. tri.triangle(0)'s own 3 edges form an
    // unknotted loop (the boundary of a single 2-simplex is trivially
    // unknotted) -- our "known link" to search for.
    regina::Triangulation<4> fourBall;
    fourBall.newSimplex();
    auto tri = fourBall.boundaryComponent(0)->build();
    regina::Triangle<3> *seed3 = tri.triangle(0);

    CobordismBuilder<3> cob(tri);
    auto &coned = cob.cone();

    // cone() builds one new pentachoron per tetrahedron of tri, in the same
    // order, preserving local vertex labels 0..3 (vertex 4 is the new
    // apex) -- so the boundary component left behind should reconstruct
    // *identically* (not just isomorphically) to the original tri.
    auto rebuiltBoundary = coned.boundaryComponent(0)->build();
    EXPECT_EQ(rebuiltBoundary.isoSig(), tri.isoSig(),
              "cone()'s remaining boundary matches the original S^3 exactly, "
              "not just up to isomorphism");

    // Find seed3's counterpart in `coned`: the pentachoron built from
    // seed3's own tetrahedron, then whichever of its 10 triangles uses
    // seed3's 3 local vertices and neither the excluded local vertex nor
    // the apex (vertex 4).
    regina::Simplex<4> *conePent =
        coned.simplex(seed3->front().simplex()->index());
    int excludedVertex = seed3->front().face();
    regina::Triangle<4> *seed4 = nullptr;
    for (int t = 0; t < 10 && !seed4; ++t) {
        regina::Perm<5> vmap = conePent->triangleMapping(t);
        bool usesApexOrExcluded = false;
        for (int i = 0; i < 3; ++i) {
            if (vmap[i] == 4 || vmap[i] == excludedVertex)
                usesApexOrExcluded = true;
        }
        if (!usesApexOrExcluded)
            seed4 = conePent->triangle(t);
    }
    EXPECT_EQ(seed4 != nullptr, true,
              "found seed3's counterpart triangle in coned");

    std::unordered_set<regina::Triangle<4> *> base{seed4};
    SurfaceFinder<4> g(coned, SurfaceCondition::boundary);
    auto &surfaces = g.findSurfaces(base);

    EXPECT_GE((int)surfaces.size(), 1,
              "seeded search finds at least one surface containing seed4");

    bool foundBoundingDisc = false;
    for (const auto &surf : surfaces) {
        if (surf.isClosed())
            continue;
        Link link = surf.boundary();
        if (link.countComponents() != 1)
            continue;
        // A 1-component boundary here is necessarily seed3's own 3-edge
        // loop: every surface in `surfaces` contains seed4, whose only
        // proper (non-improper) edges, if any survive to the final
        // boundary, are exactly those 3.
        foundBoundingDisc = true;
        break;
    }
    EXPECT_EQ(foundBoundingDisc, true,
              "found a proper surface whose boundary is a single-component "
              "(unknotted) link");
}

// ────────────────────────────────────────────────────────────────────
// S^4 (closed 4-manifold, built from 2 pentachora). Mainly a structural
// smoke-test (graph builds without error, findSurfaces() terminates), but
// the exact closed-surface count is pinned as a regression check -- not
// hand-derived (the 2-pentachoron gluing graph is too dense to enumerate by
// hand), just confirmed by the exhaustive search.
// ────────────────────────────────────────────────────────────────────
void test_four_sphere() {
    std::cout << "\n--- S^4 (closed 4-manifold) ---\n";

    auto tri = regina::Example<4>::fourSphere();
    EXPECT_EQ(tri.isClosed(), true, "S^4 is closed");
    std::cout << "    " << tri.size() << " pentachora, " << tri.countTriangles()
              << " triangles\n";

    SurfaceFinder<4> g(tri, SurfaceCondition::closed);
    std::cout << "    Gluing graph: " << g.countNodes() << " nodes, "
              << g.countEdges() << " edges\n";

    auto &s = g.findSurfaces();
    std::cout << "    Closed surfaces found: " << s.size() << "\n";
    EXPECT_EQ((int)s.size(), 15, "S^4 (2 pentachora) has 15 closed surfaces");
}

int main() {
    test_single_tetrahedron();
    test_starting_triangles_overload();
    test_starting_triangles_disconnected_base();
    test_accessor_and_stream_output();
    test_two_tetrahedra();
    test_three_sphere();
    test_cone_pipeline_finds_disc_bounding_known_unknot();
    test_four_sphere();

    std::cout << "\n"
              << bold << (failed_count > 0 ? red : green) << "=== " << passed
              << " passed, " << failed_count << " failed ===" << resetColor
              << "\n";
    return failed_count > 0 ? 1 : 0;
}
