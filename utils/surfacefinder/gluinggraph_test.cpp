// gluinggraph_test.cpp
// Sanity checks for GluingGraph: structure (nodes/edges) and surface counts
// on small triangulations where the answer is known by hand.

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <triangulation/example3.h>
#include <triangulation/example4.h>
#include <iostream>

#include "gluinggraph.h"

static int passed = 0, failed_count = 0;

#define EXPECT_EQ(actual, expected, desc)                                  \
    do {                                                                    \
        auto _a = (actual);                                                 \
        auto _e = (expected);                                               \
        if (_a == _e) {                                                     \
            std::cout << "  PASS: " << (desc) << "\n";                     \
            ++passed;                                                       \
        } else {                                                            \
            std::cout << "  FAIL: " << (desc) << "\n"                      \
                      << "        expected " << _e << ", got " << _a << "\n";\
            ++failed_count;                                                 \
        }                                                                   \
    } while (0)

#define EXPECT_GE(actual, expected, desc)                                  \
    do {                                                                    \
        auto _a = (actual);                                                 \
        auto _e = (expected);                                               \
        if (_a >= _e) {                                                     \
            std::cout << "  PASS: " << (desc) << "\n";                     \
            ++passed;                                                       \
        } else {                                                            \
            std::cout << "  FAIL: " << (desc) << "\n"                      \
                      << "        expected >= " << _e << ", got " << _a << "\n";\
            ++failed_count;                                                 \
        }                                                                   \
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
    EXPECT_EQ(tri.isClosed(), false,   "triangulation has boundary");

    // Graph structure
    {
        GluingGraph<3> g(tri, SurfaceCondition::all);
        EXPECT_EQ(g.countNodes(), 4, "K4: 4 nodes");
        // K4: each of 4 nodes has 3 neighbours → 12 directed edges total
        EXPECT_EQ(g.countEdges(), 12, "K4: 12 directed adjacency entries");
    }

    // Surface counts
    {
        GluingGraph<3> g(tri, SurfaceCondition::all);
        auto& s = g.findSurfaces();
        EXPECT_EQ((int)s.size(), 15, "--all: 15 surfaces");
    }
    {
        GluingGraph<3> g(tri, SurfaceCondition::boundary);
        auto& s = g.findSurfaces();
        EXPECT_EQ((int)s.size(), 15, "--boundary: 15 proper surfaces");
    }
    {
        GluingGraph<3> g(tri, SurfaceCondition::closed);
        auto& s = g.findSurfaces();
        EXPECT_EQ((int)s.size(), 1, "--closed: 1 closed surface (S^2)");
    }
}

// ────────────────────────────────────────────────────────────────────
// Two tetrahedra sharing one face → solid that is still B^3.
//
// Triangles: 4 + 4 − 1 = 7 total (6 boundary, 1 interior).
// The interior triangle is improper: its three edges are interior manifold
// edges, so any surface whose boundary includes one of them fails isProper().
//
// Expected:
//   --all      : ≥ 7   (at minimum one surface per triangle)
//   --boundary : ≥ 6   (at minimum one disc per boundary triangle)
//   --closed   : ≥ 1   (the boundary S^2 built from 6 triangles)
// ────────────────────────────────────────────────────────────────────
void test_two_tetrahedra() {
    std::cout << "\n--- two tetrahedra sharing one face ---\n";

    regina::Triangulation<3> tri;
    auto* t0 = tri.newSimplex();
    auto* t1 = tri.newSimplex();
    t0->join(3, t1, regina::Perm<4>());

    EXPECT_EQ(tri.countTriangles(), 7, "7 distinct triangles");

    {
        GluingGraph<3> g(tri, SurfaceCondition::all);
        EXPECT_EQ(g.countNodes(), 7, "7 nodes in gluing graph");
        auto& s = g.findSurfaces();
        EXPECT_GE((int)s.size(), 7, "--all: ≥ 7 surfaces");
    }
    {
        GluingGraph<3> g(tri, SurfaceCondition::boundary);
        auto& s = g.findSurfaces();
        EXPECT_GE((int)s.size(), 6, "--boundary: ≥ 6 proper surfaces");
    }
    {
        GluingGraph<3> g(tri, SurfaceCondition::closed);
        auto& s = g.findSurfaces();
        EXPECT_GE((int)s.size(), 1, "--closed: ≥ 1 closed surface");
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
            GluingGraph<3> g(tri, SurfaceCondition::closed);
            closed_count = g.findSurfaces().size();
            std::cout << "    --closed:   " << closed_count << "\n";
        }
        {
            GluingGraph<3> g(tri, SurfaceCondition::boundary);
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
            GluingGraph<3> g(tri, SurfaceCondition::closed);
            closed_count = g.findSurfaces().size();
            std::cout << "    --closed:   " << closed_count << "\n";
        }
        {
            GluingGraph<3> g(tri, SurfaceCondition::boundary);
            boundary_count = g.findSurfaces().size();
            std::cout << "    --boundary: " << boundary_count << "\n";
        }
        EXPECT_EQ(closed_count, boundary_count,
                  "--closed == --boundary in a closed manifold");
        EXPECT_GE((int)closed_count, 1,
                  "∂Δ^4 triangulation has at least one embedded S^2");
    }
}

// ────────────────────────────────────────────────────────────────────
// S^4 (closed 4-manifold).  Structural smoke-test: graph builds
// without error and findSurfaces() terminates.
// ────────────────────────────────────────────────────────────────────
void test_four_sphere() {
    std::cout << "\n--- S^4 (closed 4-manifold) ---\n";

    auto tri = regina::Example<4>::fourSphere();
    EXPECT_EQ(tri.isClosed(), true, "S^4 is closed");
    std::cout << "    " << tri.size() << " pentachora, "
              << tri.countTriangles() << " triangles\n";

    GluingGraph<4> g(tri, SurfaceCondition::closed);
    std::cout << "    Gluing graph: " << g.countNodes() << " nodes, "
              << g.countEdges() << " edges\n";

    auto& s = g.findSurfaces();
    std::cout << "    Closed surfaces found: " << s.size() << "\n";
    EXPECT_GE((int)s.size(), 0, "findSurfaces() completes without error");
}

int main() {
    test_single_tetrahedron();
    test_two_tetrahedra();
    test_three_sphere();
    test_four_sphere();

    std::cout << "\n=== " << passed << " passed, " << failed_count
              << " failed ===\n";
    return failed_count > 0 ? 1 : 0;
}
