// surfacefinder_test.cpp
// Sanity checks for SurfaceFinder: structure (nodes/edges) and surface counts
// on small triangulations where the answer is known by hand.

#include <iostream>
#include <sstream>
#include <string>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <triangulation/example3.h>
#include <triangulation/example4.h>
#include <unistd.h>
#include <unordered_set>
#include <vector>

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
// KnottedSurface::removeTriangle must clear improperEdges_ entries that
// belonged solely to the triangle being removed.
//
// Bug: addTriangle() flags a triangle's own un-glued facet as improper
// whenever the underlying ambient edge isn't on the manifold's boundary.
// removeTriangle() previously never undid this -- only the "this facet
// just became un-glued because its neighbour left" branch touched
// improperEdges_, not the "this facet was already un-glued and its
// owning triangle is now gone" branch. So once *any* triangle with a
// genuinely-interior un-glued edge was added and later backtracked out
// (an everyday occurrence in the DFS), isProper() would incorrectly and
// permanently report false for the rest of that surface_'s lifetime,
// silently dropping otherwise-valid surfaces under --boundary/--links.
//
// Tests this directly against KnottedSurface, bypassing SurfaceFinder's
// DFS, so the add/remove sequence is fully controlled: a boundary
// triangle (Y, from a free-standing, unglued tetrahedron -- every edge
// of every one of its faces is on the boundary) is added and kept
// throughout, so the surface never becomes trivially "closed" (which
// would mask the bug via isProper()'s isClosed() short-circuit); an
// interior triangle (X, taken from a separate, closed ∂Delta^4 S^3
// component -- every edge there is non-boundary, since the component
// itself has no boundary at all) is then added alone -- with an empty
// adjacency list, so none of its facets get glued -- making it improper,
// and then removed. isProper() must return to true afterward.
//
// (The minimal 1-tetrahedron S^3 used elsewhere in this file isn't
// suitable here: its own two triangles already have coincident local
// vertices, so adding either one alone is rejected outright as
// self-intersecting -- see test_three_sphere's comment. ∂Delta^4 is fine
// enough that its triangles are individually well-behaved.)
// ────────────────────────────────────────────────────────────────────
void test_removeTriangle_clears_improper_edges() {
    std::cout
        << "\n--- KnottedSurface::removeTriangle clears improperEdges_ ---\n";

    // Component A: the closed S^3 boundary of a single 4-simplex (5
    // tetrahedra) -- every edge is non-boundary, since the component has
    // no boundary, and its triangles are well-behaved (no coincident
    // local vertices).
    regina::Triangulation<4> fourBall;
    fourBall.newSimplex();
    regina::Triangulation<3> tri = fourBall.boundaryComponent(0)->build();
    // Component B: a second, free-standing tetrahedron, left completely
    // unglued -- every edge of every one of its faces is on the boundary.
    tri.newSimplex();

    regina::Triangle<3> *interior = nullptr;
    regina::Triangle<3> *bdryTriangle = nullptr;
    for (regina::Triangle<3> *t : tri.triangles()) {
        bool anyNonBoundary = false;
        bool allBoundary = true;
        for (int i = 0; i < 3; ++i) {
            if (t->edge(i)->isBoundary())
                ; // still eligible as a boundary edge
            else {
                anyNonBoundary = true;
                allBoundary = false;
            }
        }
        if (anyNonBoundary && interior == nullptr)
            interior = t;
        if (allBoundary && bdryTriangle == nullptr)
            bdryTriangle = t;
    }

    if (interior == nullptr || bdryTriangle == nullptr) {
        std::cout << red
                   << "  FAIL: couldn't locate suitable triangles to "
                      "set up this test\n"
                   << resetColor;
        ++failed_count;
        return;
    }

    KnottedSurface<3> ks(&tri);
    GluingNode<3>::AdjList empty;

    EXPECT_EQ(ks.addTriangle(bdryTriangle, empty), true,
              "boundary triangle adds with no neighbours present");
    EXPECT_EQ(ks.isProper(), true,
              "surface containing only the boundary triangle is proper");
    EXPECT_EQ(ks.isClosed(), false,
              "surface containing only the boundary triangle isn't closed");

    EXPECT_EQ(ks.addTriangle(interior, empty), true,
              "interior triangle adds with no neighbours present");
    EXPECT_EQ(
        ks.isProper(), false,
        "surface is improper while the interior triangle's edges are unglued");

    ks.removeTriangle(interior);

    EXPECT_EQ(ks.isClosed(), false,
              "boundary triangle is still present after removing the "
              "interior one");
    EXPECT_EQ(ks.isProper(), true,
              "removing the interior triangle must clear its improper "
              "edges (bug #1)");
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

// ────────────────────────────────────────────────────────────────────
// Every surface findSurfaces() returns must be connected and embedded (an
// injective vertex map into the ambient triangulation), and must actually
// satisfy the SurfaceCondition it was searched under -- checked here from
// scratch using regina's own primitives (isConnected(),
// countBoundaryFacets(), boundary-edge images), never KnottedSurface's own
// incremental flags, so a bug in the incremental bookkeeping can't hide
// from the very check meant to catch it.
// ────────────────────────────────────────────────────────────────────
void checkFoundSurfacesSatisfyInvariants(const std::string &label,
                                         regina::Triangulation<3> &tri,
                                         SurfaceCondition cond) {
    SurfaceFinder<3> g(tri, cond);
    auto &surfaces = g.findSurfaces();

    bool allConnected = true;
    bool allEmbedded = true;
    bool conditionHolds = true;
    bool closedFlagConsistent = true;

    for (const auto &s : surfaces) {
        if (!s.surface().isConnected())
            allConnected = false;

        std::unordered_set<const regina::Vertex<3> *> seen;
        for (regina::Vertex<2> *v : s.surface().vertices()) {
            if (!seen.insert(s.image(v)).second)
                allEmbedded = false;
        }

        if (s.isClosed() != s.surface().isClosed())
            closedFlagConsistent = false;

        if (cond == SurfaceCondition::closed) {
            if (s.surface().countBoundaryFacets() != 0)
                conditionHolds = false;
        } else if (cond == SurfaceCondition::boundary) {
            for (const regina::BoundaryComponent<2> *comp :
                s.surface().boundaryComponents()) {
                for (const regina::Edge<2> *edge : comp->edges()) {
                    if (edge->isBoundary() && !s.image(edge)->isBoundary())
                        conditionHolds = false;
                }
            }
        }
    }

    EXPECT_GE((int)surfaces.size(), 1, label + ": found at least one surface");
    EXPECT_EQ(allConnected, true, label + ": every found surface is connected");
    EXPECT_EQ(allEmbedded, true,
             label + ": every found surface has an injective vertex map");
    EXPECT_EQ(closedFlagConsistent, true,
             label + ": isClosed() agrees with surface().isClosed()");
    if (cond != SurfaceCondition::all)
        EXPECT_EQ(conditionHolds, true,
                 label +
                     ": every found surface satisfies its search condition");
}

void test_found_surfaces_satisfy_invariants() {
    std::cout << "\n--- found surfaces satisfy their invariants ---\n";

    {
        regina::Triangulation<3> tri;
        tri.newSimplex();
        checkFoundSurfacesSatisfyInvariants("single tetrahedron --all", tri,
                                            SurfaceCondition::all);
        checkFoundSurfacesSatisfyInvariants("single tetrahedron --boundary",
                                            tri, SurfaceCondition::boundary);
        checkFoundSurfacesSatisfyInvariants("single tetrahedron --closed", tri,
                                            SurfaceCondition::closed);
    }
    {
        regina::Triangulation<3> tri;
        auto *t0 = tri.newSimplex();
        auto *t1 = tri.newSimplex();
        t0->join(3, t1, regina::Perm<4>());
        checkFoundSurfacesSatisfyInvariants("two tetrahedra --all", tri,
                                            SurfaceCondition::all);
        checkFoundSurfacesSatisfyInvariants("two tetrahedra --boundary", tri,
                                            SurfaceCondition::boundary);
        checkFoundSurfacesSatisfyInvariants("two tetrahedra --closed", tri,
                                            SurfaceCondition::closed);
    }
    {
        regina::Triangulation<4> fourBall;
        fourBall.newSimplex();
        regina::Triangulation<3> tri = fourBall.boundaryComponent(0)->build();
        checkFoundSurfacesSatisfyInvariants("boundary of Delta^4 (S^3) --closed",
                                            tri, SurfaceCondition::closed);
    }
}

// ────────────────────────────────────────────────────────────────────
// KnottedSurface equality/ordering is defined purely by the *set* of
// triangle indices (see the class-level comment on operator==), so
// building the same two triangles in a different order must still
// compare equal.
// ────────────────────────────────────────────────────────────────────
void test_surface_equality_order_independent() {
    std::cout << "\n--- KnottedSurface equality is order-independent ---\n";

    regina::Triangulation<3> tri;
    tri.newSimplex();
    SurfaceFinder<3> g(tri, SurfaceCondition::all);
    regina::Triangle<3> *t0 = tri.triangle(0);
    regina::Triangle<3> *t1 = tri.triangle(1);

    KnottedSurface<3> ksAB(&tri);
    EXPECT_EQ(ksAB.addTriangle(t0, GluingNode<3>::AdjList{}), true,
             "t0 adds alone (A-then-B order)");
    EXPECT_EQ(ksAB.addTriangle(t1, g.adjacencyOf(t1)), true,
             "t1 adds and glues to t0 (A-then-B order)");

    KnottedSurface<3> ksBA(&tri);
    EXPECT_EQ(ksBA.addTriangle(t1, GluingNode<3>::AdjList{}), true,
             "t1 adds alone (B-then-A order)");
    EXPECT_EQ(ksBA.addTriangle(t0, g.adjacencyOf(t0)), true,
             "t0 adds and glues to t1 (B-then-A order)");

    EXPECT_EQ(ksAB == ksBA, true,
             "the same two triangles added in different orders compare equal");
    EXPECT_EQ((ksAB < ksBA) || (ksBA < ksAB), false,
             "neither ordering is considered less than the other");
}

// ────────────────────────────────────────────────────────────────────
// addTriangle()'s "erase on add" branch: gluing to an already-present
// neighbour along a shared, previously-unglued interior edge must clear
// that edge from improperEdges_ (the mirror image of the removeTriangle
// bug fixed above), and must NOT be mistaken for a self-intersection --
// identifying the two ends of a genuinely shared edge is exactly what a
// valid gluing is supposed to do.
// ────────────────────────────────────────────────────────────────────
void test_gluing_clears_improper_edge_and_avoids_false_self_intersection() {
    std::cout << "\n--- gluing clears improperEdges_ without a false "
                "self-intersection ---\n";

    regina::Triangulation<4> fourBall;
    fourBall.newSimplex();
    regina::Triangulation<3> tri = fourBall.boundaryComponent(0)->build();
    SurfaceFinder<3> g(tri, SurfaceCondition::all);

    regina::Triangle<3> *X = tri.triangle(0);
    const auto &adjX = g.adjacencyOf(X);
    if (adjX.empty()) {
        std::cout << red << "  FAIL: triangle 0 has no neighbours to test with\n"
                  << resetColor;
        ++failed_count;
        return;
    }
    auto firstNeighbour = adjX.begin();
    regina::Triangle<3> *Y = firstNeighbour->first->f;
    regina::Edge<3> *shared = X->edge(firstNeighbour->second[0].srcFacet);

    KnottedSurface<3> ks(&tri);
    EXPECT_EQ(ks.addTriangle(X, GluingNode<3>::AdjList{}), true,
             "X adds alone with no neighbours present");
    EXPECT_EQ((int)ks.improperEdges_.count(shared), 1,
             "X's shared, un-glued, non-boundary edge is flagged improper");

    EXPECT_EQ(ks.addTriangle(Y, g.adjacencyOf(Y)), true,
             "Y adds and glues to X along the shared edge");
    EXPECT_EQ((int)ks.improperEdges_.count(shared), 0,
             "gluing clears the shared edge from improperEdges_");
    EXPECT_EQ(ks.hasSelfIntersection(), false,
             "gluing along a shared edge is not mistaken for a self-intersection");
}

// ────────────────────────────────────────────────────────────────────
// Two triangles sharing an ambient vertex without a glued edge between
// them must be rejected as self-intersecting (two distinct components of
// the abstract surface would land on the same ambient point), and the
// rejected addTriangle() must fully roll back: no growth in surface_, no
// leftover preimage mapping, and no lingering self-intersection flag.
// ────────────────────────────────────────────────────────────────────
void test_self_intersection_detected_and_rolled_back() {
    std::cout
        << "\n--- self-intersection is detected and fully rolled back ---\n";

    regina::Triangulation<4> fourBall;
    fourBall.newSimplex();
    regina::Triangulation<3> tri = fourBall.boundaryComponent(0)->build();

    regina::Triangle<3> *A = nullptr, *B = nullptr;
    for (regina::Triangle<3> *a : tri.triangles()) {
        for (regina::Triangle<3> *b : tri.triangles()) {
            if (a == b)
                continue;
            bool sharesVertex = false;
            for (int i = 0; i < 3 && !sharesVertex; ++i)
                for (int j = 0; j < 3 && !sharesVertex; ++j)
                    if (a->vertex(i) == b->vertex(j))
                        sharesVertex = true;
            if (!sharesVertex)
                continue;
            bool sharesEdge = false;
            for (int i = 0; i < 3 && !sharesEdge; ++i)
                for (int j = 0; j < 3 && !sharesEdge; ++j)
                    if (a->edge(i) == b->edge(j))
                        sharesEdge = true;
            if (sharesEdge)
                continue;
            A = a;
            B = b;
            break;
        }
        if (A)
            break;
    }

    if (!A || !B) {
        std::cout << red
                  << "  FAIL: couldn't find two triangles sharing a vertex "
                     "but no edge\n"
                  << resetColor;
        ++failed_count;
        return;
    }

    KnottedSurface<3> ks(&tri);
    EXPECT_EQ(ks.addTriangle(A, GluingNode<3>::AdjList{}), true, "A adds alone");
    EXPECT_EQ(ks.hasSelfIntersection(), false,
             "a single triangle never self-intersects");

    EXPECT_EQ(ks.addTriangle(B, GluingNode<3>::AdjList{}), false,
             "B, sharing only a vertex with A, is rejected as self-intersecting");
    EXPECT_EQ((int)ks.surface().countTriangles(), 1,
             "the rejected triangle's abstract simplex was rolled back");
    EXPECT_EQ(ks.preimage(B) == nullptr, true,
             "the rejected triangle has no leftover preimage mapping");
    EXPECT_EQ(ks.hasSelfIntersection(), false,
             "surface_ is back to just A, with no lingering self-intersection");
}

// ────────────────────────────────────────────────────────────────────
// An ambient edge shared by three distinct triangles (the shared face
// between two tetrahedra, plus one triangle from each side -- see the
// two-tetrahedra test above) can only ever be glued to *one* of them.
// Adding the third once the other two are already present must be
// rejected as a double-glued facet, and must fully roll back.
// ────────────────────────────────────────────────────────────────────
void test_double_glued_facet_rejected_and_rolled_back() {
    std::cout << "\n--- a double-glued facet is rejected and fully rolled "
                "back ---\n";

    regina::Triangulation<3> tri;
    auto *t0 = tri.newSimplex();
    auto *t1 = tri.newSimplex();
    t0->join(3, t1, regina::Perm<4>());
    SurfaceFinder<3> g(tri, SurfaceCondition::all);

    regina::Triangle<3> *M = nullptr;
    for (regina::Triangle<3> *t : tri.triangles()) {
        if (t->degree() == 2) {
            M = t;
            break;
        }
    }
    if (!M) {
        std::cout << red
                  << "  FAIL: couldn't find the shared interior triangle\n"
                  << resetColor;
        ++failed_count;
        return;
    }

    std::vector<regina::Triangle<3> *> byFacet[3];
    for (auto &[node, gluings] : g.adjacencyOf(M))
        for (auto &gluing : gluings)
            byFacet[gluing.srcFacet].push_back(node->f);

    regina::Triangle<3> *X = nullptr, *Y = nullptr;
    for (int i = 0; i < 3; ++i) {
        if (byFacet[i].size() == 2) {
            X = byFacet[i][0];
            Y = byFacet[i][1];
            break;
        }
    }
    if (!X || !Y) {
        std::cout << red
                  << "  FAIL: couldn't find a facet of M shared by two other "
                     "triangles\n"
                  << resetColor;
        ++failed_count;
        return;
    }

    KnottedSurface<3> ks(&tri);
    EXPECT_EQ(ks.addTriangle(X, GluingNode<3>::AdjList{}), true, "X adds alone");
    EXPECT_EQ(ks.addTriangle(Y, g.adjacencyOf(Y)), true,
             "Y adds and glues to X (only X present so far)");

    EXPECT_EQ(ks.addTriangle(M, g.adjacencyOf(M)), false,
             "M can't glue the same facet to both X and Y");
    EXPECT_EQ((int)ks.surface().countTriangles(), 2,
             "the rejected triangle's abstract simplex was rolled back");
    EXPECT_EQ(ks.preimage(M) == nullptr, true,
             "the rejected triangle has no leftover preimage mapping");
}

// ────────────────────────────────────────────────────────────────────
// findSurfaces(startingTriangles) defers the self-intersection check
// while assembling its base (checkSelfIntersection=false) specifically so
// a base with an internal branch point can be built at all -- see
// addTriangle's documentation. This searches for a genuine branch point
// (two triangles A, B sharing an ambient vertex with no glued edge between
// them, plus a connector C genuinely adjacent to both) and empirically
// verifies -- by actually trying candidates and checking
// hasSelfIntersection() -- that: (a) eager per-triangle checking rejects
// B before C is present, and (b) deferred checking lets the whole base
// assemble, with the self-intersection genuinely resolving once C unites
// A and B's shared vertex.
// ────────────────────────────────────────────────────────────────────
void test_deferred_self_intersection_check_resolves_branch_point() {
    std::cout << "\n--- deferred self-intersection check resolves a branch "
                "point ---\n";

    regina::Triangulation<4> fourBall;
    fourBall.newSimplex();
    regina::Triangulation<3> tri = fourBall.boundaryComponent(0)->build();
    SurfaceFinder<3> g(tri, SurfaceCondition::all);

    auto sharesVertex = [](regina::Triangle<3> *a, regina::Triangle<3> *b) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                if (a->vertex(i) == b->vertex(j))
                    return true;
        return false;
    };
    auto sharesEdge = [](regina::Triangle<3> *a, regina::Triangle<3> *b) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                if (a->edge(i) == b->edge(j))
                    return true;
        return false;
    };

    regina::Triangle<3> *foundA = nullptr, *foundB = nullptr, *foundC = nullptr;

    for (regina::Triangle<3> *a : tri.triangles()) {
        for (regina::Triangle<3> *b : tri.triangles()) {
            if (a == b || !sharesVertex(a, b) || sharesEdge(a, b))
                continue;

            KnottedSurface<3> probe(&tri);
            if (!probe.addTriangle(a, g.adjacencyOf(a), false))
                continue;
            if (!probe.addTriangle(b, g.adjacencyOf(b), false))
                continue;
            if (!probe.hasSelfIntersection())
                continue; // not actually a branch point -- already fine

            for (auto &[node, gluing] : g.adjacencyOf(a)) {
                regina::Triangle<3> *c = node->f;
                if (c == b)
                    continue;
                KnottedSurface<3> attempt = probe;
                if (!attempt.addTriangle(c, g.adjacencyOf(c), false))
                    continue;
                if (!attempt.hasSelfIntersection()) {
                    foundA = a;
                    foundB = b;
                    foundC = c;
                    break;
                }
            }
            if (foundC)
                break;
        }
        if (foundC)
            break;
    }

    if (!foundA || !foundB || !foundC) {
        std::cout << red
                  << "  FAIL: couldn't find a branch-point triple in this "
                     "triangulation\n"
                  << resetColor;
        ++failed_count;
        return;
    }

    // Eager per-triangle checking rejects the branch point immediately:
    {
        KnottedSurface<3> eager(&tri);
        EXPECT_EQ(eager.addTriangle(foundA, g.adjacencyOf(foundA)), true,
                 "A adds alone under eager checking");
        EXPECT_EQ(
            eager.addTriangle(foundB, g.adjacencyOf(foundB)), false,
            "eager per-triangle self-intersection checking rejects B before "
            "the connecting triangle is present");
    }

    // Deferred checking lets the whole branch-point base assemble, and the
    // self-intersection genuinely resolves once the connector is in:
    KnottedSurface<3> deferred(&tri);
    EXPECT_EQ(deferred.addTriangle(foundA, g.adjacencyOf(foundA), false), true,
             "A adds alone (deferred)");
    EXPECT_EQ(deferred.addTriangle(foundB, g.adjacencyOf(foundB), false), true,
             "B adds alongside A (deferred)");
    EXPECT_EQ(deferred.hasSelfIntersection(), true,
             "the branch point genuinely isn't resolved until the connector "
             "is added");
    EXPECT_EQ(deferred.addTriangle(foundC, g.adjacencyOf(foundC), false), true,
             "connecting triangle C adds");
    EXPECT_EQ(deferred.hasSelfIntersection(), false,
             "once C is in, the full base has no self-intersection after all");
}

// ────────────────────────────────────────────────────────────────────
// Walks a single tetrahedron's boundary from one triangle up to the full
// closed S^2 and back down to three triangles, cross-checking at every
// step against regina's own countBoundaryFacets() (never trusting
// numBoundaryFacets_ in isolation), and checking detail()'s naming and
// cache invalidation across an add/remove cycle. Also exercises
// image()/preimage() as a round trip once all four triangles are in.
// ────────────────────────────────────────────────────────────────────
void test_tetrahedron_boundary_lifecycle() {
    std::cout << "\n--- single tetrahedron boundary: lifecycle of "
                "numBoundaryFacets_/detail() ---\n";

    regina::Triangulation<3> tri;
    tri.newSimplex();
    SurfaceFinder<3> g(tri, SurfaceCondition::all);
    regina::Triangle<3> *t0 = tri.triangle(0);
    regina::Triangle<3> *t1 = tri.triangle(1);
    regina::Triangle<3> *t2 = tri.triangle(2);
    regina::Triangle<3> *t3 = tri.triangle(3);

    KnottedSurface<3> ks(&tri);

    EXPECT_EQ(ks.addTriangle(t0, GluingNode<3>::AdjList{}), true,
             "t0 adds alone");
    EXPECT_EQ((int)ks.surface().countBoundaryFacets(), 3,
             "a single triangle has 3 boundary facets");
    EXPECT_EQ(ks.isClosed(), false, "a single triangle isn't closed");
    EXPECT_EQ(ks.detail(), "Disc", "a single triangle is a disc");

    EXPECT_EQ(ks.addTriangle(t1, g.adjacencyOf(t1)), true,
             "t1 adds and glues to t0");
    EXPECT_EQ(ks.detail(), "Disc",
             "two triangles glued along one edge is still a disc");

    EXPECT_EQ(ks.addTriangle(t2, g.adjacencyOf(t2)), true, "t2 adds");
    EXPECT_EQ(ks.detail(), "Disc", "three triangles (a 'cap') is still a disc");

    EXPECT_EQ(ks.addTriangle(t3, g.adjacencyOf(t3)), true,
             "t3 adds, closing up the boundary");
    EXPECT_EQ(ks.isClosed(), true, "all four triangles close up into a sphere");
    EXPECT_EQ((int)ks.surface().countBoundaryFacets(), 0,
             "a closed surface has no boundary facets");
    EXPECT_EQ(ks.detail(), "Sphere", "all four triangles form a sphere");

    // image()/preimage() round trip while all four are present.
    bool roundTripOk = true;
    for (regina::Triangle<3> *t : {t0, t1, t2, t3}) {
        regina::Triangle<2> *pre = ks.preimage(t);
        if (pre == nullptr || ks.image(pre) != t)
            roundTripOk = false;
    }
    EXPECT_EQ(roundTripOk, true,
             "preimage()/image() round-trip correctly for all four triangles");

    ks.removeTriangle(t3);
    EXPECT_EQ(ks.isClosed(), false, "removing one triangle re-opens the sphere");
    EXPECT_EQ((int)ks.surface().countBoundaryFacets(), 3,
             "removing one triangle exposes exactly its own three edges");
    EXPECT_EQ(ks.detail(), "Disc",
             "detail() recomputes correctly after removal, not left stale "
             "from being a sphere");
}

// ────────────────────────────────────────────────────────────────────
// Copying a KnottedSurface mid-construction and then mutating the copy
// must not affect the original -- targets the copy constructor's/
// operator='s requirement to rebuild inv_ from the copy's own surface_
// rather than blitting it from the source (see the class-level comment on
// inv_).
// ────────────────────────────────────────────────────────────────────
void test_copy_independence() {
    std::cout << "\n--- copying mid-construction doesn't alias state ---\n";

    regina::Triangulation<3> tri;
    tri.newSimplex();
    SurfaceFinder<3> g(tri, SurfaceCondition::all);
    regina::Triangle<3> *t0 = tri.triangle(0);
    regina::Triangle<3> *t1 = tri.triangle(1);
    regina::Triangle<3> *t2 = tri.triangle(2);

    KnottedSurface<3> original(&tri);
    original.addTriangle(t0, GluingNode<3>::AdjList{});
    original.addTriangle(t1, g.adjacencyOf(t1));

    KnottedSurface<3> copiedByCtor = original;
    copiedByCtor.addTriangle(t2, g.adjacencyOf(t2));

    EXPECT_EQ((int)original.surface().countTriangles(), 2,
             "original untouched by mutating a copy-constructed copy");
    EXPECT_EQ((int)copiedByCtor.surface().countTriangles(), 3,
             "the copy-constructed copy reflects its own mutation");
    EXPECT_EQ(original.preimage(t0) != nullptr &&
                 original.image(original.preimage(t0)) == t0,
             true, "original's own mapping is intact after the copy diverges");
    EXPECT_EQ(original.preimage(t2) == nullptr, true,
             "original has no mapping for a triangle only added to the copy");
    EXPECT_EQ(copiedByCtor.preimage(t2) != nullptr, true,
             "the copy's mapping includes the triangle only it received");

    KnottedSurface<3> copiedByAssign(&tri);
    copiedByAssign = original;
    copiedByAssign.addTriangle(t2, g.adjacencyOf(t2));

    EXPECT_EQ((int)original.surface().countTriangles(), 2,
             "original untouched by mutating an operator=-copied copy either");
    EXPECT_EQ(copiedByAssign.preimage(t2) != nullptr, true,
             "the operator=-copied copy's mapping includes the triangle only "
             "it received");
}

int main() {
    test_single_tetrahedron();
    test_starting_triangles_overload();
    test_starting_triangles_disconnected_base();
    test_accessor_and_stream_output();
    test_two_tetrahedra();
    test_removeTriangle_clears_improper_edges();
    test_three_sphere();
    test_cone_pipeline_finds_disc_bounding_known_unknot();
    test_four_sphere();

    test_found_surfaces_satisfy_invariants();
    test_surface_equality_order_independent();
    test_gluing_clears_improper_edge_and_avoids_false_self_intersection();
    test_self_intersection_detected_and_rolled_back();
    test_double_glued_facet_rejected_and_rolled_back();
    test_deferred_self_intersection_check_resolves_branch_point();
    test_tetrahedron_boundary_lifecycle();
    test_copy_independence();

    std::cout << "\n"
              << bold << (failed_count > 0 ? red : green) << "=== " << passed
              << " passed, " << failed_count << " failed ===" << resetColor
              << "\n";
    return failed_count > 0 ? 1 : 0;
}
