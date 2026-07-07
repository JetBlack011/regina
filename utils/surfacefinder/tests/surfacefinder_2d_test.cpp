// surfacefinder_2d_test.cpp
// Ground-truth surface-search tests: build a *known* 2-manifold directly as
// a regina::Triangulation<2> (so its own topological type -- orientability,
// genus, punctures -- is independently verifiable via Regina's own
// primitives, not by hand-derivation), thicken it one dimension up into an
// ambient 3-manifold via CobordismBuilder<2>::thicken(), and confirm
// SurfaceFinder<3>'s seeded search actually finds a surface that exactly
// reconstructs it (isoSig match), not just a surface of the right *count*.
//
// This complements surfacefinder_test.cpp, whose expected surface counts on
// hand-built 3-manifolds get increasingly hard to verify by hand as the
// ambient triangulation grows (see e.g. the two-tetrahedra test's own
// comment on this) -- pinning against an independently-constructed known
// surface, rather than a hand-counted total, is a strictly stronger check.
//
// Why cone() one dimension up rather than feeding a Triangulation<2>
// straight into SurfaceFinder<2>: KnottedSurface<dim> needs a
// Triangulation<dim - 1> to represent tri_'s boundary (for boundary()'s Link
// extraction), and Regina does not support Triangulation<1>
// (regina-core.h's supportedDim() starts at 2) -- so
// SurfaceFinder<2>/KnottedSurface<2> cannot compile. Separately, even if
// that were worked around, Triangle<2> aliases to Simplex<2> (see
// SafeFaceHelper in triangulation/forward.h), so a same-dimension search
// would only ever explore the ambient triangulation's own simplices under
// the identity embedding -- it could never exercise the self-intersection
// logic that only arises when the ambient dimension exceeds the surface's.
// Thickening one dimension up keeps the search in the regime the tool is
// actually meant for.
//
// Why thicken() rather than cone(): cone()'s new apex vertex's link is
// always (topologically) `expected` itself, and Regina's
// Triangulation<3>::isValid() only accepts a vertex link that's a sphere
// (ordinary interior point), a disc (ordinary boundary point), or any other
// *closed* surface (an "ideal" vertex, a deliberately-supported feature for
// cusped triangulations -- this is why coning a torus or Klein bottle is
// fine despite neither being a sphere). There is no valid vertex-link
// category for a surface that has boundary and isn't a disc, so cone()
// necessarily produces an invalid triangulation for e.g. the annulus (2
// boundary loops) or the Mobius band (non-orientable) -- confirmed directly
// by inspecting the apex vertex's link in each case. This isn't a bug, just
// an inherent limitation of coning, so it can't be worked around by fixing
// cone() -- thicken() (the literal product `expected x I`) has no apex at
// all, so this obstruction doesn't apply, and it turns out to reconstruct
// every case below (including the closed ones) just as well as cone() did.
//
// CobordismBuilder<dim>'s automatic vertex-reordering (needed by
// thicken()/cone() -- see isOrdered()) is only implemented for dim == 3, so
// every Triangulation<2> fed in here must already satisfy isOrdered() as
// built, or be relabelled into that form first (see
// tryOrderTriangulation() below, which handles this for
// regina::Example<2>::orientable()/nonOrientable() at every genus/puncture
// combination tried -- these are not ordered as built).
//
// torus() and annulus() (both isOrdered()) originally triggered a false
// "self-intersects" rejection here, regardless of whether cone() or
// thicken() built the ambient manifold. Root cause: GluingNode::AdjList
// (gluing.h) used to map straight to a single Gluing per neighbour, so two
// triangles sharing more than one edge with each other (true of both the
// minimal one-vertex 2-triangle torus and the 2-triangle annulus -- e.g. the
// torus's triangle 0 has all three of its edges independently matching
// triangle 1, on three different local-edge pairs) had all but one of those
// gluings silently overwritten, starving addTriangle()'s vertex-union step
// and producing a false positive from the self-intersection tracker on an
// otherwise perfectly valid, embeddable surface. Fixed by having AdjList
// hold a vector of gluings per neighbour (see gluing.h's class comment).
// Both cases now pass below.
//
// tryOrderTriangulation() is a test-only dimension-2 analogue of
// Triangulation<3>::order() (engine/triangulation/dim3/reorder.cpp), needed
// because CobordismBuilder's own automatic reordering is dim == 3 only (see
// above) and regina::Example<2>::orientable()/nonOrientable() -- needed for
// genus/puncture combinations beyond what the dedicated small constructions
// used elsewhere in this file provide -- are never ordered as built.
// Same strategy as the dim-3 code, one dimension down: rather than directly
// searching over each triangle's 6 possible local vertex relabellings (a
// much larger, harder-to-reconcile-between-neighbours search space),
// backtrack over a +1/-1 *orientation* for each ambient *edge*, then derive
// each triangle's local vertex ranking from in-degree counting -- exactly
// perm_from_edges()'s trick, minus the extra bookkeeping dim 3 needs for a
// tetrahedron's 6 edges and 4 sub-faces (a triangle has 3 edges and, being
// 2-dimensional, is its own only face to check). A local ranking is
// well-defined exactly when the triangle's 3 local edge-directions are
// acyclic (no directed 3-cycle) -- which also happens to be exactly the
// condition under which every gluing incident to that triangle can satisfy
// isOrdered(). Finds a solution in well under a millisecond for every case
// tried here, up to 12 triangles; scoped to this test file rather than
// folded into CobordismBuilder since dim == 3 already has its own
// (differently-implemented) solution to the same problem -- unifying the
// two is future work, not needed just to unblock these tests.

#include <array>
#include <iostream>
#include <optional>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/example2.h>
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
std::ostream &yellow(std::ostream &os) {
    return colorEnabled() ? os << "\033[33m" : os;
}
std::ostream &bold(std::ostream &os) {
    return colorEnabled() ? os << "\033[1m" : os;
}
std::ostream &resetColor(std::ostream &os) {
    return colorEnabled() ? os << "\033[0m" : os;
}
} // namespace

#define EXPECT_EQ(actual, expected, desc)                                     \
    do {                                                                     \
        auto _a = (actual);                                                  \
        auto _e = (expected);                                                \
        if (_a == _e) {                                                      \
            std::cout << green << "  PASS: " << resetColor << (desc) << "\n"; \
            ++passed;                                                        \
        } else {                                                             \
            std::cout << red << "  FAIL: " << (desc) << "\n"                 \
                      << "        expected " << _e << ", got " << _a         \
                      << resetColor << "\n";                                 \
            ++failed_count;                                                  \
        }                                                                    \
    } while (0)

#define EXPECT_GE(actual, expected, desc)                                     \
    do {                                                                     \
        auto _a = (actual);                                                  \
        auto _e = (expected);                                                \
        if (_a >= _e) {                                                      \
            std::cout << green << "  PASS: " << resetColor << (desc) << "\n"; \
            ++passed;                                                        \
        } else {                                                             \
            std::cout << red << "  FAIL: " << (desc) << "\n"                 \
                      << "        expected >= " << _e << ", got " << _a      \
                      << resetColor << "\n";                                 \
            ++failed_count;                                                  \
        }                                                                    \
    } while (0)

// ────────────────────────────────────────────────────────────────────
// See the file-level comment for why this exists and the strategy (an
// edge-orientation backtracking search, mirroring
// Triangulation<3>::order()). Returns a triangulation isomorphic to `tri`
// (same isoSig) with each triangle's own 3 vertices possibly relabelled so
// that CobordismBuilder<2>::isOrdered() holds, or nullopt if the search
// exhausts every edge-orientation assignment without finding one.
// ────────────────────────────────────────────────────────────────────
std::optional<regina::Triangulation<2>>
tryOrderTriangulation(const regina::Triangulation<2> &tri) {
    // t's local edge `localEdge` connects t's own local vertices (a, b) --
    // the two local vertices other than localEdge itself; orient[]'s +1
    // means "points from a to b" by this same convention, for whichever
    // ambient edge realises localEdge.
    auto localDir = [](regina::Triangle<2> *t, int localEdge, int &a, int &b) {
        regina::Perm<3> em = t->edgeMapping(localEdge);
        a = em[0];
        b = em[1];
    };

    int nEdges = tri.countEdges();
    std::vector<int> orient(nEdges, 0); // 0 unset, else +1/-1

    // Whether t's own 3 local edges, under the *current* (possibly still
    // partial) orient[], are acyclic. An unassigned edge makes this
    // trivially pass (a 3-cycle needs all 3 edges assigned to exist) --
    // mirrors check_consistency_on_tet()'s same early-exit for a partial
    // assignment.
    auto consistent = [&](regina::Triangle<2> *t) {
        int from[3]; // from[le]: local vertex local edge le points from
        for (int le = 0; le < 3; ++le) {
            int o = orient[t->edge(le)->index()];
            if (o == 0) {
                from[le] = -1;
                continue;
            }
            int a, b;
            localDir(t, le, a, b);
            from[le] = (o == +1) ? a : b;
        }
        if (from[0] == -1 || from[1] == -1 || from[2] == -1)
            return true;

        // Local edge le connects the two local vertices other than le.
        bool zeroToOne = (from[2] == 0);
        bool oneToTwo = (from[0] == 1);
        bool zeroToTwo = (from[1] == 0);
        bool cycleForward = zeroToOne && oneToTwo && !zeroToTwo;
        bool cycleBackward = !zeroToOne && !oneToTwo && zeroToTwo;
        return !(cycleForward || cycleBackward);
    };

    auto checkAroundEdge = [&](regina::Edge<2> *e) {
        for (const auto &emb : e->embeddings())
            if (!consistent(emb.triangle()))
                return false;
        return true;
    };

    // Same try-(+1)-then-(-1)-then-backtrack loop as ordering_iso().
    int i = 0;
    while (i < nEdges) {
        if (i < 0)
            return std::nullopt;
        regina::Edge<2> *e = tri.edge(i);
        if (orient[i] == 0) {
            orient[i] = +1;
            if (checkAroundEdge(e))
                ++i;
        } else if (orient[i] == +1) {
            orient[i] = -1;
            if (checkAroundEdge(e))
                ++i;
        } else {
            orient[i] = 0;
            --i;
        }
    }

    // Derive each triangle's local vertex ranking from in-degree counting
    // (perm_from_edges()'s trick) and rebuild with that ranking applied.
    int n = tri.countTriangles();
    std::vector<regina::Perm<3>> relabel(n); // new-local -> old-local
    for (regina::Triangle<2> *t : tri.triangles()) {
        int indeg[3] = {0, 0, 0};
        for (int le = 0; le < 3; ++le) {
            int a, b;
            localDir(t, le, a, b);
            if (orient[t->edge(le)->index()] == -1)
                std::swap(a, b);
            ++indeg[b];
        }
        // Perm<3>(indeg[0..2]) maps old-local -> new-rank; invert for the
        // new-local -> old-local convention used below.
        relabel[t->index()] =
            regina::Perm<3>(indeg[0], indeg[1], indeg[2]).inverse();
    }

    regina::Triangulation<2> rebuilt;
    rebuilt.newTriangles(n);
    std::vector<std::array<bool, 3>> done(n, {false, false, false});
    for (regina::Triangle<2> *t : tri.triangles()) {
        int i2 = t->index();
        for (int fOld = 0; fOld < 3; ++fOld) {
            regina::Triangle<2> *o = t->adjacentSimplex(fOld);
            if (!o)
                continue;
            int j = o->index();
            int fNew = relabel[i2].inverse()[fOld];
            if (done[i2][fNew])
                continue; // already joined from the other side
            regina::Perm<3> gOld = t->adjacentGluing(fOld);
            regina::Perm<3> gNew = relabel[j].inverse() * gOld * relabel[i2];
            rebuilt.triangle(i2)->join(fNew, rebuilt.triangle(j), gNew);
            done[i2][fNew] = true;
            int otherFNew = relabel[j].inverse()[gOld[fOld]];
            done[j][otherFNew] = true;
        }
    }
    return rebuilt;
}

// ────────────────────────────────────────────────────────────────────
// Thickens `expected` into an ambient 3-manifold, recovers the counterpart
// of every one of `expected`'s own triangles in the thickening's bottom
// layer, seeds a search with exactly that set, and checks that one of the
// resulting surfaces has an isoSig identical to `expected`'s -- i.e. the
// search doesn't just find *some* surface of a plausible size, it finds the
// *exact* known surface.
// ────────────────────────────────────────────────────────────────────
void test_thicken_pipeline_reconstructs_known_surface(
    const std::string &name, const regina::Triangulation<2> &input) {
    std::cout << "\n--- thicken() pipeline: " << name << " ---\n";
    std::cout << "    " << input.countTriangles() << " triangles, closed="
              << input.isClosed() << ", orientable="
              << input.isOrientable() << "\n";

    regina::Triangulation<2> expected = input;
    if (!CobordismBuilder<2>::isOrdered(expected)) {
        auto reordered = tryOrderTriangulation(expected);
        if (!reordered) {
            std::cout << yellow
                       << "    SKIPPED: no relabelling found that satisfies "
                          "CobordismBuilder-ordered form\n"
                       << resetColor;
            return;
        }
        expected = *reordered;
    }

    std::string expectedSig = expected.isoSig();

    CobordismBuilder<2> cob(expected);
    auto &thickened = cob.thicken();

    // thicken() builds, for each triangle of expected, a 3-tetrahedron
    // prism over it (see SimplicialPrism); the prism's *bottom* layer (the
    // triangle x {0} slice) is entirely contained in sub-simplex 0, as the
    // face opposite local vertex 3 -- see SimplicialPrism::capTop's class
    // comment for the same fact about the analogous *top* face. Looking
    // triangles up via cob.baseTriangulation() (not `expected` itself)
    // matters in general, since the constructor may reorder its input --
    // here isOrdered() is already known to hold, so no reordering actually
    // happens and indices line up either way, but going through
    // baseTriangulation() keeps this correct even if that ever changes.
    std::unordered_set<regina::Triangle<3> *> base;
    for (regina::Triangle<2> *t : expected.triangles()) {
        const regina::Simplex<2> *baseSimplex =
            cob.baseTriangulation().simplex(t->index());
        regina::Simplex<3> *bottomPrismSimplex = cob.currentTopSimplex(baseSimplex, 0);
        base.insert(bottomPrismSimplex->triangle(3));
    }
    EXPECT_EQ((int)base.size(), (int)expected.countTriangles(),
             name + ": found one bottom-layer counterpart triangle per "
                    "original triangle");

    SurfaceFinder<3> g(thickened, SurfaceCondition::boundary);
    auto &surfaces = g.findSurfaces(base);

    EXPECT_GE((int)surfaces.size(), 1,
              name + ": seeded search finds at least one surface");

    bool foundMatch = false;
    for (const auto &surf : surfaces) {
        if (surf.surface().isoSig() == expectedSig) {
            foundMatch = true;
            break;
        }
    }
    EXPECT_EQ(foundMatch, true,
             name + ": search finds a surface whose isoSig exactly matches "
                    "the known input triangulation");
}

int main() {
    // Small, dedicated constructions (each implemented via the
    // dimension-generic ExampleBase<dim> routines: simplicialSphere, ball,
    // ballBundle, twistedBallBundle, sphereBundle, twistedSphereBundle) --
    // already ordered as built, so tryOrderTriangulation() isn't needed.
    // Cover every combination of orientable/non-orientable and
    // closed/with-boundary, with torus and annulus alongside sphere/Klein
    // bottle and disc/Mobius band as a second closed and a second
    // with-boundary example respectively.
    test_thicken_pipeline_reconstructs_known_surface(
        "sphere", regina::Example<2>::sphereTetrahedron());
    test_thicken_pipeline_reconstructs_known_surface(
        "torus", regina::Example<2>::torus());
    test_thicken_pipeline_reconstructs_known_surface(
        "Klein bottle", regina::Example<2>::kb());
    test_thicken_pipeline_reconstructs_known_surface(
        "disc", regina::Example<2>::disc());
    test_thicken_pipeline_reconstructs_known_surface(
        "annulus", regina::Example<2>::annulus());
    test_thicken_pipeline_reconstructs_known_surface(
        "Mobius band", regina::Example<2>::mobius());

    // Higher genus / more boundary components / both -- regina::Example<2>
    // has no dedicated small constructions for these, so these go through
    // orientable()/nonOrientable() (not ordered as built) and rely on
    // tryOrderTriangulation() to fix that up.
    test_thicken_pipeline_reconstructs_known_surface(
        "genus-2 closed orientable surface", regina::Example<2>::orientable(2, 0));
    test_thicken_pipeline_reconstructs_known_surface(
        "genus-3 closed non-orientable surface",
        regina::Example<2>::nonOrientable(3, 0));
    test_thicken_pipeline_reconstructs_known_surface(
        "3-punctured sphere (pair of pants)",
        regina::Example<2>::orientable(0, 3));
    test_thicken_pipeline_reconstructs_known_surface(
        "genus-2 orientable surface with 2 punctures",
        regina::Example<2>::orientable(2, 2));

    std::cout << "\n"
              << bold << (failed_count > 0 ? red : green) << "=== " << passed
              << " passed, " << failed_count << " failed ===" << resetColor
              << "\n";
    return failed_count > 0 ? 1 : 0;
}
