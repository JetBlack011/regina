// embeddingsearch_test.cpp
// Tests for the BoundaryCondition-based output filtering added to
// EmbeddedSubmanifold/EmbeddingSearch (isClosed(), isProper(),
// boundaryComponentsMapInjectively(), satisfies()) -- see embeddingsearch.h.
// Cases exercise: a trivial "everything is boundary" ball (K4 gluing graph),
// a "proper but not connected" thickened annulus (two of its own boundary
// components map into a single ambient one), an engineered ambient-interior
// edge that makes a lone face fail isProper(), and a triangle whose three
// edges are all identified together (irreparably self-folded, excluded from
// the DFS graph entirely -- see hasIrreparableSelfFold() in
// embeddingsearch.h).

#include <iostream>
#include <maths/perm.h>
#include <sstream>
#include <stdexcept>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/example2.h>
#include <unistd.h>

#include "cobordismbuilder.h"
#include "embeddingsearch.h"

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

namespace {
// search() only ever reports its result via std::cerr (it's still a pure
// counting pass -- see embeddingsearch.h); capture that here and pull the
// "Total embedded submanifolds (<cond>): N" line back out as a count.
long long runFilteredCount(EmbeddingSearch<3, 2> &e, unsigned numThreads,
                           BoundaryCondition cond) {
    std::ostringstream captured;
    std::streambuf *oldBuf = std::cerr.rdbuf(captured.rdbuf());
    e.search(numThreads, cond);
    std::cerr.rdbuf(oldBuf);

    std::string out = captured.str();
    auto tPos = out.rfind("Total embedded submanifolds (");
    if (tPos == std::string::npos)
        throw std::runtime_error(
            "runFilteredCount(): search() did not print a total count line");
    auto colonPos = out.find("): ", tPos);
    return std::stoll(out.substr(colonPos + 3));
}
} // namespace

// ─────────────────────────────────────────────────────────────────────────────
// Single tetrahedron (B³): its 4 boundary triangles form a K4 gluing graph,
// so every non-empty subset of them is connected -- 2^4 - 1 = 15 connected
// embedded submanifolds in total, all embeddable, and every ambient edge is
// boundary (there are no gluings anywhere in this triangulation). So:
//   all = 15 (every subset)
//   closed = 1 (only the full 4-face set closes into S²)
//   proper = 15 (trivial -- every ambient edge is boundary)
//   connected = 15 (at most one submanifold boundary component in every
//                   case, one ambient boundary component total)
// ─────────────────────────────────────────────────────────────────────────────
void test_tetrahedron_boundary_conditions() {
    std::cout << "\n--- B³ (single tetrahedron): BoundaryCondition filtering "
                 "---\n";

    regina::Triangulation<3> ball;
    ball.newTetrahedron();

    EmbeddingSearch<3, 2> all(ball);
    EXPECT_EQ(runFilteredCount(all, 1, BoundaryCondition::all), 15LL,
              "--all finds all 15 connected embedded subsets");

    EmbeddingSearch<3, 2> closed(ball);
    EXPECT_EQ(runFilteredCount(closed, 1, BoundaryCondition::closed), 1LL,
              "--closed finds only the full 4-face set");

    EmbeddingSearch<3, 2> proper(ball);
    EXPECT_EQ(runFilteredCount(proper, 1, BoundaryCondition::proper), 15LL,
              "--proper finds all 15 (every ambient edge is boundary)");

    EmbeddingSearch<3, 2> connected(ball);
    EXPECT_EQ(runFilteredCount(connected, 1, BoundaryCondition::connected),
              15LL,
              "--connected finds all 15 (single ambient boundary component)");

    // Direct-construction cross-check, independent of search()/DFS: pin the
    // three predicates directly against the full embedding and a single-face
    // embedding.
    Skeleton<3, 2> skel(ball);

    EmbeddedSubmanifold<3, 2> full(skel);
    for (auto *t : ball.triangles())
        full.addFace(t->index());
    EXPECT_EQ(full.isClosed(), true, "full 4-face embedding is closed");
    EXPECT_EQ(full.isProper(), true, "full 4-face embedding is proper");
    EXPECT_EQ(full.boundaryComponentsMapInjectively(), true,
              "full 4-face embedding's boundary map is injective");

    EmbeddedSubmanifold<3, 2> single(skel);
    single.addFace(ball.triangle(0)->index());
    EXPECT_EQ(single.isClosed(), false, "single-face embedding is not closed");
    EXPECT_EQ(single.isProper(), true, "single-face embedding is proper");
    EXPECT_EQ(single.boundaryComponentsMapInjectively(), true,
              "single-face embedding's boundary map is injective");
}

// ─────────────────────────────────────────────────────────────────────────────
// Thickened annulus: doubling a *connected* annulus across both of its
// boundary circles yields a single connected closed surface (a torus), so
// the thickened cobordism has exactly one ambient boundary component even
// though the annulus itself has two boundary components. Embedding the
// front-copy layer (which reconstructs the annulus and sits entirely on the
// ambient boundary) is therefore proper, but NOT connected: its 2 boundary
// components both map into that single ambient boundary component.
// ─────────────────────────────────────────────────────────────────────────────
void test_thickened_annulus_proper_not_connected() {
    std::cout << "\n--- Thickened annulus: proper but not connected ---\n";

    auto annulus = regina::Example<2>::annulus();
    EXPECT_EQ((int)annulus.countBoundaryComponents(), 2,
              "annulus has 2 boundary components");

    CobordismBuilder<2> cob(annulus);
    auto &thickened = cob.thicken();
    EXPECT_EQ((int)thickened.countBoundaryComponents(), 1,
              "doubled annulus has a single ambient boundary component");

    std::vector<int> frontCopyIndices;
    for (regina::Triangle<2> *t : annulus.triangles()) {
        const regina::Simplex<2> *baseSimplex =
            cob.baseTriangulation().simplex(t->index());
        regina::Simplex<3> *bottomPrismSimplex =
            cob.currentTopSimplex(baseSimplex, 0);
        frontCopyIndices.push_back(bottomPrismSimplex->triangle(3)->index());
    }
    EXPECT_EQ((int)frontCopyIndices.size(), (int)annulus.countTriangles(),
              "found one front-copy triangle per original annulus triangle");

    Skeleton<3, 2> skel(thickened);
    EmbeddedSubmanifold<3, 2> emb(skel);
    for (int idx : frontCopyIndices)
        EXPECT_EQ(emb.addFace(idx), true, "front-copy face embeds");

    EXPECT_EQ(emb.isClosed(), false, "front copy of the annulus is not closed");
    EXPECT_EQ(emb.isProper(), true,
              "front copy lies entirely on the ambient boundary -- proper");
    EXPECT_EQ(emb.boundaryComponentsMapInjectively(), false,
              "both of its boundary components map into the same ambient "
              "boundary component -- not connected");
    EXPECT_EQ(emb.satisfies(BoundaryCondition::proper), true,
              "satisfies(proper) agrees with isProper()");
    EXPECT_EQ(emb.satisfies(BoundaryCondition::connected), false,
              "satisfies(connected) agrees with boundaryComponentsMapInjectively()");
}

// ─────────────────────────────────────────────────────────────────────────────
// Three tetrahedra glued cyclically around one shared edge (vertices {2,3}
// of each), so that edge is surrounded entirely by internal faces and is
// never boundary, while each tetrahedron's faces 2 and 3 stay unglued. A
// lone face containing that edge is embeddable on its own, but not proper:
// one of its 3 unglued subtri_ edges is the engineered ambient-interior one.
// ─────────────────────────────────────────────────────────────────────────────
void test_engineered_interior_edge_not_proper() {
    std::cout << "\n--- Engineered ambient-interior edge: not proper ---\n";

    regina::Triangulation<3> tri;
    auto *t0 = tri.newTetrahedron();
    auto *t1 = tri.newTetrahedron();
    auto *t2 = tri.newTetrahedron();
    regina::Perm<4> g(1, 0, 2, 3);
    t0->join(1, t1, g);
    t1->join(1, t2, g);
    t2->join(1, t0, g);

    EXPECT_EQ(tri.isValid(), true, "construction is a valid triangulation");

    // Precondition, verified directly against Regina's own skeleton (not
    // just hand-derived): the shared edge is genuinely ambient-interior.
    const regina::Edge<3> *sharedEdge = t0->edge(2, 3);
    EXPECT_EQ(sharedEdge->isBoundary(), false,
              "central edge is genuinely ambient-interior");

    Skeleton<3, 2> skel(tri);
    EmbeddedSubmanifold<3, 2> emb(skel);
    int idx = t0->triangle(1)->index(); // vertices {0,2,3}, contains the edge
    EXPECT_EQ(emb.addFace(idx), true, "single face embeds alone");
    EXPECT_EQ(emb.isProper(), false,
              "a lone face exposing the engineered interior edge is not "
              "proper");
    EXPECT_EQ(emb.satisfies(BoundaryCondition::proper), false,
              "satisfies(proper) agrees with isProper()");
}

// ─────────────────────────────────────────────────────────────────────────────
// A single tetrahedron self-glued on two of its own faces (isoSig "bkaaid":
// 1 vertex, 1 edge, 2 triangles, 1 tetrahedron) so that BOTH of its triangles
// end up with all three of their local edges identified to that one ambient
// edge -- a genuine 3-way self-collapse, not the ordinary 2-facet fold this
// code already supports (e.g. a Möbius band from one triangle). A subdim-cell
// can touch any one ambient facet at most twice, so neither triangle can ever
// be part of a valid embedded submanifold; hasIrreparableSelfFold() should
// flag both, EmbeddingSearch should exclude both from its DFS graph, and a
// full search() over this triangulation should complete instantly and safely
// instead of crashing (the original bug: Regina's own join() throwing
// "cannot join facets... already joined" -- see embeddingsearch.h).
// ─────────────────────────────────────────────────────────────────────────────
void test_triple_self_fold_excluded() {
    std::cout << "\n--- Triple self-fold (all 3 edges of a triangle "
                 "identified): excluded from the DFS graph ---\n";

    regina::Triangulation<3> tri;
    auto *t = tri.newTetrahedron();
    t->join(0, t, regina::Perm<4>(1, 2, 0, 3));
    t->join(2, t, regina::Perm<4>(0, 2, 3, 1));

    EXPECT_EQ(tri.isValid(), true, "construction is a valid triangulation");
    EXPECT_EQ((int)tri.countTriangles(), 2, "two distinct triangle faces");
    for (int i = 0; i < 2; ++i)
        EXPECT_EQ(tri.triangle(i)->edge(0) == tri.triangle(i)->edge(1) &&
                       tri.triangle(i)->edge(1) == tri.triangle(i)->edge(2),
                   true,
                   "triangle's 3 local edges are all the same ambient edge");

    Skeleton<3, 2> skel(tri);
    EXPECT_EQ((hasIrreparableSelfFold<3, 2>(skel.getNodes()[0].gluings)), true,
              "triangle 0 is flagged as irreparably self-folded");
    EXPECT_EQ((hasIrreparableSelfFold<3, 2>(skel.getNodes()[1].gluings)), true,
              "triangle 1 is flagged as irreparably self-folded");

    // Negative control: the ordinary 2-facet self-fold engineered in
    // test_engineered_interior_edge_not_proper() must NOT be flagged.
    {
        regina::Triangulation<3> foldTri;
        auto *t0 = foldTri.newTetrahedron();
        auto *t1 = foldTri.newTetrahedron();
        auto *t2 = foldTri.newTetrahedron();
        regina::Perm<4> g(1, 0, 2, 3);
        t0->join(1, t1, g);
        t1->join(1, t2, g);
        t2->join(1, t0, g);
        Skeleton<3, 2> foldSkel(foldTri);
        bool anyFlagged = false;
        for (size_t i = 0; i < foldSkel.numFaces(); ++i)
            anyFlagged = anyFlagged ||
                         hasIrreparableSelfFold<3, 2>(
                             foldSkel.getNodes()[i].gluings);
        EXPECT_EQ(anyFlagged, false,
                  "an ordinary self-fold triangulation flags no faces");
    }

    EmbeddingSearch<3, 2> search(tri);
    EXPECT_EQ(search.numEmbeddableFaces(), 0ULL,
              "both triangles are excluded from the DFS graph");

    // The original crash: EmbeddingSearch::search() over this triangulation
    // used to throw partway through addFace(). It should now complete
    // instantly (there's nothing left to search) instead.
    EXPECT_EQ(runFilteredCount(search, 1, BoundaryCondition::all), 0LL,
              "search() completes safely and finds nothing");
}

// ─────────────────────────────────────────────────────────────────────────────
// Batch boundary-link recognition (EmbeddingSearch<4, 2>::processSurfaceBoundaries()
// / linkTally()): a single pentachoron (B^4, dim = 4) whose boundary is
// ∂Δ^4 = S^3, exactly the "single tetrahedron" scenario from
// test_tetrahedron_boundary_conditions one dimension up -- every ambient
// edge is boundary, so every connected subset of the 10 triangles is
// proper/connected. In particular a single triangle T is a proper surface
// (a disc) whose boundary is T's own 3 edges -- trivially an unknotted loop,
// since it's the boundary of an embedded disc already sitting inside this
// S^3. search() with BoundaryCondition::connected should find it, batch it,
// and processSurfaceBoundaries() (called internally as part of search()'s
// final flush) should recognize its complement as the unknot and record it
// bounding a Disc.
// ─────────────────────────────────────────────────────────────────────────────
void test_boundary_link_batch_recognizes_unknot() {
    std::cout << "\n--- EmbeddingSearch<4,2>: batch boundary-link "
                 "recognition finds the unknot ---\n";

    regina::Triangulation<4> fourBall;
    fourBall.newSimplex();
    EXPECT_EQ((int)fourBall.countBoundaryComponents(), 1,
              "single pentachoron has one boundary component (S^3)");

    EmbeddingSearch<4, 2> e(fourBall);
    e.search(1, BoundaryCondition::connected);

    std::string summary = e.linkTally().summary();
    EXPECT_EQ(summary.find("Unknot") != std::string::npos, true,
              "batch processing recognized the unknot among the boundary "
              "links found");
    EXPECT_EQ(summary.find("Disc") != std::string::npos, true,
              "the unknot was recorded as bounding (at least) a Disc");
}

template <typename F> void run(const char *name, F fn) {
    std::cout << "\nRunning " << name << "...\n";
    try {
        fn();
    } catch (const std::exception &e) {
        std::cout << red << "  EXCEPTION: " << e.what() << resetColor << "\n";
        ++failed_count;
    }
}

int main() {
    run("test_tetrahedron_boundary_conditions",
        test_tetrahedron_boundary_conditions);
    run("test_thickened_annulus_proper_not_connected",
        test_thickened_annulus_proper_not_connected);
    run("test_engineered_interior_edge_not_proper",
        test_engineered_interior_edge_not_proper);
    run("test_triple_self_fold_excluded", test_triple_self_fold_excluded);
    run("test_boundary_link_batch_recognizes_unknot",
        test_boundary_link_batch_recognizes_unknot);

    std::cout << "\n"
              << bold << (failed_count > 0 ? red : green) << "=== " << passed
              << " passed, " << failed_count << " failed ===" << resetColor
              << "\n";
    return failed_count > 0 ? 1 : 0;
}
