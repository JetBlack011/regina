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

#include <algorithm>
#include <exception>
#include <iostream>
#include <sstream>
#include <maths/perm.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/example2.h>
#include <unistd.h>

#include "cobordismbuilder.h"
#include "embeddingsearch.h"
#include "surfacesearch.h"

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
long long runFilteredCount(EmbeddingSearch<3, 2> &e, unsigned numThreads,
                           BoundaryCondition cond) {
    return e.search(numThreads, cond).satisfyingCount;
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
    EXPECT_EQ((EmbeddedSubmanifold<3, 2>::hasIrreparableSelfGluing(
                  skel.getNodes()[0].gluings)), true,
              "triangle 0 is flagged as irreparably self-folded");
    EXPECT_EQ((EmbeddedSubmanifold<3, 2>::hasIrreparableSelfGluing(
                  skel.getNodes()[1].gluings)), true,
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
                         EmbeddedSubmanifold<3, 2>::hasIrreparableSelfGluing(
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
// Batch boundary-link recognition (SurfaceSearch::processSurfaceBoundaries()
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
    std::cout << "\n--- SurfaceSearch: batch boundary-link "
                 "recognition finds the unknot ---\n";

    regina::Triangulation<4> fourBall;
    fourBall.newSimplex();
    EXPECT_EQ((int)fourBall.countBoundaryComponents(), 1,
              "single pentachoron has one boundary component (S^3)");

    SurfaceSearch e(fourBall);
    e.search(1, BoundaryCondition::connected);

    std::string summary = e.linkTally().summary();
    EXPECT_EQ(summary.find("Unknot") != std::string::npos, true,
              "batch processing recognized the unknot among the boundary "
              "links found");
    EXPECT_EQ(summary.find("Disc") != std::string::npos, true,
              "the unknot was recorded as bounding (at least) a Disc");
}

// ─────────────────────────────────────────────────────────────────────────────
// Memory-bounding regression: SurfaceSearchLimits::pendingSurfaceCap adds
// backpressure (a DFS worker thread helps drain pendingSurfaces_ itself once
// it's over cap -- see SurfaceSearch::ThreadHook::onFlush()) purely to bound
// memory on a long search. This must never change *what* is found, only
// when/by whom it gets processed: running the same search with a
// deliberately tiny cap (forcing backpressure on nearly every find) must
// produce byte-identical final tallies to running it with a large cap.
// ─────────────────────────────────────────────────────────────────────────────
void test_backpressure_does_not_drop_or_double_count_surfaces() {
    std::cout << "\n--- SurfaceSearch: pendingSurfaceCap backpressure changes "
                 "timing only, never drops/double-counts a surface ---\n";

    // Same fixture as test_boundary_link_batch_recognizes_unknot: a single
    // pentachoron's boundary (S^3), which under --connected yields several
    // distinct satisfying surfaces -- enough to exceed a cap of 1.
    auto buildFourBall = [] {
        regina::Triangulation<4> fourBall;
        fourBall.newSimplex();
        return fourBall;
    };

    regina::Triangulation<4> triLarge = buildFourBall();
    SurfaceSearch large(triLarge);
    large.configureLimits(SurfaceSearchLimits{}); // defaults (cap 20000)
    SearchStats largeStats = large.search(1, BoundaryCondition::connected);

    regina::Triangulation<4> triTiny = buildFourBall();
    SurfaceSearch tiny(triTiny);
    SurfaceSearchLimits tinyLimits;
    tinyLimits.pendingSurfaceCap = 1; // forces backpressure on nearly every find
    tiny.configureLimits(tinyLimits);
    SearchStats tinyStats = tiny.search(1, BoundaryCondition::connected);

    EXPECT_EQ(tinyStats.satisfyingCount, largeStats.satisfyingCount,
              "a tiny pendingSurfaceCap finds the same number of satisfying "
              "surfaces as a large one");
    EXPECT_EQ(tinyStats.embeddedCount, largeStats.embeddedCount,
              "...and the same total embedded-submanifold count");
    EXPECT_EQ(tiny.linkTally().summary(), large.linkTally().summary(),
              "...and an identical final boundary-link tally");
    EXPECT_EQ(tiny.surfaceTypeTally().summary(),
              large.surfaceTypeTally().summary(),
              "...and an identical final surface-type tally");

    EXPECT_EQ(tiny.pendingSurfaceQueueStats().producerDrainEvents > 0, true,
              "the tiny cap actually triggered producer-side draining at "
              "least once (otherwise this test isn't exercising "
              "backpressure at all)");
}

// ─────────────────────────────────────────────────────────────────────────────
// ConnectedInducedSubgraphEnumerator's seeded mode, directly, on a small
// hand-built graph: 1 (seed) -- 3 -- 2, i.e. adj[1]={3}, adj[3]={1,2},
// adj[2]={3}. This is exactly the shape that would defeat a naive
// "just call enumerateFromRoot(w, visit) on each sibling" implementation:
// vertex 2 is only reachable through sibling root w=3, and 2 < 3, so an
// implementation that (incorrectly) re-anchors at w=3 would floor new
// candidates at "> 3" and never discover 2 at all. The correct
// implementation keeps the seed (vertex 1) as the anchor throughout, so 2
// must still be found.
// ─────────────────────────────────────────────────────────────────────────────
void test_seeded_enumerator_preserves_anchor() {
    std::cout << "\n--- ConnectedInducedSubgraphEnumerator: seeded mode "
                 "keeps the seed as anchor, not the sibling ---\n";

    std::vector<std::vector<int>> adj(4);
    adj[1] = {3};
    adj[3] = {1, 2};
    adj[2] = {3};

    ConnectedInducedSubgraphEnumerator enumerator(3, adj, /*isSeeded=*/true);
    EXPECT_EQ(enumerator.getRoots().size(), 1ULL,
              "vertex 1's only neighbor (3) is its only root");
    EXPECT_EQ(enumerator.getRoots()[0], 3,
              "the single root is vertex 3");

    std::vector<std::vector<int>> results;
    enumerator.enumerate(
        [&](const std::vector<int> &U) { results.push_back(U); });

    bool foundSeedPlusThree = false, foundAll = false;
    for (auto &U : results) {
        std::vector<int> sorted = U;
        std::sort(sorted.begin(), sorted.end());
        if (sorted == std::vector<int>{1, 3})
            foundSeedPlusThree = true;
        if (sorted == std::vector<int>{1, 2, 3})
            foundAll = true;
    }
    EXPECT_EQ(results.size(), 2ULL,
              "exactly 2 connected supersets of the seed are found");
    EXPECT_EQ(foundSeedPlusThree, true, "{1,3} is found");
    EXPECT_EQ(foundAll, true,
              "{1,2,3} is found -- vertex 2 is reachable via sibling 3 "
              "despite 2 < 3, because the seed (not the sibling) stays the "
              "anchor throughout");
}

// ─────────────────────────────────────────────────────────────────────────────
// EmbeddingSearch's seeded constructor, on the same B³ single tetrahedron as
// test_tetrahedron_boundary_conditions (K4 gluing graph among its 4
// triangle faces, so every non-empty subset is connected). Seeding with one
// face restricts results to connected subsets CONTAINING that face: exactly
// the 2^3 = 8 subsets of the other 3 faces union the seed. Of those, only
// the full 4-face set closes into S². This also exercises the "seed alone"
// special case in search() (the seed-only 1-face result can only be
// produced by that code path, never by any root's subtree).
// ─────────────────────────────────────────────────────────────────────────────
void test_seeded_search_tetrahedron() {
    std::cout << "\n--- B³ (single tetrahedron): seeded search from one "
                 "face ---\n";

    regina::Triangulation<3> ball;
    ball.newTetrahedron();
    std::vector<int> seed = {static_cast<int>(ball.triangle(0)->index())};

    EmbeddingSearch<3, 2> all(ball, seed);
    EXPECT_EQ(runFilteredCount(all, 1, BoundaryCondition::all), 8LL,
              "--all finds all 8 connected supersets of the seed face");

    EmbeddingSearch<3, 2> closed(ball, seed);
    EXPECT_EQ(runFilteredCount(closed, 1, BoundaryCondition::closed), 1LL,
              "--closed finds only the full 4-face set");

    EmbeddingSearch<3, 2> proper(ball, seed);
    EXPECT_EQ(runFilteredCount(proper, 1, BoundaryCondition::proper), 8LL,
              "--proper finds all 8 (every ambient edge is boundary)");

    EmbeddingSearch<3, 2> connected(ball, seed);
    EXPECT_EQ(runFilteredCount(connected, 2, BoundaryCondition::connected),
              8LL,
              "--connected finds all 8, using multiple threads over the "
              "seed's siblings");
}

// ─────────────────────────────────────────────────────────────────────────────
// A seed face excluded by buildGraph_'s pre-filter (the triple self-fold
// from test_triple_self_fold_excluded) must be rejected immediately by
// EmbeddingSearch's seeded constructor, not discovered later inside
// search().
// ─────────────────────────────────────────────────────────────────────────────
void test_seeded_search_rejects_invalid_seed() {
    std::cout << "\n--- Seeded EmbeddingSearch: invalid seed face throws at "
                 "construction ---\n";

    regina::Triangulation<3> tri;
    auto *t = tri.newTetrahedron();
    t->join(0, t, regina::Perm<4>(1, 2, 0, 3));
    t->join(2, t, regina::Perm<4>(0, 2, 3, 1));

    bool threw = false;
    try {
        EmbeddingSearch<3, 2> search(tri, {0});
    } catch (const regina::InvalidArgument &) {
        threw = true;
    }
    EXPECT_EQ(threw, true,
              "constructing with an excluded face as the seed throws "
              "regina::InvalidArgument");
}

// ─────────────────────────────────────────────────────────────────────────────
// DepthCappedPredicate, directly, on a hand-built path graph 1-2-3-4-5
// (adj[1]={2}, adj[2]={1,3}, adj[3]={2,4}, adj[4]={3,5}, adj[5]={4}).
// Wrapped around a trivial always-succeeding predicate, a cap of maxDepth
// should permit exactly the connected induced subgraphs of size <=
// maxDepth and reject any attempt to grow one larger -- this is the
// mechanism EmbeddingSearch::runSearch_'s iterative-deepening round loop
// relies on to guarantee every capped pass finishes quickly regardless of
// thread count (see enumerate_cis.h).
// ─────────────────────────────────────────────────────────────────────────────
namespace {
class AlwaysTruePredicate : public ConditionalPredicate {
  public:
    bool tryAdd(int) override { return true; }
    void undo(int) override {}
};
} // namespace

void test_depth_capped_predicate_bounds_depth() {
    std::cout
        << "\n--- DepthCappedPredicate: bounds enumerated subgraph size ---\n";

    std::vector<std::vector<int>> adj(6);
    adj[1] = {2};
    adj[2] = {1, 3};
    adj[3] = {2, 4};
    adj[4] = {3, 5};
    adj[5] = {4};

    ConnectedInducedSubgraphEnumerator enumerator(5, adj);

    AlwaysTruePredicate inner;
    DepthCappedPredicate capped(inner, 3);

    size_t maxSize = 0;
    size_t countAtCap = 0;
    enumerator.enumerateFiltered(
        [&](const std::vector<int> &U) {
            maxSize = std::max(maxSize, U.size());
            if (U.size() == 3)
                ++countAtCap;
        },
        capped);

    EXPECT_EQ(maxSize, 3ULL, "no reported subgraph exceeds the cap");
    // Connected induced subgraphs of size exactly 3 in this path: {1,2,3},
    // {2,3,4}, {3,4,5} -- 3 of them.
    EXPECT_EQ(countAtCap, 3ULL,
              "every size-3 connected subgraph is still found");
}

void test_depth_capped_predicate_clamps_to_one() {
    std::cout
        << "\n--- DepthCappedPredicate: clamps maxDepth to at least 1 ---\n";

    AlwaysTruePredicate inner;
    DepthCappedPredicate zeroCapped(inner, 0);
    EXPECT_EQ(zeroCapped.tryAdd(1), true,
              "a maxDepth <= 0 is clamped to 1, so the first tryAdd still "
              "succeeds -- required for seeded searches, whose seed commit "
              "must always be allowed through (see seedFastForward_)");
    EXPECT_EQ(zeroCapped.tryAdd(2), false,
              "the second tryAdd is rejected once the (clamped) cap of 1 is "
              "reached");
}

// ─────────────────────────────────────────────────────────────────────────────
// End-to-end iterative deepening (search()'s iddfsIterations/iddfsStep
// parameters, i.e. surfer's --iddfs-iterations/--iddfs-step): running the
// same B³ tetrahedron search (see test_tetrahedron_boundary_conditions,
// whose K4 gluing graph gives 15 connected embedded subsets) both without
// and with iterative deepening enabled, across several (iterations, step)
// combinations -- including a step so large the very first capped pass
// already covers everything, and a cap that lands exactly on the total
// face count -- must produce identical SearchStats in every case. The
// capped passes plus the suppression guard in runSearch_'s round loop must
// neither miss nor double-count anything, regardless of how many rounds
// run or where the round boundaries fall.
// ─────────────────────────────────────────────────────────────────────────────
void test_iddfs_matches_single_pass_unseeded() {
    std::cout << "\n--- IDDFS: unseeded results match a single unbounded "
                 "pass, for several (iterations, step) combinations ---\n";

    regina::Triangulation<3> ball;
    ball.newTetrahedron();

    EmbeddingSearch<3, 2> plain(ball);
    SearchStats plainStats = plain.search(1, BoundaryCondition::all);
    EXPECT_EQ(plainStats.satisfyingCount, 15LL,
              "sanity check: a single pass finds all 15 connected subsets");

    struct Combo {
        unsigned iterations;
        long long step;
        std::optional<long long> start = std::nullopt;
    };
    std::vector<Combo> combos = {
        {1, 1},        // cap=1: only single-face results in the capped pass
        {2, 1},        // caps 1, 2 (iddfsStart unset -- defaults to step)
        {1, 100},      // cap far exceeds every possible result -- the
                       // capped pass alone finds everything; the final
                       // pass must then find nothing new
        {3, 2},        // caps 2, 4, 6 -- the second cap (4) lands exactly
                       // on this graph's total face count
        {2, 3, 1LL},   // explicit iddfsStart=1, step=3: caps 1, 4 -- a
                       // different sequence than {2,3} (unset start) would
                       // give (caps 3, 6), exercising --iddfs-start
                       // specifically
        {2, 1, 1LL},   // explicit iddfsStart == iddfsStep: must match the
                       // unset-start case {2, 1} above exactly
    };

    for (const auto &combo : combos) {
        EmbeddingSearch<3, 2> iddfs(ball);
        SearchStats iddfsStats =
            iddfs.search(1, BoundaryCondition::all, {}, combo.iterations,
                        combo.step, combo.start);
        std::ostringstream desc;
        desc << "iterations=" << combo.iterations << " step=" << combo.step
             << " start=" << (combo.start ? std::to_string(*combo.start)
                                          : std::string("unset"));
        EXPECT_EQ(iddfsStats.satisfyingCount, plainStats.satisfyingCount,
                  "satisfyingCount matches (" + desc.str() + ")");
        EXPECT_EQ(iddfsStats.foundCount, plainStats.foundCount,
                  "foundCount matches (" + desc.str() + ")");
        EXPECT_EQ(iddfsStats.embeddedCount, plainStats.embeddedCount,
                  "embeddedCount matches (" + desc.str() + ")");
        EXPECT_EQ(iddfsStats.satisfyingFaceSum, plainStats.satisfyingFaceSum,
                  "satisfyingFaceSum matches (" + desc.str() + ")");
        EXPECT_EQ(iddfsStats.largestSatisfying, plainStats.largestSatisfying,
                  "largestSatisfying matches (" + desc.str() + ")");
    }
}

// As above, but seeded with 2 of the tetrahedron's 4 faces (a seed with more
// than one face), so the seeded added-face accounting in runSearch_'s
// worker (addedFaceCount = U.size() - 1, rather than just U.size()) is
// actually exercised, not just the single-face-seed case
// test_seeded_search_tetrahedron already covers.
void test_iddfs_matches_single_pass_seeded() {
    std::cout << "\n--- IDDFS: seeded results (seedFaceCount > 1) match a "
                 "single unbounded pass ---\n";

    regina::Triangulation<3> ball;
    ball.newTetrahedron();
    std::vector<int> seed = {static_cast<int>(ball.triangle(0)->index()),
                             static_cast<int>(ball.triangle(1)->index())};

    EmbeddingSearch<3, 2> plain(ball, seed);
    SearchStats plainStats = plain.search(1, BoundaryCondition::all);
    EXPECT_EQ(plainStats.satisfyingCount, 4LL,
              "sanity check: a 2-face seed out of 4 finds 2^(4-2) = 4 "
              "connected supersets");

    struct Combo {
        unsigned iterations;
        long long step;
    };
    std::vector<Combo> combos = {{1, 1}, {2, 1}, {1, 100}};

    for (const auto &combo : combos) {
        EmbeddingSearch<3, 2> iddfs(ball, seed);
        SearchStats iddfsStats = iddfs.search(
            2, BoundaryCondition::all, {}, combo.iterations, combo.step);
        std::ostringstream desc;
        desc << "iterations=" << combo.iterations << " step=" << combo.step;
        EXPECT_EQ(iddfsStats.satisfyingCount, plainStats.satisfyingCount,
                  "satisfyingCount matches (" + desc.str() + ")");
        EXPECT_EQ(iddfsStats.foundCount, plainStats.foundCount,
                  "foundCount matches (" + desc.str() + ")");
        EXPECT_EQ(iddfsStats.satisfyingFaceSum, plainStats.satisfyingFaceSum,
                  "satisfyingFaceSum matches (" + desc.str() + ")");
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SearchStats::iddfsRound/iddfsTotalRounds/iddfsCapped: the default
// (iddfsIterations == 0) must report exactly 1 (unbounded) round, matching
// this method's behavior from before iterative deepening existed; with
// iddfsIterations == N, the final SearchStats (as seen by
// onSearchComplete) must report round N+1 of N+1, uncapped -- the final
// pass, not one of the capped ones.
// ─────────────────────────────────────────────────────────────────────────────
void test_iddfs_search_stats_fields() {
    std::cout << "\n--- IDDFS: SearchStats round/cap fields ---\n";

    regina::Triangulation<3> ball;
    ball.newTetrahedron();

    {
        EmbeddingSearch<3, 2> plain(ball);
        SearchCallbacks callbacks;
        SearchStats finalStats;
        callbacks.onSearchComplete = [&](const SearchStats &s) {
            finalStats = s;
        };
        plain.search(1, BoundaryCondition::all, callbacks);
        EXPECT_EQ(finalStats.iddfsTotalRounds, 1U,
                  "default (iddfsIterations=0) reports exactly 1 round");
        EXPECT_EQ(finalStats.iddfsRound, 1U,
                  "default's final round is round 1");
        EXPECT_EQ(finalStats.iddfsCapped, false,
                  "default's only round is the unbounded one");
    }

    {
        EmbeddingSearch<3, 2> iddfs(ball);
        SearchCallbacks callbacks;
        SearchStats finalStats;
        callbacks.onSearchComplete = [&](const SearchStats &s) {
            finalStats = s;
        };
        iddfs.search(1, BoundaryCondition::all, callbacks,
                    /*iddfsIterations=*/3, /*iddfsStep=*/1);
        EXPECT_EQ(finalStats.iddfsTotalRounds, 4U,
                  "3 capped passes + 1 final pass = 4 total rounds");
        EXPECT_EQ(finalStats.iddfsRound, 4U,
                  "final stats report the last round (the final pass)");
        EXPECT_EQ(finalStats.iddfsCapped, false,
                  "final stats' round is the unbounded one, not capped");
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// iddfsCapForRound(): the pure formula runSearch_'s round loop uses to turn
// (iddfsStart, iddfsStep) into each capped pass's face-count cap. Tested
// directly (the exact function runSearch_ calls, not a re-derivation of
// it) since the aggregate SearchStats a full search produces are, by
// design, invariant to exactly which cap sequence was used -- correctness
// there doesn't depend on iddfsStart being wired through correctly, only
// on the cap sequence being monotonically increasing. This is what
// actually pins down that --iddfs-start controls the first pass's cap
// (not --iddfs-step, which only sets the increment after it).
// ─────────────────────────────────────────────────────────────────────────────
void test_iddfs_cap_for_round_formula() {
    std::cout << "\n--- iddfsCapForRound(): start + (iter - 1) * step ---\n";

    EXPECT_EQ(iddfsCapForRound(1, 5, 3), 5LL,
              "pass 1 is capped at exactly iddfsStart");
    EXPECT_EQ(iddfsCapForRound(2, 5, 3), 8LL,
              "pass 2 is capped at iddfsStart + 1 * iddfsStep");
    EXPECT_EQ(iddfsCapForRound(3, 5, 3), 11LL,
              "pass 3 is capped at iddfsStart + 2 * iddfsStep");
    EXPECT_EQ(iddfsCapForRound(1, 1, 100), 1LL,
              "a small iddfsStart caps pass 1 far below iddfsStep, "
              "distinguishing --iddfs-start from --iddfs-step");
}

// ─────────────────────────────────────────────────────────────────────────────
// iddfsMaxDepth(): the pure formula runSearch_'s worker uses to turn a
// round's face-count cap into DepthCappedPredicate's maxDepth. Tested
// directly for the same reason as iddfsCapForRound() above -- this pins
// down that the seeded case adds exactly +1 (for the seed's one-time
// tryAdd(), regardless of how many faces the seed itself contains), rather
// than an offset derived from the seed's size (the earlier, buggy formula
// this replaced -- see runSearch_'s worker before this fix, which computed
// `capFaces - seedFaceCount + 1` and so wasted every round whose capFaces
// was below the seed's face count).
// ─────────────────────────────────────────────────────────────────────────────
void test_iddfs_max_depth_formula() {
    std::cout << "\n--- iddfsMaxDepth(): +1 for a seed's one-time commit, "
                 "regardless of seed size ---\n";

    EXPECT_EQ(iddfsMaxDepth(5, false), 5LL,
              "unseeded: the cap applies directly, no offset");
    EXPECT_EQ(iddfsMaxDepth(5, true), 6LL,
              "seeded: +1 for the seed's single tryAdd(), independent of "
              "the seed's own face count");
    EXPECT_EQ(iddfsMaxDepth(1, true), 2LL,
              "even a tiny cap (1) still gets exactly +1 when seeded -- "
              "not clamped or scaled down by a large seed");
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
    run("test_backpressure_does_not_drop_or_double_count_surfaces",
        test_backpressure_does_not_drop_or_double_count_surfaces);
    run("test_boundary_link_batch_recognizes_unknot",
        test_boundary_link_batch_recognizes_unknot);
    run("test_seeded_enumerator_preserves_anchor",
        test_seeded_enumerator_preserves_anchor);
    run("test_seeded_search_tetrahedron", test_seeded_search_tetrahedron);
    run("test_seeded_search_rejects_invalid_seed",
        test_seeded_search_rejects_invalid_seed);
    run("test_depth_capped_predicate_bounds_depth",
        test_depth_capped_predicate_bounds_depth);
    run("test_depth_capped_predicate_clamps_to_one",
        test_depth_capped_predicate_clamps_to_one);
    run("test_iddfs_matches_single_pass_unseeded",
        test_iddfs_matches_single_pass_unseeded);
    run("test_iddfs_matches_single_pass_seeded",
        test_iddfs_matches_single_pass_seeded);
    run("test_iddfs_search_stats_fields", test_iddfs_search_stats_fields);
    run("test_iddfs_cap_for_round_formula", test_iddfs_cap_for_round_formula);
    run("test_iddfs_max_depth_formula", test_iddfs_max_depth_formula);

    std::cout << "\n"
              << bold << (failed_count > 0 ? red : green) << "=== " << passed
              << " passed, " << failed_count << " failed ===" << resetColor
              << "\n";
    return failed_count > 0 ? 1 : 0;
}
