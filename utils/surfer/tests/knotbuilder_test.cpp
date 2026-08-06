// knotbuilder_test.cpp
// Tests for knotbuilder's PD-code -> triangulated-S³ pipeline (Block,
// Block::glue(), buildLink()), and its integration with CobordismBuilder's
// thicken()/cone(). See cobordismbuilder_test.cpp for CobordismBuilder's own
// tests independent of knotbuilder.

#include <iostream>
#include <link/link.h>
#include <string>
#include <triangulation/dim3.h>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>

#include "cobordismbuilder.h"
#include "knotbuilder.h"
#include "linkcomplement.h"

static int passed = 0, failed_count = 0;

// ANSI color manipulators for PASS/FAIL output. Auto-disabled when stdout
// isn't a terminal (e.g. piped to a file), so logs don't fill up with
// escape codes.
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

// ═════════════════════════════════════════════════════════════════════════
// Knotbuilder integration.
//
// PD codes below use the standard KnotInfo/Knot Atlas convention: one
// X(a,b,c,d) per crossing (1-indexed), flattened into a single list.
// They are not just trusted: each is checked to actually produce a valid,
// closed triangulation of S³ (isSphere()), and Link is checked to find the
// expected number of components from the returned link edges, before being
// used in any CobordismBuilder pipeline test.
// ═════════════════════════════════════════════════════════════════════════
const char *TREFOIL_PD = "1 4 2 5 3 6 4 1 5 2 6 3"; // 3_1
const char *HOPF_LINK_PD = "1 4 2 3 3 2 4 1";       // 2-component link
const char *FIGURE_EIGHT_PD = "4 2 5 1 8 6 1 5 6 3 7 4 2 7 3 8"; // 4_1

regina::Triangulation<3>
buildFromPD(const char *pd, std::vector<const regina::Edge<3> *> &edges) {
    knotbuilder::PDCode pdcode = knotbuilder::parsePDCode(pd);
    // Build into a named local and std::move() it out explicitly (rather
    // than destructuring and returning the destructured triangulation),
    // since `edges` must keep pointing into whichever Triangulation<3>
    // object ends up owning the simplices -- an explicit move guarantees
    // that, whereas returning a structured binding by value is not
    // guaranteed to elide/move rather than copy.
    auto result = knotbuilder::buildLink(pdcode);
    edges = std::move(result.edges);
    return std::move(result.tri);
}

// ─────────────────────────────────────────────────────────────────────────────
// Independent structural check on the raw edge list returned by buildLink(),
// deliberately *not* reusing Link's own tip-to-tail grouping logic (that's
// exactly what this is meant to double-check). Verifies:
//   1. no duplicate edges;
//   2. every vertex touched by the edge set has degree exactly 2 within that
//      set (necessary for the edges to be a disjoint union of simple closed
//      curves, rather than dangling ends or a branch point);
//   3. tracing tip-to-tail from any edge always returns to its own starting
//      vertex (each component is a genuine closed loop, not an open path
//      that Link's greedy grouping would otherwise silently accept); and
//   4. the edges partition into exactly the expected number of such loops,
//      consuming every edge exactly once.
// This says nothing about which knot/link it is (that's what the homology/
// census/isomorphism checks are for) -- only that getLinkEdges() produced
// something combinatorially sane for Link to work with in the first place.
// ─────────────────────────────────────────────────────────────────────────────
void checkEdgesFormSimpleClosedLoops(
    const std::vector<const regina::Edge<3> *> &edges, int expectedComponents,
    const std::string &name) {
    std::unordered_set<const regina::Edge<3> *> edgeSet(edges.begin(),
                                                        edges.end());
    EXPECT_EQ(edgeSet.size(), edges.size(),
              name + ": no duplicate edges in the link edge list");

    std::unordered_map<const regina::Vertex<3> *, int> degree;
    for (const regina::Edge<3> *e : edges) {
        ++degree[e->vertex(0)];
        ++degree[e->vertex(1)];
    }
    bool allDegreeTwo = true;
    for (auto &[v, d] : degree) {
        if (d != 2)
            allDegreeTwo = false;
    }
    EXPECT_EQ(allDegreeTwo, true,
              name + ": every vertex touched by the link edges has degree "
                     "exactly 2 (no dangling ends or branch points)");

    std::unordered_set<const regina::Edge<3> *> remaining(edges.begin(),
                                                          edges.end());
    int componentsFound = 0;
    bool allClosed = true;
    while (!remaining.empty()) {
        const regina::Edge<3> *start = *remaining.begin();
        const regina::Vertex<3> *startVertex = start->vertex(0);
        const regina::Vertex<3> *currVertex = start->vertex(1);
        remaining.erase(start);

        while (currVertex != startVertex) {
            const regina::Edge<3> *next = nullptr;
            for (const regina::Edge<3> *e : remaining) {
                if (e->vertex(0) == currVertex || e->vertex(1) == currVertex) {
                    next = e;
                    break;
                }
            }
            if (next == nullptr)
                break; // dangling: this component never closes up

            currVertex = (next->vertex(0) == currVertex) ? next->vertex(1)
                                                         : next->vertex(0);
            remaining.erase(next);
        }

        if (currVertex != startVertex)
            allClosed = false;
        ++componentsFound;
    }

    EXPECT_EQ(allClosed, true,
              name + ": every component closes up into a single loop "
                     "(tip-to-tail trace returns to its own start)");
    EXPECT_EQ(componentsFound, expectedComponents,
              name + ": edges trace exactly the expected number of simple "
                     "closed loops");
}

// ─────────────────────────────────────────────────────────────────────────────
// Directed analogue of checkEdgesFormSimpleClosedLoops(), driven by
// buildLink()'s new `reversed` data: `edges[i]` runs vertex(1)->vertex(0)
// iff `reversed[i]`. This is the diagnostic §2 of the orientation-tracking
// plan calls for -- confirming the tagged directions produce one single,
// genuinely closed, *consistently-signed* traversal per component, not
// just that the code compiles and runs. If the wall/vertex correspondences
// baked into Block::getLinkEdgeDirections(), or the PD-code walk in
// computeStrandDirections(), were wrong -- even for just the understrand or
// just the overstrand -- this is expected to fail sharply: a vertex would
// end up with in-degree/out-degree other than exactly 1, since a genuine
// per-crossing direction error breaks the head-to-tail chaining that a
// purely undirected check (checkEdgesFormSimpleClosedLoops) can't see at
// all.
// ─────────────────────────────────────────────────────────────────────────────
void checkEdgesFormDirectedClosedLoops(
    const std::vector<const regina::Edge<3> *> &edges,
    const std::vector<bool> &reversed, int expectedComponents,
    const std::string &name) {
    EXPECT_EQ(reversed.size(), edges.size(),
              name + ": reversed has one entry per edge");

    std::unordered_map<const regina::Edge<3> *, bool> reversedOf;
    std::unordered_map<const regina::Vertex<3> *, const regina::Edge<3> *>
        outEdgeFrom;
    std::unordered_map<const regina::Vertex<3> *, int> inDegree, outDegree;
    for (size_t i = 0; i < edges.size(); ++i) {
        reversedOf[edges[i]] = reversed[i];
        const regina::Vertex<3> *tail =
            reversed[i] ? edges[i]->vertex(1) : edges[i]->vertex(0);
        const regina::Vertex<3> *head =
            reversed[i] ? edges[i]->vertex(0) : edges[i]->vertex(1);
        outEdgeFrom[tail] = edges[i];
        ++outDegree[tail];
        ++inDegree[head];
    }

    bool allDegreeOne = true;
    for (auto &[v, d] : outDegree)
        if (d != 1)
            allDegreeOne = false;
    for (auto &[v, d] : inDegree)
        if (d != 1)
            allDegreeOne = false;
    EXPECT_EQ(allDegreeOne, true,
              name + ": every vertex has exactly one outgoing and one "
                     "incoming directed edge (a per-crossing direction "
                     "error would break this)");

    std::unordered_set<const regina::Edge<3> *> visited;
    int componentsFound = 0;
    bool allClosed = true;
    for (const regina::Edge<3> *startEdge : edges) {
        if (visited.contains(startEdge))
            continue;

        const regina::Vertex<3> *startVertex =
            reversedOf[startEdge] ? startEdge->vertex(1) : startEdge->vertex(0);
        const regina::Vertex<3> *currVertex = startVertex;
        bool closed = false;
        for (size_t step = 0; step <= edges.size(); ++step) {
            auto it = outEdgeFrom.find(currVertex);
            if (it == outEdgeFrom.end())
                break; // dangling: shouldn't happen given the degree check
                       // above, but don't trust that blindly here.
            const regina::Edge<3> *next = it->second;
            visited.insert(next);
            currVertex =
                reversedOf[next] ? next->vertex(0) : next->vertex(1);
            if (currVertex == startVertex) {
                closed = true;
                break;
            }
        }
        if (!closed)
            allClosed = false;
        ++componentsFound;
    }

    EXPECT_EQ(allClosed, true,
              name + ": every component closes up into a single directed "
                     "loop (head-to-tail trace returns to its own start)");
    EXPECT_EQ(componentsFound, expectedComponents,
              name + ": directed edges trace exactly the expected number "
                     "of components");
}

// ─────────────────────────────────────────────────────────────────────────────
// Layer 1: knotbuilder's raw output, on its own, independent of
// CobordismBuilder. If buildLink() ever stops producing a valid closed S³
// (e.g. a knotbuilder change breaks the block gluing pattern), this should
// fail here rather than surfacing as a confusing CobordismBuilder/glue()
// failure downstream.
// ─────────────────────────────────────────────────────────────────────────────
void test_knotbuilder_trefoil_is_valid_s3() {
    std::cout << "\n--- knotbuilder: trefoil PD code -> closed S³ ---\n";

    std::vector<const regina::Edge<3> *> edges;
    auto tri = buildFromPD(TREFOIL_PD, edges);

    EXPECT_EQ(tri.isValid(), true, "trefoil triangulation is valid");
    EXPECT_EQ(tri.isClosed(), true, "trefoil triangulation is closed");
    EXPECT_EQ(tri.isIdeal(), false,
              "trefoil triangulation has no ideal vertices");
    EXPECT_EQ(tri.isSphere(), true, "trefoil triangulation is S³");
    EXPECT_EQ((int)edges.size(), 9, "3 crossings × 3 link edges per block = 9");

    checkEdgesFormSimpleClosedLoops(edges, 1, "trefoil");

    Link link(tri, edges);
    EXPECT_EQ((int)link.comps_.size(), 1, "trefoil is a 1-component knot");
}

void test_knotbuilder_hopf_link_two_components() {
    std::cout
        << "\n--- knotbuilder: Hopf link PD code -> 2-component link ---\n";

    std::vector<const regina::Edge<3> *> edges;
    auto tri = buildFromPD(HOPF_LINK_PD, edges);

    EXPECT_EQ(tri.isValid(), true, "Hopf link triangulation is valid");
    EXPECT_EQ(tri.isSphere(), true, "Hopf link triangulation is S³");

    checkEdgesFormSimpleClosedLoops(edges, 2, "Hopf link");

    Link link(tri, edges);
    EXPECT_EQ((int)link.comps_.size(), 2, "Hopf link has 2 components");
}

// ─────────────────────────────────────────────────────────────────────────────
// The orientation-tracking plan's own dedicated diagnostic: confirms
// buildLink()'s new `reversed` data produces one single, genuinely closed,
// consistently-signed directed traversal per component -- across a battery
// of known knots/links, including L6a3{0} and L6a3{1} (the exact pair the
// fix exists for). PD strings taken directly from
// links_4d_smooth_slice_genus_11_crossings_pd_codes.csv; parsePDCode()
// strips all non-digit characters, so the "PD[X[...]; ...]" KnotInfo
// notation can be fed in as-is.
// ─────────────────────────────────────────────────────────────────────────────
const char *L6A3_0_PD = "PD[X[8;1;9;2];X[2;9;3;10];X[10;3;11;4];X[12;5;7;6];"
                       "X[6;7;1;8];X[4;11;5;12]]";
const char *L6A3_1_PD = "PD[X[10;2;11;1];X[2;10;3;9];X[8;4;9;3];X[12;6;7;5];"
                       "X[6;12;1;11];X[4;8;5;7]]";

void checkDirectedTraversal(const char *pd, int expectedComponents,
                            const std::string &name) {
    knotbuilder::PDCode pdcode = knotbuilder::parsePDCode(pd);
    auto result = knotbuilder::buildLink(pdcode);

    EXPECT_EQ(result.tri.isValid(), true, name + ": triangulation is valid");
    EXPECT_EQ(result.tri.isSphere(), true, name + ": triangulation is S³");

    checkEdgesFormDirectedClosedLoops(result.edges, result.reversed,
                                      expectedComponents, name);
}

void test_knotbuilder_reversed_directions_are_consistent() {
    std::cout << "\n--- knotbuilder: reversed[] forms one consistently-"
                 "signed directed traversal per component ---\n";

    checkDirectedTraversal(TREFOIL_PD, 1, "trefoil");
    checkDirectedTraversal(HOPF_LINK_PD, 2, "Hopf link");
    checkDirectedTraversal(FIGURE_EIGHT_PD, 1, "figure-8");
    checkDirectedTraversal(L6A3_0_PD, 2, "L6a3{0}");
    checkDirectedTraversal(L6A3_1_PD, 2, "L6a3{1}");
}

// ─────────────────────────────────────────────────────────────────────────────
// L6a3{0} and L6a3{1} are the same diagram with one component's orientation
// reversed -- the exact pair that motivated this feature (see the
// orientation-tracking plan's Context). Their complements are identical
// (identify() names them the same, "L206001"), so a complement-based check
// can never distinguish them; this confirms buildLink()'s new `reversed`
// tagging, driven purely by each PD code's own strand directions, actually
// does.
// ─────────────────────────────────────────────────────────────────────────────
void test_knotbuilder_l6a3_variants_tag_differently() {
    std::cout << "\n--- knotbuilder: L6a3{0} and L6a3{1} produce genuinely "
                 "different reversed[] tags ---\n";

    auto result0 = knotbuilder::buildLink(knotbuilder::parsePDCode(L6A3_0_PD));
    auto result1 = knotbuilder::buildLink(knotbuilder::parsePDCode(L6A3_1_PD));

    EXPECT_EQ(result0.reversed.size(), result1.reversed.size(),
              "both PD codes have the same crossing count, so the same "
              "number of tagged edges");
    EXPECT_EQ(result0.reversed != result1.reversed, true,
              "L6a3{0}'s and L6a3{1}'s reversed[] sequences differ -- "
              "reversing one component's orientation in the PD code "
              "produces a genuinely different tagging under this one fixed "
              "convention, exactly as the underlying links themselves "
              "differ despite sharing a complement");
}

// ─────────────────────────────────────────────────────────────────────────────
// knotbuilder::reduceVertices(): pinches away every internal edge except the
// preserved knot/link edges and loop edges, so a caller can shrink
// buildLink()'s (many-tetrahedra-per-crossing) output down while continuing
// to work with the same knot/link edges. Checks (for both a knot and a
// link, so a mistake specific to component count isn't missed):
//   1. the result is still a valid, closed, non-ideal S³;
//   2. it actually has no more vertices than the input (the stated purpose);
//   3. the returned edge list still corresponds 1-1 to the input, and still
//      traces the same knot/link diagram (checkEdgesFormSimpleClosedLoops +
//      Link's own component count);
//   4. every edge NOT in that returned list is a loop (the pinch loop
//      reached a genuine fixpoint, not stopping early); and
//   5. calling it again on its own output is a safe, idempotent no-op.
// ─────────────────────────────────────────────────────────────────────────────
void checkReduceVertices(const char *name, const char *pd,
                          int expectedComponents) {
    std::vector<const regina::Edge<3> *> edges;
    auto tri = buildFromPD(pd, edges);

    auto reduced = knotbuilder::reduceVertices(tri, edges);

    EXPECT_EQ(reduced.tri.isValid(), true,
              std::string(name) + ": reduced triangulation is valid");
    EXPECT_EQ(reduced.tri.isClosed(), true,
              std::string(name) + ": reduced triangulation is closed");
    EXPECT_EQ(reduced.tri.isIdeal(), false,
              std::string(name) +
                  ": reduced triangulation has no ideal vertices");
    EXPECT_EQ(reduced.tri.isSphere(), true,
              std::string(name) + ": reduced triangulation is still S³");
    EXPECT_EQ(reduced.tri.countVertices() <= tri.countVertices(), true,
              std::string(name) +
                  ": reduceVertices() does not increase vertex count");

    EXPECT_EQ(reduced.edges.size(), edges.size(),
              std::string(name) +
                  ": reduced edge list is the same size as the input");
    checkEdgesFormSimpleClosedLoops(reduced.edges, expectedComponents,
                                    std::string(name) + " (reduced)");

    Link reducedLink(reduced.tri, reduced.edges);
    EXPECT_EQ((int)reducedLink.comps_.size(), expectedComponents,
              std::string(name) +
                  ": reduced triangulation's edges still trace the same "
                  "number of components");

    // Fixpoint check: reduceVertices() must never pinch an edge incident to
    // a preserved edge's endpoint (doing so could turn the preserved edge
    // into a loop, or otherwise disturb the knot/link it traces -- see the
    // "protected vertices" logic in knotbuilder.cpp). So the loop is only
    // guaranteed to reach a fixpoint where every remaining non-preserved
    // edge is *either* a loop *or* touches a preserved edge's vertex --
    // not necessarily a loop on its own.
    std::unordered_set<const regina::Edge<3> *> preserved(
        reduced.edges.begin(), reduced.edges.end());
    std::unordered_set<const regina::Vertex<3> *> preservedVertices;
    for (const regina::Edge<3> *e : reduced.edges) {
        preservedVertices.insert(e->vertex(0));
        preservedVertices.insert(e->vertex(1));
    }
    bool fixpointReached = true;
    for (const regina::Edge<3> *e : reduced.tri.edges()) {
        if (preserved.contains(e))
            continue;
        bool isLoop = e->vertex(0) == e->vertex(1);
        bool touchesPreservedVertex =
            preservedVertices.contains(e->vertex(0)) ||
            preservedVertices.contains(e->vertex(1));
        if (!isLoop && !touchesPreservedVertex)
            fixpointReached = false;
    }
    EXPECT_EQ(fixpointReached, true,
              std::string(name) +
                  ": every edge other than the preserved ones is a loop or "
                  "touches a preserved vertex (the pinch loop reached a "
                  "fixpoint)");

    auto reducedAgain =
        knotbuilder::reduceVertices(reduced.tri, reduced.edges);
    EXPECT_EQ(reducedAgain.tri.countVertices(), reduced.tri.countVertices(),
              std::string(name) +
                  ": reduceVertices() is idempotent on its own output");
}

void test_knotbuilder_reduce_vertices() {
    std::cout << "\n--- knotbuilder: reduceVertices() shrinks vertex count, "
                 "preserves link ---\n";
    checkReduceVertices("trefoil", TREFOIL_PD, 1);
    checkReduceVertices("Hopf link", HOPF_LINK_PD, 2);
}

// ─────────────────────────────────────────────────────────────────────────────
// Regression test: this specific diagram is what originally exposed the bug
// in Block::glue(). It's the Hopf link's shadow with crossing 1's tuple
// rotated by one position, which flips that crossing's over/under role and
// makes the diagram non-alternating (both transitions become same-parity
// instead of alternating). The rotated diagram is the 2-component unlink,
// not the Hopf link.
//
// The old parity-based glue() produced a triangulation that looked fine
// (isValid() == true, no non-manifold points) but was structurally wrong:
// the two S²×I boundary components came out as tori instead of spheres, so
// finiteToIdeal() turned them into genuine cusps (isIdeal() == true) instead
// of closing up to S³. isSphere() == false was the tell. The current
// glue() (matching walls by role -- own/anchor/private -- rather than by
// wall-index parity) has no such branch to get wrong, but this is kept as a
// dedicated, minimal, always-run check specifically for this failure mode.
// ─────────────────────────────────────────────────────────────────────────────
void test_knotbuilder_nonalternating_regression() {
    std::cout << "\n--- knotbuilder: non-alternating regression (flipped "
                 "Hopf shadow -> unlink) ---\n";

    // Same 4-valent shadow as the Hopf link (crossing 0 unchanged), but
    // crossing 1's tuple [2,1,3,0] is cyclically rotated by one position to
    // [1,3,0,2], flipping which pair of arms is the under-strand.
    knotbuilder::PDCode pd = {{0, 3, 1, 2}, {1, 3, 0, 2}};

    auto [tri, edges, reversed] = knotbuilder::buildLink(pd);

    EXPECT_EQ(tri.isValid(), true,
              "non-alternating shadow: triangulation is valid");
    EXPECT_EQ(tri.isClosed(), true,
              "non-alternating shadow: triangulation is closed");
    EXPECT_EQ(tri.isIdeal(), false,
              "non-alternating shadow: no ideal vertices (an ideal vertex "
              "here would mean a torus boundary wasn't closed into a "
              "sphere before capping)");
    EXPECT_EQ(tri.isSphere(), true,
              "non-alternating shadow: triangulation is S³");

    checkEdgesFormSimpleClosedLoops(edges, 2, "non-alternating shadow");

    Link link(tri, edges);
    EXPECT_EQ((int)link.comps_.size(), 2,
              "non-alternating shadow: 2 components (the unlink)");
}

// ─────────────────────────────────────────────────────────────────────────────
// Broader robustness sweep: real PD codes for a range of named knots (3 to
// 13 crossings, both alternating and non-alternating -- the "a"/"n" suffix
// in the 11+ crossing names is KnotInfo's own alternating/non-alternating
// classification), from KnotInfo (https://knotinfo.org/). Each is checked
// for a valid, closed S³ triangulation with the correct single-component
// count, and its complement's homology is cross-checked against Regina's
// own independent Link::fromPD() construction of the same PD code -- a
// cheap (no exponential search needed) but genuine correctness signal
// beyond "didn't crash", especially for the larger examples where an exact
// isomorphism/retriangulation check would be too slow to run routinely.
//
// (This is a curated sample. utils/surfer/pd_codes.csv has PD codes
// for every knot up to 13 crossings -- all 12,467 of them pass the same
// valid/closed/sphere checks, verified separately as a one-off sweep; that
// full run is too slow (~13 minutes) to bake into the routine suite.)
// ─────────────────────────────────────────────────────────────────────────────
struct NamedKnot {
    const char *name;
    const char *pd;
};

const NamedKnot NAMED_KNOTS[] = {
    {"3_1", "1,5,2,4 3,1,4,6 5,3,6,2"},
    {"4_1", "4,2,5,1 8,6,1,5 6,3,7,4 2,7,3,8"},
    {"5_1", "2,8,3,7 4,10,5,9 6,2,7,1 8,4,9,3 10,6,1,5"},
    {"5_2", "1,5,2,4 3,9,4,8 5,1,6,10 7,3,8,2 9,7,10,6"},
    {"6_1", "1,7,2,6 3,10,4,11 5,3,6,2 7,1,8,12 9,4,10,5 11,9,12,8"},
    {"6_2", "1,8,2,9 3,11,4,10 5,1,6,12 7,2,8,3 9,7,10,6 11,5,12,4"},
    {"6_3", "4,2,5,1 8,4,9,3 12,9,1,10 10,5,11,6 6,11,7,12 2,8,3,7"},
    {"7_4",
     "2,10,3,9 4,12,5,11 6,14,7,13 8,4,9,3 10,2,11,1 12,8,13,7 14,6,1,5"},
    {"8_19 (non-alternating)",
     "2,14,3,13 5,11,6,10 7,15,8,14 9,5,10,4 11,7,12,6 12,2,13,1 15,9,16,8 "
     "16,4,1,3"},
    {"8_20 (non-alternating)",
     "1,7,2,6 4,13,5,14 5,9,6,8 7,3,8,2 10,15,11,16 12,9,13,10 14,3,15,4 "
     "16,11,1,12"},
    {"8_21 (non-alternating)",
     "1,7,2,6 4,13,5,14 5,9,6,8 7,3,8,2 9,13,10,12 11,1,12,16 14,3,15,4 "
     "15,11,16,10"},
    {"10_165",
     "1,8,2,9 3,12,4,13 4,17,5,18 7,2,8,3 9,14,10,15 11,17,12,16 13,6,14,7 "
     "15,20,16,1 18,5,19,6 19,11,20,10"},
    {"11a_1",
     "4,2,5,1 10,6,11,5 8,3,9,4 2,9,3,10 16,12,17,11 14,7,15,8 6,15,7,16 "
     "20,14,21,13 22,18,1,17 18,22,19,21 12,20,13,19"},
    {"11n_1 (non-alternating)",
     "4,2,5,1 10,6,11,5 8,3,9,4 2,9,3,10 11,16,12,17 14,7,15,8 6,15,7,16 "
     "13,20,14,21 17,22,18,1 21,18,22,19 19,12,20,13"},
    {"12a_1",
     "1,5,2,4 3,8,4,9 5,11,6,10 7,14,8,15 9,2,10,3 11,17,12,16 13,21,14,20 "
     "15,6,16,7 17,22,18,23 19,13,20,12 21,24,22,1 23,18,24,19"},
    {"12n_1 (non-alternating)",
     "1,5,2,4 3,8,4,9 5,11,6,10 14,8,15,7 9,2,10,3 16,11,17,12 20,13,21,14 "
     "6,16,7,15 22,18,23,17 12,19,13,20 24,22,1,21 18,24,19,23"},
    {"13a_1",
     "3,1,4,26 1,8,2,9 7,2,8,3 9,5,10,4 5,14,6,15 13,6,14,7 15,11,16,10 "
     "11,20,12,21 19,12,20,13 21,17,22,16 17,24,18,25 23,18,24,19 "
     "25,23,26,22"},
    {"13n_1 (non-alternating)",
     "3,1,4,26 1,8,2,9 7,2,8,3 9,5,10,4 14,6,15,5 6,14,7,13 10,15,11,16 "
     "20,12,21,11 12,20,13,19 16,21,17,22 24,18,25,17 18,24,19,23 "
     "22,25,23,26"},
};

void test_knotbuilder_many_named_knots() {
    std::cout << "\n--- knotbuilder: sweep over named knots (3-13 "
                 "crossings, alternating and non-alternating), from "
                 "KnotInfo ---\n";

    for (const auto &knot : NAMED_KNOTS) {
        std::vector<const regina::Edge<3> *> edges;
        auto tri = buildFromPD(knot.pd, edges);

        EXPECT_EQ(tri.isValid(), true,
                  std::string(knot.name) + ": triangulation is valid");
        EXPECT_EQ(tri.isClosed(), true,
                  std::string(knot.name) + ": triangulation is closed");
        EXPECT_EQ(tri.isSphere(), true,
                  std::string(knot.name) + ": triangulation is S³");

        checkEdgesFormSimpleClosedLoops(edges, 1, knot.name);

        Link link(tri, edges);
        EXPECT_EQ((int)link.comps_.size(), 1,
                  std::string(knot.name) + ": is a 1-component knot");

        // Cross-check against Regina's own independent Link::fromPD()
        // pipeline for the same PD code: the two complements' homology
        // (H_1 of a knot complement is always Z, so this mainly catches
        // gross errors, but it's cheap enough to run on every case,
        // including the larger ones where an exact isomorphism check via
        // retriangulation would be too slow).
        auto mine = link.buildComplement();
        regina::Link reginaLink = regina::Link::fromPD(knot.pd);
        auto reginaComp = reginaLink.complement();
        EXPECT_EQ(mine.homology() == reginaComp.homology(), true,
                  std::string(knot.name) +
                      ": complement homology matches Regina's independent "
                      "Link::fromPD() construction");
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Layer 2: the specific property the user asked to check. knotbuilder's raw
// output is *not* ordered (block construction doesn't aim for that), but it
// must be orderable, since that's what CobordismBuilder<3>'s constructor
// relies on (it calls Triangulation<3>::order() itself when needed — see
// cobordismbuilder.h). Recording the raw isOrdered() value directly here,
// rather than assuming it, means a future knotbuilder change that breaks
// orderability entirely (as opposed to just not ordering by default) gets
// caught right here instead of as a mysterious CobordismBuilder exception.
// ─────────────────────────────────────────────────────────────────────────────
void test_knotbuilder_output_is_orderable() {
    std::cout << "\n--- knotbuilder output: raw ordering + orderability ---\n";

    const std::pair<const char *, const char *> knots[] = {
        {"trefoil", TREFOIL_PD},
        {"Hopf link", HOPF_LINK_PD},
        {"figure-8", FIGURE_EIGHT_PD},
    };

    for (const auto &[name, pd] : knots) {
        std::vector<const regina::Edge<3> *> edges;
        auto tri = buildFromPD(pd, edges);

        std::cout << "  " << name << ": raw isOrdered() = " << tri.isOrdered()
                  << "\n";

        regina::Triangulation<3> ordered(tri);
        EXPECT_EQ(ordered.order(), true,
                  std::string(name) + ": knotbuilder's output is orderable");
        EXPECT_EQ(ordered.isOrdered(), true,
                  std::string(name) +
                      ": order() actually produces an ordered triangulation");

        // The property that matters end to end: CobordismBuilder<3>'s
        // constructor requires isOrdered() and calls order() itself, so it
        // must accept knotbuilder's raw (unordered) output without throwing.
        bool constructedOk = true;
        try {
            CobordismBuilder<3> cob(tri);
        } catch (const std::exception &) {
            constructedOk = false;
        }
        EXPECT_EQ(constructedOk, true,
                  std::string(name) +
                      ": CobordismBuilder<3> accepts knotbuilder's raw output");
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Layer 3: the full pipeline requested — knotbuilder builds S³ with the
// link as an edge circuit, CobordismBuilder thickens it (giving room for
// surfer to later search in), then cone() caps the top into a
// genuine triangulated 4-ball. The one boundary component left afterward
// must be combinatorially identical to knotbuilder's original S³ — that's
// what will let a future surfer test locate the link's edges inside
// the boundary of this B⁴.
// ─────────────────────────────────────────────────────────────────────────────
void checkKnotToBallPipeline(const char *name, const char *pd, int layers) {
    std::vector<const regina::Edge<3> *> edges;
    auto tri = buildFromPD(pd, edges);

    CobordismBuilder<3> cob(tri);
    if (layers > 0)
        cob.thicken(layers);
    auto &result = cob.cone();

    EXPECT_EQ(result.isValid(), true,
              std::string(name) + ": thicken(" + std::to_string(layers) +
                  ")+cone() triangulation is valid");
    EXPECT_EQ((int)result.countBoundaryComponents(), 1,
              std::string(name) + ": thicken(" + std::to_string(layers) +
                  ")+cone() has exactly 1 boundary component");
    EXPECT_EQ(result.eulerCharManifold(), 1L,
              std::string(name) + ": thicken(" + std::to_string(layers) +
                  ")+cone() is contractible (Euler characteristic 1)");
    EXPECT_EQ(
        result.boundaryComponent(0)->build().isIsomorphicTo(tri).has_value(),
        true,
        std::string(name) + ": remaining boundary ≅ knotbuilder's original S³");
}

void test_knotbuilder_trefoil_to_ball_pipeline() {
    std::cout << "\n--- knotbuilder trefoil -> thicken(n)+cone() -> B⁴ ---\n";
    checkKnotToBallPipeline("trefoil (0 layers)", TREFOIL_PD, 0);
    checkKnotToBallPipeline("trefoil (1 layer)", TREFOIL_PD, 1);
    checkKnotToBallPipeline("trefoil (3 layers)", TREFOIL_PD, 3);
}

void test_knotbuilder_hopf_link_to_ball_pipeline() {
    std::cout << "\n--- knotbuilder Hopf link -> thicken(n)+cone() -> B⁴ ---\n";
    checkKnotToBallPipeline("Hopf link (1 layer)", HOPF_LINK_PD, 1);
    checkKnotToBallPipeline("Hopf link (2 layers)", HOPF_LINK_PD, 2);
}

void test_knotbuilder_figure_eight_to_ball_pipeline() {
    std::cout << "\n--- knotbuilder figure-8 -> thicken(n)+cone() -> B⁴ ---\n";
    checkKnotToBallPipeline("figure-8 (1 layer)", FIGURE_EIGHT_PD, 1);
    checkKnotToBallPipeline("figure-8 (2 layers)", FIGURE_EIGHT_PD, 2);
}

template <typename F>
void run(const char *name, F fn) {
    std::cout << "\nRunning " << name << "...\n";
    try {
        fn();
    } catch (const std::exception &e) {
        std::cout << red << "  EXCEPTION: " << e.what() << resetColor << "\n";
        ++failed_count;
    }
}

int main() {
    run("test_knotbuilder_trefoil_is_valid_s3",
        test_knotbuilder_trefoil_is_valid_s3);
    run("test_knotbuilder_hopf_link_two_components",
        test_knotbuilder_hopf_link_two_components);
    run("test_knotbuilder_reversed_directions_are_consistent",
        test_knotbuilder_reversed_directions_are_consistent);
    run("test_knotbuilder_l6a3_variants_tag_differently",
        test_knotbuilder_l6a3_variants_tag_differently);
    run("test_knotbuilder_reduce_vertices", test_knotbuilder_reduce_vertices);
    run("test_knotbuilder_nonalternating_regression",
        test_knotbuilder_nonalternating_regression);
    run("test_knotbuilder_many_named_knots", test_knotbuilder_many_named_knots);
    run("test_knotbuilder_output_is_orderable",
        test_knotbuilder_output_is_orderable);
    run("test_knotbuilder_trefoil_to_ball_pipeline",
        test_knotbuilder_trefoil_to_ball_pipeline);
    run("test_knotbuilder_hopf_link_to_ball_pipeline",
        test_knotbuilder_hopf_link_to_ball_pipeline);
    run("test_knotbuilder_figure_eight_to_ball_pipeline",
        test_knotbuilder_figure_eight_to_ball_pipeline);

    std::cout << "\n"
              << bold << (failed_count > 0 ? red : green) << "=== " << passed
              << " passed, " << failed_count << " failed ===" << resetColor
              << "\n";
    return failed_count > 0 ? 1 : 0;
}
