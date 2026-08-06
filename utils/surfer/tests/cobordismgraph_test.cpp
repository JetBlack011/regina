// cobordismgraph_test.cpp
//
// Tests for ../cobordismgraph.h/.cpp: the name/genus resolution graph
// verifyslicegenus.cpp builds up across its --input rows (addEdge(),
// resolvable(), propagateGraph(), recordWitness()) and the boundary-side
// safety classification (splitBoundary()) that feeds it. All pure logic on
// plain data types -- no triangulations, no search, no census -- so these
// tests run instantly and exercise the previously-untested (only via full
// CLI runs) core of verifyslicegenus's genus deductions directly.

#include <iostream>
#include <string>
#include <unistd.h>

#include <triangulation/dim3.h>

#include "../cobordismgraph.h"

using namespace cobordismgraph;

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

#define EXPECT_EQ(actual, expected, desc)                                     \
    do {                                                                      \
        auto _a = (actual);                                                   \
        auto _e = (expected);                                                 \
        if (_a == _e) {                                                       \
            std::cout << green << "  PASS: " << resetColor << (desc) << "\n"; \
            ++passed;                                                         \
        } else {                                                              \
            std::cout << red << "  FAIL: " << (desc) << "\n"                  \
                      << "        expected " << _e << ", got " << _a          \
                      << resetColor << "\n";                                  \
            ++failed_count;                                                   \
        }                                                                     \
    } while (0)

namespace {

std::function<std::string()> constPairSig(std::string sig) {
    return [sig] { return sig; };
}

void test_direct_witness_resolves() {
    Graph graph;
    std::unordered_map<std::string, int> knownGenus;
    std::unordered_map<std::string, OutputRow> outputRows;

    auto outcome = recordWitness("K", 1, 1, 1, "", 1, constPairSig("sig1"),
                                 "direct", graph, knownGenus, outputRows);

    EXPECT_EQ(outcome.resolved, true,
              "a direct witness at exactly target resolves the subject");
    EXPECT_EQ(outcome.fatalBug, false, "not a fatal bug");
    EXPECT_EQ(knownGenus.contains("K"), true, "knownGenus[K] is now set");
    EXPECT_EQ(knownGenus["K"], 1, "knownGenus[K] == the resolved genus");
    EXPECT_EQ(outputRows["K"].status, std::string("resolved"),
              "lo == hi -> status \"resolved\"");
    EXPECT_EQ(outputRows["K"].witnessKind, std::string("direct"),
              "witnessKind is passed through verbatim");
    EXPECT_EQ(outputRows["K"].witnessPairSig, std::string("sig1"),
              "witnessPairSig captured from the lazy callback");
}

void test_direct_witness_range_status() {
    Graph graph;
    std::unordered_map<std::string, int> knownGenus;
    std::unordered_map<std::string, OutputRow> outputRows;

    recordWitness("K", 1, 3, 3, "", 3, constPairSig(""), "direct", graph,
                 knownGenus, outputRows);

    EXPECT_EQ(outputRows["K"].status, std::string("range"),
              "lo != hi -> status \"range\" even though this run's target "
              "(hi) was reached exactly");
}

void test_direct_witness_fatal_bug() {
    Graph graph;
    std::unordered_map<std::string, int> knownGenus;
    std::unordered_map<std::string, OutputRow> outputRows;

    auto outcome = recordWitness("K", 2, 3, 3, "", 1, constPairSig(""),
                                 "direct", graph, knownGenus, outputRows);

    EXPECT_EQ(outcome.fatalBug, true,
              "a direct witness genus below subjectLo is a mathematical "
              "impossibility");
    EXPECT_EQ(outcome.resolved, false, "not resolved");
    EXPECT_EQ(knownGenus.contains("K"), false,
              "no knownGenus entry written on a fatal bug");
    EXPECT_EQ(outcome.fatalBugMessage.find("K") != std::string::npos, true,
              "fatal bug message names the subject");
}

void test_direct_witness_improvement_does_not_resolve() {
    Graph graph;
    std::unordered_map<std::string, int> knownGenus;
    std::unordered_map<std::string, OutputRow> outputRows;

    // subjectLo=1, target=3 (subjectHi), witnessGenus=2: >= lo, < target --
    // a genuine improvement over the known upper bound, but this run's
    // specific target (3) was not witnessed, so nothing resolves.
    auto outcome = recordWitness("K", 1, 3, 3, "", 2, constPairSig(""),
                                 "direct", graph, knownGenus, outputRows);

    EXPECT_EQ(outcome.resolved, false,
              "a genuine improvement below target still isn't a resolution "
              "of THIS run's target");
    EXPECT_EQ(outcome.fatalBug, false, "not a fatal bug");
    EXPECT_EQ(knownGenus.contains("K"), false, "no knownGenus entry written");
}

void test_cobordism_resolves_forward_direction() {
    Graph graph;
    std::unordered_map<std::string, int> knownGenus;
    std::unordered_map<std::string, OutputRow> outputRows;
    knownGenus["V"] = 5;

    // target == witnessGenus + h (8 == 3 + 5).
    auto outcome = recordWitness("K", 6, 8, 8, "V", 3, constPairSig("sig2"),
                                 "cobordism", graph, knownGenus, outputRows);

    EXPECT_EQ(outcome.resolved, true,
              "target == witnessGenus + h resolves via the constructive "
              "direction");
    EXPECT_EQ(knownGenus["K"], 8, "knownGenus[K] == target");
    EXPECT_EQ(outputRows["K"].viaKnot, std::string("V"),
              "viaKnot records which name this was resolved through");
    EXPECT_EQ(outputRows["K"].viaEdgeGenus, 3, "viaEdgeGenus == witnessGenus");
    EXPECT_EQ(hasEdgeWithGenus(graph, "K", "V", 3), true,
              "the witnessed edge is recorded");
    EXPECT_EQ(hasEdgeWithGenus(graph, "V", "K", 3), true,
              "addEdge() records both directions");
}

void test_cobordism_resolves_reverse_direction() {
    Graph graph;
    std::unordered_map<std::string, int> knownGenus;
    std::unordered_map<std::string, OutputRow> outputRows;
    knownGenus["V"] = 10;

    // target == h - witnessGenus (7 == 10 - 3).
    auto outcome = recordWitness("K", 3, 7, 7, "V", 3, constPairSig(""),
                                 "cobordism", graph, knownGenus, outputRows);

    EXPECT_EQ(outcome.resolved, true,
              "target == h - witnessGenus resolves via the contrapositive "
              "direction");
    EXPECT_EQ(knownGenus["K"], 7, "knownGenus[K] == target");
}

void test_cobordism_fatal_bug() {
    Graph graph;
    std::unordered_map<std::string, int> knownGenus;
    std::unordered_map<std::string, OutputRow> outputRows;
    knownGenus["V"] = 2;

    // Implied upper bound for K is 2 + 1 = 3, below subjectLo=5.
    auto outcome = recordWitness("K", 5, 9, 9, "V", 1, constPairSig(""),
                                 "cobordism", graph, knownGenus, outputRows);

    EXPECT_EQ(outcome.fatalBug, true,
              "an implied upper bound below subjectLo is a mathematical "
              "impossibility");
    EXPECT_EQ(outcome.fatalBugMessage.find("K") != std::string::npos &&
                  outcome.fatalBugMessage.find("V") != std::string::npos,
              true, "fatal bug message names both the subject and viaName");
}

void test_cobordism_edge_recorded_without_resolving_then_propagates() {
    Graph graph;
    std::unordered_map<std::string, int> knownGenus;
    std::unordered_map<std::string, OutputRow> outputRows;
    // Neither K nor V is known yet.

    auto outcome = recordWitness("K", 1, 10, 10, "V", 4, constPairSig("sig3"),
                                 "cobordism", graph, knownGenus, outputRows);

    EXPECT_EQ(outcome.resolved, false,
              "neither endpoint known yet -- nothing resolves immediately");
    EXPECT_EQ(outcome.fatalBug, false, "not a fatal bug");
    EXPECT_EQ(hasEdgeWithGenus(graph, "K", "V", 4), true,
              "the edge is still recorded for later use");

    // V becomes known later (e.g. its own row resolves independently).
    knownGenus["V"] = 6;
    auto resolved = resolvable("K", 10, graph, knownGenus);
    EXPECT_EQ(resolved.has_value(), true,
              "resolvable() finds the earlier-recorded edge once V is known "
              "(4 + 6 == 10)");
    EXPECT_EQ(resolved->viaKnot, std::string("V"), "resolved via V");

    // The same scenario via propagateGraph(), matching how
    // verifyslicegenus.cpp actually drives this after every row.
    std::vector<InputRow> pending = {InputRow{"K", "", 1, 10, 0}};
    std::vector<std::string> newlyResolved =
        propagateGraph(pending, graph, knownGenus, outputRows);
    EXPECT_EQ(knownGenus["K"], 10,
              "propagateGraph() resolves K purely from the graph, with no "
              "search run for it");
    EXPECT_EQ(outputRows["K"].witnessKind, std::string("propagated"),
              "propagateGraph()'s own witnessKind");
    EXPECT_EQ(newlyResolved.size(), static_cast<size_t>(1),
              "propagateGraph() reports exactly the one name it newly "
              "resolved, for the caller to announce");
    EXPECT_EQ(newlyResolved[0], std::string("K"),
              "the reported name is K");
}

void test_dedup_skips_pairsig_capture_on_repeat_non_resolving_edge() {
    Graph graph;
    std::unordered_map<std::string, int> knownGenus;
    std::unordered_map<std::string, OutputRow> outputRows;
    addEdge(graph, "K", "V", 4, "existing-sig");

    int captureCalls = 0;
    auto countingCapture = [&] {
        ++captureCalls;
        return std::string("should-not-be-used");
    };

    // V still unknown, so this can't resolve -- and an edge at this exact
    // genus already exists, so recordWitness() should bail out before ever
    // invoking the (expensive) pairSig capture.
    recordWitness("K", 1, 100, 100, "V", 4, countingCapture, "cobordism",
                 graph, knownGenus, outputRows);

    EXPECT_EQ(captureCalls, 0,
              "a duplicate non-resolving edge never captures pairSig -- "
              "capturing it is expensive and this is the overwhelmingly "
              "common case in a real search");
}

void test_split_boundary_single_curve_is_safe() {
    std::vector<BoundaryComponentNames> components = {
        BoundaryComponentNames{0, {"3_1"}, std::nullopt},
        BoundaryComponentNames{1, {"4_1"}, std::nullopt},
    };
    BoundarySplit split = splitBoundary(components, 0, "3_1");

    EXPECT_EQ(split.searchCurveCount, static_cast<size_t>(1),
              "component 0 (== searchSideBC) is the search side");
    EXPECT_EQ(split.otherSides.size(), static_cast<size_t>(1),
              "exactly one other side");
    EXPECT_EQ(split.otherSides[0].name, std::string("4_1"),
              "a single curve is named directly");
    EXPECT_EQ(split.otherSides[0].safe, true,
              "a single curve has no orientation ambiguity -- always safe");
}

void test_split_boundary_unlink_is_safe() {
    std::vector<BoundaryComponentNames> components = {
        BoundaryComponentNames{0, {"3_1"}, std::nullopt},
        BoundaryComponentNames{
            1, {"a", "b"}, std::optional<std::string>("2-component unlink")},
    };
    BoundarySplit split = splitBoundary(components, 0, "3_1");

    EXPECT_EQ(split.otherSides.size(), static_cast<size_t>(1), "one other side");
    EXPECT_EQ(split.otherSides[0].name, std::string("2-component unlink"),
              "named via linkName");
    EXPECT_EQ(split.otherSides[0].safe, true,
              "a split unlink's components don't interact -- orientation "
              "doesn't matter, safe for a deduction");
}

void test_split_boundary_linked_multicomponent_is_unsafe() {
    std::vector<BoundaryComponentNames> components = {
        BoundaryComponentNames{0, {"3_1"}, std::nullopt},
        BoundaryComponentNames{1, {"a", "b"},
                               std::optional<std::string>("L6a3")},
    };
    BoundarySplit split = splitBoundary(components, 0, "3_1");

    EXPECT_EQ(split.otherSides[0].name, std::string("L6a3"), "named");
    EXPECT_EQ(split.otherSides[0].safe, false,
              "a genuinely linked multi-component name is NOT safe -- "
              "different orientations of the same link can have different "
              "true slice genus while sharing one complement");
}

void test_split_boundary_multiple_other_sides_not_collapsed() {
    std::vector<BoundaryComponentNames> components = {
        BoundaryComponentNames{0, {"3_1"}, std::nullopt}, // search side
        BoundaryComponentNames{1, {"4_1"}, std::nullopt}, // other #1
        BoundaryComponentNames{2, {"5_1"}, std::nullopt}, // other #2
    };
    BoundarySplit split = splitBoundary(components, 0, "3_1");

    EXPECT_EQ(split.otherSides.size(), static_cast<size_t>(2),
              "both other sides are kept -- neither silently overwrites the "
              "other (the old BoundarySplit collapsed multiple \"other\" "
              "components into a single farName, discarding all but the "
              "last)");
}

void test_split_boundary_search_side_name_mismatch_is_not_search_side() {
    // Same shape as the real L6a3{0} fatal-bug repro: component 0 ==
    // searchSideBC holds exactly as many curves as the row's own component
    // count (2), but they don't actually identify as this row's own link
    // -- the DFS wandered onto an unrelated 2-component link that just
    // happens to have the same curve count. splitBoundary() must not
    // trust the geometric position alone.
    std::vector<BoundaryComponentNames> components = {
        BoundaryComponentNames{0, {"Unknot", "Unknot"},
                               std::optional<std::string>("L206001")},
        BoundaryComponentNames{1, {"Unknot"}, std::nullopt},
    };
    BoundarySplit split = splitBoundary(components, 0, "L6a3{0}");

    EXPECT_EQ(split.searchCurveCount, static_cast<size_t>(0),
              "component 0's identified name (L206001) doesn't match this "
              "row's own name (L6a3{0}), so it is NOT treated as the search "
              "side even though it's geometrically on searchSideBC and even "
              "though its curve count matches this row's component count");
    EXPECT_EQ(split.otherSides.size(), static_cast<size_t>(2),
              "both components are treated as \"other\" sides instead");
    EXPECT_EQ(split.otherSides[0].name, std::string("L206001"),
              "the mismatched component is named via its own linkName");
    EXPECT_EQ(split.otherSides[0].safe, false,
              "L206001 is a genuine (unproven-split) multi-component "
              "linkName -- not safe for a deduction");
    EXPECT_EQ(split.otherSides[1].name, std::string("Unknot"), "");
    EXPECT_EQ(split.otherSides[1].safe, true, "a single curve is always safe");
}

void test_split_boundary_search_side_name_match_is_search_side() {
    // Sanity check paired with the mismatch test above: when the name
    // DOES match, component == searchSideBC is accepted as the search side
    // exactly as before, even with more than one curve.
    std::vector<BoundaryComponentNames> components = {
        BoundaryComponentNames{0, {"Unknot", "Unknot"},
                               std::optional<std::string>("L6a3{0}")},
        BoundaryComponentNames{1, {"Unknot"}, std::nullopt},
    };
    BoundarySplit split = splitBoundary(components, 0, "L6a3{0}");

    EXPECT_EQ(split.searchCurveCount, static_cast<size_t>(2),
              "component 0's name matches this row's own name, so it is "
              "the search side despite holding more than one curve");
    EXPECT_EQ(split.otherSides.size(), static_cast<size_t>(1),
              "only the genuinely-other component remains");
}

// ─────────────────────────────────────────────────────────────────────────────
// matchesRowOrientation()'s own decision logic, isolated from
// buildRowOrientation()'s isomorphism/geometry machinery (validated
// separately -- see this feature's own diagnostic and the end-to-end
// L6a3{0}/L6a3{1} repro) by hand-constructing RowOrientation/OrientedCurve
// directly against a single tetrahedron's own edges, rather than going
// through a real knotbuilder+CobordismBuilder pipeline. Two of its own
// edges stand in for two independent "components".
// ─────────────────────────────────────────────────────────────────────────────
void test_matches_row_orientation_logic() {
    regina::Triangulation<3> tri;
    tri.newTetrahedron(); // tetrahedron 0: e0, e1 below, kept vertex-disjoint
    tri.newTetrahedron(); // tetrahedron 1, ungled to the first: fully
                          // disjoint from it, for e2 below
    // Regina's standard tetrahedron edge ordering is
    // {(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)}, so edge 0 == (0,1) and edge 5
    // == (2,3) of the same tetrahedron are its one pair of opposite
    // (vertex-disjoint) edges -- needed so the two "components" below
    // don't silently clobber each other's entry in
    // RowOrientation::headOf (keyed by tail vertex index, so two edges
    // sharing an endpoint would collide).
    regina::Edge<3> *e0 = tri.tetrahedron(0)->edge(0);
    regina::Edge<3> *e1 = tri.tetrahedron(0)->edge(5);
    // A third edge from the other, entirely disjoint tetrahedron -- not
    // one of the row's own tagged edges at all.
    regina::Edge<3> *e2 = tri.tetrahedron(1)->edge(0);

    RowOrientation row;
    row.headOf[e0->vertex(0)->index()] = e0->vertex(1)->index();
    row.headOf[e1->vertex(0)->index()] = e1->vertex(1)->index();

    std::vector<OrientedCurve> allMatch = {{{e0, false}}, {{e1, false}}};
    EXPECT_EQ(matchesRowOrientation(row, allMatch), true,
              "every curve's induced direction agrees with the row's own "
              "tag -- accepted");

    std::vector<OrientedCurve> allFlipped = {{{e0, true}}, {{e1, true}}};
    EXPECT_EQ(matchesRowOrientation(row, allFlipped), true,
              "every curve's induced direction disagrees with the row's "
              "own tag, but uniformly -- a global flip, always allowed "
              "(it's just the surface's other orientation choice)");

    std::vector<OrientedCurve> mixed = {{{e0, false}}, {{e1, true}}};
    EXPECT_EQ(matchesRowOrientation(row, mixed), false,
              "one component agrees, the other doesn't -- exactly the "
              "L6a3{0}/L6a3{1} misattribution signature, rejected");

    std::vector<OrientedCurve> unknownEdge = {{{e2, false}}};
    EXPECT_EQ(matchesRowOrientation(row, unknownEdge), false,
              "a curve edge that isn't one of the row's own tagged edges "
              "at all is rejected, not silently ignored");
}

} // namespace

void run(const std::string &name, void (*fn)()) {
    std::cout << bold << "\n=== " << name << " ===" << resetColor << "\n";
    fn();
}

int main() {
    run("direct_witness_resolves", test_direct_witness_resolves);
    run("direct_witness_range_status", test_direct_witness_range_status);
    run("direct_witness_fatal_bug", test_direct_witness_fatal_bug);
    run("direct_witness_improvement_does_not_resolve",
        test_direct_witness_improvement_does_not_resolve);
    run("cobordism_resolves_forward_direction",
        test_cobordism_resolves_forward_direction);
    run("cobordism_resolves_reverse_direction",
        test_cobordism_resolves_reverse_direction);
    run("cobordism_fatal_bug", test_cobordism_fatal_bug);
    run("cobordism_edge_recorded_without_resolving_then_propagates",
        test_cobordism_edge_recorded_without_resolving_then_propagates);
    run("dedup_skips_pairsig_capture_on_repeat_non_resolving_edge",
        test_dedup_skips_pairsig_capture_on_repeat_non_resolving_edge);
    run("split_boundary_single_curve_is_safe",
        test_split_boundary_single_curve_is_safe);
    run("split_boundary_unlink_is_safe", test_split_boundary_unlink_is_safe);
    run("split_boundary_linked_multicomponent_is_unsafe",
        test_split_boundary_linked_multicomponent_is_unsafe);
    run("split_boundary_multiple_other_sides_not_collapsed",
        test_split_boundary_multiple_other_sides_not_collapsed);
    run("split_boundary_search_side_name_mismatch_is_not_search_side",
        test_split_boundary_search_side_name_mismatch_is_not_search_side);
    run("split_boundary_search_side_name_match_is_search_side",
        test_split_boundary_search_side_name_match_is_search_side);
    run("matches_row_orientation_logic", test_matches_row_orientation_logic);

    std::cout << bold << "\n=== Summary: " << passed << " passed, "
              << failed_count << " failed ===" << resetColor << "\n";
    return failed_count > 0 ? 1 : 0;
}
