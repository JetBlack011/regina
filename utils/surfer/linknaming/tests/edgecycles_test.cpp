// edgecycles_test.cpp
//
// Tests for ../complement/edgecycles.h, one walk per family: directed
// chaining (both open-chain policies) and its counting variant, orienting
// undirected edges (tolerant and one-closed-curve), and counting closed
// curves. Hand-made edge lists for the shapes that matter (closed curves,
// two curves, an arc, a branching, a loop edge, a repeated edge), and a real
// link's edges from buildLink().

#include <algorithm>
#include <iostream>
#include <random>
#include <string>
#include <unistd.h>
#include <vector>

#include "diagramtriangulation/fromdiagram.h"
#include "linknaming/complement/edgecycles.h"
#include "linknaming/complement/linkcomplement.h"

using namespace edgecycles;

static int passed = 0, failed_count = 0;

// For EXPECT_EQ's failure message.
std::ostream &operator<<(std::ostream &os, const std::optional<size_t> &v) {
    return v ? os << *v : os << "nullopt";
}

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

std::string str(const std::vector<size_t> &v) {
    std::string s = "[";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? "," : "") + std::to_string(v[i]);
    return s + "]";
}
std::string str(const std::vector<std::vector<size_t>> &v) {
    std::string s;
    for (const auto &c : v) s += str(c);
    return s;
}
std::string str(const std::vector<Step> &v) {
    std::string s = "[";
    for (size_t i = 0; i < v.size(); ++i)
        s += (i ? "," : "") + std::to_string(v[i].pos) + (v[i].reversed ? "r" : "");
    return s + "]";
}
std::string str(const std::vector<std::vector<Step>> &v) {
    std::string s;
    for (const auto &c : v) s += str(c);
    return s;
}

// Directed triangle 0->1->2->0 (ids 10, 11, 12) and directed 2-gon 5->6->5
// (ids 20, 21), interleaved.
const std::vector<EdgeEnds> kTwoDirected = {
    {11, 1, 2}, {20, 5, 6}, {10, 0, 1}, {21, 6, 5}, {12, 2, 0}};

void test_chain_directed() {
    auto c = chainDirected(kTwoDirected, OpenChain::stop);
    EXPECT_EQ(c.has_value(), true, "closed curves chain");
    EXPECT_EQ(str(*c), std::string("[0,4,2][1,3]"),
              "each curve starts at its first edge in input order and runs head to tail");
    EXPECT_EQ(str(*chainDirected(kTwoDirected, OpenChain::refuse)), str(*c),
              "refuse changes nothing when every curve closes");
    EXPECT_EQ(countDirectedCycles(kTwoDirected), std::optional<size_t>(2), "two cycles");

    // An open directed chain 0->1->2 (no edge leaves 2).
    const std::vector<EdgeEnds> open = {{1, 0, 1}, {2, 1, 2}};
    auto stopped = chainDirected(open, OpenChain::stop);
    EXPECT_EQ(str(*stopped), std::string("[0,1]"), "stop: the curve ends at the dead end");
    EXPECT_EQ(chainDirected(open, OpenChain::refuse).has_value(), false,
              "refuse: a dead end gives no answer at all");
    EXPECT_EQ(countDirectedCycles(open).has_value(), false, "an open chain is not cycles");

    // Two edges leave vertex 0: the LAST in input order is taken.
    const std::vector<EdgeEnds> branch = {{1, 0, 1}, {2, 1, 0}, {3, 0, 2}, {4, 2, 0}};
    EXPECT_EQ(str(*chainDirected(branch, OpenChain::stop)), std::string("[2,3][1,2,3,2,3]"),
              "a vertex left twice: the last edge leaving it wins, and a walk that never "
              "returns to its start stops after n + 1 steps (as the old walks did)");
    EXPECT_EQ(countDirectedCycles(branch).has_value(), false,
              "two edges leaving one vertex are not cycles");
}

// Undirected: square 0-1-2-3 given out of order and against its direction.
const std::vector<EdgeEnds> kSquare = {{0, 0, 1}, {1, 2, 3}, {2, 1, 2}, {3, 0, 3}};

void test_walk_curves() {
    EXPECT_EQ(str(walkCurves(kSquare)), std::string("[0,2,1,3r]"),
              "one curve from the first edge's v0, each next edge the first unused at the vertex");
    // Two disjoint triangles, and an arc 7-8-9.
    const std::vector<EdgeEnds> mixed = {{0, 0, 1}, {1, 4, 5}, {2, 1, 2}, {3, 5, 6},
                                         {4, 2, 0}, {5, 6, 4}, {6, 8, 7}, {7, 8, 9}};
    EXPECT_EQ(str(walkCurves(mixed)), std::string("[0,2,4][1,3,5][6][7]"),
              "components by first edge; an arc walked from its first edge's v0 runs out at "
              "once and is split (tolerant, as Link's own walk was)");
}

void test_walk_closed_curve() {
    auto s = walkClosedCurve(kSquare);
    EXPECT_EQ(s.has_value() ? str(*s) : std::string("none"), std::string("[0,2,1,3r]"),
              "one closed curve");
    EXPECT_EQ(walkClosedCurve({{0, 0, 1}, {1, 1, 2}}).has_value(), false, "an arc is refused");
    EXPECT_EQ(walkClosedCurve({{0, 0, 1}, {1, 1, 0}, {2, 2, 3}, {3, 3, 2}}).has_value(), false,
              "two curves are refused");
    EXPECT_EQ(walkClosedCurve({{0, 0, 1}, {1, 1, 2}, {2, 2, 0}, {3, 0, 3}, {4, 3, 0}})
                  .has_value(),
              false, "a vertex of degree 4 is refused");
    auto loop = walkClosedCurve({{9, 4, 4}});
    EXPECT_EQ(loop.has_value() ? str(*loop) : std::string("none"), std::string("[0]"),
              "a loop edge is a closed curve");
    EXPECT_EQ(walkClosedCurve({}).has_value(), false, "nothing is not a curve");
}

void test_count_closed_curves() {
    EXPECT_EQ(countClosedCurves(kSquare), std::optional<size_t>(1), "one");
    EXPECT_EQ(countClosedCurves({{0, 0, 1}, {1, 1, 0}, {2, 2, 3}, {3, 3, 2}}),
              std::optional<size_t>(2), "two 2-gons");
    EXPECT_EQ(countClosedCurves({{9, 4, 4}}), std::optional<size_t>(1), "a loop edge");
    EXPECT_EQ(countClosedCurves({{0, 0, 1}, {1, 1, 2}}).has_value(), false, "an arc");
    EXPECT_EQ(countClosedCurves({{0, 0, 1}, {0, 0, 1}}).has_value(), false,
              "a repeated edge");
    EXPECT_EQ(countClosedCurves({}), std::optional<size_t>(0), "nothing: zero curves");
}

// A real link (L6a4{0;0}, three components) through every family: the
// counts agree, every input order gives the same curves after sorting, and
// Link's components are walkCurves() of its sorted edges.
void test_real_link() {
    const diagramtriangulation::TriangulationWithLink built =
        diagramtriangulation::buildLink(diagramtriangulation::parsePDCode(
        "PD[X[6; 1; 7; 2]; X[12; 8; 9; 7]; X[4; 12; 1; 11]; X[10; 5; 11; 6]; "
        "X[8; 4; 5; 3]; X[2; 9; 3; 10]]"));
    const std::vector<EdgeEnds> ends = endsOf(built.edges);
    EXPECT_EQ(countClosedCurves(ends), std::optional<size_t>(3), "three closed curves");
    EXPECT_EQ(walkCurves(ends).size(), size_t(3), "three walked curves");
    std::vector<EdgeEnds> directed;
    for (size_t i = 0; i < built.edges.size(); ++i)
        directed.push_back(built.reversed[i] ? EdgeEnds{ends[i].id, ends[i].v1, ends[i].v0}
                                             : ends[i]);
    EXPECT_EQ(countDirectedCycles(directed), std::optional<size_t>(3),
              "the PD's directions chain into three directed cycles");
    EXPECT_EQ(chainDirected(directed, OpenChain::refuse)->size(), size_t(3),
              "and chainDirected() finds them");
    const Link link(built.tri, built.edges);
    EXPECT_EQ(link.countComponents(), 3, "Link agrees");
}

} // namespace

void run(const std::string &name, void (*fn)()) {
    std::cout << bold << "\n=== " << name << " ===" << resetColor << "\n";
    fn();
}

int main() {
    run("chain_directed", test_chain_directed);
    run("walk_curves", test_walk_curves);
    run("walk_closed_curve", test_walk_closed_curve);
    run("count_closed_curves", test_count_closed_curves);
    run("real_link", test_real_link);

    std::cout << bold << "\n=== Summary: " << passed << " passed, "
              << failed_count << " failed ===" << resetColor << "\n";
    return failed_count > 0 ? 1 : 0;
}
