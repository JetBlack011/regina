// linkcomplement_test.cpp
//
// Tests for EdgeComplement::edgeIndices() and Link's split into components
// (see ../complement/linkcomplement.h): pure edge-set representation,
// independent of any naming/census logic.
//
// See surfer's namecache_test.cpp for BoundarySignatureCache and
// censusnaming_test.cpp for the complement cache -- those exercise
// census::nameComplement() and its caches, not this file's
// representation-only concern.

#include <algorithm>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include "diagramtriangulation/fromdiagram.h"
#include "linknaming/complement/linkcomplement.h"

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

// Formats a vector<size_t> as "[a, b, c]", purely so EXPECT_EQ has
// something streamable to print on failure.
std::string toString(const std::vector<size_t> &v) {
    std::ostringstream out;
    out << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        if (i)
            out << ", ";
        out << v[i];
    }
    out << "]";
    return out.str();
}

// A single, unglued pentachoron's boundary: 5 tetrahedra triangulating S^3,
// built the exact same way SurfaceSearch builds each ambient boundary
// component (BoundaryComponent<4>::build()).
regina::Triangulation<3> testBoundary() {
    regina::Triangulation<4> pent;
    pent.newSimplex();
    return pent.boundaryComponent(0)->build();
}

void test_edge_indices() {
    regina::Triangulation<3> boundary = testBoundary();
    const regina::Edge<3> *e0 = boundary.edge(0);
    const regina::Edge<3> *e2 = boundary.edge(2);

    EdgeComplement ec(boundary, {e2, e0});
    EXPECT_EQ(toString(ec.edgeIndices()), toString({0, 2}),
              "edgeIndices() returns tracked edges sorted by index, "
              "regardless of insertion order");
}

// A Link's components as their edges' indices, in the Link's own order:
// "[[a, b, ...], [c, ...]]".
std::string componentsOf(const Link &link) {
    std::ostringstream out;
    out << "[";
    for (size_t c = 0; c < link.comps_.size(); ++c) {
        std::vector<size_t> indices;
        for (const regina::Edge<3> *e : link.comps_[c].edges())
            indices.push_back(e->index());
        out << (c ? ", " : "") << toString(indices);
    }
    out << "]";
    return out.str();
}

// Link::Link() splits an edge set into components by a walk whose every
// choice is made by edge index (phase 3.0): the components, their order and
// each one's edge sequence depend on the set alone -- not on the order the
// edges were supplied in, and never on their addresses, which once made
// goal runs explore different links under different heap layouts.
void test_link_components_independent_of_input_order() {
    const char *pds[] = {
        "PD[X[4; 1; 3; 2]; X[2; 3; 1; 4]]",                        // L2a1{0}
        "PD[X[6; 1; 7; 2]; X[8; 3; 5; 4]; X[2; 5; 3; 6]; X[4; 7; 1; 8]]", // L4a1{0}
        "PD[X[6; 1; 7; 2]; X[12; 8; 9; 7]; X[4; 12; 1; 11]; X[10; 5; 11; 6]; "
        "X[8; 4; 5; 3]; X[2; 9; 3; 10]]",                          // L6a4{0;0}
        "[[1;5;2;4];[3;1;4;6];[5;3;6;2]]",                         // 3_1
    };
    const size_t expectedComponents[] = {2, 2, 3, 1};
    std::mt19937 rng(20261001);
    for (size_t p = 0; p < std::size(pds); ++p) {
        const diagramtriangulation::TriangulationWithLink built =
            diagramtriangulation::buildLink(diagramtriangulation::parsePDCode(pds[p]));
        std::vector<const regina::Edge<3> *> edges = built.edges;
        std::ranges::sort(edges, {}, [](const regina::Edge<3> *e) { return e->index(); });

        const Link sorted(built.tri, edges);
        const std::string reference = componentsOf(sorted);
        EXPECT_EQ(sorted.comps_.size(), expectedComponents[p],
                  std::string(pds[p]) + ": component count");

        // Each component starts from its lowest edge, and the components
        // come in the order of those lowest edges.
        bool startsLowest = true, ordered = true;
        size_t previousFirst = 0;
        for (size_t c = 0; c < sorted.comps_.size(); ++c) {
            const auto &es = sorted.comps_[c].edges();
            const size_t first = es.front()->index();
            for (const regina::Edge<3> *e : es)
                if (e->index() < first) startsLowest = false;
            if (c > 0 && first < previousFirst) ordered = false;
            previousFirst = first;
        }
        EXPECT_EQ(startsLowest, true,
                  std::string(pds[p]) + ": each component starts from its lowest edge");
        EXPECT_EQ(ordered, true,
                  std::string(pds[p]) + ": components ordered by their lowest edge");

        // The same set in other orders, and copied to fresh allocations
        // (so the edges' addresses relative to each other differ too).
        bool same = true;
        std::vector<const regina::Edge<3> *> shuffled = edges;
        std::ranges::reverse(shuffled);
        if (componentsOf(Link(built.tri, shuffled)) != reference) same = false;
        for (int k = 0; k < 20; ++k) {
            std::ranges::shuffle(shuffled, rng);
            if (componentsOf(Link(built.tri, shuffled)) != reference) same = false;
            auto copy = std::make_unique<regina::Triangulation<3>>(built.tri);
            std::vector<const regina::Edge<3> *> copied;
            for (const regina::Edge<3> *e : shuffled) copied.push_back(copy->edge(e->index()));
            if (componentsOf(Link(*copy, copied)) != reference) same = false;
        }
        EXPECT_EQ(same, true,
                  std::string(pds[p]) + ": the same components and edge sequences "
                  "from every input order, in every copy");
    }
}

} // namespace

void run(const std::string &name, void (*fn)()) {
    std::cout << bold << "\n=== " << name << " ===" << resetColor << "\n";
    fn();
}

int main() {
    run("edge_indices", test_edge_indices);
    run("link_components_independent_of_input_order",
        test_link_components_independent_of_input_order);

    std::cout << bold << "\n=== Summary: " << passed << " passed, "
              << failed_count << " failed ===" << resetColor << "\n";
    return failed_count > 0 ? 1 : 0;
}
