// linkcomplement_test.cpp
//
// Tests for EdgeComplement::edgeIndices() (see ../linkcomplement.h): pure
// edge-set representation, independent of any recognition/census logic.
//
// See identifycomplement_test.cpp for BoundarySignatureCache and
// recognition-cache tests -- those exercise identify::identify() and its
// caches, not this file's representation-only concern.

#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include "../linkcomplement.h"

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

} // namespace

void run(const std::string &name, void (*fn)()) {
    std::cout << bold << "\n=== " << name << " ===" << resetColor << "\n";
    fn();
}

int main() {
    run("edge_indices", test_edge_indices);

    std::cout << bold << "\n=== Summary: " << passed << " passed, "
              << failed_count << " failed ===" << resetColor << "\n";
    return failed_count > 0 ? 1 : 0;
}
