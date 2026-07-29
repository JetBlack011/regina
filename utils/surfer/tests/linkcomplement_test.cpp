// linkcomplement_test.cpp
//
// Tests for EdgeComplement::edgeIndices() and BoundarySignatureCache (see
// ../linkcomplement.h): the pre-triangulation dedup layer that
// canonicalizes a marked edge set against its ambient (fixed) boundary
// triangulation's own automorphism group, so a boundary curve already seen
// -- exactly, or up to a symmetry of that triangulation -- short-circuits
// before EdgeComplement::identify() ever runs
// buildComplement()/simplify()/isoSig().
//
// The property that matters most in the second test below is automorphism
// invariance: two marked edge sets related by an automorphism of the fixed
// ambient triangulation must be treated as the same cache entry. That test
// re-derives the automorphism's action on edges independently of
// BoundarySignatureCache's own internals (which reuse pairsig.h's
// faceDescriptor()/applyIsomorphism()/resolveFaceIndex()) -- mirroring
// pairsig_test.cpp's own independent mapFace() helper, for the same reason:
// a test of automorphism-invariance shouldn't be circular with the
// machinery it's checking.

#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

#include <maths/perm.h>
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

// Independently re-derives edge `e`'s image under automorphism `iso` of
// `t` -- the same "front().vertices() + FaceNumbering::faceNumber"
// technique BoundarySignatureCache itself uses internally (via pairsig.h),
// but written standalone here so this test doesn't depend on that
// machinery to construct a *known* automorphic edge set.
size_t mapEdge(const regina::Triangulation<3> &t,
               const regina::Isomorphism<3> &iso, int e) {
    const auto &emb = t.edge(e)->front();
    size_t srcSimplex = emb.simplex()->index();
    regina::Perm<4> p = emb.vertices();
    auto destSimplex = static_cast<size_t>(iso.simpImage(srcSimplex));
    regina::Perm<4> q = iso.facetPerm(srcSimplex) * p;
    int localEdge = regina::FaceNumbering<3, 1>::faceNumber(q);
    return t.simplex(destSimplex)->edge(localEdge)->index();
}

// A single, unglued pentachoron's boundary: 5 tetrahedra triangulating S^3
// with the full S5 symmetry group (120 automorphisms) acting on its 10
// edges -- built the exact same way SurfaceSearch builds each ambient
// boundary component (BoundaryComponent<4>::build()), so this exercises
// BoundarySignatureCache against a realistic, richly symmetric example
// without needing the rest of the search pipeline.
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

void test_memoization() {
    regina::Triangulation<3> boundary = testBoundary();
    BoundarySignatureCache cache(boundary);

    int computeCalls = 0;
    auto compute = [&] {
        ++computeCalls;
        return std::string("result-A");
    };

    std::string r1 = cache.identifyCached({0}, compute);
    std::string r2 = cache.identifyCached({0}, compute);

    EXPECT_EQ(r1, std::string("result-A"),
              "first call for a new edge set returns compute()'s result");
    EXPECT_EQ(r2, std::string("result-A"),
              "second call for the same edge set returns the cached result");
    EXPECT_EQ(computeCalls, 1,
              "compute() is only invoked once across both calls");
    EXPECT_EQ(cache.stats().checks, 2LL, "both calls count as checks");
    EXPECT_EQ(cache.stats().hits, 1LL,
              "exactly the second (repeat) call is a cache hit");
    EXPECT_EQ(cache.size(), static_cast<size_t>(1),
              "one distinct canonical signature seen so far");
}

void test_automorphism_invariance() {
    regina::Triangulation<3> boundary = testBoundary();

    // Find a nontrivial automorphism of `boundary` that moves edge 0 to a
    // different edge -- guaranteed to exist for a single pentachoron's
    // boundary (S5 acts transitively on its 10 edges).
    std::optional<size_t> movedEdge;
    boundary.findAllIsomorphisms(boundary,
        [&](const regina::Isomorphism<3> &iso) {
            size_t image = mapEdge(boundary, iso, 0);
            if (image != 0) {
                movedEdge = image;
                return true; // stop -- found one
            }
            return false;
        });

    if (!movedEdge) {
        std::cout << red
                  << "  FAIL: could not find a nontrivial automorphism "
                     "moving edge 0 -- test triangulation assumption "
                     "violated"
                  << resetColor << "\n";
        ++failed_count;
        return;
    }

    BoundarySignatureCache cache(boundary);
    EdgeComplement real(boundary, {boundary.edge(0)});
    std::string realResult = real.identify();

    std::string r1 = cache.identifyCached({0}, [&] { return realResult; });

    bool sentinelCalled = false;
    std::string r2 = cache.identifyCached({*movedEdge}, [&] {
        sentinelCalled = true;
        return std::string("SENTINEL_SHOULD_NOT_BE_CALLED");
    });

    EXPECT_EQ(r1, realResult,
              "identifyCached({0}, ...) returns EdgeComplement::identify()'s "
              "real result");
    EXPECT_EQ(r2, realResult,
              "identifyCached() for edge 0's image under a boundary "
              "automorphism returns the SAME result as edge 0 itself");
    EXPECT_EQ(sentinelCalled, false,
              "...without ever invoking its own compute() callback (a "
              "genuine cache hit, not a coincidental equal result)");
    EXPECT_EQ(cache.stats().hits, 1LL,
              "exactly one hit recorded (the automorphism-equivalent "
              "lookup)");
}

} // namespace

void run(const std::string &name, void (*fn)()) {
    std::cout << bold << "\n=== " << name << " ===" << resetColor << "\n";
    fn();
}

int main() {
    run("edge_indices", test_edge_indices);
    run("memoization", test_memoization);
    run("automorphism_invariance", test_automorphism_invariance);

    std::cout << bold << "\n=== Summary: " << passed << " passed, "
              << failed_count << " failed ===" << resetColor << "\n";
    return failed_count > 0 ? 1 : 0;
}
