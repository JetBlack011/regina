// identifycomplement_test.cpp
//
// Tests for identify::BoundarySignatureCache and the recognition cache (see
// ../identifycomplement.h): the pre-triangulation dedup layer that
// canonicalizes a marked edge set against its ambient (fixed) boundary
// triangulation's own automorphism group, so a boundary curve already seen
// -- exactly, or up to a symmetry of that triangulation -- short-circuits
// before identify::identify() ever runs
// buildComplement()/simplify()/isoSig(); and the isoSig-keyed recognition
// cache behind identify::identify() itself.
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
#include <triangulation/example3.h>

#include "../identifycomplement.h"
#include "../knotbuilder.h"
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

void test_memoization() {
    regina::Triangulation<3> boundary = testBoundary();
    identify::BoundarySignatureCache cache(boundary);

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

    identify::BoundarySignatureCache cache(boundary);
    EdgeComplement real(boundary, {boundary.edge(0)});
    std::string realResult = identify::identify(real);

    std::string r1 = cache.identifyCached({0}, [&] { return realResult; });

    bool sentinelCalled = false;
    std::string r2 = cache.identifyCached({*movedEdge}, [&] {
        sentinelCalled = true;
        return std::string("SENTINEL_SHOULD_NOT_BE_CALLED");
    });

    EXPECT_EQ(r1, realResult,
              "identifyCached({0}, ...) returns identify::identify()'s "
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

void test_boundary_signature_cache_clear_threshold() {
    regina::Triangulation<3> boundary = testBoundary();
    identify::BoundarySignatureCache cache(boundary, /*clearThreshold=*/2);

    int computeCalls = 0;
    auto compute = [&] {
        ++computeCalls;
        return "r" + std::to_string(computeCalls);
    };

    // Edge sets of different sizes can never canonicalize to the same key
    // (automorphisms are bijections, so they preserve cardinality) -- a
    // simple, topology-agnostic way to force several genuinely distinct
    // cache entries regardless of this boundary's own automorphism group.
    cache.identifyCached({0}, compute);       // 1 edge -- distinct key
    cache.identifyCached({0, 1}, compute);    // 2 edges -- distinct key; cache now at its threshold
    EXPECT_EQ(cache.stats().cacheResets, 0LL,
              "no reset yet -- the threshold is only checked before "
              "admitting the NEXT new key");

    cache.identifyCached({0, 1, 2}, compute); // 3 edges -- triggers the clear before inserting
    EXPECT_EQ(cache.stats().cacheResets, 1LL,
              "exactly one reset after exceeding the threshold");
    EXPECT_EQ(cache.size(), static_cast<size_t>(1),
              "post-reset, only the entry that triggered the reset is "
              "present");

    int callsBefore = computeCalls;
    cache.identifyCached({0}, compute); // was evicted by the reset above
    EXPECT_EQ(computeCalls, callsBefore + 1,
              "a pre-reset key is looked up as a fresh miss after the "
              "cache was cleared -- a clean miss, not a crash or a wrong "
              "hit against an unrelated key (string-keyed, no id-reuse "
              "risk)");
}

void test_recognition_cache_clear_threshold() {
    // Unlike BoundarySignatureCache's pre-triangulation combinatorial key,
    // recognitionCache is keyed by the drilled complement's POST-simplify
    // isoSig -- a topological invariant of the resulting manifold, not of
    // how many edges were drilled. So this needs two genuinely
    // topologically distinct complements, not just different edge counts.
    // Drilling no edges at all from a fixed triangulation just recognizes
    // that triangulation itself (buildComplement()'s pinch loop is a no-op
    // on an empty edge set) -- so pairing the trivial pentachoron-boundary
    // case (an unknot complement, genus 1) with the figure-eight knot
    // complement (genuinely hyperbolic, not any handlebody) guarantees two
    // distinct isoSigs.
    regina::Triangulation<3> boundary = testBoundary();
    regina::Triangulation<3> figureEight = regina::Example<3>::figureEight();

    identify::resetRecognitionCacheForTesting();
    size_t defaultLimit = identify::recognitionCacheLimit.load();
    identify::recognitionCacheLimit.store(1);

    identify::identify(EdgeComplement(boundary, {boundary.edge(0)}));
    identify::identify(EdgeComplement(figureEight, {}));

    EXPECT_EQ(identify::recognitionCacheStats().cacheResets >= 1, true,
              "recognitionCache reset at least once after exceeding its "
              "(deliberately tiny) limit");

    identify::recognitionCacheLimit.store(defaultLimit);
    identify::resetRecognitionCacheForTesting();
}

// identify::identify(const Link&): the split-unlink fast path
// (groupProvesUnlink(), identifycomplement.cpp), generalizing
// identify(const EdgeComplement&)'s genus-1/"Unknot" check from n == 1 to
// any n via free-group recognition rather than a handlebody genus check
// (a split n-component unlink's complement has n SEPARATE torus boundary
// components, so -- unlike a solid torus -- it is never itself a
// handlebody; recogniseHandlebody() correctly returns -1 for it, not n,
// which is what makes a genus-based check wrong here).
void test_identify_link_unlink() {
    // The 2-component unlink: the Hopf link's shadow with crossing 1's
    // tuple cyclically rotated by one position, flipping that crossing's
    // over/under role -- the exact fixture knotbuilder_test.cpp's own
    // "non-alternating regression" test uses, already independently
    // checked there (isSphere(), 2 components) to be two split, unknotted
    // loops.
    knotbuilder::PDCode pd = {{0, 3, 1, 2}, {1, 3, 0, 2}};
    auto [tri, edges, reversed] = knotbuilder::buildLink(pd);
    Link link(tri, edges);

    EXPECT_EQ(link.countComponents(), 2,
              "fixture sanity: the flipped Hopf shadow has 2 components");
    EXPECT_EQ(identify::identify(link), std::string("2-component unlink"),
              "identify(const Link&) names a split 2-component unlink via "
              "the free-fundamental-group fast path, instead of falling "
              "back to a bare isoSig");
}

// The flip side of the above: identify() must not mistake a genuinely
// LINKED multi-component complement for a handlebody either. The
// (unflipped) Hopf link has the same component count as the unlink fixture
// above, but its complement (T^2 x I) is not a handlebody at all.
void test_identify_link_hopf_not_unknot() {
    knotbuilder::PDCode pd = knotbuilder::parsePDCode("1 4 2 3 3 2 4 1");
    auto [tri, edges, reversed] = knotbuilder::buildLink(pd);
    Link link(tri, edges);

    EXPECT_EQ(link.countComponents(), 2,
              "fixture sanity: the Hopf link has 2 components");
    EXPECT_EQ(identify::identify(link) != "Unknot", true,
              "the (linked) Hopf link's complement is not a handlebody, "
              "so identify() must not call it \"Unknot\"");
}

} // namespace

void run(const std::string &name, void (*fn)()) {
    std::cout << bold << "\n=== " << name << " ===" << resetColor << "\n";
    fn();
}

int main() {
    run("memoization", test_memoization);
    run("automorphism_invariance", test_automorphism_invariance);
    run("boundary_signature_cache_clear_threshold",
        test_boundary_signature_cache_clear_threshold);
    run("recognition_cache_clear_threshold",
        test_recognition_cache_clear_threshold);
    run("identify_link_unlink", test_identify_link_unlink);
    run("identify_link_hopf_not_unknot", test_identify_link_hopf_not_unknot);

    std::cout << bold << "\n=== Summary: " << passed << " passed, "
              << failed_count << " failed ===" << resetColor << "\n";
    return failed_count > 0 ? 1 : 0;
}
