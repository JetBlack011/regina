// censusnaming_test.cpp
//
// Tests for census::identify() and the recognition cache behind it (see
// ../census/censusnaming.h and ../complement/complementcache.h): the cache's
// limit and its behaviour under concurrent clears, the split-unlink fast
// path, and the Pachner-search policy.
//
// Split from identifycomplement_test.cpp (refactor phase 2), whose
// BoundarySignatureCache tests went to surfer's namecache_test.cpp and whose
// certifiesUnlink()/capInCone() tests went to unlinknaming_test.cpp.

#include <atomic>
#include <chrono>
#include <iostream>
#include <string>
#include <thread>
#include <unistd.h>
#include <vector>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <triangulation/example3.h>

#include "linknaming/census/censusnaming.h"
#include "linknaming/complement/complementcache.h"
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

void test_complement_cache_clear_threshold() {
    // Unlike BoundarySignatureCache's pre-triangulation combinatorial key,
    // recognitionCache is keyed by the drilled complement's POST-simplify
    // isoSig -- a topological invariant of the resulting manifold, not of
    // how many edges were drilled. So this needs two genuinely
    // topologically distinct complements, not just different edge counts.
    // Drilling no edges at all from a fixed triangulation just recognizes
    // that triangulation itself (buildComplement()'s pinch loop is a no-op
    // on an empty edge set) -- so pairing the pentachoron-boundary case (one
    // edge with distinct ends, which pinching just collapses, leaving S^3)
    // with the figure-eight knot complement (genuinely hyperbolic, not any
    // handlebody) guarantees two distinct isoSigs.
    regina::Triangulation<3> boundary = testBoundary();
    regina::Triangulation<3> figureEight = regina::Example<3>::figureEight();

    complement::resetCacheForTesting();
    size_t defaultLimit = complement::cacheLimit.load();
    complement::cacheLimit.store(1);

    census::nameComplement(EdgeComplement(boundary, {boundary.edge(0)}));
    census::nameComplement(EdgeComplement(figureEight, {}));

    EXPECT_EQ(complement::cacheStats().cacheResets >= 1, true,
              "recognitionCache reset at least once after exceeding its "
              "(deliberately tiny) limit");

    complement::cacheLimit.store(defaultLimit);
    complement::resetCacheForTesting();
}

void test_complement_cache_clear_race() {
    // resolveRecognition() takes the genus from cachedGenus(), which stores
    // it and releases recognitionCacheMutex, and then dereferences
    // lookupRecognition(sig) under a second lock. storeRecognition() clears
    // the whole cache once it is at recognitionCacheLimit, so another
    // thread's store can land between the two and leave an empty optional.
    // With the limit at 1, every store of a different signature clears it:
    // four threads name one complement while four others churn another.
    // Each answer must be the one a quiet, single-threaded call gives. The
    // churners also clear the cache directly, as a store at the limit does,
    // so the narrow window between the two lookups is actually reached.
    // The three edges of one triangle: an unknotted circle, whose complement
    // is a solid torus, so resolveRecognition() takes the genus != -1 branch
    // that dereferences the second lookup. (A single edge with distinct ends
    // just collapses, leaving S^3, which takes the census branch instead.)
    auto unknot = [](const regina::Triangulation<3> &t) {
        const regina::Triangle<3> *f = t.triangle(0);
        return EdgeComplement(t, {f->edge(0), f->edge(1), f->edge(2)});
    };
    regina::Triangulation<3> boundary = testBoundary();
    complement::resetCacheForTesting();
    const std::string expected = census::nameComplement(unknot(boundary));
    EXPECT_EQ(expected, std::string("Unknot"),
              "the quiet call takes the genus != -1 branch");

    const size_t defaultLimit = complement::cacheLimit.load();
    complement::cacheLimit.store(1);
    complement::resetCacheForTesting();

    std::atomic<long> calls{0}, wrong{0};
    std::atomic<bool> stop{false};
    const auto deadline =
        std::chrono::steady_clock::now() + std::chrono::seconds(20);
    auto namer = [&] {
        regina::Triangulation<3> own = testBoundary();
        for (int i = 0; i < 4000 && !stop.load(); ++i) {
            std::string name = census::nameComplement(unknot(own));
            ++calls;
            if (name != expected && ++wrong <= 5)
                std::cout << "  wrong name: '" << name << "'\n";
            if (std::chrono::steady_clock::now() > deadline)
                stop.store(true);
        }
    };
    auto churner = [&] {
        regina::Triangulation<3> own = regina::Example<3>::figureEight();
        while (!stop.load()) {
            census::nameComplement(EdgeComplement(own, {}));
            for (int k = 0; k < 1000 && !stop.load(); ++k)
                complement::resetCacheForTesting();
            if (std::chrono::steady_clock::now() > deadline)
                stop.store(true);
        }
    };
    std::vector<std::thread> threads;
    for (int t = 0; t < 4; ++t)
        threads.emplace_back(namer);
    for (int t = 0; t < 4; ++t)
        threads.emplace_back(churner);
    for (int t = 0; t < 4; ++t)
        threads[t].join();
    stop.store(true);
    for (size_t t = 4; t < threads.size(); ++t)
        threads[t].join();

    std::cout << "  " << calls.load() << " calls naming '" << expected
              << "'\n";
    EXPECT_EQ(wrong.load(), 0L,
              "every concurrent identification agrees with the quiet one, "
              "however often the cache is cleared underneath it");

    complement::cacheLimit.store(defaultLimit);
    complement::resetCacheForTesting();
}

// census::identify(const Link&): the split-unlink fast path
// (groupProvesUnlink(), identifycomplement.cpp), generalizing
// identify(const EdgeComplement&)'s genus-1/"Unknot" check from n == 1 to
// any n via free-group recognition rather than a handlebody genus check
// (a split n-component unlink's complement has n SEPARATE torus boundary
// components, so -- unlike a solid torus -- it is never itself a
// handlebody; recogniseHandlebody() correctly returns -1 for it, not n,
// which is what makes a genus-based check wrong here).
void test_name_link_unlink() {
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
    EXPECT_EQ(census::nameComplement(link), std::string("2-component unlink"),
              "identify(const Link&) names a split 2-component unlink via "
              "the free-fundamental-group fast path, instead of falling "
              "back to a bare isoSig");
}

// The flip side of the above: identify() must not mistake a genuinely
// LINKED multi-component complement for a handlebody either. The
// (unflipped) Hopf link has the same component count as the unlink fixture
// above, but its complement (T^2 x I) is not a handlebody at all.
void test_name_link_hopf_not_unknot() {
    knotbuilder::PDCode pd = knotbuilder::parsePDCode("1 4 2 3 3 2 4 1");
    auto [tri, edges, reversed] = knotbuilder::buildLink(pd);
    Link link(tri, edges);

    EXPECT_EQ(link.countComponents(), 2,
              "fixture sanity: the Hopf link has 2 components");
    EXPECT_EQ(census::nameComplement(link) != "Unknot", true,
              "the (linked) Hopf link's complement is not a handlebody, "
              "so identify() must not call it \"Unknot\"");
}

// The Pachner-search policy. With no local census, every non-trivial
// complement misses the cheap rungs, so whether retriangulateAndLookup() runs
// is decided by the policy alone: never for a link complement (its name bears
// no bound during a search) unless retriangulateLinks is set, and still for
// a knot's. Regina's census holds hyperbolic manifolds only, so the trefoil's
// Seifert-fibred complement misses it too.
void test_pachner_search_policy() {
    census::resetCensusForTesting();
    complement::resetCacheForTesting();
    const bool oldOnMiss = census::retriangulateOnMiss.load();
    const bool oldLinks = census::retriangulateLinks.load();
    const long long oldBudget = census::retriangulateTimeBudgetSeconds.load();
    census::retriangulateOnMiss.store(true);
    census::retriangulateLinks.store(false);
    census::retriangulateTimeBudgetSeconds.store(1);

    {
        knotbuilder::PDCode pd = knotbuilder::parsePDCode("1 4 2 3 3 2 4 1");
        auto [tri, edges, reversed] = knotbuilder::buildLink(pd);
        Link hopf(tri, edges);
        census::nameComplement(hopf);
        census::nameComplement(hopf); // a repeat must not retry either
    }
    auto afterLink = complement::cacheStats();
    EXPECT_EQ(afterLink.pachnerLinks.attempts, 0LL,
              "a link complement that misses the census is NOT sent to the "
              "Pachner search");

    {
        knotbuilder::PDCode pd =
            knotbuilder::parsePDCode("[[1;5;2;4];[3;1;4;6];[5;3;6;2]]");
        auto [tri, edges, reversed] = knotbuilder::buildLink(pd);
        Link trefoil(tri, edges);
        census::nameComplement(trefoil);
    }
    auto afterKnot = complement::cacheStats();
    EXPECT_EQ(afterKnot.pachnerKnots.attempts >= 1, true,
              "a knot complement that misses the census still is: a knot's "
              "name can bear a bound");
    EXPECT_EQ(afterKnot.pachnerLinks.attempts, 0LL,
              "and no link attempt crept in");

    census::retriangulateOnMiss.store(oldOnMiss);
    census::retriangulateLinks.store(oldLinks);
    census::retriangulateTimeBudgetSeconds.store(oldBudget);
    complement::resetCacheForTesting();
}

} // namespace

void run(const std::string &name, void (*fn)()) {
    std::cout << bold << "\n=== " << name << " ===" << resetColor << "\n";
    fn();
}

int main() {
    run("recognition_cache_clear_threshold",
        test_complement_cache_clear_threshold);
    run("recognition_cache_clear_race", test_complement_cache_clear_race);
    run("identify_link_unlink", test_name_link_unlink);
    run("identify_link_hopf_not_unknot", test_name_link_hopf_not_unknot);
    run("pachner_search_policy", test_pachner_search_policy);


    std::cout << bold << "\n=== Summary: " << passed << " passed, "
              << failed_count << " failed ===" << resetColor << "\n";
    return failed_count > 0 ? 1 : 0;
}
