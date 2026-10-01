// vertexlinks_test.cpp
//
// Direct unit tests of PetalCache, independent of KnottedSurface: petal
// identity is order-independent (a corner set interns to the same id
// regardless of insertion order), distinct corner sets get distinct ids,
// unknot/linking lookups correctly report miss-then-hit, and the pairwise
// linking cache is order-independent in its two petal ids.

#include <iostream>
#include <string>
#include <unistd.h>

#include "../vertexlinks.h"

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

#define EXPECT_TRUE(actual, desc) EXPECT_EQ((actual), true, desc)

void test_intern_same_corners_different_order_same_id() {
    std::cout << "\n--- internPetal(): order-independent identity ---\n";

    PetalCache cache;
    auto idA = cache.internPetal({{3, 0}, {1, 1}, {2, 2}});
    auto idB = cache.internPetal({{1, 1}, {2, 2}, {3, 0}});
    EXPECT_EQ(idA, idB,
              "the same corner set interns to the same id regardless of "
              "insertion order");
}

void test_intern_distinct_corners_distinct_ids() {
    std::cout << "\n--- internPetal(): distinct petals get distinct ids ---\n";

    PetalCache cache;
    auto idA = cache.internPetal({{1, 0}, {2, 0}});
    auto idB = cache.internPetal({{1, 0}, {2, 1}}); // differs only in local vertex
    auto idC = cache.internPetal({{1, 0}, {3, 0}}); // differs only in face
    EXPECT_EQ(idA == idB, false, "differing local vertex gives a distinct id");
    EXPECT_EQ(idA == idC, false, "differing face gives a distinct id");
}

void test_unknot_cache_miss_then_hit() {
    std::cout << "\n--- lookupUnknot()/recordUnknot(): miss then hit ---\n";

    PetalCache cache;
    auto id = cache.internPetal({{1, 0}, {2, 1}, {3, 2}});

    EXPECT_EQ(cache.lookupUnknot(id).has_value(), false,
              "a freshly-interned petal has no cached isUnknot() result yet");
    cache.recordUnknot(id, true);
    auto cached = cache.lookupUnknot(id);
    EXPECT_TRUE(cached.has_value(), "isUnknot() result is cached after recordUnknot()");
    EXPECT_EQ(*cached, true, "the cached value matches what was recorded");
}

void test_linking_cache_order_independent() {
    std::cout << "\n--- lookupLinksNonzero()/recordLinksNonzero(): order-independent pair key ---\n";

    PetalCache cache;
    auto idA = cache.internPetal({{1, 0}});
    auto idB = cache.internPetal({{2, 0}});

    EXPECT_EQ(cache.lookupLinksNonzero(idA, idB).has_value(), false,
              "no cached linking result before recording one");
    cache.recordLinksNonzero(idA, idB, true);

    auto viaAB = cache.lookupLinksNonzero(idA, idB);
    auto viaBA = cache.lookupLinksNonzero(idB, idA);
    EXPECT_TRUE(viaAB.has_value(), "(a, b) hits after recording");
    EXPECT_TRUE(viaBA.has_value(), "(b, a) hits the same entry as (a, b)");
    EXPECT_EQ(*viaAB, true, "(a, b) reports the recorded value");
    EXPECT_EQ(*viaBA, true, "(b, a) reports the same recorded value");
}

void test_stats_track_hits_and_misses() {
    std::cout << "\n--- Stats: unknot/linking checks and cache hits are counted ---\n";

    PetalCache cache;
    auto idA = cache.internPetal({{1, 0}});
    auto idB = cache.internPetal({{2, 0}});

    cache.lookupUnknot(idA);   // miss
    cache.recordUnknot(idA, false);
    cache.lookupUnknot(idA);   // hit

    cache.lookupLinksNonzero(idA, idB); // miss
    cache.recordLinksNonzero(idA, idB, false);
    cache.lookupLinksNonzero(idA, idB); // hit

    const auto &stats = cache.stats();
    EXPECT_EQ(stats.unknotChecks, 2LL, "two lookupUnknot() calls were made");
    EXPECT_EQ(stats.unknotCacheHits, 1LL, "exactly one of them was a cache hit");
    EXPECT_EQ(stats.linkingChecks, 2LL, "two lookupLinksNonzero() calls were made");
    EXPECT_EQ(stats.linkingCacheHits, 1LL, "exactly one of them was a cache hit");
}

void test_clear_threshold_triggers_reset() {
    std::cout << "\n--- setClearThreshold(): a full clear triggers once over the threshold ---\n";

    PetalCache cache;
    cache.setClearThreshold(3);
    cache.internPetal({{1, 0}});
    cache.internPetal({{2, 0}});
    cache.internPetal({{3, 0}}); // fills the cache to its threshold
    EXPECT_EQ(cache.stats().cacheResets, 0LL,
              "no reset yet -- the threshold is only checked before minting the NEXT new id");
    cache.internPetal({{4, 0}}); // this one triggers the clear before minting id 3
    EXPECT_EQ(cache.stats().cacheResets, 1LL,
              "exactly one reset after exceeding the threshold");
    EXPECT_EQ(cache.size(), size_t(1),
              "post-clear, only the petal that triggered the clear is present");
}

void test_stale_id_lookup_is_clean_miss_not_wrong_hit() {
    std::cout << "\n--- lookupUnknot(): a pre-clear id is a clean miss, never a wrong hit ---\n";

    PetalCache cache;
    cache.setClearThreshold(1);
    auto idBeforeClear = cache.internPetal({{1, 0}});
    cache.recordUnknot(idBeforeClear, true);

    // The next internPetal() call is guaranteed to trigger a clear (threshold
    // is 1, and one entry is already present), then mint a fresh id -- which,
    // pre-epoch-tagging, would numerically collide with idBeforeClear.
    auto idAfterClear = cache.internPetal({{2, 0}});
    EXPECT_EQ(cache.stats().cacheResets, 1LL, "the second intern triggered the clear");

    EXPECT_EQ(cache.lookupUnknot(idBeforeClear).has_value(), false,
              "the pre-clear id is recognized as stale (clean miss), not "
              "reinterpreted against the post-clear table");
    EXPECT_EQ(cache.lookupUnknot(idAfterClear).has_value(), false,
              "the post-clear id correctly has no recorded value yet either");
}

void test_stale_id_record_is_safe_noop() {
    std::cout << "\n--- recordUnknot()/recordLinksNonzero(): a pre-clear id is a safe no-op ---\n";

    PetalCache cache;
    cache.setClearThreshold(1);
    auto idBeforeClear = cache.internPetal({{1, 0}});
    auto idAfterClear = cache.internPetal({{2, 0}}); // triggers the clear

    // Recording against the stale id must not corrupt whatever now lives at
    // that same numeric local-id slot post-clear.
    cache.recordUnknot(idBeforeClear, false);
    EXPECT_EQ(cache.lookupUnknot(idAfterClear).has_value(), false,
              "recording against a stale id left the post-clear id's entry untouched");

    cache.recordLinksNonzero(idBeforeClear, idAfterClear, true);
    EXPECT_EQ(cache.lookupLinksNonzero(idAfterClear, idAfterClear).has_value(),
              false,
              "recording a linking result with a stale id was a safe no-op");
}

void run(const std::string &name, void (*fn)()) {
    std::cout << bold << "\n=== " << name << " ===" << resetColor << "\n";
    fn();
}

int main() {
    run("intern_same_corners_different_order_same_id",
        test_intern_same_corners_different_order_same_id);
    run("intern_distinct_corners_distinct_ids",
        test_intern_distinct_corners_distinct_ids);
    run("unknot_cache_miss_then_hit", test_unknot_cache_miss_then_hit);
    run("linking_cache_order_independent",
        test_linking_cache_order_independent);
    run("stats_track_hits_and_misses", test_stats_track_hits_and_misses);
    run("clear_threshold_triggers_reset", test_clear_threshold_triggers_reset);
    run("stale_id_lookup_is_clean_miss_not_wrong_hit",
        test_stale_id_lookup_is_clean_miss_not_wrong_hit);
    run("stale_id_record_is_safe_noop", test_stale_id_record_is_safe_noop);

    std::cout << bold << "\n=== Summary: " << passed << " passed, "
              << failed_count << " failed ===" << resetColor << "\n";
    return failed_count > 0 ? 1 : 0;
}
