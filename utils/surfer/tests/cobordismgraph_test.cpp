// cobordismgraph_test.cpp
//
// Tests for ../cobordismgraph.h/.cpp: the interval solver verifyslicegenus.cpp
// runs over its accumulated witness set (componentsFromName(), NameTable,
// propagate(), judge()) and the boundary classification (splitBoundary())
// that feeds it. All pure logic on plain data types -- no triangulations, no
// search, no census -- so these tests run instantly and exercise the core of
// verifyslicegenus's genus deductions directly.
//
// The cases below deliberately pin down the two things that are easy to get
// silently wrong: the component-count terms in the cobordism inequality
// (which vanish for knots, so a knots-only test suite would never notice
// them missing), and the max/min over an orientation-ambiguous far side's
// candidate set.

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

// Builds a cobordism witness between `subject` (with `nA` components) and
// `other` (with `nB`), at genus `g`. `candidates` defaults to just `other`.
Witness cobordism(const std::string &subject, int nA, const std::string &other,
                  int nB, int g,
                  std::vector<std::string> candidates = {}) {
    Witness w;
    w.kind = WitnessKind::cobordism;
    w.subject = subject;
    w.subjectComponents = nA;
    w.other = other;
    w.otherComponents = nB;
    w.genus = g;
    w.otherCandidates = candidates.empty()
                            ? std::vector<std::string>{other}
                            : std::move(candidates);
    return w;
}

Witness direct(const std::string &subject, int nA, int g, bool tubed = false) {
    Witness w;
    w.kind = WitnessKind::direct;
    w.subject = subject;
    w.subjectComponents = nA;
    w.genus = g;
    w.tubed = tubed;
    return w;
}

// ─────────────────────────────────────────────────────────────────────────
// Component counts read off a name
// ─────────────────────────────────────────────────────────────────────────

void test_components_from_name() {
    EXPECT_EQ(componentsFromName("3_1"), 1, "a knot name is one component");
    EXPECT_EQ(componentsFromName("Unknot"), 1, "the unknot is one component");
    EXPECT_EQ(componentsFromName("cPcbbbadu"), 1,
              "a bare isoSig fallback is treated as one component");
    // LinkInfo's tag carries an orientation choice per component AFTER the
    // first, so the count is one more than the number of entries.
    EXPECT_EQ(componentsFromName("L2a1{0}"), 2, "L2a1 is the Hopf link");
    EXPECT_EQ(componentsFromName("L6a5{0;1}"), 3,
              "L6a5 (Borromean) is three components");
    EXPECT_EQ(componentsFromName("L11n459{1;0;0}"), 4, "four components");
    EXPECT_EQ(componentsFromName("2-component unlink"), 2,
              "an unlink states its own component count");
    EXPECT_EQ(componentsFromName("5-component unlink"), 5, "");
}

void test_base_name() {
    EXPECT_EQ(baseName("L6a3{0}"), std::string("L6a3"), "tag stripped");
    EXPECT_EQ(baseName("L6a5{0;1}"), std::string("L6a5"), "multi-entry tag");
    EXPECT_EQ(baseName("10_132"), std::string("10_132"), "no tag, unchanged");
}

// ─────────────────────────────────────────────────────────────────────────
// NameTable: turning an orientation-blind name into candidates
// ─────────────────────────────────────────────────────────────────────────

void test_candidates_expand_orientation_variants() {
    NameTable names;
    names.addLiterature("L4a1{0}", 0, 0);
    names.addLiterature("L4a1{1}", 1, 1);
    names.addLiterature("4_1", 1, 1);

    auto cands = names.candidates("L4a1");
    EXPECT_EQ(cands.size(), static_cast<size_t>(2),
              "a complement-derived name expands to every oriented variant "
              "it could be -- the complement cannot tell them apart");

    auto knotCands = names.candidates("4_1");
    EXPECT_EQ(knotCands.size(), static_cast<size_t>(1),
              "a knot name has no oriented variants to disambiguate");
    EXPECT_EQ(knotCands[0], std::string("4_1"), "");

    auto unknown = names.candidates("cPcbbbadu");
    EXPECT_EQ(unknown.size(), static_cast<size_t>(1),
              "an unregistered name (a bare isoSig) is its own only "
              "candidate");
}

void test_candidates_filtered_by_observed_component_count() {
    NameTable names;
    names.addLiterature("L4a1{0}", 0, 0);  // 2 components
    names.addLiterature("L4a1{1}", 1, 1);  // 2 components
    names.addLiterature("L4a1{0;0}", 3, 3); // contrived 3-component variant

    auto cands = names.candidates("L4a1", 2);
    EXPECT_EQ(cands.size(), static_cast<size_t>(2),
              "a variant whose component count disagrees with what was "
              "actually observed on the boundary simply isn't what was "
              "found, so it is filtered out");
}

// ─────────────────────────────────────────────────────────────────────────
// propagate(): the cobordism inequality
// ─────────────────────────────────────────────────────────────────────────

void test_direct_witness_gives_upper_bound() {
    NameTable names;
    names.addLiterature("K", 1, 1);
    auto bounds = propagate({direct("K", 1, 1)}, names);

    EXPECT_EQ(bounds["K"].haveUpper(), true, "a direct witness bounds K");
    EXPECT_EQ(bounds["K"].hi, 1, "at exactly the surface's own genus");
    EXPECT_EQ(bounds["K"].basis == Basis::constructive, true,
              "a direct witness depends on nothing else, so the bound is "
              "constructive -- nothing taken on faith");
}

void test_knot_cobordism_reduces_to_the_classic_rule() {
    // With n_0 = n_1 = 1 the component terms vanish and the inequality is
    // the familiar |g_4(K_0) - g_4(K_1)| <= g. This is the regression guard
    // for the whole knot-only code path.
    NameTable names;
    names.addLiterature("K", 0, 5);
    names.addLiterature("J", 2, 2);
    auto bounds = propagate({direct("J", 1, 2), cobordism("K", 1, "J", 1, 3)},
                            names);

    EXPECT_EQ(bounds["K"].hi, 5, "g_4(K) <= g_4(J) + g = 2 + 3, the n=1 rule");
    EXPECT_EQ(bounds["K"].haveLower(), false,
              "g_4(K) >= g_4(J) - g = 2 - 3 = -1 is no constraint at all "
              "(a genus is never negative), so nothing is recorded");
}

void test_component_correction_on_the_upper_bound() {
    // THE case the old knot-only formula got wrong. A genus-0 cobordism
    // from a knot to a genuinely linked 2-component far side does NOT make
    // the knot slice: that far side must be capped with a CONNECTED surface,
    // which costs + n_1 - 1 = 1 in genus.
    //
    // Deliberately not an unlink. An unlink far side caps with disjoint
    // discs instead and carries no penalty at all -- see
    // test_unlink_far_side_carries_no_component_penalty, which is the case
    // this test originally (and wrongly) used.
    NameTable names;
    names.addLiterature("K", 1, 1);
    names.addLiterature("L4a1{1}", 0, 0);
    auto bounds = propagate({cobordism("K", 1, "L4a1{1}", 2, 0)}, names);

    EXPECT_EQ(bounds["K"].hi, 1,
              "g_4(K) <= 0 + 0 + (2 - 1) = 1. The uncorrected rule would "
              "give 0 here, i.e. would 'prove' K slice on evidence that "
              "says no such thing");
}

void test_component_correction_on_the_lower_bound() {
    // The mirror case: deducing a lower bound for a MULTI-component
    // subject costs - n_0 + 1.
    NameTable names;
    names.addLiterature("L", 0, 9); // 2 components, per the witness below
    names.addLiterature("K", 3, 3);
    auto bounds = propagate({cobordism("L", 2, "K", 1, 0)}, names);

    EXPECT_EQ(bounds["L"].lo, 2,
              "g_4(L) >= g_4(K) - g - n_0 + 1 = 3 - 0 - 2 + 1 = 2. The "
              "uncorrected rule would claim 3, a bound that is simply not "
              "justified");
}

void test_unlink_far_side_carries_no_component_penalty() {
    // The `+ n_far - 1` term assumes the far side is capped with its minimal
    // CONNECTED surface. An n-component unlink instead bounds n DISJOINT
    // discs, and gluing those onto a connected cobordism still yields a
    // connected surface -- so chi* = n, not 2 - n, and the penalty vanishes.
    //
    // Worked case, found live: a connected genus-0 cobordism from 6_1 to the
    // 2-component unlink has chi(Sigma) = 2 - 0 - 1 - 2 = -1. Capping with
    // two discs gives chi = 1, i.e. 2 - 2G - 1 = 1, i.e. G = 0 -- exactly
    // 6_1's slice disc. Scoring it as genus 1 (the connected cap) made it
    // exceed 6_1's literature value of 0, so it was thrown away.
    NameTable names;
    names.addLiterature("6_1", 0, 0);
    auto bounds =
        propagate({cobordism("6_1", 1, "2-component unlink", 2, 0)}, names);

    EXPECT_EQ(bounds["6_1"].hi, 0,
              "a genus-0 cobordism to the 2-component unlink certifies that "
              "6_1 is slice: the penalty is 0, not n - 1");
    EXPECT_EQ(bounds["6_1"].basis == Basis::constructive, true,
              "the unlink is an axiom, so nothing is taken on faith");

    Verdict v = judge("6_1", bounds["6_1"], names);
    EXPECT_EQ(v.status == Status::verified, true, "and the row verifies");
}

void test_non_unlink_far_side_keeps_the_penalty() {
    // The correction is specific to unlinks. A genuinely linked far side has
    // no disjoint-disc cap available, so the connected-surface penalty stands.
    NameTable names;
    names.addLiterature("K", 0, 9);
    names.addLiterature("L4a1{0}", 0, 0);
    auto bounds = propagate({cobordism("K", 1, "L4a1{0}", 2, 0)}, names);
    EXPECT_EQ(bounds["K"].hi, 1,
              "0 + 0 + (2 - 1) = 1 -- unchanged for a non-unlink far side");
}

void test_unlink_axiom_is_constructive() {
    NameTable names;
    names.addLiterature("K", 0, 9); // wide enough that the derived 2 is news
    auto bounds =
        propagate({cobordism("K", 1, "3-component unlink", 3, 0)}, names);

    EXPECT_EQ(bounds["3-component unlink"].hi, 0,
              "an n-component unlink bounds n discs, which tube into a "
              "connected planar surface of genus 0");
    EXPECT_EQ(bounds["K"].basis == Basis::constructive, true,
              "the unlink axiom is a fact, not a literature value, so a "
              "bound resting on it stays constructive");
    EXPECT_EQ(bounds["K"].hi, 0,
              "0 + 0 + 0 = 0: an unlink far side caps with n DISJOINT discs, "
              "so no component penalty applies -- a genus-0 cobordism to any "
              "unlink certifies the subject is slice");
}

void test_slice_composite_axiom_is_constructive() {
    NameTable names;
    names.addLiterature("K", 0, 9); // wide enough that the derived 0 is news
    auto bounds = propagate({cobordism("K", 1, "3_1#m3_1", 1, 0)}, names);

    EXPECT_EQ(bounds["3_1#m3_1"].hi, 0,
              "K # m(K^r) is the identity of the concordance group and bounds "
              "an explicit ribbon disc, so it grounds a chain like the unknot");
    EXPECT_EQ(bounds["K"].basis == Basis::constructive, true,
              "the ribbon disc is a theorem, not a literature value, so a "
              "bound resting on it stays constructive");
    EXPECT_EQ(bounds["K"].hi, 0,
              "0 + 0 + (1 - 1) = 0: a genus-0 cobordism to a slice knot "
              "certifies the subject is slice");
}

void test_slice_composite_allowlist_is_not_a_pattern() {
    // The rule is K # m(K^r), and which SPELLING satisfies it depends on the
    // summand's symmetry. 8_17 is the first non-invertible knot, so
    // 8_17#m8_17 is NOT the concordance inverse and is not known to be slice.
    // Anything matching on the shape "A#mA" would wrongly axiom it.
    NameTable names;
    names.addLiterature("K", 0, 9);
    auto bounds = propagate({cobordism("K", 1, "8_17#m8_17", 1, 0)}, names);

    EXPECT_EQ(bounds["8_17#m8_17"].hi != 0, true,
              "8_17 is non-invertible, so 8_17#m8_17 is not the ribbon case "
              "and must not be axiomed -- the allowlist is not a pattern");
    EXPECT_EQ(bounds.contains("K") && bounds["K"].haveUpper(), false,
              "with no bound on the far side there is nothing to propagate");
}

// ─────────────────────────────────────────────────────────────────────────
// propagate(): orientation-ambiguous far sides
// ─────────────────────────────────────────────────────────────────────────

void test_candidate_set_takes_the_worst_case_for_an_upper_bound() {
    // L4a1{0} has slice genus 0 and L4a1{1} has 1, and they share one
    // complement -- so a far side identified only as "L4a1" must be bounded
    // as though it were the worse of the two.
    NameTable names;
    names.addLiterature("K", 0, 9);
    names.addLiterature("L4a1{0}", 0, 0);
    names.addLiterature("L4a1{1}", 1, 1);

    auto bounds = propagate(
        {cobordism("K", 1, "L4a1", 2, 0, {"L4a1{0}", "L4a1{1}"})}, names);

    EXPECT_EQ(bounds["K"].hi, 2,
              "max(0, 1) + 0 + (2 - 1) = 2 -- sound whichever variant it "
              "actually was, where taking the min would be wrong and "
              "discarding the cobordism would waste it");
    EXPECT_EQ(bounds["K"].basis == Basis::literatureAssisted, true,
              "the candidates' own bounds came from the literature, so the "
              "conclusion is conditional on it and says so");
}

void test_candidate_set_takes_the_best_case_for_a_lower_bound() {
    NameTable names;
    names.addLiterature("K", 0, 9);
    names.addLiterature("L4a1{0}", 2, 2);
    names.addLiterature("L4a1{1}", 5, 5);

    auto bounds = propagate(
        {cobordism("K", 1, "L4a1", 2, 0, {"L4a1{0}", "L4a1{1}"})}, names);

    EXPECT_EQ(bounds["K"].lo, 2,
              "min(2, 5) - 0 - (1 - 1) = 2 -- the weaker of the two, since "
              "we don't know which variant we found");
}

void test_unambiguous_base_costs_nothing() {
    // About 46% of the base names in the links table have every variant at
    // the same genus, so the max/min collapse and the ambiguity is free.
    NameTable names;
    names.addLiterature("K", 0, 9);
    names.addLiterature("L2a1{0}", 0, 0);
    names.addLiterature("L2a1{1}", 0, 0);

    auto bounds = propagate(
        {cobordism("K", 1, "L2a1", 2, 0, {"L2a1{0}", "L2a1{1}"})}, names);

    EXPECT_EQ(bounds["K"].hi, 1,
              "both variants agree at 0, so the bound is exactly what an "
              "unambiguous identification would have given: 0 + 0 + 1");
}

void test_one_unbounded_candidate_blocks_the_upper_bound() {
    NameTable names;
    // K's own literature is deliberately left unregistered here: with it
    // registered, K's literature upper bound would be a perfectly valid
    // (if useless) bound, and the point of this test is the DERIVED one.
    names.addLiterature("Kother", 0, 9);
    names.addLiterature("L4a1{0}", 0, 0);
    // L4a1{1} deliberately unregistered: nothing bounds it, so nothing can
    // bound K through a far side that might BE it.
    auto bounds = propagate(
        {cobordism("K", 1, "L4a1", 2, 0, {"L4a1{0}", "L4a1{1}"})}, names);

    EXPECT_EQ(bounds["K"].haveUpper(), false,
              "an upper bound must hold for EVERY candidate; one with no "
              "known bound leaves the whole deduction unavailable");
}

// ─────────────────────────────────────────────────────────────────────────
// propagate(): chaining
// ─────────────────────────────────────────────────────────────────────────

void test_chains_through_an_unnamed_isosig_node() {
    // identify() falls back to a bare isoSig when the census misses. Such a
    // node has no bounds of its own, but it still CONNECTS: two rows that
    // each cobound with it are thereby related to each other. Keeping these
    // nodes in the graph is free information.
    NameTable names;
    names.addLiterature("X", 0, 9);
    names.addLiterature("Y", 0, 0);

    std::vector<Witness> witnesses = {
        direct("Y", 1, 0),
        cobordism("Y", 1, "gLLMQacdefeffhhnkxk", 1, 0),
        cobordism("X", 1, "gLLMQacdefeffhhnkxk", 1, 1),
    };
    auto bounds = propagate(witnesses, names);

    EXPECT_EQ(bounds["gLLMQacdefeffhhnkxk"].hi, 0,
              "the unnamed node picks up Y's own bound through their "
              "shared cobordism");
    EXPECT_EQ(bounds["X"].hi, 1,
              "and passes it on to X, which never met Y directly");
    EXPECT_EQ(bounds["X"].basis == Basis::constructive, true,
              "every step rested on a surface we actually found");
}

void test_tubed_witness_genus_is_taken_at_face_value() {
    // A disconnected find is recorded with its TUBED genus, so the solver
    // needs no special case: two disjoint discs bounding a 2-component
    // unlink tube into a connected planar surface of genus 0.
    NameTable names;
    names.addLiterature("L", 0, 1);
    auto bounds = propagate({direct("L", 2, 0, /*tubed=*/true)}, names);

    EXPECT_EQ(bounds["L"].hi, 0, "genus 0, as tubed");
    EXPECT_EQ(bounds["L"].tubed, true,
              "recorded as tubed, so a reader knows the pair signature "
              "names a disconnected complex");
}

void test_propagation_terminates_on_a_cycle() {
    // Both relaxations are monotone and every step subtracts a
    // non-negative amount, so a cycle cannot pump a bound indefinitely.
    // If this ever regresses it hangs rather than failing, which is
    // exactly why it is worth pinning down.
    NameTable names;
    names.addLiterature("A", 0, 9);
    names.addLiterature("B", 0, 9);
    names.addLiterature("C", 0, 9);
    std::vector<Witness> witnesses = {
        direct("A", 1, 1),
        cobordism("A", 1, "B", 1, 0),
        cobordism("B", 1, "C", 1, 0),
        cobordism("C", 1, "A", 1, 0),
    };
    auto bounds = propagate(witnesses, names);
    EXPECT_EQ(bounds["B"].hi, 1, "the cycle settles at A's own bound");
    EXPECT_EQ(bounds["C"].hi, 1, "");
}

void test_self_cobordism_cannot_confirm_the_literature() {
    // Regression for a real bug. The search routinely finds a trivial
    // self-cobordism (a knot to its own collar image, genus 0). Bounding a
    // name through itself is vacuous -- but upperOf() falls back to a
    // name's LITERATURE upper bound when nothing is derived yet, so a
    // genus-0 self-loop would "derive" exactly the literature value and
    // then get reported as a verification of it. The literature confirming
    // itself, with no surface of ours involved.
    NameTable names;
    names.addLiterature("6_1", 0, 0);
    auto bounds = propagate({cobordism("6_1", 1, "6_1", 1, 0)}, names);

    EXPECT_EQ(bounds["6_1"].haveUpper(), false,
              "a self-cobordism bounds nothing -- it must NOT launder the "
              "literature value into a derived one");

    Verdict v = judge("6_1", bounds["6_1"], names);
    EXPECT_EQ(v.status == Status::unresolved, true,
              "and so the row stays honestly unresolved");
}

void test_self_cobordism_still_allows_a_real_bound_from_elsewhere() {
    // The self-loop is skipped, not poisonous: a genuine witness on the
    // same name still lands.
    NameTable names;
    names.addLiterature("K", 0, 5);
    auto bounds = propagate(
        {cobordism("K", 1, "K", 1, 0), direct("K", 1, 2)}, names);
    EXPECT_EQ(bounds["K"].hi, 2, "the direct witness still bounds K");
}

void test_ambiguous_far_side_gets_no_reverse_bound() {
    // Regression for a real contradiction caught by a live run. A genus-0
    // cobordism from L4a1{0} to a far side identified only as "L7n1"
    // bounds whichever variant the surface actually witnesses -- but not
    // the other one. Pushing the bound onto every candidate gave
    // L7n1{0} (true genus 2) a bound of 1.
    NameTable names;
    names.addLiterature("L4a1{0}", 0, 0);
    names.addLiterature("L7n1{0}", 2, 2);
    names.addLiterature("L7n1{1}", 0, 0);

    auto bounds = propagate(
        {direct("L4a1{0}", 2, 0),
         cobordism("L4a1{0}", 2, "L7n1", 2, 0, {"L7n1{0}", "L7n1{1}"})},
        names);

    EXPECT_EQ(bounds["L4a1{0}"].hi, 0, "the subject is still bounded");
    EXPECT_EQ(bounds["L7n1{0}"].haveUpper(), false,
              "an ambiguous far side receives NO bound: the surface "
              "witnesses exactly one of the candidates and we cannot tell "
              "which, so neither may be bounded individually");
    EXPECT_EQ(bounds["L7n1{1}"].haveUpper(), false, "likewise the other");
}

void test_unambiguous_far_side_does_get_a_reverse_bound() {
    // The restriction above is about ambiguity, not about direction: with
    // a single candidate the reverse bound is perfectly sound and must
    // still fire, or link-to-knot cobordisms would stop teaching us
    // anything about the knot.
    NameTable names;
    names.addLiterature("L", 0, 0);
    names.addLiterature("K", 0, 9);

    auto bounds = propagate(
        {direct("L", 2, 0), cobordism("L", 2, "K", 1, 0)}, names);

    EXPECT_EQ(bounds["K"].hi, 1,
              "g_4(K) <= g_4(L) + g + n_L - 1 = 0 + 0 + 2 - 1 = 1");
}

void test_two_step_cycle_through_an_alias_is_refused() {
    // Regression for a bug spotted in real output. The census names the
    // same knot twice -- 8_14 is also L108014, an untranslated Christy name
    // for a one-component complement -- so a trivial product annulus shows
    // up as a cobordism between two DIFFERENT-looking nodes. 8_14's own
    // literature value then flowed out to L108014 and came back as a
    // "derived" bound on 8_14, reported as having verified the literature.
    //
    // The one-step guard (never bound a name from itself) cannot see this;
    // the support set can.
    NameTable names;
    names.addLiterature("8_14", 1, 1);
    // L108014 deliberately absent from the tables: it is the same knot
    // under a second census name, so nothing bounds it independently.
    auto bounds = propagate({cobordism("8_14", 1, "L108014", 1, 0)}, names);

    EXPECT_EQ(bounds["8_14"].haveUpper(), false,
              "a bound on 8_14 that rests on 8_14's own literature value, "
              "however many hops away, is circular and must be refused");

    Verdict v = judge("8_14", bounds["8_14"], names);
    EXPECT_EQ(v.status == Status::unresolved, true,
              "so the row stays honestly unresolved rather than claiming a "
              "verification earned by a trivial annulus");
}

void test_longer_cycle_is_refused() {
    // Same principle at length three: X -> Y -> Z -> X.
    NameTable names;
    names.addLiterature("X", 1, 1);
    auto bounds = propagate({cobordism("X", 1, "Y", 1, 0),
                             cobordism("Y", 1, "Z", 1, 0),
                             cobordism("Z", 1, "X", 1, 0)},
                            names);
    EXPECT_EQ(bounds["X"].haveUpper(), false,
              "no cycle length launders a literature value back onto its "
              "own name");
}

void test_support_set_records_what_a_bound_rests_on() {
    NameTable names;
    names.addLiterature("K", 0, 9);
    names.addLiterature("J", 2, 2); // asserted, never witnessed
    auto bounds = propagate({cobordism("K", 1, "J", 1, 1)}, names);

    EXPECT_EQ(bounds["K"].hi, 3, "g_4(K) <= 2 + 1 + 0");
    EXPECT_EQ(bounds["K"].support.size(), static_cast<size_t>(1),
              "the bound rests on exactly one literature value");
    EXPECT_EQ(bounds["K"].support[0], std::string("J"), "namely J's");
    EXPECT_EQ(bounds["K"].basis == Basis::literatureAssisted, true,
              "and is therefore assisted, not constructive");

    // Witness J for real and the support empties out.
    auto bounds2 = propagate(
        {cobordism("K", 1, "J", 1, 1), direct("J", 1, 2)}, names);
    EXPECT_EQ(bounds2["K"].support.empty(), true,
              "with J actually witnessed, nothing is taken on faith");
    EXPECT_EQ(bounds2["K"].basis == Basis::constructive, true, "");
}

void test_unregistered_multicomponent_far_side_gets_no_reverse_bound() {
    // Regression for four live contradictions that killed three phases of a
    // sweep. candidates() falls back to `{name}` for any name absent from
    // the literature tables, so a far side known only by its COMPLEMENT --
    // here a Christy census name for a 2-component link -- arrives looking
    // like an unambiguous singleton. It is the opposite: a complement says
    // nothing about how its components are oriented.
    //
    // What actually happened: 6_1 -g0-> L204001 (2 curves) pushed a bound
    // onto L204001, which then bounded L6a3{0} at 1 against a literature
    // value of 2.
    NameTable names;
    names.addLiterature("6_1", 0, 0);
    names.addLiterature("L6a3{0}", 2, 2);
    names.addLiterature("L6a3{1}", 0, 0);
    // L204001 deliberately unregistered -- that is the whole point.

    std::vector<Witness> witnesses = {
        cobordism("6_1", 1, "L204001", 2, 0),
        cobordism("L6a3{0}", 2, "L204001", 2, 0),
    };
    auto bounds = propagate(witnesses, names);

    EXPECT_EQ(bounds["L204001"].haveUpper(), false,
              "an unregistered multi-component far side is maximally "
              "orientation-ambiguous, so no bound may be pushed onto it");
    EXPECT_EQ(bounds["L6a3{0}"].haveUpper(), false,
              "and so nothing propagates back out of it");

    Verdict v = judge("L6a3{0}", bounds["L6a3{0}"], names);
    EXPECT_EQ(v.status == Status::contradiction, false,
              "no contradiction, where the unguarded version produced one");
}

void test_unregistered_SINGLE_component_far_side_still_chains() {
    // The guard keys on component count, not on being unregistered: a
    // single curve has no orientation freedom, so a bare isoSig knot
    // complement must still chain. That path carries most of the graph's
    // connectivity and would be expensive to lose.
    NameTable names;
    names.addLiterature("Y", 0, 0);
    names.addLiterature("X", 0, 9);
    std::vector<Witness> witnesses = {
        direct("Y", 1, 0),
        cobordism("Y", 1, "gLLMQacdefeffhhnkxk", 1, 0),
        cobordism("X", 1, "gLLMQacdefeffhhnkxk", 1, 1),
    };
    auto bounds = propagate(witnesses, names);
    EXPECT_EQ(bounds["gLLMQacdefeffhhnkxk"].hi, 0,
              "a one-component unnamed node still receives a bound");
    EXPECT_EQ(bounds["X"].hi, 1, "and still passes it on");
}

void test_normalize_identified_name() {
    // identify() decorates a translated census hit; the --input tables do
    // not. Left undecorated, "4_1 (m004 : #1)" would be a DIFFERENT graph
    // node from the "4_1" row for the very same knot, and nothing would
    // ever chain through it.
    EXPECT_EQ(normalizeIdentifiedName("4_1 (m004 : #1)"), std::string("4_1"),
              "the census annotation is stripped for node identity");
    EXPECT_EQ(normalizeIdentifiedName("L6a1 (s780 : #6)"),
              std::string("L6a1"), "same for links");
    EXPECT_EQ(normalizeIdentifiedName("Unknot"), std::string("Unknot"),
              "an undecorated name is untouched");
    EXPECT_EQ(normalizeIdentifiedName("3-component unlink"),
              std::string("3-component unlink"), "unlinks are untouched");
    EXPECT_EQ(normalizeIdentifiedName("cPcbbbadu"), std::string("cPcbbbadu"),
              "a bare isoSig is untouched");
}

// ─────────────────────────────────────────────────────────────────────────
// judge()
// ─────────────────────────────────────────────────────────────────────────

void test_judge_distinguishes_assisted_verification() {
    // Reaching the literature's lower bound only by trusting ANOTHER
    // name's literature value is a real deduction but not an independent
    // verification, and the paper's tables depend on the difference.
    NameTable names;
    names.addLiterature("K", 1, 1);
    names.addLiterature("J", 0, 0); // J's genus is asserted, never witnessed
    auto bounds = propagate({cobordism("K", 1, "J", 1, 1)}, names);
    Verdict v = judge("K", bounds["K"], names);

    EXPECT_EQ(v.status == Status::verifiedAssisted, true,
              "K's bound rests on J's literature value, so it is reported "
              "as assisted rather than as a verification");

    // Now witness J for real; K's own bound becomes independent.
    auto bounds2 = propagate(
        {cobordism("K", 1, "J", 1, 1), direct("J", 1, 0)}, names);
    Verdict v2 = judge("K", bounds2["K"], names);
    EXPECT_EQ(v2.status == Status::verified, true,
              "with J actually witnessed, the whole chain is constructive");
}

void test_judge_verified() {
    NameTable names;
    names.addLiterature("K", 1, 1);
    auto bounds = propagate({direct("K", 1, 1)}, names);
    Verdict v = judge("K", bounds["K"], names);

    EXPECT_EQ(v.status == Status::verified, true,
              "we built a surface meeting the literature's lower bound, so "
              "the genus is exactly that");
    EXPECT_EQ(v.value, 1, "");
}

void test_judge_improved() {
    NameTable names;
    names.addLiterature("K", 0, 3); // literature range, upper bound 3
    auto bounds = propagate({direct("K", 1, 2)}, names);
    Verdict v = judge("K", bounds["K"], names);

    EXPECT_EQ(v.status == Status::improved, true,
              "genus 2 beats the literature's upper bound of 3 without "
              "reaching its lower bound of 0 -- a genuinely new bound");
    EXPECT_EQ(v.value, 2, "");
}

void test_judge_contradiction() {
    NameTable names;
    names.addLiterature("K", 2, 2);
    auto bounds = propagate({direct("K", 1, 1)}, names);
    Verdict v = judge("K", bounds["K"], names);

    EXPECT_EQ(v.status == Status::contradiction, true,
              "a surface below the literature's proven lower bound is a "
              "mathematical impossibility -- it means a bug in the search, "
              "not a new result");
}

void test_judge_unresolved() {
    NameTable names;
    names.addLiterature("K", 1, 1);
    auto bounds = propagate({}, names);
    Verdict v = judge("K", bounds["K"], names);
    EXPECT_EQ(v.status == Status::unresolved, true, "nothing derived");
}

// ─────────────────────────────────────────────────────────────────────────
// Bookkeeping
// ─────────────────────────────────────────────────────────────────────────

void test_have_witness_dedup() {
    std::vector<Witness> witnesses = {cobordism("K", 1, "J", 1, 2)};
    EXPECT_EQ(haveWitness(witnesses, cobordism("K", 1, "J", 1, 2)), true,
              "an identical witness is a duplicate -- this is what keeps a "
              "harvest run from capturing thousands of pair signatures for "
              "the same fact");
    EXPECT_EQ(haveWitness(witnesses, cobordism("K", 1, "J", 1, 3)), false,
              "a different genus is a different fact");
    EXPECT_EQ(haveWitness(witnesses, cobordism("K", 1, "L", 1, 2)), false,
              "a different far side is a different fact");
}

void test_build_depends_on_chain() {
    NameTable names;
    names.addLiterature("X", 0, 9);
    names.addLiterature("Y", 0, 9);
    std::vector<Witness> witnesses = {
        cobordism("Y", 1, "Unknot", 1, 0),
        cobordism("X", 1, "Y", 1, 0),
    };
    auto bounds = propagate(witnesses, names);
    EXPECT_EQ(buildDependsOn("Y", bounds), std::string("Y;Unknot"),
              "the chain bottoms out at the Unknot");
}

void test_build_depends_on_is_cycle_safe() {
    std::unordered_map<std::string, Bounds> bounds;
    bounds["A"].viaName = "B";
    bounds["A"].hi = 1;
    bounds["B"].viaName = "A";
    bounds["B"].hi = 1;
    EXPECT_EQ(buildDependsOn("A", bounds), std::string("A;B"),
              "a cycle terminates instead of looping forever");
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
    EXPECT_EQ(split.otherSides[0].components, 1,
              "a single curve is a one-component far side");
}

void test_split_boundary_unlink_components() {
    std::vector<BoundaryComponentNames> components = {
        BoundaryComponentNames{0, {"3_1"}, std::nullopt},
        BoundaryComponentNames{
            1, {"a", "b"}, std::optional<std::string>("2-component unlink")},
    };
    BoundarySplit split = splitBoundary(components, 0, "3_1");

    EXPECT_EQ(split.otherSides.size(), static_cast<size_t>(1), "one other side");
    EXPECT_EQ(split.otherSides[0].name, std::string("2-component unlink"),
              "named via linkName");
    EXPECT_EQ(split.otherSides[0].components, 2,
              "the component count is OBSERVED from the curve count, not "
              "inferred from the name -- this is the n_1 the cobordism "
              "inequality's + n_1 - 1 term refers to");
}

void test_split_boundary_linked_multicomponent_is_kept() {
    std::vector<BoundaryComponentNames> components = {
        BoundaryComponentNames{0, {"3_1"}, std::nullopt},
        BoundaryComponentNames{1, {"a", "b"},
                               std::optional<std::string>("L6a3")},
    };
    BoundarySplit split = splitBoundary(components, 0, "3_1");

    EXPECT_EQ(split.otherSides[0].name, std::string("L6a3"), "named");
    EXPECT_EQ(split.otherSides[0].components, 2,
              "a genuinely linked multi-component far side is KEPT, not "
              "refused: its orientation ambiguity is represented downstream "
              "as a candidate set propagate() bounds over, where the old "
              "design discarded the cobordism outright");
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
    EXPECT_EQ(split.otherSides[0].components, 2,
              "L206001 is a genuine multi-component linkName, kept with "
              "its observed curve count");
    EXPECT_EQ(split.otherSides[1].name, std::string("Unknot"), "");
    EXPECT_EQ(split.otherSides[1].components, 1, "a single curve");
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
    run("components_from_name", test_components_from_name);
    run("base_name", test_base_name);
    run("candidates_expand_orientation_variants",
        test_candidates_expand_orientation_variants);
    run("candidates_filtered_by_observed_component_count",
        test_candidates_filtered_by_observed_component_count);
    run("direct_witness_gives_upper_bound",
        test_direct_witness_gives_upper_bound);
    run("knot_cobordism_reduces_to_the_classic_rule",
        test_knot_cobordism_reduces_to_the_classic_rule);
    run("component_correction_on_the_upper_bound",
        test_component_correction_on_the_upper_bound);
    run("component_correction_on_the_lower_bound",
        test_component_correction_on_the_lower_bound);
    run("unlink_far_side_carries_no_component_penalty",
        test_unlink_far_side_carries_no_component_penalty);
    run("non_unlink_far_side_keeps_the_penalty",
        test_non_unlink_far_side_keeps_the_penalty);
    run("unlink_axiom_is_constructive", test_unlink_axiom_is_constructive);
    run("slice_composite_axiom_is_constructive",
        test_slice_composite_axiom_is_constructive);
    run("slice_composite_allowlist_is_not_a_pattern",
        test_slice_composite_allowlist_is_not_a_pattern);
    run("candidate_set_takes_the_worst_case_for_an_upper_bound",
        test_candidate_set_takes_the_worst_case_for_an_upper_bound);
    run("candidate_set_takes_the_best_case_for_a_lower_bound",
        test_candidate_set_takes_the_best_case_for_a_lower_bound);
    run("unambiguous_base_costs_nothing", test_unambiguous_base_costs_nothing);
    run("one_unbounded_candidate_blocks_the_upper_bound",
        test_one_unbounded_candidate_blocks_the_upper_bound);
    run("chains_through_an_unnamed_isosig_node",
        test_chains_through_an_unnamed_isosig_node);
    run("tubed_witness_genus_is_taken_at_face_value",
        test_tubed_witness_genus_is_taken_at_face_value);
    run("propagation_terminates_on_a_cycle",
        test_propagation_terminates_on_a_cycle);
    run("self_cobordism_cannot_confirm_the_literature",
        test_self_cobordism_cannot_confirm_the_literature);
    run("self_cobordism_still_allows_a_real_bound_from_elsewhere",
        test_self_cobordism_still_allows_a_real_bound_from_elsewhere);
    run("ambiguous_far_side_gets_no_reverse_bound",
        test_ambiguous_far_side_gets_no_reverse_bound);
    run("unambiguous_far_side_does_get_a_reverse_bound",
        test_unambiguous_far_side_does_get_a_reverse_bound);
    run("two_step_cycle_through_an_alias_is_refused",
        test_two_step_cycle_through_an_alias_is_refused);
    run("longer_cycle_is_refused", test_longer_cycle_is_refused);
    run("support_set_records_what_a_bound_rests_on",
        test_support_set_records_what_a_bound_rests_on);
    run("unregistered_multicomponent_far_side_gets_no_reverse_bound",
        test_unregistered_multicomponent_far_side_gets_no_reverse_bound);
    run("unregistered_SINGLE_component_far_side_still_chains",
        test_unregistered_SINGLE_component_far_side_still_chains);
    run("normalize_identified_name", test_normalize_identified_name);
    run("judge_verified", test_judge_verified);
    run("judge_distinguishes_assisted_verification",
        test_judge_distinguishes_assisted_verification);
    run("judge_improved", test_judge_improved);
    run("judge_contradiction", test_judge_contradiction);
    run("judge_unresolved", test_judge_unresolved);
    run("have_witness_dedup", test_have_witness_dedup);
    run("build_depends_on_chain", test_build_depends_on_chain);
    run("build_depends_on_is_cycle_safe", test_build_depends_on_is_cycle_safe);
    run("split_boundary_single_curve_is_safe",
        test_split_boundary_single_curve_is_safe);
    run("split_boundary_unlink_components",
        test_split_boundary_unlink_components);
    run("split_boundary_linked_multicomponent_is_kept",
        test_split_boundary_linked_multicomponent_is_kept);
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
