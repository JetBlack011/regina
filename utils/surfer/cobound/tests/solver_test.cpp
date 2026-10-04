// solver_test.cpp
//
// Tests for ../solver/solver.h/.cpp: the interval solver `cobound solve`
// runs over the database's cobordisms (componentsFromName(), NameTable,
// propagate(), judge()) and the boundary classification (splitBoundary())
// that feeds it. All pure logic on plain data types -- no triangulations, no
// search, no census -- so these tests run instantly and exercise the core of
// the solver's genus deductions directly.
//
// The cases below deliberately pin down the two things that are easy to get
// silently wrong: the component-count terms in the cobordism inequality
// (which vanish for knots, so a knots-only test suite would never notice
// them missing), and WHICH outgoing links may carry a bound at all. An
// outgoing link is named from its complement, and a complement determines a knot
// (Gordon-Luecke) but not a link (Rolfsen twisting: infinitely many links
// per exterior). So a multi-component outgoing link bounds nothing unless it
// is a structurally proven unlink -- see outgoingBearsBound() and solver.h
// \ref cg_outgoing. For months the solver took a max/min over a link's
// ORIENTATION variants as if that were the whole candidate set; 178 of 333
// "verified" rows rested on it. The tests here are what stop that coming
// back.

#include <iostream>
#include <map>
#include <string>
#include <unistd.h>

#include <triangulation/dim3.h>

#include "cobound/cobordisms/cobordism.h"
#include "cobound/search/incoming.h"
#include "cobound/search/preconditions.h"
#include "cobound/solver/literature.h"
#include "cobound/solver/solver.h"
#include "linknaming/names.h"

using namespace cobordisms;
using namespace search;
using namespace solver;
using namespace linknaming;

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

// Builds a cobordism between `subject` (with `nA` components) and
// `other` (with `nB`), at genus `g`. `candidates` defaults to just `other`.
Cobordism cobordism(const std::string &subject, int nA, const std::string &other,
                  int nB, int g,
                  std::vector<std::string> candidates = {}) {
    Cobordism w;
    w.kind = CobordismKind::cobordism;
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

Cobordism direct(const std::string &subject, int nA, int g, bool tubed = false) {
    Cobordism w;
    w.kind = CobordismKind::direct;
    w.subject = subject;
    w.subjectComponents = nA;
    w.genus = g;
    w.tubed = tubed;
    return w;
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

void test_direct_cobordism_gives_upper_bound() {
    NameTable names;
    names.addLiterature("K", 1, 1);
    auto bounds = propagate({direct("K", 1, 1)}, names);

    EXPECT_EQ(bounds["K"].haveUpper(), true, "a direct witness bounds K");
    EXPECT_EQ(bounds["K"].hi, 1, "at exactly the surface's own genus");
    EXPECT_EQ(bounds["K"].basis == Basis::constructive, true,
              "a direct witness depends on nothing else, so the bound is "
              "constructive -- nothing taken on faith");
}

void test_external_proofs() {
    // certified_bounds: a certified proof is a base case like a direct
    // cobordism, constructive only when it names no literature value, and it
    // grounds chains through the cobordisms.
    NameTable names;
    names.addLiterature("K", 1, 1);
    names.addLiterature("J", 0, 1);
    names.addLiterature("M", 2, 2);
    names.addLiterature("S", 0, 0);
    {
        auto bounds = propagate({cobordism("J", 1, "K", 1, 0)}, names,
                                {{"K", 1, {}, "cascade:run/K"}});
        EXPECT_EQ(bounds["K"].hi, 1, "a proof bounds its target");
        EXPECT_EQ(bounds["K"].basis == Basis::constructive, true,
                  "with no literature leaves it is constructive");
        EXPECT_EQ(bounds["K"].viaName, std::string("cascade:run/K"),
                  "the report names the proof");
        EXPECT_EQ(bounds["J"].hi, 1, "and grounds a chain: J <= K + 0");
        EXPECT_EQ(bounds["J"].basis == Basis::constructive, true,
                  "still constructive along the chain");
    }
    {
        auto bounds = propagate({}, names, {{"K", 1, {"S", "M"}, "cascade:run/K"}});
        EXPECT_EQ(bounds["K"].basis == Basis::literatureAssisted, true,
                  "literature leaves make it assisted");
        EXPECT_EQ(bounds["K"].support.size(), static_cast<size_t>(2),
                  "the leaves are its support");
    }
    {
        auto bounds = propagate({}, names, {{"K", 1, {"K"}, "cascade:run/K"}});
        EXPECT_EQ(bounds["K"].haveUpper(), false,
                  "a proof resting on the target's own literature value is circular");
    }
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
    // between a knot and a genuinely linked 3-component link does NOT make
    // the knot slice: the link must be capped with a CONNECTED surface,
    // which costs + n - 1 = 2 in genus.
    //
    // The link is the SUBJECT here, not the outgoing link. A linked outgoing
    // link is named from its complement and so bounds nothing at all (see
    // test_linked_outgoing_bounds_nothing); the only linked endpoint the
    // solver may reason from is one known by construction, i.e. a searched link, and
    // the bound then flows in the reverse direction onto the knot.
    NameTable names;
    names.addLiterature("L", 0, 0);
    names.addLiterature("K", 0, 9);
    auto bounds = propagate(
        {direct("L", 3, 0), cobordism("L", 3, "K", 1, 0)}, names);

    EXPECT_EQ(bounds["K"].hi, 2,
              "g_4(K) <= 0 + 0 + (3 - 1) = 2. The uncorrected rule would "
              "give 0 here, i.e. would 'prove' K slice on evidence that "
              "says no such thing");
}

void test_component_correction_on_the_lower_bound() {
    // The mirror case: deducing a lower bound for a MULTI-component
    // subject costs - n_0 + 1.
    NameTable names;
    names.addLiterature("L", 0, 9); // 2 components, per the cobordism below
    names.addLiterature("K", 3, 3);
    auto bounds = propagate({cobordism("L", 2, "K", 1, 0)}, names);

    EXPECT_EQ(bounds["L"].lo, 2,
              "g_4(L) >= g_4(K) - g - n_0 + 1 = 3 - 0 - 2 + 1 = 2. The "
              "uncorrected rule would claim 3, a bound that is simply not "
              "justified");
}

void test_unlink_outgoing_carries_no_component_penalty() {
    // The `+ n_b - 1` term assumes the outgoing link is capped with its minimal
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

void test_derived_lower_above_derived_upper_is_a_contradiction() {
    // K bounds a disc (a direct genus-0 cobordism: hi 0, constructive), and a
    // genus-1 cobordism runs from K to a knot outgoing link whose every candidate
    // has literature g4 = 3, which transports lo(K) >= 3 - 1 - 1 + 1 = 2.
    // Two candidates keep the reverse rule out, so nothing else can notice.
    // Both bounds sit inside K's literature interval [0, 2], but together
    // they say 2 <= g4(K) <= 0: a cobordism, a name or a table is wrong, and
    // the status is the computation's only consistency check (paper,
    // def:status). It must not come out "verified".
    NameTable names;
    names.addLiterature("K", 0, 2);
    names.addLiterature("J", 3, 3);
    names.addLiterature("J2", 3, 3);
    auto bounds = propagate({direct("K", 1, 0),
                             cobordism("K", 1, "J", 1, 1, {"J", "J2"})},
                            names);

    EXPECT_EQ(bounds["K"].hi, 0, "the disc gives hi(K) = 0");
    EXPECT_EQ(bounds["K"].lo, 2, "the cobordism gives lo(K) = 2");
    Verdict v = judge("K", bounds["K"], names);
    EXPECT_EQ(v.status == Status::contradiction, true,
              "derived lo 2 above derived hi 0 is a contradiction, not a "
              "verification");
}

void test_knot_outgoing_named_as_a_link_bounds_nothing() {
    // One curve was observed on the outgoing link, but the name it was given
    // is a link's base, so no registered variant has the observed count: the
    // name and the geometry disagree. candidates() then hands back
    // every link variant, and a one-curve outgoing link bears a bound
    // (outgoingBearsBound), so the subject would be bounded by a LINK's g4 as
    // if it were this knot's. The cobordism is built as the driver builds it
    // (cobordisms/pending.cpp: otherCandidates = candidates(name, observed)).
    NameTable names;
    names.addLiterature("K", 0, 3);
    names.addLiterature("L2a1{0}", 0, 0);
    names.addLiterature("L2a1{1}", 0, 0);
    Cobordism w = cobordism("K", 1, "L2a1", 1, 0,
                          names.candidates("L2a1", 1));
    auto bounds = propagate({w}, names);

    EXPECT_EQ(bounds.contains("K") && bounds.at("K").haveUpper(), false,
              "a knot far side whose name only matches links carries no "
              "bound: no candidate is a knot");
}

void test_linked_outgoing_bounds_nothing() {
    // The unlink exemption is the ONLY way a multi-component outgoing link gets
    // to carry a bound. Every other multi-component name is a statement
    // about a complement, and a link complement belongs to infinitely many
    // links (Rolfsen twisting), so "L4a1" here does not mean L4a1 -- it
    // means "some link whose exterior is L4a1's", which has no slice genus.
    // Both variants are registered and agree, and it still bounds nothing:
    // agreement across the ORIENTATION variants is not agreement across the
    // actual candidate set, which cannot be enumerated.
    NameTable names;
    names.addLiterature("K", 0, 9);
    names.addLiterature("L4a1{0}", 0, 0);
    names.addLiterature("L4a1{1}", 0, 0);
    auto bounds = propagate(
        {cobordism("K", 1, "L4a1", 2, 0, {"L4a1{0}", "L4a1{1}"})}, names);
    EXPECT_EQ(bounds["K"].haveUpper(), false,
              "a complement-named 2-component far side carries no upper "
              "bound, however its oriented variants are tabulated");
    EXPECT_EQ(bounds["K"].haveLower(), false, "nor a lower one");
    EXPECT_EQ(bounds["L4a1{0}"].haveUpper(), false,
              "and receives none in the reverse direction either");
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
// propagate(): orientation-ambiguous outgoing links
// ─────────────────────────────────────────────────────────────────────────

void test_orientation_variants_are_not_a_candidate_set() {
    // L4a1{0} has slice genus 0 and L4a1{1} has 1, and they share one
    // complement. The solver used to bound K by max(0, 1) + 0 + (2 - 1) = 2,
    // "sound whichever variant it actually was". It is not: the outgoing link
    // need not be EITHER variant. Our own peripheral tests found the census
    // name m129 standing for 18 different links, g_4 from 0 to 2. So the
    // max over {L4a1{0}, L4a1{1}} is a max over the wrong set.
    NameTable names;
    names.addLiterature("K", 0, 9);
    names.addLiterature("L4a1{0}", 0, 0);
    names.addLiterature("L4a1{1}", 1, 1);

    auto bounds = propagate(
        {cobordism("K", 1, "L4a1", 2, 0, {"L4a1{0}", "L4a1{1}"})}, names);

    EXPECT_EQ(bounds["K"].haveUpper(), false,
              "no upper bound from a max over orientation variants");
    EXPECT_EQ(bounds["K"].haveLower(), false,
              "no lower bound from a min over them");
}

void test_outgoing_bears_bound() {
    // The predicate the solver gates on, stated directly.
    EXPECT_EQ(outgoingBearsBound(cobordism("K", 1, "4_1", 1, 0)), true,
              "a knot far side bounds (Gordon-Luecke)");
    EXPECT_EQ(outgoingBearsBound(
                  cobordism("K", 1, "gLLMQacdefeffhhnkxk", 1, 0)),
              true, "so does a one-component bare isoSig");
    EXPECT_EQ(outgoingBearsBound(cobordism("K", 1, "Unknot", 1, 0)), true,
              "and the unknot");
    EXPECT_EQ(outgoingBearsBound(
                  cobordism("K", 1, "2-component unlink", 2, 0)),
              true, "an unlink is a structural proof of the link itself");
    EXPECT_EQ(outgoingBearsBound(cobordism("K", 1, "L4a1{0}", 2, 0)), false,
              "a Thistlethwaite name with 2 curves does not");
    EXPECT_EQ(outgoingBearsBound(cobordism("K", 1, "L204001", 2, 0)), false,
              "nor a Christy census name");
    EXPECT_EQ(outgoingBearsBound(cobordism("K", 1, "m129", 2, 0)), false,
              "nor a SnapPea census name");
    EXPECT_EQ(outgoingBearsBound(
                  cobordism("K", 1, "gLLMQacdefeffhhnkxk", 3, 0)),
              false, "nor a bare isoSig with 3 curves");
    // The gate is on the OBSERVED count, so a name that looks like a knot
    // cannot smuggle a two-curve outgoing link through.
    EXPECT_EQ(outgoingBearsBound(cobordism("K", 1, "4_1", 2, 0)), false,
              "a knot-shaped name on a 2-curve far side is refused");
}

// ─────────────────────────────────────────────────────────────────────────
// propagate(): chaining
// ─────────────────────────────────────────────────────────────────────────

void test_chains_through_an_unnamed_isosig_link() {
    // census::nameComplement() falls back to a bare isoSig when the census
    // misses. Such a name has no bounds of its own, but it still CONNECTS:
    // two subjects that each cobound with it are thereby related to each other.
    // Keeping these names in the graph is free information.
    NameTable names;
    names.addLiterature("X", 0, 9);
    names.addLiterature("Y", 0, 0);

    std::vector<Cobordism> cobordisms = {
        direct("Y", 1, 0),
        cobordism("Y", 1, "gLLMQacdefeffhhnkxk", 1, 0),
        cobordism("X", 1, "gLLMQacdefeffhhnkxk", 1, 1),
    };
    auto bounds = propagate(cobordisms, names);

    EXPECT_EQ(bounds["gLLMQacdefeffhhnkxk"].hi, 0,
              "the unnamed node picks up Y's own bound through their "
              "shared cobordism");
    EXPECT_EQ(bounds["X"].hi, 1,
              "and passes it on to X, which never met Y directly");
    EXPECT_EQ(bounds["X"].basis == Basis::constructive, true,
              "every step rested on a surface we actually found");
}

void test_tubed_cobordism_genus_is_taken_at_face_value() {
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
    std::vector<Cobordism> cobordisms = {
        direct("A", 1, 1),
        cobordism("A", 1, "B", 1, 0),
        cobordism("B", 1, "C", 1, 0),
        cobordism("C", 1, "A", 1, 0),
    };
    auto bounds = propagate(cobordisms, names);
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
    // The self-loop is skipped, not poisonous: a genuine cobordism on the
    // same name still lands.
    NameTable names;
    names.addLiterature("K", 0, 5);
    auto bounds = propagate(
        {cobordism("K", 1, "K", 1, 0), direct("K", 1, 2)}, names);
    EXPECT_EQ(bounds["K"].hi, 2, "the direct witness still bounds K");
}

void test_ambiguous_outgoing_gets_no_reverse_bound() {
    // Regression for a real contradiction caught by a live run. A genus-0
    // cobordism from L4a1{0} to an outgoing link named only as "L7n1" was
    // once pushed onto every candidate, giving L7n1{0} (true genus 2) a
    // bound of 1. The first fix stopped at "we cannot tell which variant";
    // the real reason is stronger -- the outgoing link need not be any L7n1 at
    // all, since L7n1 shares its exterior with L5a1, L8n2, L9n3, L10n9.
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

void test_unambiguous_outgoing_does_get_a_reverse_bound() {
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
    // up as a cobordism between two DIFFERENT-looking names. 8_14's own
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
    names.addLiterature("J", 2, 2); // asserted, never constructed
    auto bounds = propagate({cobordism("K", 1, "J", 1, 1)}, names);

    EXPECT_EQ(bounds["K"].hi, 3, "g_4(K) <= 2 + 1 + 0");
    EXPECT_EQ(bounds["K"].support.size(), static_cast<size_t>(1),
              "the bound rests on exactly one literature value");
    EXPECT_EQ(bounds["K"].support[0], std::string("J"), "namely J's");
    EXPECT_EQ(bounds["K"].basis == Basis::literatureAssisted, true,
              "and is therefore assisted, not constructive");

    // A cobordism for J, for real, and the support empties out.
    auto bounds2 = propagate(
        {cobordism("K", 1, "J", 1, 1), direct("J", 1, 2)}, names);
    EXPECT_EQ(bounds2["K"].support.empty(), true,
              "with J actually witnessed, nothing is taken on faith");
    EXPECT_EQ(bounds2["K"].basis == Basis::constructive, true, "");
}

void test_unregistered_multicomponent_outgoing_gets_no_reverse_bound() {
    // Regression for four live contradictions that killed three phases of a
    // sweep. candidates() falls back to `{name}` for any name absent from
    // the literature tables, so an outgoing link known only by its COMPLEMENT --
    // here a Christy census name for a 2-component link -- arrives looking
    // like an unambiguous singleton. It is the opposite: a complement says
    // nothing about WHICH 2-component link this is, let alone how it is
    // oriented -- see outgoingBearsBound().
    //
    // What actually happened: 6_1 -g0-> L204001 (2 curves) pushed a bound
    // onto L204001, which then bounded L6a3{0} at 1 against a literature
    // value of 2.
    NameTable names;
    names.addLiterature("6_1", 0, 0);
    names.addLiterature("L6a3{0}", 2, 2);
    names.addLiterature("L6a3{1}", 0, 0);
    // L204001 deliberately unregistered -- that is the whole point.

    std::vector<Cobordism> cobordisms = {
        cobordism("6_1", 1, "L204001", 2, 0),
        cobordism("L6a3{0}", 2, "L204001", 2, 0),
    };
    auto bounds = propagate(cobordisms, names);

    EXPECT_EQ(bounds["L204001"].haveUpper(), false,
              "an unregistered multi-component far side is a complement, "
              "not a link, so no bound may be pushed onto it");
    EXPECT_EQ(bounds["L6a3{0}"].haveUpper(), false,
              "and so nothing propagates back out of it");

    Verdict v = judge("L6a3{0}", bounds["L6a3{0}"], names);
    EXPECT_EQ(v.status == Status::contradiction, false,
              "no contradiction, where the unguarded version produced one");
}

void test_unregistered_SINGLE_component_outgoing_still_chains() {
    // The guard keys on component count, not on being unregistered: a
    // single curve has no orientation freedom, so a bare isoSig knot
    // complement must still chain. That path carries most of the graph's
    // connectivity and would be expensive to lose.
    NameTable names;
    names.addLiterature("Y", 0, 0);
    names.addLiterature("X", 0, 9);
    std::vector<Cobordism> cobordisms = {
        direct("Y", 1, 0),
        cobordism("Y", 1, "gLLMQacdefeffhhnkxk", 1, 0),
        cobordism("X", 1, "gLLMQacdefeffhhnkxk", 1, 1),
    };
    auto bounds = propagate(cobordisms, names);
    EXPECT_EQ(bounds["gLLMQacdefeffhhnkxk"].hi, 0,
              "a one-component unnamed node still receives a bound");
    EXPECT_EQ(bounds["X"].hi, 1, "and still passes it on");
}

void test_two_component_isosig_link_does_not_chain() {
    // The mirror of the test above. An unnamed TWO-curve name is a shared
    // exterior, and two cobordisms landing on the same exterior may have
    // found two different links -- so the name must not connect them. This
    // is the m129 case from results/far_side_identification_census.csv,
    // where one name was 18 links, and it is the regression the golden-
    // truffle plan asks for: a relapse into complement-only reasoning.
    NameTable names;
    names.addLiterature("Y", 0, 0);
    names.addLiterature("X", 0, 9);
    std::vector<Cobordism> cobordisms = {
        direct("Y", 1, 0),
        cobordism("Y", 1, "m129", 2, 0),
        cobordism("X", 1, "m129", 2, 1),
    };
    auto bounds = propagate(cobordisms, names);
    EXPECT_EQ(bounds["m129"].haveUpper(), false,
              "a two-curve exterior receives no bound from Y");
    EXPECT_EQ(bounds["X"].haveUpper(), false,
              "and so cannot pass one on to X: the two surfaces may have "
              "found different links with that exterior");
}

void test_unlink_outgoing_still_chains_and_is_penalty_free() {
    // The one multi-component outgoing link that DOES bear a bound, checked
    // alongside the ones that do not so the gate is seen to be selective
    // rather than a blanket refusal of links.
    NameTable names;
    names.addLiterature("K", 0, 9);
    auto bounds = propagate(
        {cobordism("K", 1, "3-component unlink", 3, 1)}, names);
    EXPECT_EQ(bounds["K"].hi, 1,
              "g_4(K) <= 0 + 1 + 0: the unlink is an axiom and carries no "
              "component penalty");
    EXPECT_EQ(bounds["K"].basis == Basis::constructive, true,
              "and rests on nothing from the literature");
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
    names.addLiterature("J", 0, 0); // J's genus is asserted, never constructed
    auto bounds = propagate({cobordism("K", 1, "J", 1, 1)}, names);
    Verdict v = judge("K", bounds["K"], names);

    EXPECT_EQ(v.status == Status::verifiedAssisted, true,
              "K's bound rests on J's literature value, so it is reported "
              "as assisted rather than as a verification");

    // Now a cobordism for J, for real; K's own bound becomes independent.
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

void test_have_cobordism_dedup() {
    std::vector<Cobordism> cobordisms = {cobordism("K", 1, "J", 1, 2)};
    EXPECT_EQ(haveCobordism(cobordisms, cobordism("K", 1, "J", 1, 2)), true,
              "an identical witness is a duplicate -- this is what keeps a "
              "harvest run from capturing thousands of pair signatures for "
              "the same fact");
    EXPECT_EQ(haveCobordism(cobordisms, cobordism("K", 1, "J", 1, 3)), false,
              "a different genus is a different fact");
    EXPECT_EQ(haveCobordism(cobordisms, cobordism("K", 1, "L", 1, 2)), false,
              "a different far side is a different fact");
}

void test_build_depends_on_chain() {
    NameTable names;
    names.addLiterature("X", 0, 9);
    names.addLiterature("Y", 0, 9);
    std::vector<Cobordism> cobordisms = {
        cobordism("Y", 1, "Unknot", 1, 0),
        cobordism("X", 1, "Y", 1, 0),
    };
    auto bounds = propagate(cobordisms, names);
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
        BoundaryComponentNames{0, {"3_1"}, std::nullopt, {1, 2, 3}},
        BoundaryComponentNames{1, {"4_1"}, std::nullopt, {4, 5, 6}},
    };
    BoundarySplit split = splitBoundary(components, 0);

    EXPECT_EQ(split.searchCurveCount, static_cast<size_t>(1),
              "component 0 (== incomingBC) is the incoming side");
    EXPECT_EQ(split.otherSides.size(), static_cast<size_t>(1),
              "exactly one other side");
    EXPECT_EQ(split.otherSides[0].name, std::string("4_1"),
              "a single curve is named directly");
    EXPECT_EQ(split.otherSides[0].components, 1,
              "a single curve is a one-component far side");
}

void test_split_boundary_unlink_components() {
    std::vector<BoundaryComponentNames> components = {
        BoundaryComponentNames{0, {"3_1"}, std::nullopt, {1, 2, 3}},
        BoundaryComponentNames{
            1, {"?", "?"}, std::optional<std::string>("2-component unlink"),
            {4, 5, 6, 7, 8, 9}},
    };
    BoundarySplit split = splitBoundary(components, 0);

    EXPECT_EQ(split.otherSides.size(), static_cast<size_t>(1), "one other side");
    EXPECT_EQ(split.otherSides[0].name, std::string("2-component unlink"),
              "named via linkName, whatever the per-curve placeholders say");
    EXPECT_EQ(split.otherSides[0].components, 2,
              "the component count is OBSERVED from the curve count, not "
              "inferred from the name -- this is the n_1 the cobordism "
              "inequality's + n_1 - 1 term refers to");
}

void test_split_boundary_linked_multicomponent_is_kept() {
    std::vector<BoundaryComponentNames> components = {
        BoundaryComponentNames{0, {"3_1"}, std::nullopt, {1, 2, 3}},
        BoundaryComponentNames{1, {"?", "?"},
                               std::optional<std::string>("L6a3"),
                               {4, 5, 6, 7}},
    };
    BoundarySplit split = splitBoundary(components, 0);

    EXPECT_EQ(split.otherSides[0].name, std::string("L6a3"), "named");
    EXPECT_EQ(split.otherSides[0].components, 2,
              "a genuinely linked multi-component outgoing link is KEPT, not "
              "refused: it is recorded, and outgoingBearsBound() is what "
              "declines to bound anything by it");
}

void test_split_boundary_multiple_other_sides_not_collapsed() {
    std::vector<BoundaryComponentNames> components = {
        BoundaryComponentNames{0, {"3_1"}, std::nullopt, {1}},
        BoundaryComponentNames{1, {"4_1"}, std::nullopt, {2}},
        BoundaryComponentNames{2, {"5_1"}, std::nullopt, {3}},
    };
    BoundarySplit split = splitBoundary(components, 0);

    EXPECT_EQ(split.otherSides.size(), static_cast<size_t>(2),
              "both other sides are kept -- neither silently overwrites the "
              "other");
}

void test_split_boundary_seeded_ignores_names() {
    // The D1 regression (2026-09-26). Seeded, the incoming side is L by
    // construction, and its name must never be consulted: gdb on 8_8 caught
    // the searched link's own name as "8_8 (o9_37770 : #17)" and the incoming side as
    // "8_8 (o9_37770 : #6)" -- the same manifold, a different census entry
    // number -- and the old name comparison then discarded every surface of
    // the search (270 of the atlas's rows, with nothing logged).
    std::vector<BoundaryComponentNames> components = {
        BoundaryComponentNames{0, {"8_8 (o9_37770 : #6)"}, std::nullopt,
                               {7, 8, 9}},
        BoundaryComponentNames{1, {"Unknot"}, std::nullopt, {1, 2}},
    };
    BoundarySplit split = splitBoundary(components, 0);

    EXPECT_EQ(split.searchCurveCount, static_cast<size_t>(1),
              "the incoming side is component incomingBC, whatever it was "
              "named");
    EXPECT_EQ(split.otherSides.size(), static_cast<size_t>(1),
              "the far side is still classified");
}

void test_split_boundary_unnamed_side_flagged() {
    std::vector<BoundaryComponentNames> components = {
        BoundaryComponentNames{0, {"3_1"}, std::nullopt, {1}},
        BoundaryComponentNames{1, {"?", "?"}, std::nullopt, {2, 3}},
    };
    BoundarySplit split = splitBoundary(components, 0);

    EXPECT_EQ(split.unnamedSide, true,
              "a multi-curve far side with no link name is reported, never "
              "silently dropped (dropping it would turn a cobordism into a "
              "'direct' witness)");
    EXPECT_EQ(split.otherSides.empty(), true, "and not classified");
}

// ─────────────────────────────────────────────────────────────────────────────
// classifyIncomingOrientation()'s decision logic, isolated from
// buildIncomingOrientation()'s geometry by hand-constructing IncomingOrientation and
// OrientedCurve against a triangulation's own edges. Two vertex-disjoint
// edges stand in for two components of the incoming link.
// ─────────────────────────────────────────────────────────────────────────────
void test_classify_incoming_orientation() {
    regina::Triangulation<3> tri;
    tri.newTetrahedron();
    tri.newTetrahedron(); // unglued: its edges are disjoint from the first's
    // Regina's tetrahedron edge ordering is {(0,1),(0,2),(0,3),(1,2),(1,3),
    // (2,3)}: edges 0 and 5 of one tetrahedron are vertex-disjoint.
    regina::Edge<3> *e0 = tri.tetrahedron(0)->edge(0);
    regina::Edge<3> *e1 = tri.tetrahedron(0)->edge(5);
    regina::Edge<3> *e2 = tri.tetrahedron(1)->edge(0); // not the incoming link's

    IncomingOrientation incoming;
    incoming.tailOf[e0->index()] = e0->vertex(0)->index();
    incoming.tailOf[e1->index()] = e1->vertex(0)->index();

    const std::map<const regina::Edge<3> *, size_t> oneComponent = {
        {e0, 0}, {e1, 0}};
    const std::map<const regina::Edge<3> *, size_t> twoComponents = {
        {e0, 0}, {e1, 1}};
    using V = OrientationVerdict;

    std::vector<OrientedCurve> allMatch = {{{e0, false}}, {{e1, false}}};
    EXPECT_EQ(classifyIncomingOrientation(incoming, allMatch, oneComponent) == V::match,
              true, "both curves agree with the row's tag: match");

    std::vector<OrientedCurve> allFlipped = {{{e0, true}}, {{e1, true}}};
    EXPECT_EQ(classifyIncomingOrientation(incoming, allFlipped, oneComponent) ==
                  V::match,
              true,
              "both disagree, uniformly: the surface's other orientation");

    std::vector<OrientedCurve> mixed = {{{e0, false}}, {{e1, true}}};
    EXPECT_EQ(classifyIncomingOrientation(incoming, mixed, oneComponent) ==
                  V::mismatch,
              true,
              "ONE surface component inducing a mixed pattern witnesses a "
              "different oriented variant -- the L6a3{0}/L6a3{1} "
              "misattribution -- and is rejected");

    EXPECT_EQ(classifyIncomingOrientation(incoming, mixed, twoComponents) == V::match,
              true,
              "the D2 regression (2026-09-26): the same pattern split across "
              "TWO surface components is a match, since each component is "
              "oriented independently and they tube together either way. "
              "Treating it as a mismatch discarded every disconnected "
              "surface of whole rows");

    std::vector<OrientedCurve> foreign = {{{e2, false}}};
    const std::map<const regina::Edge<3> *, size_t> foreignComponent = {
        {e2, 0}};
    EXPECT_EQ(classifyIncomingOrientation(incoming, foreign, foreignComponent) ==
                  V::foreignEdge,
              true, "an edge that is not the row's own is reported as such");

    std::vector<OrientedCurve> incoherent = {{{e0, false}, {e1, true}}};
    EXPECT_EQ(classifyIncomingOrientation(incoming, incoherent, oneComponent) ==
                  V::incoherentCurve,
              true, "one curve whose edges disagree in direction");

    std::vector<OrientedCurve> single = {{{e0, false}}};
    EXPECT_EQ(classifyIncomingOrientation(incoming, single, {}) == V::incoherentCurve,
              true, "a curve with no known surface component is not guessed");

    // The judgement's flips (phase 3: the one walk behind both
    // classifyIncomingOrientation() and outgoing::incomingFlips()).
    using Flips = std::map<size_t, int>;
    const Flips keep = {{0, 1}}, reverse = {{0, -1}}, each = {{0, 1}, {1, -1}};
    EXPECT_EQ(judgeIncomingOrientation(incoming, allMatch, oneComponent).flips == keep, true,
              "matching curves keep their component (+1)");
    EXPECT_EQ(judgeIncomingOrientation(incoming, allFlipped, oneComponent).flips == reverse, true,
              "uniformly reversed curves reverse it (-1)");
    EXPECT_EQ(judgeIncomingOrientation(incoming, mixed, twoComponents).flips == each, true,
              "two components, each its own flip");
    const IncomingOrientationJudgement none = judgeIncomingOrientation(incoming, {}, oneComponent);
    EXPECT_EQ(none.verdict == V::mismatch && none.noCurves, true,
              "no curves: a mismatch (nothing witnesses the row)");
    EXPECT_EQ(none.consistentFlips() == std::optional<Flips>(Flips()), true,
              "...whose flips are the empty map, as incomingFlips() always gave");
    EXPECT_EQ(!judgeIncomingOrientation(incoming, mixed, oneComponent).consistentFlips(), true,
              "a mismatch has no flips");
    EXPECT_EQ(!judgeIncomingOrientation(incoming, foreign, foreignComponent).consistentFlips(), true,
              "nor has a foreign edge");
    EXPECT_EQ(judgeIncomingOrientation(incoming, mixed, twoComponents).consistentFlips() ==
                  std::optional<Flips>(each),
              true, "a match's flips are the judgement's");
}

void test_cobordism_identity() {
    Cobordism a = cobordism("K", 2, "L6a3", 2, 0);
    Cobordism b = a;
    EXPECT_EQ(cobordismIdentity(a) == cobordismIdentity(b), true,
              "identical witnesses share an identity");

    b.resolvedVertices = 1;
    EXPECT_EQ(cobordismIdentity(a) == cobordismIdentity(b), false,
              "the D5 fix: an embedded witness is never a duplicate of a "
              "resolved one");

    Cobordism c = a;
    c.genus = 1;
    EXPECT_EQ(cobordismIdentity(a) == cobordismIdentity(c), false,
              "a different genus is a different witness");

    std::vector<Cobordism> cobordisms = {a};
    EXPECT_EQ(haveCobordism(cobordisms, b), false, "haveCobordism() agrees");
    EXPECT_EQ(haveCobordism(cobordisms, c), false, "haveCobordism() agrees");
    cobordisms.push_back(b);
    EXPECT_EQ(haveCobordism(cobordisms, b), true, "haveCobordism() agrees");
}

} // namespace

// ─────────────────────────────────────────────────────────────────────────
// Per-cobordism proved outgoing links (outgoing_resolutions)
// ─────────────────────────────────────────────────────────────────────────

void test_proved_link_outgoing_bears_bound() {
    // The counterpart of test_linked_outgoing_bounds_nothing. There the
    // outgoing link was named from its complement, which does not determine a link.
    // Here applyOutgoingResolutions() has marked it proved PER COBORDISM -- an
    // isometry carrying this cobordism's own meridians -- and complement plus
    // meridians does determine the link, so its oriented variants ARE the
    // complete candidate set and the max over them is sound.
    NameTable names;
    names.addLiterature("K", 0, 9);
    names.addLiterature("L4a1{0}", 0, 0);
    names.addLiterature("L4a1{1}", 0, 0);
    Cobordism w = cobordism("K", 1, "L4a1", 2, 0, {"L4a1{0}", "L4a1{1}"});
    EXPECT_EQ(outgoingBearsBound(w), false, "unproved: the gate refuses it");
    w.outgoingProved = true;
    EXPECT_EQ(outgoingBearsBound(w), true, "proved: the gate accepts it");
    auto bounds = propagate({w}, names);
    EXPECT_EQ(bounds["K"].hi, 1,
              "g_4(K) <= max(0, 0) + 0 + (2 - 1) = 1, by hand");
}

// ─────────────────────────────────────────────────────────────────────────
// Split (disjoint-union) outgoing names: "A u B"
// ─────────────────────────────────────────────────────────────────────────

void test_split_outgoing_upper_bound() {
    // g_4(A u B) <= g_4(A) + g_4(B): tube the factors' minimal surfaces
    // together (constructive). With 3_1 at 1 and the Unknot axiom at 0, a
    // genus-0 cobordism K -> (3_1 u Unknot) gives
    //     g_4(K) <= (1 + 0) + 0 + (2 - 1) = 2.
    NameTable names;
    names.addLiterature("K", 0, 9);
    names.addLiterature("3_1", 1, 1);
    names.addLiterature("m3_1", 1, 1);
    Cobordism w = cobordism("K", 1, "3_1 u Unknot", 2, 0);
    w.outgoingProved = true;
    auto bounds = propagate({w}, names);
    EXPECT_EQ(bounds["K"].hi, 2, "g_4(K) <= 1 + 0 + 0 + 1 = 2, by hand");
    EXPECT_EQ(bounds["K"].basis == Basis::literatureAssisted, true,
              "it rests on 3_1's literature value, so it is assisted");

    Cobordism m = cobordism("K", 1, "3_1|m3_1 u Unknot", 2, 0);
    m.outgoingProved = true;
    auto mb = propagate({m}, names);
    EXPECT_EQ(mb["K"].hi, 2,
              "a mirror alternation takes the worse of 3_1 and m3_1 -- the "
              "same 1 -- so the same bound");
}

// The split rule observed through the public solver: a proved genus-0
// cobordism from a knot K (literature [0, 9], so it constrains nothing) to a
// outgoing link F gives, from the cobordism inequalities with n_0 = 1,
//     lo(K) = lo(F) - 0 - 1 + 1 = lo(F),    hi(K) = hi(F) + 0 + (n_F - 1).
Bounds probeOutgoing(const std::string &outgoing, int nOutgoing, NameTable names) {
    names.addLiterature("K", 0, 9);
    Cobordism w = cobordism("K", 1, outgoing, nOutgoing, 0);
    w.outgoingProved = true;
    auto bounds = propagate({w}, names);
    return bounds["K"];
}

void test_split_unknot_factor_changes_nothing() {
    // g_4(L u U) = g_4(L): a split unknot is tubed on (<=) or capped off by a
    // disc in a collar (>=). So 3_1 u Unknot has g_4 exactly 1.
    NameTable names;
    names.addLiterature("3_1", 1, 1);
    Bounds k = probeOutgoing("3_1 u Unknot", 2, names);
    EXPECT_EQ(k.lo, 1, "lo(3_1 u Unknot) = g_4(3_1) = 1: the unknot drops out");
    EXPECT_EQ(k.hi, 2, "hi = (1 + 0) + (2 - 1) = 2, by hand");
    Bounds k2 = probeOutgoing("3_1 u Unknot u Unknot", 3, names);
    EXPECT_EQ(k2.lo, 1, "any number of split unknots drop out");
}

void test_split_lower_bound_is_not_additive() {
    // THE counterexample to the old additive rule. For every knot K, K u -K
    // bounds an annulus in B^4, so g_4(4_1 u 4_1) = 0 although each factor
    // has g_4 = 1 (4_1 is amphicheiral, so -4_1 = 4_1). The additive rule
    // claimed 2 and, on 2026-09-24, "proved" g_4(L10a91{0}) >= 1 against a
    // literature value of 0. The sound bound,
    //     g_4(4_1 # 4_1) - 1 >= (1 - 1) - 1 = -1,
    // says nothing.
    NameTable names;
    names.addLiterature("4_1", 1, 1);
    Bounds k = probeOutgoing("4_1 u 4_1", 2, names);
    EXPECT_EQ(k.haveLower() && k.lo > 0, false,
              "no positive lower bound through 4_1 u 4_1, whose g_4 is 0");
    EXPECT_EQ(k.hi, 3, "upper still the constructive (1 + 1) + (2 - 1) = 3");

    // The regression, end to end: a slice two-component subject with a
    // genus-0 cobordism to a proved 4_1 u 4_1 is not pushed above its
    // literature value, and nothing is judged a contradiction.
    names.addLiterature("S{0}", 0, 0);
    Cobordism w = cobordism("S{0}", 2, "4_1 u 4_1", 2, 0);
    w.outgoingProved = true;
    auto bounds = propagate({w}, names);
    EXPECT_EQ(bounds["S{0}"].haveLower() && bounds["S{0}"].lo > 0, false,
              "the slice subject keeps a lower bound of at most 0");
    EXPECT_EQ(judge("S{0}", bounds["S{0}"], names).status ==
                  Status::contradiction,
              false, "no contradiction");
}

void test_split_lower_bound_three_factors() {
    // With f nontrivial knot factors, g_4(u) >= g_4(#) - (f - 1), and the sum
    // is bounded by the composite-knot rule. For 3_1 u 5_1 u Unknot the
    // unknot drops out (f = 2):
    //     g_4(3_1 # 5_1) >= g_4(5_1) - g_4(3_1) = 2 - 1 = 1,
    //     g_4(3_1 u 5_1 u Unknot) >= 1 - (2 - 1) = 0,
    // which says nothing (a derived 0 never improves the floor of 0).
    NameTable names;
    names.addLiterature("3_1", 1, 1);
    names.addLiterature("5_1", 2, 2);
    Bounds k = probeOutgoing("3_1 u 5_1 u Unknot", 3, names);
    EXPECT_EQ(k.haveLower() && k.lo > 0, false,
              "(2 - 1) - (2 - 1) = 0: no positive lower bound, by hand");
    EXPECT_EQ(k.hi, 5, "hi = (1 + 2 + 0) + (3 - 1) = 5, by hand");

    // A gap big enough to be positive. The unknot must DROP OUT: counted as
    // a factor it would make f = 3 and the bound (3 - 1) - 2 = 0.
    names.addLiterature("7_1", 3, 3);
    EXPECT_EQ(probeOutgoing("3_1 u 7_1 u Unknot", 3, names).lo, 1,
              "(3 - 1) - (2 - 1) = 1, by hand: the unknot is not a factor");

    // Three nontrivial factors: 11a_367 is T(2,11), g_4 = 5.
    //     g_4(3_1 # 3_1 # 11a_367) >= 5 - 1 - 1 = 3,
    //     g_4(3_1 u 3_1 u 11a_367) >= 3 - (3 - 1) = 1.
    // (The additive rule claimed 7; dropping the -(f - 1) would claim 3.)
    names.addLiterature("11a_367", 5, 5);
    EXPECT_EQ(probeOutgoing("3_1 u 3_1 u 11a_367", 3, names).lo, 1,
              "5 - 1 - 1 - (3 - 1) = 1, by hand");
}

void test_split_mirror_alternatives_use_the_unmirrored_name() {
    // A mirror image has the same g_4, and the table holds only 3_1. A
    // factor written "3_1|m3_1" must still be bounded, both ways, with
    // nothing registered for m3_1.
    NameTable names;
    names.addLiterature("3_1", 1, 1);
    Bounds k = probeOutgoing("3_1|m3_1 u Unknot", 2, names);
    EXPECT_EQ(k.hi, 2, "upper found through the unmirrored name: 1 + 0 + 1");
    EXPECT_EQ(k.lo, 1, "lower likewise: g_4(3_1) = 1");
}

void test_split_composite_alternative_keeps_its_mirror() {
    // Only a PRIME knot's m is stripped. "m3_1#3_1" is the square knot
    // (slice); stripping its leading m would give "3_1#3_1", the granny
    // (g_4 = 2), and borrow the granny's bound. Registered here as literature
    // to make the mix-up visible: the square knot's outgoing link must get only
    // the composite rule's 1 - 1 = 0.
    NameTable names;
    names.addLiterature("3_1", 1, 1);
    names.addLiterature("3_1#3_1", 2, 2);
    Bounds k = probeOutgoing("m3_1#3_1 u Unknot", 2, names);
    EXPECT_EQ(k.haveLower() && k.lo > 0, false,
              "no lower bound borrowed from the granny: 1 - 1 = 0");
    EXPECT_EQ(k.hi, 3, "upper from the summands: (1 + 1) + (2 - 1) = 3");
    // The granny itself, by name, does carry its value.
    EXPECT_EQ(probeOutgoing("3_1#3_1 u Unknot", 2, names).lo, 2,
              "3_1#3_1 u Unknot: g_4(3_1#3_1) = 2");
}

void test_split_with_a_link_factor_has_no_lower_bound() {
    // The band inequality needs a lower bound on g_4 of the connected sum,
    // which only the composite-KNOT rule gives; a link factor leaves it open.
    NameTable names;
    names.addLiterature("L2a1{0}", 0, 0);
    Bounds k = probeOutgoing("L2a1{0} u Unknot", 3, names);
    EXPECT_EQ(k.haveLower(), false, "no lower bound claimed");
}

void test_unproved_split_outgoing_bounds_nothing() {
    // A split name reaching the solver WITHOUT a per-cobordism proof is still
    // a multi-component outgoing link, and the gate refuses it.
    NameTable names;
    names.addLiterature("K", 0, 9);
    names.addLiterature("3_1", 1, 1);
    auto bounds = propagate({cobordism("K", 1, "3_1 u Unknot", 2, 0)}, names);
    EXPECT_EQ(bounds["K"].haveUpper(), false, "no proof, no bound");
}

// ─────────────────────────────────────────────────────────────────────────
// Composite outgoing names: "K #_c L"
// ─────────────────────────────────────────────────────────────────────────

// A named outgoing link bounds by ITS variant, not the worst of its base's, and
// receives a bound from the subject whatever its component count.
void test_named_outgoing() {
    NameTable names;
    names.addLiterature("S", 0, 9);
    names.addLiterature("L7n1{0}", 2, 2);
    names.addLiterature("L7n1{1}", 0, 0);
    Cobordism w = cobordism("S", 1, "L7n1{1}", 2, 0, nameCandidates("L7n1{1}"));
    w.outgoingProved = true;
    w.outgoingNamed = true;
    auto bounds = propagate({w, direct("S", 1, 0)}, names);
    EXPECT_EQ(bounds["S"].hi, 0, "the subject is bounded by its own genus-0 witness");
    EXPECT_EQ(bounds["L7n1{1}"].haveUpper(), true, "an exact far side receives a bound");
    EXPECT_EQ(bounds["L7n1{1}"].hi, 0, "g4(L1) <= g4(S) + g + n(S) - 1 = 0");

    NameTable names2;
    names2.addLiterature("S", 0, 9);
    names2.addLiterature("L7n1{0}", 2, 2);
    names2.addLiterature("L7n1{1}", 0, 0);
    Cobordism v = cobordism("S", 1, "L7n1{1}", 2, 0, nameCandidates("L7n1{1}"));
    v.outgoingProved = true; // proved but not an identity
    auto b2 = propagate({v}, names2);
    EXPECT_EQ(b2["S"].hi, 1, "forward: g4(L7n1{1}) + 0 + 2 - 1, not L7n1{0}'s 2 + 1");
    EXPECT_EQ(b2["L7n1{1}"].haveUpper(), false, "a description receives nothing");
}

// --sum-rules: a sum along components and a split with a link factor bound
// the subject from their pieces; without the flag they bound nothing.
void test_sum_rules() {
    for (bool on : {false, true}) {
        NameTable names;
        names.setSumRules(on);
        names.addLiterature("S", 0, 9);
        names.addLiterature("T", 0, 9);
        names.addLiterature("L7n1{0}", 2, 2);
        names.addLiterature("L2a1{0}", 0, 0);
        Cobordism sum = cobordism("S", 1, "#{L2a1{0}[?] # L7n1{0}[?]}", 3, 0,
                                nameCandidates("#{L2a1{0}[?] # L7n1{0}[?]}"));
        sum.outgoingProved = true;
        Cobordism split = cobordism("T", 1, "L7n1{0} u L2a1{0}", 4, 0,
                                  nameCandidates("L7n1{0} u L2a1{0}"));
        split.outgoingProved = true;
        auto bounds = propagate({sum, split}, names);
        const std::string tag = on ? " (sum rules on)" : " (sum rules off)";
        // Sum: g4 <= 2 + 0; >= 2 - (0 + 2 - 1) = 1. Subject: hi 2 + 0 + 3 - 1
        // = 4, lo 1 - 0 - 1 + 1 = 1.
        EXPECT_EQ(bounds["S"].haveUpper(), on, "sum: an upper bound only with the rule" + tag);
        if (on) EXPECT_EQ(bounds["S"].hi, 4, "sum: g4(S) <= (2 + 0) + 0 + 3 - 1");
        EXPECT_EQ(bounds["S"].haveLower() && bounds["S"].lo == 1, on,
                  "sum: g4(S) >= (2 - (0 + 2 - 1)) - 0 - 1 + 1 = 1 only with the rule" + tag);
        // Split: g4 >= 2 - (0 + 2 - 1) = 1, capping the Hopf link off.
        EXPECT_EQ(bounds["T"].haveLower() && bounds["T"].lo == 1, on,
                  "split with a link factor: lower bound 1 only with the rule" + tag);
    }
}

void test_composite_outgoing_upper_bound() {
    // g_4(K #_c L) <= g_4(K) + g_4(L). With 3_1 at 1 and both orientations of
    // L2a1 at 0, a proved genus-0 cobordism K -> (m3_1 #_0 L2a1) gives
    //     g_4(K) <= (1 + max(0, 0)) + 0 + (2 - 1) = 2.
    NameTable names;
    names.addLiterature("K", 0, 9);
    names.addLiterature("3_1", 1, 1);
    names.addLiterature("L2a1{0}", 0, 0);
    names.addLiterature("L2a1{1}", 0, 0);
    Cobordism w = cobordism("K", 1, "m3_1 #_0 L2a1", 2, 0);
    w.outgoingProved = true;
    auto bounds = propagate({w}, names);
    EXPECT_EQ(bounds["K"].hi, 2, "g_4(K) <= 1 + 0 + 0 + 1 = 2, by hand");
    EXPECT_EQ(bounds["K"].haveLower(), false,
              "the only lower bound is g_4(L) - g_4(K) = 0 - 1 < 0: nothing");
}

void test_composite_takes_the_worst_orientation() {
    // L proved up to orientation, so the bound must hold for either variant.
    NameTable names;
    names.addLiterature("K", 0, 9);
    names.addLiterature("3_1", 1, 1);
    names.addLiterature("L4a1{0}", 0, 0);
    names.addLiterature("L4a1{1}", 1, 1);
    Cobordism w = cobordism("K", 1, "3_1 #_0 L4a1", 2, 0);
    w.outgoingProved = true;
    auto bounds = propagate({w}, names);
    EXPECT_EQ(bounds["K"].hi, 3, "g_4(K) <= 1 + max(0, 1) + 0 + 1 = 3");
}

void test_unproved_composite_bounds_nothing() {
    NameTable names;
    names.addLiterature("K", 0, 9);
    names.addLiterature("3_1", 1, 1);
    names.addLiterature("L2a1{0}", 0, 0);
    names.addLiterature("L2a1{1}", 0, 0);
    auto bounds = propagate({cobordism("K", 1, "3_1 #_0 L2a1", 2, 0)}, names);
    EXPECT_EQ(bounds["K"].haveUpper(), false, "no proof, no bound");
}

// ─────────────────────────────────────────────────────────────────────────
// Composite knots, and sliceness from the concordance group
// ─────────────────────────────────────────────────────────────────────────

namespace {
NameTable symmetryTable() {
    NameTable names;
    names.setSymmetry("3_1", SymmetryType::reversible);
    names.setSymmetry("5_2", SymmetryType::reversible);
    names.setSymmetry("4_1", SymmetryType::fullyAmphicheiral);
    names.setSymmetry("8_17", SymmetryType::negativeAmphicheiral);
    names.setSymmetry("9_32", SymmetryType::chiral);
    return names;
}
} // namespace

void test_elementary_slice_is_a_constructive_anchor() {
    NameTable names = symmetryTable();
    names.addLiterature("K", 0, 9);
    names.addLiterature("5_2", 1, 1);
    auto bounds = propagate({cobordism("K", 1, "5_2#m5_2", 1, 0)}, names);
    EXPECT_EQ(bounds["K"].hi, 0, "g_4(K) <= 0 + 0 + 0 = 0 through the anchor");
    EXPECT_EQ(bounds["K"].basis == Basis::constructive, true,
              "5_2 # m5_2 bounds a ribbon disc: no literature is used, even "
              "though 5_2's own value (1) is in the table");
}

void test_composite_knot_bounds() {
    // g_4 is subadditive, and K_i ~ A # -(rest) gives the lower bound
    //     g_4(A) >= g_4(K_i) - sum_{j != i} g_4(K_j).
    // For 3_1 # 6_1 (6_1 slice): upper 1 + 0 = 1, lower 1 - 0 = 1, so a
    // genus-0 cobordism K -> 3_1#6_1 pins g_4(K) = 1.
    NameTable names;
    names.addLiterature("K", 0, 9);
    names.addLiterature("3_1", 1, 1);
    names.addLiterature("6_1", 0, 0);
    auto bounds = propagate({cobordism("K", 1, "3_1#6_1", 1, 0)}, names);
    EXPECT_EQ(bounds["K"].hi, 1, "g_4(K) <= (1 + 0) + 0 = 1, by hand");
    EXPECT_EQ(bounds["K"].lo, 1, "g_4(K) >= (1 - 0) - 0 - 1 + 1 = 1, by hand");

    NameTable n2;
    n2.addLiterature("K", 0, 9);
    n2.addLiterature("3_1", 1, 1);
    n2.addLiterature("4_1", 1, 1);
    auto b2 = propagate({cobordism("K", 1, "3_1#4_1", 1, 0)}, n2);
    EXPECT_EQ(b2["K"].hi, 2, "3_1 # 4_1: g_4(K) <= 1 + 1 = 2");
    EXPECT_EQ(b2["K"].haveLower(), false,
              "and only |1 - 1| = 0 below: nothing");
}

void test_composite_link_lower_bound() {
    // g_4(K #_c L) >= g_4(L) - g_4(K): summing -K into the same component
    // undoes K up to concordance. With L's variants at 2 and K = 3_1 at 1,
    //     lower(3_1 #_0 L) = 2 - 1 = 1,
    // and a genus-0 cobordism S -> it (proved) gives g_4(S) >= 1 - 0 - 1 + 1.
    NameTable names;
    names.addLiterature("S", 0, 9);
    names.addLiterature("3_1", 1, 1);
    names.addLiterature("L7a1{0}", 2, 2);
    names.addLiterature("L7a1{1}", 2, 2);
    Cobordism w = cobordism("S", 1, "3_1 #_0 L7a1", 2, 0);
    w.outgoingProved = true;
    auto bounds = propagate({w}, names);
    EXPECT_EQ(bounds["S"].lo, 1, "g_4(S) >= (2 - 1) - 0 - 1 + 1 = 1, by hand");
    EXPECT_EQ(bounds["S"].hi, 4, "and g_4(S) <= (1 + 2) + 0 + (2 - 1) = 4");
}

void run(const std::string &name, void (*fn)()) {
    std::cout << bold << "\n=== " << name << " ===" << resetColor << "\n";
    fn();
}

int main() {
    run("external_proofs", test_external_proofs);
    run("derived_lower_above_derived_upper_is_a_contradiction",
        test_derived_lower_above_derived_upper_is_a_contradiction);
    run("knot_outgoing_named_as_a_link_bounds_nothing",
        test_knot_outgoing_named_as_a_link_bounds_nothing);
    run("proved_link_outgoing_bears_bound",
        test_proved_link_outgoing_bears_bound);
    run("split_outgoing_upper_bound", test_split_outgoing_upper_bound);
    run("split_unknot_factor_changes_nothing",
        test_split_unknot_factor_changes_nothing);
    run("split_lower_bound_is_not_additive",
        test_split_lower_bound_is_not_additive);
    run("split_lower_bound_three_factors",
        test_split_lower_bound_three_factors);
    run("split_mirror_alternatives_use_the_unmirrored_name",
        test_split_mirror_alternatives_use_the_unmirrored_name);
    run("split_composite_alternative_keeps_its_mirror",
        test_split_composite_alternative_keeps_its_mirror);
    run("split_with_a_link_factor_has_no_lower_bound",
        test_split_with_a_link_factor_has_no_lower_bound);
    run("unproved_split_outgoing_bounds_nothing",
        test_unproved_split_outgoing_bounds_nothing);
    run("composite_outgoing_upper_bound", test_composite_outgoing_upper_bound);
    run("composite_takes_the_worst_orientation",
        test_composite_takes_the_worst_orientation);
    run("unproved_composite_bounds_nothing",
        test_unproved_composite_bounds_nothing);
    run("elementary_slice_is_a_constructive_anchor",
        test_elementary_slice_is_a_constructive_anchor);
    run("composite_knot_bounds", test_composite_knot_bounds);
    run("composite_link_lower_bound", test_composite_link_lower_bound);
    run("candidates_expand_orientation_variants",
        test_candidates_expand_orientation_variants);
    run("candidates_filtered_by_observed_component_count",
        test_candidates_filtered_by_observed_component_count);
    run("direct_cobordism_gives_upper_bound",
        test_direct_cobordism_gives_upper_bound);
    run("knot_cobordism_reduces_to_the_classic_rule",
        test_knot_cobordism_reduces_to_the_classic_rule);
    run("component_correction_on_the_upper_bound",
        test_component_correction_on_the_upper_bound);
    run("component_correction_on_the_lower_bound",
        test_component_correction_on_the_lower_bound);
    run("unlink_outgoing_carries_no_component_penalty",
        test_unlink_outgoing_carries_no_component_penalty);
    run("linked_far_side_bounds_nothing", test_linked_outgoing_bounds_nothing);
    run("unlink_axiom_is_constructive", test_unlink_axiom_is_constructive);
    run("slice_composite_axiom_is_constructive",
        test_slice_composite_axiom_is_constructive);
    run("slice_composite_allowlist_is_not_a_pattern",
        test_slice_composite_allowlist_is_not_a_pattern);
    run("orientation_variants_are_not_a_candidate_set",
        test_orientation_variants_are_not_a_candidate_set);
    run("outgoing_bears_bound", test_outgoing_bears_bound);
    run("chains_through_an_unnamed_isosig_link",
        test_chains_through_an_unnamed_isosig_link);
    run("tubed_cobordism_genus_is_taken_at_face_value",
        test_tubed_cobordism_genus_is_taken_at_face_value);
    run("propagation_terminates_on_a_cycle",
        test_propagation_terminates_on_a_cycle);
    run("self_cobordism_cannot_confirm_the_literature",
        test_self_cobordism_cannot_confirm_the_literature);
    run("self_cobordism_still_allows_a_real_bound_from_elsewhere",
        test_self_cobordism_still_allows_a_real_bound_from_elsewhere);
    run("ambiguous_outgoing_gets_no_reverse_bound",
        test_ambiguous_outgoing_gets_no_reverse_bound);
    run("unambiguous_outgoing_does_get_a_reverse_bound",
        test_unambiguous_outgoing_does_get_a_reverse_bound);
    run("two_step_cycle_through_an_alias_is_refused",
        test_two_step_cycle_through_an_alias_is_refused);
    run("longer_cycle_is_refused", test_longer_cycle_is_refused);
    run("support_set_records_what_a_bound_rests_on",
        test_support_set_records_what_a_bound_rests_on);
    run("unregistered_multicomponent_outgoing_gets_no_reverse_bound",
        test_unregistered_multicomponent_outgoing_gets_no_reverse_bound);
    run("unregistered_SINGLE_component_outgoing_still_chains",
        test_unregistered_SINGLE_component_outgoing_still_chains);
    run("two_component_isosig_node_does_not_chain",
        test_two_component_isosig_link_does_not_chain);
    run("unlink_outgoing_still_chains_and_is_penalty_free",
        test_unlink_outgoing_still_chains_and_is_penalty_free);
    run("judge_verified", test_judge_verified);
    run("judge_distinguishes_assisted_verification",
        test_judge_distinguishes_assisted_verification);
    run("judge_improved", test_judge_improved);
    run("judge_contradiction", test_judge_contradiction);
    run("judge_unresolved", test_judge_unresolved);
    run("have_cobordism_dedup", test_have_cobordism_dedup);
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
    run("split_boundary_seeded_ignores_names",
        test_split_boundary_seeded_ignores_names);
    run("split_boundary_unnamed_side_flagged",
        test_split_boundary_unnamed_side_flagged);
    run("classify_incoming_orientation", test_classify_incoming_orientation);
    run("cobordism_identity", test_cobordism_identity);
    run("named_outgoing", test_named_outgoing);
    run("sum_rules", test_sum_rules);

    std::cout << bold << "\n=== Summary: " << passed << " passed, "
              << failed_count << " failed ===" << resetColor << "\n";
    return failed_count > 0 ? 1 : 0;
}
