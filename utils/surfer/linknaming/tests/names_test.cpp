// names_test.cpp
//
// Tests for ../names.h: the grammar of the names an outgoing link is recorded
// under -- component counts, bases and tags, census suffixes, splits,
// composites, sums and knot marks. Pure string functions. The cases from
// componentsFromName() to knotSummands() moved here, verbatim, from the
// solver's test (solver_test.cpp) when the name helpers became one module
// (phase 3).

#include <iostream>
#include <string>
#include <unistd.h>
#include <utility>

#include "linknaming/names.h"

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

void test_normalize_complement_name() {
    // census::nameComplement() decorates a translated census hit; the
    // tables do not. Left undecorated, "4_1 (m004 : #1)" would be a
    // DIFFERENT name to the solver from the "4_1" row for the very same knot, and nothing would
    // ever chain through it.
    EXPECT_EQ(normalizeComplementName("4_1 (m004 : #1)"), std::string("4_1"),
              "the census annotation is stripped for node identity");
    EXPECT_EQ(normalizeComplementName("L6a1 (s780 : #6)"),
              std::string("L6a1"), "same for links");
    EXPECT_EQ(normalizeComplementName("Unknot"), std::string("Unknot"),
              "an undecorated name is untouched");
    EXPECT_EQ(normalizeComplementName("3-component unlink"),
              std::string("3-component unlink"), "unlinks are untouched");
    EXPECT_EQ(normalizeComplementName("cPcbbbadu"), std::string("cPcbbbadu"),
              "a bare isoSig is untouched");
}

void test_split_names() {
    EXPECT_EQ(componentsFromName("3_1 u Unknot"), 2, "knot u knot: 1 + 1");
    EXPECT_EQ(componentsFromName("Unknot u L2a1{0}"), 3,
              "a tagged link factor states its own count: 1 + 2");
    EXPECT_EQ(componentsFromName("3_1|m3_1 u Unknot"), 2,
              "a mirror alternation is still one component");
    EXPECT_EQ(baseName("L2a1{0} u Unknot"),
              std::string("L2a1{0} u Unknot"),
              "a split name has no base: stripping at '{' would read it as "
              "the 2-component link L2a1 and hand back that link's "
              "orientation variants");
}

void test_composite_parts() {
    auto a = compositeParts("m3_1 #_0 L2a1");
    EXPECT_EQ(a.has_value(), true, "parses");
    EXPECT_EQ(a->knot, std::string("m3_1"), "knot keeps its chirality");
    EXPECT_EQ(a->component, 0, "component");
    EXPECT_EQ(a->link, std::string("L2a1"), "link base name");
    EXPECT_EQ(compositeParts("3_1 #_1 L10n77")->component, 1, "component 1");
    EXPECT_EQ(compositeParts("3_1 u Unknot").has_value(), false,
              "a split name is not composite");
    EXPECT_EQ(compositeParts("3_1#m3_1").has_value(), false,
              "the slice composite (a KNOT sum) is not a knot-into-link sum");
    EXPECT_EQ(compositeParts("3_1 #_ L2a1").has_value(), false,
              "no component index: refused, since K # L is not well defined "
              "without one");
    // A name writes L with its orientation tag, and the component as
    // "?" when the namer does not compute it.
    EXPECT_EQ(compositeParts("3_1 #_0 L2a1{0}").has_value(), true,
              "a TAGGED link part (exact names) is accepted");
    EXPECT_EQ(compositeParts("3_1 #_0 L2a1{0}")->link, std::string("L2a1{0}"),
              "and kept with its tag");
    EXPECT_EQ(compositeParts("3_1 #_? L2a1{0}")->component, -1,
              "an unindexed component is -1");
    EXPECT_EQ(compositeParts("mr8_17 #_? L2a1{0}")->knot, std::string("mr8_17"),
              "a marked knot keeps its marks");
    EXPECT_EQ(compositeParts("3_1 #_0 L2a1{x}").has_value(), false,
              "a malformed tag is refused");
}

void test_knot_marks() {
    EXPECT_EQ(stripKnotMarks("mr8_17"), std::string("8_17"), "m and r stripped");
    EXPECT_EQ(stripKnotMarks("r8_17"), std::string("8_17"), "r stripped");
    EXPECT_EQ(stripKnotMarks("m3_1"), std::string("3_1"), "m stripped");
    EXPECT_EQ(stripKnotMarks("11n_34"), std::string("11n_34"), "no marks, unchanged");
    EXPECT_EQ(stripKnotMarks("mL2a1{0}"), std::string("mL2a1{0}"), "not a knot: unchanged");
    EXPECT_EQ(knotSummands("3_1#mr8_17").size(), size_t(2), "a marked summand is a summand");
}

void test_sum_pieces() {
    auto p = sumPieces("5_1 #_? 3_1 #_? L4a1{1}");
    EXPECT_EQ(p.has_value() && p->size() == 3, true, "knots summed into a link: three pieces");
    EXPECT_EQ(p && (*p)[2] == std::make_pair(std::string("L4a1{1}"), 2), true, "the link, 2 components");
    auto q = sumPieces("#{L2a1{0}[?] # L2a1{0}[?] ; L2a1{0}[?] # L6a3{0}[?]}");
    EXPECT_EQ(q.has_value() && q->size() == 4, true,
              "every written occurrence is a piece (over-counts, which only weakens)");
    EXPECT_EQ(sumPieces("L2a1{0}#3_1").has_value(), false, "not a sum along components");
    EXPECT_EQ(sumPieces("3_1 #_? L2a1").has_value(), false, "an untagged link states no count");
    EXPECT_EQ(nameCandidates("L8n2{0}|L8n2{1}").size(), size_t(2), "proved alternatives");
    EXPECT_EQ(nameCandidates("3_1#3_1|3_1#m3_1 u Unknot").size(), size_t(1),
              "a split is one candidate");
}

void test_knot_summands() {
    EXPECT_EQ(knotSummands("3_1#m5_2").size(), size_t(2), "two summands");
    EXPECT_EQ(knotSummands("3_1#m5_2")[1], std::string("m5_2"), "chirality kept");
    EXPECT_EQ(knotSummands("11a_367#m11a_367").size(), size_t(2),
              "11+ crossing names (with a/n) parse");
    EXPECT_EQ(knotSummands("5_2").size(), size_t(0), "a prime is not composite");
    EXPECT_EQ(knotSummands("m129 : #2").size(), size_t(0),
              "a census hit suffix is not a connected sum");
    EXPECT_EQ(knotSummands("3_1 #_0 L2a1").size(), size_t(0),
              "a knot-into-link sum is not a composite KNOT");
    EXPECT_EQ(knotSummands("L2a1{0}#3_1").size(), size_t(0), "links refused");
}

// The plain cuts the name helpers share (phase 3: each was written out by
// hand at several call sites).
void test_strip_orientation_tag() {
    EXPECT_EQ(stripOrientationTag("L6a3{0}"), std::string("L6a3"), "tag cut");
    EXPECT_EQ(stripOrientationTag("L11n459{0;1;0}"), std::string("L11n459"),
              "multi-entry tag cut");
    EXPECT_EQ(stripOrientationTag("10_132"), std::string("10_132"),
              "no tag, unchanged");
    EXPECT_EQ(stripOrientationTag("L2a1{0} u Unknot"), std::string("L2a1"),
              "cut at the first brace, split or not (baseName() is the "
              "split-aware one)");
    EXPECT_EQ(baseName("L2a1{0} u Unknot"), std::string("L2a1{0} u Unknot"),
              "and baseName() still leaves a split alone");
}

void test_strip_census_suffix() {
    EXPECT_EQ(stripCensusSuffix("m004 : #1"), std::string("m004"),
              "the varying #N is cut");
    EXPECT_EQ(stripCensusSuffix("L104001"), std::string("L104001"),
              "a Christy name has no suffix");
    EXPECT_EQ(stripCensusSuffix("4_1"), std::string("4_1"), "a table name");
    EXPECT_EQ(stripCensusSuffix("s780 : #6 extra"), std::string("s780"),
              "everything from the first \" : \" goes");
}

} // namespace

void run(const std::string &name, void (*fn)()) {
    std::cout << bold << "\n=== " << name << " ===" << resetColor << "\n";
    fn();
}

int main() {
    run("components_from_name", test_components_from_name);
    run("base_name", test_base_name);
    run("normalize_complement_name", test_normalize_complement_name);
    run("split_names", test_split_names);
    run("composite_parts", test_composite_parts);
    run("knot_marks", test_knot_marks);
    run("sum_pieces", test_sum_pieces);
    run("knot_summands", test_knot_summands);
    run("strip_orientation_tag", test_strip_orientation_tag);
    run("strip_census_suffix", test_strip_census_suffix);

    std::cout << bold << "\n=== Summary: " << passed << " passed, "
              << failed_count << " failed ===" << resetColor << "\n";
    return failed_count > 0 ? 1 : 0;
}
