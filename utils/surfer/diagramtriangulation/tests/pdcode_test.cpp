// pdcode_test.cpp
//
// Tests for ../pdcode.h: the one parser of PD text into T's crossings
// (parsePDCode(), labels renumbered from 0; pdLabels(), labels as written) and
// the one formatter (formatPDCode(), in each stored spelling, also respelling a
// given code). The round trip over all 17,153
// table PD codes is in linknaming/tests/tables_test.cpp, beside the tables'
// reader.

#include <array>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

#include "diagramtriangulation/pdcode.h"

using namespace diagramtriangulation;

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

const char *const kTrefoil = "[[1;5;2;4];[3;1;4;6];[5;3;6;2]]";

void test_parse() {
    const PDCode a = parsePDCode(kTrefoil);
    EXPECT_EQ(a.size(), size_t(3), "three crossings");
    EXPECT_EQ(a[0] == (std::array<int, 4>{0, 4, 1, 3}), true,
              "1-based labels are renumbered from 0");
    const PDCode b = parsePDCode("PD[X[0; 4; 1; 3]; X[2; 0; 3; 5]; X[4; 2; 5; 1]]");
    EXPECT_EQ(a == b, true, "a code holding a 0 is read as 0-based, any punctuation");
    EXPECT_EQ(parsePDCode("1 5 2 4 3 1 4 6 5 3 6 2") == a, true, "bare integers");
}

void test_format() {
    const std::vector<std::array<long, 4>> pd{{1, 5, 2, 4}, {3, 1, 4, 6}, {5, 3, 6, 2}};
    EXPECT_EQ(formatPDCode(pd, PDSpelling::semicolons), std::string(kTrefoil),
              "semicolons: the knot table's spelling through 10 crossings, byte for byte");
    EXPECT_EQ(formatPDCode(pd, PDSpelling::commas),
              std::string("[[1,5,2,4],[3,1,4,6],[5,3,6,2]]"),
              "commas: farsidediagram's pd= spelling");
    EXPECT_EQ(formatPDCode(std::vector<std::array<int, 4>>{}, PDSpelling::semicolons),
              std::string("[]"), "the empty code (the unknot's)");
}

// A PD given in another spelling is recorded as formatPDCode(pdLabels(text)):
// labels and crossing order as written, so parsePDCode() -- T's input -- reads
// the record exactly as it read the text.
void test_respell() {
    const std::string linkinfo = "PD[X[4; 1; 3; 2]; X[2; 3; 1; 4]]";
    EXPECT_EQ(pdLabels(linkinfo) == (PDCode{{4, 1, 3, 2}, {2, 3, 1, 4}}), true,
              "pdLabels: LinkInfo's spelling, labels and crossing order as written");
    EXPECT_EQ(pdLabels("PD[X[0; 4; 1; 3]; X[2; 0; 3; 5]; X[4; 2; 5; 1]]")[0] ==
                  (std::array<int, 4>{0, 4, 1, 3}),
              true, "pdLabels: a 0-based code is not renumbered");
    EXPECT_EQ(formatPDCode(pdLabels(linkinfo), PDSpelling::semicolons),
              std::string("[[4;1;3;2];[2;3;1;4]]"), "LinkInfo's spelling respelt with semicolons");
    EXPECT_EQ(formatPDCode(pdLabels("[[1; 5; 2; 4]; [3; 1; 4; 6]; [5; 3; 6; 2]]"),
                           PDSpelling::semicolons),
              std::string(kTrefoil), "the knot table's spaced spelling respelt without spaces");
    for (const char *text : {"PD[X[4; 1; 3; 2]; X[2; 3; 1; 4]]",
                             "PD[X[6; 1; 7; 2]; X[10; 7; 5; 8]; X[4; 5; 1; 6]; X[2; 10; 3; 9]; X[8; 4; 9; 3]]",
                             "PD[X[0; 4; 1; 3]; X[2; 0; 3; 5]; X[4; 2; 5; 1]]",
                             "[[1; 5; 2; 4]; [3; 1; 4; 6]; [5; 3; 6; 2]]"})
        EXPECT_EQ(parsePDCode(formatPDCode(pdLabels(text), PDSpelling::semicolons)) ==
                      parsePDCode(text),
                  true, std::string(text) + ": the respelt code parses as the text does");
}

void test_round_trip() {
    for (const char *text : {kTrefoil, "[[4;2;5;1];[8;6;1;5];[6;3;7;4];[2;7;3;8]]",
                             "PD[X[6; 1; 7; 2]; X[8; 3; 5; 4]; X[2; 5; 3; 6]; X[4; 7; 1; 8]]"}) {
        const PDCode p = parsePDCode(text);
        for (PDSpelling s : {PDSpelling::semicolons, PDSpelling::commas})
            EXPECT_EQ(parsePDCode(formatPDCode(p, s)) == p, true,
                      std::string(text) + ": parse(format(parse)) == parse");
    }
}

} // namespace

void run(const std::string &name, void (*fn)()) {
    std::cout << bold << "\n=== " << name << " ===" << resetColor << "\n";
    fn();
}

int main() {
    run("parse", test_parse);
    run("format", test_format);
    run("respell", test_respell);
    run("round_trip", test_round_trip);

    std::cout << bold << "\n=== Summary: " << passed << " passed, "
              << failed_count << " failed ===" << resetColor << "\n";
    return failed_count > 0 ? 1 : 0;
}
