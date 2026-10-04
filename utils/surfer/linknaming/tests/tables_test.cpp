// tables_test.cpp
//
// Tests for ../tables.h: the one reader of the knot and link tables, the one
// parser of their 4-genus field, the one symmetry enum and its reader, the
// PD text -> regina::Link parser and the slice-composite anchors
// (isElementarySlice); and, over every table PD code, the round trip through
// diagramtriangulation/pdcode.h's parser and formatter. The two isElementarySlice cases moved here from the
// solver's test (solver_test.cpp) with the function (phase 3), their
// NameTable becoming the SymmetryTable it now takes.

#include <array>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>
#include <utility>

#include "diagramtriangulation/pdcode.h"
#include "linknaming/tables.h"

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

// A file in a fresh temporary directory, removed with it.
struct TempFile {
    std::filesystem::path dir, path;
    explicit TempFile(const std::string &text) {
        char tmpl[] = "/tmp/tables_test.XXXXXX";
        dir = mkdtemp(tmpl);
        path = dir / "table.csv";
        std::ofstream(path, std::ios::binary) << text;
    }
    ~TempFile() { std::filesystem::remove_all(dir); }
};

// The table's literature 4-genus: what parses, and that a malformed value
// never becomes a bound (moved from the cascade's leaves test, F1).
void test_parse_table_g4() {
    EXPECT_EQ(parseTableG4("2") == std::make_pair(2, 2), true, "plain value");
    EXPECT_EQ(parseTableG4("[0;1]") == std::make_pair(0, 1), true, "interval");
    EXPECT_EQ(parseTableG4("[1;2]") == std::make_pair(1, 2), true, "interval 1..2");
    for (const char *bad : {"", "x", "[1;0]", "[;1]", "[0;]", "[0,1]", "-1", "1.5", "[0;1",
                            " 1", "1 ", "1x", "[0;1] "})
        EXPECT_EQ(parseTableG4(bad).has_value(), false, std::string("malformed refused: '") + bad + "'");
}

void test_read_table_rows() {
    TempFile f("Name,PD Notation,Genus-4D\r\n"
               "3_1,[[1;5;2;4];[3;1;4;6];[5;3;6;2]],1\r\n"
               "\n"
               "no commas here\n"
               "only,one comma\n"
               "L2a1{0},PD[X[4; 1; 3; 2]; X[2; 3; 1; 4]],0\n"
               "L10n74{1;0},PD[X[1; 2; 3; 4]],[0;1]");
    const std::vector<TableRow> rows = readTableRows(f.path);
    EXPECT_EQ(rows.size(), size_t(3), "header, empty lines and lines without two commas skipped");
    EXPECT_EQ(rows[0].name, std::string("3_1"), "name");
    EXPECT_EQ(rows[0].pd, std::string("[[1;5;2;4];[3;1;4;6];[5;3;6;2]]"), "PD as written");
    EXPECT_EQ(rows[0].g4, std::string("1"), "g4 without the CR");
    EXPECT_EQ(rows[1].pd, std::string("PD[X[4; 1; 3; 2]; X[2; 3; 1; 4]]"), "a link PD, spaces kept");
    EXPECT_EQ(rows[2].name, std::string("L10n74{1;0}"), "a tag's ';' is not a separator");
    EXPECT_EQ(rows[2].g4, std::string("[0;1]"), "an interval, unterminated last line");
    bool threw = false;
    try {
        readTableRows("/nonexistent/table.csv");
    } catch (const regina::InvalidArgument &) {
        threw = true;
    }
    EXPECT_EQ(threw, true, "an unreadable table is an error, never an empty one");
}

// Every spelling KnotInfo uses, which both readers of knot_symmetry.csv
// accepted before there was one (review S3.14); nothing else.
void test_symmetry_types() {
    EXPECT_EQ(parseSymmetryType("chiral") == SymmetryType::chiral, true, "chiral");
    EXPECT_EQ(parseSymmetryType("reversible") == SymmetryType::reversible, true, "reversible");
    EXPECT_EQ(parseSymmetryType("positive amphicheiral") == SymmetryType::positiveAmphicheiral,
              true, "positive amphicheiral");
    EXPECT_EQ(parseSymmetryType("negative amphicheiral") == SymmetryType::negativeAmphicheiral,
              true, "negative amphicheiral");
    EXPECT_EQ(parseSymmetryType("fully amphicheiral") == SymmetryType::fullyAmphicheiral, true,
              "fully amphicheiral");
    for (const char *bad : {"", "Chiral", "amphicheiral", "fully amphichiral", "reversible "})
        EXPECT_EQ(parseSymmetryType(bad).has_value(), false, std::string("refused: '") + bad + "'");

    TempFile f("name,symmetry_type,basis\r\n"
               "3_1,reversible,knotinfo+snappy\r\n"
               "4_1,fully amphicheiral,knotinfo+snappy\n"
               "8_17,negative amphicheiral\n"
               "x_1,unknown,knotinfo\n"
               "nocomma\n");
    const SymmetryTable t = readSymmetryTable(f.path);
    EXPECT_EQ(t.size(), size_t(3), "the rows whose type parses");
    EXPECT_EQ(t.at("3_1") == SymmetryType::reversible, true, "3_1");
    EXPECT_EQ(t.at("4_1") == SymmetryType::fullyAmphicheiral, true, "4_1 (CRLF before it)");
    EXPECT_EQ(t.at("8_17") == SymmetryType::negativeAmphicheiral, true, "a row without basis");
}

// The two spellings of a table PD code read as the same diagram, labels as
// written.
void test_link_from_table_pd() {
    regina::Link a = linkFromTablePD("[[4;1;3;2];[2;3;1;4]]");
    regina::Link b = linkFromTablePD("PD[X[4; 1; 3; 2]; X[2; 3; 1; 4]]");
    EXPECT_EQ(a.sig<2>(false, false, true), b.sig<2>(false, false, true),
              "knot-table and link-table spellings give the same oriented diagram");
    EXPECT_EQ(a.countComponents(), size_t(2), "the Hopf link");
    bool threw = false;
    try {
        linkFromTablePD("[[1;2;3]]");
    } catch (const regina::InvalidArgument &) {
        threw = true;
    }
    EXPECT_EQ(threw, true, "a label count not divisible by 4 is refused");
}

// The atlas's tables, from SURFER_TEST_ATLAS_DATA or the usual checkouts;
// empty when there are none (the whole-table cases then say they skipped).
std::filesystem::path atlasData() {
    if (const char *d = std::getenv("SURFER_TEST_ATLAS_DATA")) return d;
    const char *home = std::getenv("HOME");
    for (const char *rel : {"/Projects/cobordism-atlas/data", "/Projects/triangles/cobordism-atlas/data"})
        if (home && std::filesystem::exists(std::string(home) + rel + "/knot_symmetry.csv"))
            return std::string(home) + rel;
    return {};
}

// Every table PD code, both ways (plan: "the PD formatter and parsers,
// round-tripped over all 17,153 table PDs"): T's parser (parsePDCode) and
// the formatter, and PD text -> regina::Link (linkFromTablePD) and back
// through the formatter; and the formatter reproduces the knot table's own
// spelling of every knot.
void test_round_trip_every_table_pd() {
    const std::filesystem::path data = atlasData();
    if (data.empty()) {
        std::cout << "  SKIPPED: no atlas tables (set SURFER_TEST_ATLAS_DATA)\n";
        return;
    }
    size_t rows = 0, tParse = 0, link = 0, knotSpelling = 0, knots = 0;
    for (const char *file : {"4d_smooth_slice_genus_13_crossings_pd_codes.csv",
                             "links_4d_smooth_slice_genus_11_crossings_pd_codes.csv"})
        for (const TableRow &row : readTableRows(data / file)) {
            ++rows;
            const knotbuilder::PDCode p = knotbuilder::parsePDCode(row.pd);
            bool ok = true;
            for (knotbuilder::PDSpelling sp :
                 {knotbuilder::PDSpelling::semicolons, knotbuilder::PDSpelling::commas})
                ok = ok && knotbuilder::parsePDCode(knotbuilder::formatPDCode(p, sp)) == p;
            tParse += ok;
            const regina::Link l = linkFromTablePD(row.pd);
            const std::string text =
                knotbuilder::formatPDCode(l.pdData(), knotbuilder::PDSpelling::semicolons);
            const regina::Link back = linkFromTablePD(text);
            link += back.sig<2>(false, false, true) == l.sig<2>(false, false, true) &&
                    knotbuilder::formatPDCode(back.pdData(),
                                              knotbuilder::PDSpelling::semicolons) == text;
            if (row.pd.rfind("[[", 0) == 0) {
                ++knots;
                std::vector<std::array<int, 4>> asWritten = p;
                for (auto &x : asWritten)
                    for (int &v : x) ++v;
                // The table writes "[[1;5;2;4];..." through 10 crossings and
                // "[[3; 1; 4; 26]; [1; ..." from 11: the same up to spaces.
                std::string compact = row.pd;
                std::erase(compact, ' ');
                knotSpelling +=
                    knotbuilder::formatPDCode(asWritten, knotbuilder::PDSpelling::semicolons) ==
                    compact;
            }
        }
    EXPECT_EQ(rows, size_t(17153), "every table row read");
    EXPECT_EQ(tParse, rows, "parsePDCode(formatPDCode(p)) == p in both spellings, every row");
    EXPECT_EQ(link, rows,
              "linkFromTablePD(format(pdData)) is the same oriented diagram, every row");
    EXPECT_EQ(knotSpelling, knots,
              "the formatter spells every knot table PD as the table does, up to spaces");
}

// The signature table diagram naming reads, from the tables the process
// already loaded (one table load, phase 5), is the one fromTables() reads
// from the files: the same signatures and the same first-wins names. On a
// small table (always), and on the atlas's whole tables when present.
void test_signature_table_from_loaded_tables() {
    char tmpl[] = "/tmp/sigtable_XXXXXX";
    const std::filesystem::path dir = mkdtemp(tmpl);
    std::ofstream(dir / "knots.csv")
        << "Name,PD Notation,Genus-4D\n"
           "3_1,[[1;5;2;4];[3;1;4;6];[5;3;6;2]],1\n"
           "4_1,[[4;2;5;1];[8;6;1;5];[6;3;7;4];[2;7;3;8]],1\n";
    std::ofstream(dir / "links.csv")
        << "Name,PD Notation,Genus-4D\n"
           "L2a1{0},PD[X[4; 1; 3; 2]; X[2; 3; 1; 4]],0\n"
           "L2a1{1},PD[X[4; 2; 3; 1]; X[2; 4; 1; 3]],0\n";
    auto compare = [](const std::filesystem::path &k, const std::filesystem::path &l,
                      const char *what) {
        const linknaming::SignatureTable fromFiles =
            linknaming::SignatureTable::fromTables(k.string(), l.string());
        const Tables tables = Tables::load(k.string(), l.string(), "");
        const linknaming::SignatureTable shared = linknaming::SignatureTable::fromTables(tables);
        EXPECT_EQ(shared.knots(), fromFiles.knots(), std::string(what) + ": knot signatures");
        EXPECT_EQ(shared.links(), fromFiles.links(), std::string(what) + ": link signatures");
        EXPECT_EQ(shared == fromFiles, true,
                  std::string(what) + ": the same table from the loaded tables as from the files");
    };
    compare(dir / "knots.csv", dir / "links.csv", "small table");
    std::filesystem::remove_all(dir);
    const std::filesystem::path data = atlasData();
    if (data.empty()) {
        std::cout << "  SKIPPED the whole tables: no atlas tables (set SURFER_TEST_ATLAS_DATA)\n";
        return;
    }
    compare(data / "4d_smooth_slice_genus_13_crossings_pd_codes.csv",
            data / "links_4d_smooth_slice_genus_11_crossings_pd_codes.csv", "the atlas tables");
}

SymmetryTable symmetryTable() {
    SymmetryTable names;
    names["3_1"] = SymmetryType::reversible;
    names["5_2"] = SymmetryType::reversible;
    names["4_1"] = SymmetryType::fullyAmphicheiral;
    names["8_17"] = SymmetryType::negativeAmphicheiral;
    names["9_32"] = SymmetryType::chiral;
    return names;
}

void test_elementary_slice() {
    SymmetryTable none;
    EXPECT_EQ(isElementarySlice("3_1#m3_1", none), true,
              "the long-standing anchors survive with no symmetry data");
    EXPECT_EQ(isElementarySlice("4_1#4_1", none), true, "both of them");
    EXPECT_EQ(isElementarySlice("5_2#m5_2", none), false,
              "anything else needs the summand's symmetry type");

    SymmetryTable n = symmetryTable();
    EXPECT_EQ(isElementarySlice("m3_1#3_1", n), true,
              "reversible: -3_1 = m3_1, in either spelling");
    EXPECT_EQ(isElementarySlice("3_1#3_1", n), false,
              "the granny knot is not slice");
    EXPECT_EQ(isElementarySlice("4_1#m4_1", n), true,
              "fully amphicheiral: m4_1 = 4_1 = -4_1");
    EXPECT_EQ(isElementarySlice("5_2#m5_2#4_1#4_1", n), true,
              "a sum of inverse pairs is slice");
    EXPECT_EQ(isElementarySlice("3_1#m3_1#4_1", n), false,
              "an unpaired summand spoils it");
    EXPECT_EQ(isElementarySlice("8_17#8_17", n), true,
              "negative amphicheiral: -8_17 = 8_17, so 8_17 # 8_17 is slice");
    EXPECT_EQ(isElementarySlice("8_17#m8_17", n), false,
              "but 8_17 # m8_17 = 8_17 # 8_17^r is not -- the classic trap");
    EXPECT_EQ(isElementarySlice("9_32#m9_32", n), false,
              "chiral non-invertible: -K = m(K^r) is not what the name says");
    EXPECT_EQ(isElementarySlice("6_1#m6_1", n), false,
              "unknown symmetry type: refused, never guessed");
}

void test_elementary_slice_with_marks() {
    SymmetryTable names;
    names["3_1"] = SymmetryType::reversible;
    names["9_32"] = SymmetryType::chiral;
    names["8_17"] = SymmetryType::negativeAmphicheiral;
    names["12a_1"] = SymmetryType::positiveAmphicheiral;
    EXPECT_EQ(isElementarySlice("9_32#mr9_32", names), true, "chiral: K # mrK = K # -K");
    EXPECT_EQ(isElementarySlice("9_32#m9_32", names), false, "chiral: K # mK is not");
    EXPECT_EQ(isElementarySlice("3_1#mr3_1", names), true, "reversible: the r is meaningless");
    EXPECT_EQ(isElementarySlice("8_17#8_17", names), true, "negative amphicheiral: -K = K");
    EXPECT_EQ(isElementarySlice("8_17#r8_17", names), false, "8_17 # r8_17 = 8_17 # m8_17 is not");
    EXPECT_EQ(isElementarySlice("12a_1#r12a_1", names), true, "positive amphicheiral: -K = rK");
    EXPECT_EQ(isElementarySlice("12a_1#12a_1", names), false, "and K # K is not");
}

} // namespace

void run(const std::string &name, void (*fn)()) {
    std::cout << bold << "\n=== " << name << " ===" << resetColor << "\n";
    fn();
}

int main() {
    run("parse_table_g4", test_parse_table_g4);
    run("read_table_rows", test_read_table_rows);
    run("symmetry_types", test_symmetry_types);
    run("link_from_table_pd", test_link_from_table_pd);
    run("elementary_slice", test_elementary_slice);
    run("elementary_slice_with_marks", test_elementary_slice_with_marks);
    run("round_trip_every_table_pd", test_round_trip_every_table_pd);
    run("signature_table_from_loaded_tables", test_signature_table_from_loaded_tables);

    std::cout << bold << "\n=== Summary: " << passed << " passed, "
              << failed_count << " failed ===" << resetColor << "\n";
    return failed_count > 0 ? 1 : 0;
}
