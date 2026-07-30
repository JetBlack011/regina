// census_test.cpp
//
// Tests for census::localCensusLookup() (see ../identifycomplement.h/.cpp):
// the SQLite-backed local census that identify::resolveRecognition() checks
// before falling back to the real, mutex-guarded regina::Census::lookup().
// These tests build a small scratch .sqlite fixture directly (rather than
// the real, multi-hundred-thousand-row census tools/gen_census.py
// produces), so they exercise the lookup/formatting/fallthrough logic in
// isolation, independent of whether the real census has been generated on
// this machine.

#include <cstdio>
#include <iostream>
#include <optional>
#include <string>
#include <unistd.h>

#include <sqlite3.h>

#include <triangulation/dim3.h>
#include <triangulation/example3.h>

#include "../identifycomplement.h"
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

const char *FIXTURE_PATH = "census_test_fixture.sqlite";

// Builds a tiny scratch census at FIXTURE_PATH: one source='regina' row
// whose raw name (with a " : #N" suffix, matching CensusHit::name()'s real
// format) has a known Rolfsen translation in rolfsentable.h, one whose raw
// name doesn't, and one source='snappy' row with an already-pretty name
// (as identify_boundaries.py's _pick_best_name() would produce).
void buildFixture() {
    std::remove(FIXTURE_PATH);

    sqlite3 *db = nullptr;
    sqlite3_open(FIXTURE_PATH, &db);
    sqlite3_exec(db,
        "CREATE TABLE census ("
        "  isosig TEXT PRIMARY KEY, name TEXT NOT NULL, source TEXT NOT NULL"
        ");"
        "INSERT INTO census VALUES"
        "  ('fake-sig-regina-rolfsen', 'L104001 : #1', 'regina'),"
        "  ('fake-sig-regina-no-rolfsen', 'not-in-rolfsen-table', 'regina'),"
        "  ('fake-sig-snappy', 'K12n124', 'snappy');",
        nullptr, nullptr, nullptr);
    sqlite3_close(db);
}

void test_regina_hit_formats_via_rolfsen() {
    buildFixture();
    census::setCensusPath(FIXTURE_PATH);

    EXPECT_EQ(
        census::localCensusLookup("fake-sig-regina-rolfsen").value_or("<MISS>"),
        std::string("4_1 (L104001 : #1)"),
        "a source='regina' hit whose base name is in rolfsentable.h is "
        "formatted the same way censusLookupName() formats a real "
        "Census::lookup() hit");
}

void test_regina_hit_raw_when_no_rolfsen_entry() {
    buildFixture();
    census::setCensusPath(FIXTURE_PATH);

    EXPECT_EQ(
        census::localCensusLookup("fake-sig-regina-no-rolfsen").value_or("<MISS>"),
        std::string("not-in-rolfsen-table"),
        "a source='regina' hit with no rolfsentable.h entry returns the "
        "raw census name unchanged");
}

void test_snappy_hit_returns_verbatim() {
    buildFixture();
    census::setCensusPath(FIXTURE_PATH);

    EXPECT_EQ(census::localCensusLookup("fake-sig-snappy").value_or("<MISS>"),
              std::string("K12n124"),
              "a source='snappy' hit is returned exactly as stored, never "
              "passed through rolfsenName() (its name isn't a raw census "
              "name, so that would only ever miss)");
}

void test_miss_returns_nullopt() {
    buildFixture();
    census::setCensusPath(FIXTURE_PATH);

    EXPECT_EQ(census::localCensusLookup("not-a-real-isosig").has_value(), false,
              "an isoSig absent from the census is a clean miss "
              "(nullopt), for callers to fall through to "
              "censusLookupName()");
}

void test_missing_file_falls_through_cleanly() {
    census::resetCensusForTesting();

    EXPECT_EQ(census::localCensusLookup("fake-sig-regina-rolfsen").has_value(),
              false,
              "with no census file present at all, every lookup misses "
              "cleanly -- no crash, same as an isoSig genuinely absent "
              "from a present census");
}

void test_reenable_after_reset() {
    // Exercises the generation-bump invalidation mechanism itself:
    // CensusConnection_'s thread-local connection, opened for the fixture
    // in an earlier test, must notice resetCensusForTesting()'s path
    // change above and NOT still be pointing at (or reusing a stale
    // negative result from) the fixture or the nonexistent path.
    buildFixture();
    census::setCensusPath(FIXTURE_PATH);

    EXPECT_EQ(census::localCensusLookup("fake-sig-snappy").value_or("<MISS>"),
              std::string("K12n124"),
              "switching setCensusPath() back to a real fixture after "
              "resetCensusForTesting() re-enables lookups");
}

// Builds a single-row fixture keyed by a caller-supplied signature -- used
// below to insert a row keyed by a REAL triangulation's REAL isoSig(),
// rather than the synthetic placeholder strings the tests above use.
void buildSingleRowFixture(const std::string &sig, const std::string &name,
                           const std::string &source) {
    std::remove(FIXTURE_PATH);

    sqlite3 *db = nullptr;
    sqlite3_open(FIXTURE_PATH, &db);
    sqlite3_exec(db,
        "CREATE TABLE census ("
        "  isosig TEXT PRIMARY KEY, name TEXT NOT NULL, source TEXT NOT NULL"
        ");",
        nullptr, nullptr, nullptr);

    sqlite3_stmt *stmt = nullptr;
    sqlite3_prepare_v2(db, "INSERT INTO census VALUES (?, ?, ?)", -1, &stmt,
                       nullptr);
    sqlite3_bind_text(stmt, 1, sig.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 2, name.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 3, source.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_step(stmt);
    sqlite3_finalize(stmt);
    sqlite3_close(db);
}

// Guards against the exact bug this census once had: tools/gen_census.py
// originally keyed regina-sourced rows by Triangulation<3>::neoSig() (what
// the real census databases -- and regina::Census::lookup() -- actually use,
// see census-impl.h), while every runtime query uses Triangulation<3>::
// isoSig() (an unrelated, older signature scheme -- see that script's
// docstring for the full story). That mismatch made every census lookup
// miss silently, even for isoSigs genuinely covered by the real census, and
// none of the tests above would have caught it: they all insert and query
// with the same synthetic placeholder string by construction, so a
// key-format mismatch in *real* usage can't show up there. This test uses a
// real triangulation and its real isoSig() throughout, going through the
// full identify::identify() path (not just census::localCensusLookup()
// directly), so a future regression reintroducing any key-format mismatch
// fails here.
void test_real_triangulation_isosig_key_matches_production_query() {
    regina::Triangulation<3> figureEight = regina::Example<3>::figureEight();
    std::string sig = figureEight.isoSig();

    buildSingleRowFixture(sig, "synthetic-name-for-fig8", "snappy");
    census::setCensusPath(FIXTURE_PATH);
    identify::resetRecognitionCacheForTesting();

    EXPECT_EQ(census::localCensusLookup(sig).value_or("<MISS>"),
              std::string("synthetic-name-for-fig8"),
              "a fixture row keyed by a REAL triangulation's isoSig() is a "
              "direct census::localCensusLookup() hit");

    EXPECT_EQ(identify::identify(EdgeComplement(figureEight, {})),
              std::string("synthetic-name-for-fig8"),
              "...and the SAME key is what identify::resolveRecognition()/"
              "identify::identify() actually queries with, end to end (not "
              "some other signature computed independently)");

    identify::resetRecognitionCacheForTesting();
}

} // namespace

void run(const std::string &name, void (*fn)()) {
    std::cout << bold << "\n=== " << name << " ===" << resetColor << "\n";
    fn();
}

int main() {
    run("regina_hit_formats_via_rolfsen", test_regina_hit_formats_via_rolfsen);
    run("regina_hit_raw_when_no_rolfsen_entry",
        test_regina_hit_raw_when_no_rolfsen_entry);
    run("snappy_hit_returns_verbatim", test_snappy_hit_returns_verbatim);
    run("miss_returns_nullopt", test_miss_returns_nullopt);
    run("missing_file_falls_through_cleanly",
        test_missing_file_falls_through_cleanly);
    run("reenable_after_reset", test_reenable_after_reset);
    run("real_triangulation_isosig_key_matches_production_query",
        test_real_triangulation_isosig_key_matches_production_query);

    std::remove(FIXTURE_PATH);

    std::cout << bold << "\n=== Summary: " << passed << " passed, "
              << failed_count << " failed ===" << resetColor << "\n";
    return failed_count > 0 ? 1 : 0;
}
