//
//  gen_knot_census_names.cpp
//
//  Standalone generation utility (not a CTest test -- see CMakeLists.txt):
//  for every knot in a KnotInfo-style PD code CSV ("Name,PD Notation" rows)
//  up to a given crossing number, and (with --links) every link in a
//  Thistlethwaite-Link-Table-style PD code CSV ("Name,PD Notation,
//  Genus-4D" rows, deduped to one row per base name -- see
//  stripOrientationSuffix() below, since component orientation doesn't
//  change the complement), builds the complement the same way
//  surfer.cpp/verifyslicegenus.cpp does at runtime (knotbuilder::buildLink()
//  -> Knot/Link -> buildComplement()) and records every
//  regina::Census::lookup() hit against it.
//
//  This produces the raw data behind utils/surfer/knot_census_names.csv
//  (and, with --links, utils/surfer/link_census_names.csv), which
//  utils/surfer/tests/gen_knot_names_header.py then turns into
//  utils/surfer/linknames.h (a census-name -> classical-name lookup table
//  covering both Rolfsen knot names and Thistlethwaite link names).
//
//  Usage:
//      gen_knot_census_names <pd_codes.csv> [maxCrossings]
//      gen_knot_census_names --links <links_pd_codes.csv> [maxCrossings]
//
//  CSV data rows (name,crossings,census_db,census_name,isosig) go to stdout
//  (one per census hit; a single row with empty census_db/census_name
//  fields for entries with no hit, e.g. torus knots or non-hyperbolic
//  links, which cannot appear in any of Regina's hyperbolic-only census
//  databases -- the isosig column is still populated on a miss, so it can
//  be fed to an offline SnapPy identification pass the same way
//  custom_analysis/identify_boundaries.py already does for knot
//  complements). Progress goes to stderr, so redirecting stdout to a file
//  works cleanly:
//
//      ./gen_knot_census_names pd_codes.csv > ../knot_census_names.csv
//      ./gen_knot_census_names --links links_pd_codes.csv > ../link_census_names.csv

#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <unordered_set>
#include <vector>

#include <census/census.h>
#include <triangulation/dim3.h>

#include "../knotbuilder.h"
#include "../linkcomplement.h"

namespace {
void usage(const char *progName, const std::string &error = std::string()) {
    if (!error.empty())
        std::cerr << error << "\n\n";

    std::cerr << "Usage:\n";
    std::cerr << "    " << progName << " <pd_codes.csv> [maxCrossings]\n";
    std::cerr << "    " << progName
              << " --links <links_pd_codes.csv> [maxCrossings]\n\n";
    std::cerr << "    <pd_codes.csv> : a CSV file with \"Name,PD Notation\" "
                 "rows (a header\n"
                 "                     row is skipped automatically)\n";
    std::cerr << "    --links <links_pd_codes.csv> : a CSV file with "
                 "\"Name,PD Notation,Genus-4D\"\n"
                 "                     rows (e.g. "
                 "links_4d_smooth_slice_genus_11_crossings_pd_codes.csv);\n"
                 "                     rows sharing a base name (differing "
                 "only by a trailing\n"
                 "                     {orientation} suffix) are deduped "
                 "to one representative,\n"
                 "                     since component orientation doesn't "
                 "change the complement\n";
    std::cerr << "    [maxCrossings] : optional cap on crossing number "
                 "(default 10 for knots,\n"
                 "                     unbounded for links)\n";
    exit(1);
}

// Writes a CSV field, quoting it if it contains a comma, quote or newline.
void writeField(std::ostream &out, const std::string &field) {
    if (field.find_first_of(",\"\n") == std::string::npos) {
        out << field;
        return;
    }
    out << '"';
    for (char c : field) {
        if (c == '"')
            out << '"';
        out << c;
    }
    out << '"';
}

void writeRow(const std::string &name, int crossings, const std::string &db,
             const std::string &censusName, const std::string &isoSig) {
    writeField(std::cout, name);
    std::cout << ',' << crossings << ',';
    writeField(std::cout, db);
    std::cout << ',';
    writeField(std::cout, censusName);
    std::cout << ',';
    writeField(std::cout, isoSig);
    std::cout << '\n';
}

// Splits a CSV line on plain commas -- safe here (not a general CSV
// parser) because neither input file ever quotes a field or embeds a
// comma inside one: PD notation uses semicolons as its internal separator
// in both the knot CSV ("[[1;5;2;4];...]") and the link CSV
// ("PD[X[4; 1; 3; 2]; ...]"), never a comma.
std::vector<std::string> splitCsvLine(const std::string &line) {
    std::vector<std::string> fields;
    size_t start = 0;
    while (true) {
        size_t comma = line.find(',', start);
        if (comma == std::string::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, comma - start));
        start = comma + 1;
    }
    return fields;
}

// The crossing number is just the leading digits of the Name column, up to
// the first '_' (KnotInfo's own convention, e.g. "10_139" -> 10). Returns
// -1 if there's no such leading digit run (e.g. a link name like "L6a1",
// which has no leading digit at all) -- callers fall back to
// crossingsFromPDCode() in that case.
int crossingsFromName(const std::string &name) {
    size_t underscore = name.find('_');
    if (underscore == std::string::npos)
        return -1;
    try {
        return std::stoi(name.substr(0, underscore));
    } catch (const std::exception &) {
        return -1;
    }
}

// Thistlethwaite link names have no leading digit run to parse a crossing
// count from (see crossingsFromName()'s doc comment) -- derive it from the
// PD code itself instead, the same pattern verifyslicegenus.cpp's
// loadInputCsv() uses for link input rows.
int crossingsFromPDCode(const knotbuilder::PDCode &pdcode) {
    return static_cast<int>(pdcode.size());
}

// Strips a trailing "{...}" orientation-variant suffix from a Thistlethwaite
// link name (e.g. "L11n459{0;1;0}" -> "L11n459"). Component orientation
// doesn't change the complement, so every row sharing a base name builds
// the exact same manifold -- processing only one representative per base
// name is both correct and (for names with many orientation variants)
// significantly cheaper.
std::string stripOrientationSuffix(const std::string &name) {
    size_t brace = name.find('{');
    if (brace == std::string::npos)
        return name;
    return name.substr(0, brace);
}

// Builds `pdcode`'s complement (as a Knot for a single-component diagram,
// a Link otherwise -- both just forward to EdgeComplement::buildComplement()
// unchanged, but going through the same types surfer.cpp/verifyslicegenus.cpp
// use at runtime keeps this generator honest about what it's actually
// mimicking), and returns both the complement's isoSig (the same key
// census::localCensusLookup()/identify::identify() query by -- written out
// so a miss row can be fed to an offline SnapPy identification pass, the
// same way custom_analysis/identify_boundaries.py already does for knot
// complements) and every regina::Census::lookup() hit against it.
std::pair<std::string, std::list<regina::CensusHit>>
censusHitsFor(const knotbuilder::PDCode &pdcode) {
    auto [tri, edges] = knotbuilder::buildLink(pdcode);
    Link link(tri, edges);
    regina::Triangulation<3> complement = link.buildComplement();
    std::string isoSig = complement.isoSig();
    return {isoSig, regina::Census::lookup(complement)};
}

// Tallies for one generation pass, threaded through processRow() and both
// runKnots()/runLinks() -- bundled into one struct (rather than several
// same-typed out-params) so counters can't be swapped by mistake at a call
// site.
struct Counters {
    long total = 0;      // rows read (oriented rows, for links)
    long considered = 0; // rows within maxCrossings
    long hits = 0;
    long misses = 0;
    long failures = 0;
};

// Processes one row (already crossing-filtered), building its complement
// and writing every census hit (or a single empty-hit row on a miss) to
// stdout. Shared by both the knot and link passes below.
void processRow(const std::string &name, int crossings,
                const std::string &pdStr, Counters &counters) {
    try {
        knotbuilder::PDCode pdcode = knotbuilder::parsePDCode(pdStr);
        auto [isoSig, censusHits] = censusHitsFor(pdcode);
        if (censusHits.empty()) {
            writeRow(name, crossings, "", "", isoSig);
            ++counters.misses;
        } else {
            for (const regina::CensusHit &hit : censusHits)
                writeRow(name, crossings, hit.db().desc(), hit.name(), isoSig);
            ++counters.hits;
        }
    } catch (const std::exception &e) {
        std::cerr << "\nFAILED: " << name << ": " << e.what() << "\n";
        ++counters.failures;
    }
}

void runKnots(std::istream &file, int maxCrossings) {
    std::string line;
    std::getline(file, line); // header row

    Counters c;

    while (std::getline(file, line)) {
        std::vector<std::string> fields = splitCsvLine(line);
        if (fields.size() < 2)
            continue;
        const std::string &name = fields[0];
        const std::string &pdStr = fields[1];
        ++c.total;

        int crossings = crossingsFromName(name);
        if (crossings < 0 || crossings > maxCrossings)
            continue;
        ++c.considered;

        std::cerr << "\r[*] " << c.considered << " considered (" << c.hits
                  << " hits, " << c.misses << " misses, " << c.failures
                  << " failures), on " << name << "               "
                  << std::flush;

        processRow(name, crossings, pdStr, c);
    }

    std::cerr << "\n=== " << c.considered << "/" << c.total
              << " knots considered (<= " << maxCrossings
              << " crossings): " << c.hits << " with a census hit, "
              << c.misses << " with none, " << c.failures
              << " failed to build ===\n";

    if (c.failures != 0)
        exit(1);
}

void runLinks(std::istream &file, int maxCrossings) {
    std::string line;
    std::getline(file, line); // header row

    Counters c;
    long deduped = 0;
    std::unordered_set<std::string> seenBaseNames;

    while (std::getline(file, line)) {
        std::vector<std::string> fields = splitCsvLine(line);
        if (fields.size() < 2)
            continue;
        const std::string &rawName = fields[0];
        const std::string &pdStr = fields[1];
        ++c.total;

        std::string baseName = stripOrientationSuffix(rawName);
        if (!seenBaseNames.insert(baseName).second)
            continue; // already processed this base link under another
                       // orientation variant
        ++deduped;

        knotbuilder::PDCode pdcode;
        try {
            pdcode = knotbuilder::parsePDCode(pdStr);
        } catch (const std::exception &e) {
            std::cerr << "\nFAILED: " << baseName << ": " << e.what() << "\n";
            ++c.failures;
            continue;
        }

        int crossings = crossingsFromPDCode(pdcode);
        if (crossings > maxCrossings)
            continue;
        ++c.considered;

        std::cerr << "\r[*] " << c.considered << " considered (" << c.hits
                  << " hits, " << c.misses << " misses, " << c.failures
                  << " failures), on " << baseName << "               "
                  << std::flush;

        processRow(baseName, crossings, pdStr, c);
    }

    std::cerr << "\n=== " << c.considered << "/" << deduped
              << " base links considered (<= " << maxCrossings
              << " crossings, from " << c.total << " oriented rows): "
              << c.hits << " with a census hit, " << c.misses
              << " with none, " << c.failures << " failed to build ===\n";

    if (c.failures != 0)
        exit(1);
}
} // namespace

int main(int argc, char *argv[]) {
    if (argc < 2 || argc > 4)
        usage(argv[0], "Please provide a PD code CSV file.");

    bool linksMode = std::string(argv[1]) == "--links";
    if (linksMode && argc < 3)
        usage(argv[0], "--links requires a CSV file path.");

    const char *csvPath = linksMode ? argv[2] : argv[1];
    std::ifstream file(csvPath);
    if (!file)
        usage(argv[0], std::string("Could not open file: ") + csvPath);

    int maxCrossingsArgIndex = linksMode ? 3 : 2;
    int defaultMaxCrossings = linksMode ? 999 : 10;
    int maxCrossings = (argc > maxCrossingsArgIndex)
                            ? std::stoi(argv[maxCrossingsArgIndex])
                            : defaultMaxCrossings;

    if (linksMode)
        runLinks(file, maxCrossings);
    else
        runKnots(file, maxCrossings);

    return 0;
}
