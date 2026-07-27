//
//  gen_knot_census_names.cpp
//
//  Standalone generation utility (not a CTest test -- see CMakeLists.txt):
//  for every knot in a KnotInfo-style PD code CSV ("Name,PD Notation" rows)
//  up to a given crossing number, builds the knot complement the same way
//  surfer.cpp does at runtime (knotbuilder::buildLink() -> Knot ->
//  buildComplement()) and records every regina::Census::lookup() hit
//  against it.
//
//  This produces the raw data behind utils/surfer/knot_census_names.csv,
//  which utils/surfer/tests/gen_knot_names_header.py then turns into
//  utils/surfer/rolfsentable.h (a census-name -> Rolfsen-name lookup table).
//
//  Usage:
//      gen_knot_census_names pd_codes.csv [maxCrossings]
//
//  CSV data rows go to stdout (one per census hit; a single row with empty
//  census_db/census_name fields for knots with no hit, e.g. torus knots,
//  which are not hyperbolic and so cannot appear in any of Regina's
//  hyperbolic-only census databases). Progress goes to stderr, so
//  redirecting stdout to a file works cleanly:
//
//      ./gen_knot_census_names pd_codes.csv > ../knot_census_names.csv

#include <fstream>
#include <iostream>
#include <list>
#include <string>

#include <census/census.h>
#include <triangulation/dim3.h>

#include "../knotbuilder.h"
#include "../linkcomplement.h"

namespace {
void usage(const char *progName, const std::string &error = std::string()) {
    if (!error.empty())
        std::cerr << error << "\n\n";

    std::cerr << "Usage:\n";
    std::cerr << "    " << progName << " <pd_codes.csv> [maxCrossings]\n\n";
    std::cerr << "    <pd_codes.csv> : a CSV file with \"Name,PD Notation\" "
                 "rows (a header\n"
                 "                     row is skipped automatically)\n";
    std::cerr << "    [maxCrossings] : optional cap on crossing number "
                 "(default 10)\n";
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
             const std::string &censusName) {
    writeField(std::cout, name);
    std::cout << ',' << crossings << ',';
    writeField(std::cout, db);
    std::cout << ',';
    writeField(std::cout, censusName);
    std::cout << '\n';
}

// The crossing number is just the leading digits of the Name column, up to
// the first '_' (KnotInfo's own convention, e.g. "10_139" -> 10).
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
} // namespace

int main(int argc, char *argv[]) {
    if (argc < 2 || argc > 3)
        usage(argv[0], "Please provide a PD code CSV file.");

    std::ifstream file(argv[1]);
    if (!file)
        usage(argv[0], std::string("Could not open file: ") + argv[1]);

    int maxCrossings = (argc == 3) ? std::stoi(argv[2]) : 10;

    std::string line;
    std::getline(file, line); // header row

    long total = 0, considered = 0, hits = 0, misses = 0, failures = 0;

    while (std::getline(file, line)) {
        size_t comma = line.find(',');
        if (comma == std::string::npos)
            continue;
        std::string name = line.substr(0, comma);
        std::string pdStr = line.substr(comma + 1);
        ++total;

        int crossings = crossingsFromName(name);
        if (crossings < 0 || crossings > maxCrossings)
            continue;
        ++considered;

        std::cerr << "\r[*] " << considered << " considered (" << hits
                  << " hits, " << misses << " misses, " << failures
                  << " failures), on " << name << "               "
                  << std::flush;

        try {
            knotbuilder::PDCode pdcode = knotbuilder::parsePDCode(pdStr);
            auto [tri, edges] = knotbuilder::buildLink(pdcode);
            Knot knot(tri, edges);
            regina::Triangulation<3> complement = knot.buildComplement();

            std::list<regina::CensusHit> censusHits =
                regina::Census::lookup(complement);
            if (censusHits.empty()) {
                writeRow(name, crossings, "", "");
                ++misses;
            } else {
                for (const regina::CensusHit &hit : censusHits)
                    writeRow(name, crossings, hit.db().desc(), hit.name());
                ++hits;
            }
        } catch (const std::exception &e) {
            std::cerr << "\nFAILED: " << name << ": " << e.what() << "\n";
            ++failures;
        }
    }

    std::cerr << "\n=== " << considered << "/" << total
              << " knots considered (<= " << maxCrossings
              << " crossings): " << hits << " with a census hit, " << misses
              << " with none, " << failures << " failed to build ===\n";

    return (failures == 0) ? 0 : 1;
}
