//
//  triangulate_knot_table.cpp
//
//  Standalone validation utility (not a CTest test -- see CMakeLists.txt):
//  sweeps knotbuilder::buildLink() over every PD code in a CSV file (one
//  "Name,PD Notation" row per line, as downloaded from KnotInfo,
//  https://knotinfo.org/) and checks that each produces a valid, closed S³
//  triangulation with a well-formed Link. Reports any failures individually
//  and a final pass/fail summary.
//
//  Usage:
//      triangulate_knot_table pd_codes.csv [limit]
//
//  With the full KnotInfo table (knots up to 13 crossings, ~12,467 rows)
//  this takes on the order of 10-15 minutes.

#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <triangulation/dim3.h>

#include "../knotbuilder.h"
#include "../knottedsurfaces.h"

namespace {
void usage(const char *progName, const std::string &error = std::string()) {
    if (!error.empty())
        std::cerr << error << "\n\n";

    std::cerr << "Usage:\n";
    std::cerr << "    " << progName << " <pd_codes.csv> [limit]\n\n";
    std::cerr << "    <pd_codes.csv> : a CSV file with \"Name,PD Notation\" "
                 "rows (a header\n"
                 "                     row is skipped automatically)\n";
    std::cerr << "    [limit]        : optional cap on the number of rows "
                 "to check\n";
    exit(1);
}
} // namespace

int main(int argc, char *argv[]) {
    if (argc < 2 || argc > 3)
        usage(argv[0], "Please provide a PD code CSV file.");

    std::ifstream file(argv[1]);
    if (!file)
        usage(argv[0], std::string("Could not open file: ") + argv[1]);

    long limit = (argc == 3) ? std::stol(argv[2])
                             : std::numeric_limits<long>::max();

    std::string line;
    std::getline(file, line); // header row

    long total = 0, passed = 0;
    long failedBuild = 0, failedValid = 0, failedClosed = 0, failedSphere = 0,
         failedLinkException = 0;

    auto start = std::chrono::steady_clock::now();

    while (total < limit && std::getline(file, line)) {
        size_t comma = line.find(',');
        if (comma == std::string::npos)
            continue;
        std::string name = line.substr(0, comma);
        std::string pdStr = line.substr(comma + 1);

        ++total;
        knotbuilder::PDCode pdcode = knotbuilder::parsePDCode(pdStr);
        if (pdcode.empty()) {
            std::cout << "EMPTY PD CODE: " << name << "\n";
            continue;
        }

        regina::Triangulation<3> tri;
        std::vector<const regina::Edge<3> *> edges;
        if (!knotbuilder::buildLink(tri, pdcode, edges)) {
            ++failedBuild;
            std::cout << "BUILD FAILED: " << name << "\n";
            continue;
        }

        bool ok = true;
        if (!tri.isValid()) {
            ++failedValid;
            ok = false;
        }
        if (!tri.isClosed()) {
            ++failedClosed;
            ok = false;
        }
        if (!tri.isSphere()) {
            ++failedSphere;
            ok = false;
        }

        try {
            Link link(tri, edges);
        } catch (const std::exception &e) {
            ++failedLinkException;
            std::cout << "LINK EXCEPTION: " << name << ": " << e.what()
                      << "\n";
            ok = false;
        }

        if (ok) {
            ++passed;
        } else {
            std::cout << "FAIL: " << name << " (" << pdcode.size()
                      << " crossings) valid=" << tri.isValid()
                      << " closed=" << tri.isClosed()
                      << " sphere=" << tri.isSphere() << "\n";
        }
    }

    double elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();

    std::cout << "\n=== " << passed << "/" << total << " passed in "
              << elapsed << "s ===\n";
    std::cout << "failedBuild=" << failedBuild << " failedValid=" << failedValid
              << " failedClosed=" << failedClosed
              << " failedSphere=" << failedSphere
              << " failedLinkException=" << failedLinkException << "\n";

    return (passed == total) ? 0 : 1;
}
