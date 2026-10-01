//
//  pdcode.h
//
//  PD codes as text: the one parser of a PD code into the crossings T is
//  built from, and the one formatter of a PD code into text.
//

#ifndef SURFER_DIAGRAMTRIANGULATION_PDCODE_H
#define SURFER_DIAGRAMTRIANGULATION_PDCODE_H

#include <array>
#include <sstream>
#include <string>
#include <vector>

/*! \file utils/surfer/diagramtriangulation/pdcode.h
 *  \brief PD (planar diagram) codes as text, both ways.
 *
 *  Reading: parsePDCode() is how T is built from a PD code: whatever the
 *  punctuation, the integers in fours, renumbered from 0 (a code that
 *  already holds a 0 is taken as 0-based). Reading a PD code into a
 *  regina::Link, labels as written, is linknaming/tables.h's
 *  (exactnaming::linkFromTablePD()).
 *
 *  Writing: formatPDCode() is the one formatter. The files that store a PD
 *  code spell it in two ways, both frozen until the atlas task picks one
 *  (PDSpelling).
 */

namespace knotbuilder {

/** A planar diagram code: one 4-tuple of strand labels per crossing. */
using PDCode = std::vector<std::array<int, 4>>;

/** Parses a PD (planar diagram) code string into a PDCode. */
PDCode parsePDCode(std::string pdcode_str);

/** How a stored PD code is spelt. */
enum class PDSpelling {
    /** `[[1;5;2;4];[3;1;4;6];...]`: every PD a cascade writes (kept.csv's and
     *  the `.rows.csv` sidecar's row PDs, nodes.csv, node_bounds.jsonl,
     *  certificates), and the knot table's but for spaces (from 11 crossings
     *  it writes `[[3; 1; 4; 26]; [1; ...`). */
    semicolons,
    /** `[[1,5,2,4],[3,1,4,6],...]`: farsidediagram's `pd=` field. */
    commas,
};

/** `pd` as text, spelt as `spelling` says. Labels are written as given. */
template <typename Int>
std::string formatPDCode(const std::vector<std::array<Int, 4>> &pd,
                         PDSpelling spelling) {
    const char sep = spelling == PDSpelling::semicolons ? ';' : ',';
    std::ostringstream o;
    o << '[';
    for (size_t i = 0; i < pd.size(); ++i) {
        if (i) o << sep;
        o << '[' << pd[i][0] << sep << pd[i][1] << sep << pd[i][2] << sep << pd[i][3]
          << ']';
    }
    o << ']';
    return o.str();
}

} // namespace knotbuilder

#endif // SURFER_DIAGRAMTRIANGULATION_PDCODE_H
