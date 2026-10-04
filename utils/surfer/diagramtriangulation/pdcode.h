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
 *  punctuation, the integers in fours (pdLabels()), renumbered from 0 (a
 *  code that already holds a 0 is taken as 0-based). Reading a PD code into
 *  a regina::Link, labels as written, is linknaming/tables.h's
 *  (linknaming::linkFromTablePD()).
 *
 *  Writing: formatPDCode() is the one formatter. The files that store a PD
 *  code spell it in two ways, both frozen until the atlas task picks one
 *  (PDSpelling). A PD given in another spelling (LinkInfo's
 *  `PD[X[4; 1; 3; 2]; ...]`, the knot table's spaced one from 11 crossings)
 *  is respelt as formatPDCode(pdLabels(text), spelling): the same integers
 *  in the same order, so the same T.
 */

namespace diagramtriangulation {

/** A planar diagram code: one 4-tuple of strand labels per crossing. */
using PDCode = std::vector<std::array<int, 4>>;

/** A PD code's text as its integers in fours, labels and crossing order
 *  exactly as written, whatever the punctuation. */
PDCode pdLabels(std::string pdcode_str);

/** Parses a PD (planar diagram) code string into a PDCode: pdLabels(),
 *  renumbered from 0. */
PDCode parsePDCode(std::string pdcode_str);

/** How a stored PD code is spelt. */
enum class PDSpelling {
    /** `[[1;5;2;4];[3;1;4;6];...]`: every PD a goal run writes (kept.csv's and
     *  the `.rows.csv` sidecar's `row_pd`s, the search's log.txt, nodes.csv,
     *  node_bounds.jsonl, certificates; its target's given PD respelt so),
     *  and the knot table's but for spaces (from 11 crossings it writes
     *  `[[3; 1; 4; 26]; [1; ...`). */
    semicolons,
    /** `[[1,5,2,4],[3,1,4,6],...]`: `cobound draw`'s `pd=` field. */
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

} // namespace diagramtriangulation

#endif // SURFER_DIAGRAMTRIANGULATION_PDCODE_H
