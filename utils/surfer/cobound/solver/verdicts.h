//
//  verdicts.h
//
//  The verdicts file: one line per table row (verify_genus_v2.csv's
//  columns, frozen).
//

#ifndef SURFER_COBOUND_VERDICTS_H
#define SURFER_COBOUND_VERDICTS_H

#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

#include "cobound/cobordisms/database.h"
#include "cobound/solver/literature.h"
#include "cobound/solver/solver.h"

/*! \file utils/surfer/cobound/solver/verdicts.h
 *  \brief The verdicts file (`--output` of verifyslicegenus; the atlas's
 *  verify_genus_v2.csv and every campaign's verify_shard.csv): its columns
 *  are a frozen format, which merge_cobordisms.py reads by name (knot,
 *  exhausted_depth, searched_faces, search_outcome).
 *
 *  `solve` writes every line from the solver's verdicts; a search run
 *  writes each searched row's search record and leaves its status and
 *  bounds as they were (plan divergence 10).
 */

namespace verdicts {

/// One row of the verdicts file: a name's current standing plus, when
/// something was derived, how.
struct OutputRow {
  std::string knot;
  int resolvedGenus = 0;
      /**< The established genus when `status` pins one; otherwise the best
           derived upper bound, or 0. */
  std::string status;
      // verified | improved | pinned | bounded | unresolved | skipped
  std::string witnessKind; // direct | cobordism | none
  std::string witnessPairSig;
  std::string viaKnot;
  int viaEdgeGenus = 0;
  std::string dependsOn;
  int literatureLo = 0;
  int literatureHi = 0;

  // Added alongside the interval solver.
  std::string derivedLo; // empty when no lower bound was derived
  std::string derivedHi; // empty when no upper bound was derived
  std::string witnessBasis; // constructive | literature-assisted | empty
  bool tubed = false;
      /**< Whether the witness surface was disconnected as found, with the
           recorded genus being its tubed genus. */

  // T5 bookkeeping: how hard this row was actually tried, so a later,
  // bigger-budget pass knows what is worth re-searching and what is
  // already settled. Without this a resume re-runs an identical search and
  // learns nothing.
  long long searchedFaces = 0; // the max_faces used; 0 means unbounded
  std::string searchOutcome;
      // exhausted | timeout | quiescent | stopped | empty (never searched)
  long long exhaustedDepth = -1;
      /**< The largest face cap at which EVERY root was enumerated to
           completion, or -1 if no round finished.

           This is the row's only exhaustive claim, and it is much stronger
           than `searchOutcome`: it says no cobordism exists for this object
           with at most this many added faces, rather than merely that we
           looked for a while. A timed-out run covers only a prefix of the
           root list, so it leaves this at -1 however long it ran.

           Kept as the best ever achieved for the row: a later, shallower
           run must not erase a deeper exhaustive result. */
};

/// The file's header line (frozen).
extern const char *const OUTPUT_HEADER;

std::string formatOutputRow(const OutputRow &r);

/// Loads a previously-written verdicts file, if present, keyed by knot name.
std::unordered_map<std::string, OutputRow>
loadOutputCsv(const std::filesystem::path &path);

/**
 * Rewrites the whole verdicts file from `outputRows` (atomicWrite), so a
 * crash mid-write never corrupts the previous, already-durable version.
 *
 * The file is a single unified table shared across every run ever pointed
 * at it, regardless of what any one run's targets cover. So this writes
 * every row in `outputRows`, not just this run's: `rows` only keeps this
 * run's own rows in their familiar order at the top of the file; every
 * other row follows, sorted by name for a stable, diffable order.
 */
void writeOutputCsv(const std::filesystem::path &path,
                    const std::vector<cobordismgraph::InputRow> &rows,
                    const std::unordered_map<std::string, OutputRow> &outputRows);

const char *statusName(cobordismgraph::Status s);

/**
 * Rebuilds `name`'s row from the solver's verdict, keeping whatever search
 * bookkeeping (searchedFaces/searchOutcome/exhaustedDepth) the row already
 * carried -- that records what we DID, which no amount of re-solving
 * changes. A bound's witness pair signature is read back through `reader`
 * when the bound holds only its offset.
 */
OutputRow rowFromVerdict(
    const std::string &name, const cobordismgraph::Verdict &v, const OutputRow *existing,
    const std::unordered_map<std::string, cobordismgraph::Bounds> &bounds,
    witnessstore::PairSigReader &reader);

} // namespace verdicts

#endif // SURFER_COBOUND_VERDICTS_H
