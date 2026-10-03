//
//  verifyslicegenus.cpp
//
//  Created by John Teague on 07/29/2026.
//

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <atomic>
#include <cassert>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include "cobound/bounds/axioms.h"
#include "cobound/bounds/searchjudge.h"
#include "cobound/cobordisms/cobordism.h"
#include "cobound/search/incoming.h"
#include "cobound/search/preconditions.h"
#include "cobound/solver/literature.h"
#include "cobound/solver/solver.h"
#include "linknaming/names.h"
#include "surfer/report/atomicwrite.h"
#include "surfer/report/progress.h"
#include "surfer/report/csvwriter.h"
#include "surfer/enumeration/submanifoldsearch.h"
#include "surfer/enumeration/surfacesearch.h"
#include "cobound/cobordisms/cobordismkey.h"
#include "cobound/cobordisms/database.h"
#include "cobound/cobordisms/pairsigner.h"
#include "cobound/cobordisms/pending.h"
#include "linknaming/tables.h"
#include "cobound/outgoing/outgoingnamer.h"
#include "diagramtriangulation/fromdiagram.h"
#include "linknaming/complement/linkcomplement.h"
#include "surfer/submanifold/linkingnumber.h"
#include "linknaming/census/censusnaming.h"
#include "linknaming/complement/complementcache.h"
#include "linknaming/complement/unlinknaming.h"
#include "cobound/search/search.h"
#include "cobound/search/searchreport.h"

using namespace cobordismgraph;

namespace {

// ─────────────────────────────────────────────────────────────────────────
// Fatal-bug detection: a found surface implying a genus BELOW an already-
// established true lower bound is a mathematical impossibility, not a data
// problem -- it means the search/embedding library itself computed
// something wrong (a bad genus, a false claim of embeddedness, etc.), and
// nothing else this run subsequently reports can be trusted either. So
// this halts the whole program, not just the current row.
// ─────────────────────────────────────────────────────────────────────────

// Set (once) from inside a search worker thread's onSurfaceBoundaryProcessed
// callback via flagFatalBug() below; only ever read/acted on back on the
// main thread, after the offending search's e.search() call has returned
// and any watchdog thread has joined -- never torn down the process from
// within a callback, which may be running on one of several concurrent
// worker threads.
std::atomic<bool> fatalBugDetected_{false};
std::mutex fatalBugMutex_;
std::string fatalBugMessage_;

// Records `message` as the reason for a fatal halt, if none has been
// recorded yet (first detector wins; later ones are redundant once the
// whole program is about to stop anyway). Safe to call from any thread.
void flagFatalBug(std::string message) {
  bool expected = false;
  if (!fatalBugDetected_.compare_exchange_strong(expected, true))
    return;
  std::lock_guard<std::mutex> lock(fatalBugMutex_);
  fatalBugMessage_ = std::move(message);
}

// Prints an unmissable banner and terminates the process (exit code 2,
// distinct from usage()'s exit code 1) if flagFatalBug() was ever called.
// Must only be called from the main thread with no search worker threads
// still running -- see fatalBugDetected_'s own comment. Call this after
// every e.search() this driver runs.
void haltIfFatalBugDetected() {
  if (!fatalBugDetected_.load())
    return;
  std::string message;
  {
    std::lock_guard<std::mutex> lock(fatalBugMutex_);
    message = fatalBugMessage_;
  }
  std::cerr
      << "\n\x1b[1;31m"
      << "################################################################"
         "################\n"
         "#####  FATAL: SEARCH LIBRARY BUG DETECTED -- HALTING NOW  #####\n"
         "################################################################"
         "################\x1b[0m\n\n"
      << message << "\n\n"
      << "This is a mathematical impossibility, not a data or literature "
         "issue -- it means\n"
         "surfacesearch.h/embeddedsubmanifold.h computed something wrong "
         "(a bad genus, a\n"
         "false claim of embeddedness, etc.). Nothing else this run could "
         "report from here\n"
         "on can be trusted, so the program is stopping immediately rather "
         "than continuing\n"
         "to the next knot/link. Whatever was already durably written to "
         "--output before\n"
         "this point is unaffected and safe to keep.\n";
  std::exit(2);
}

// ─────────────────────────────────────────────────────────────────────────
// CSV parsing helpers
// ─────────────────────────────────────────────────────────────────────────

// The witness store and the table readers live in witnessstore.h, shared
// with cascadesearch.
using witnessstore::loadNameTable;
using witnessstore::loadWitnesses;
using witnessstore::rewriteWitnessFile;


std::vector<InputRow> loadInputCsv(const std::filesystem::path &path) {
  std::vector<InputRow> rows;
  for (const exactnaming::TableRow &table : exactnaming::readTableRows(path)) {
    // A row to search must state its literature bounds: a malformed field
    // stops the run (it always did), never becomes a bound.
    const auto g4 = exactnaming::parseTableG4(table.g4);
    if (!g4)
      throw std::runtime_error("malformed 4-genus field '" + table.g4 + "' for " +
                               table.name + " in " + path.string());
    InputRow row;
    row.name = table.name;
    row.pdNotation = table.pd;
    row.lo = g4->first;
    row.hi = g4->second;
    const std::string &pd = row.pdNotation;
    // Crossing count is derived from the PD code itself (works uniformly
    // for both knot names like "13n_1109" and link names like "L10a1{0}",
    // which have no leading digit run to parse) rather than from `name`.
    row.crossings =
        static_cast<int>(knotbuilder::parsePDCode(pd).size());
    rows.push_back(std::move(row));
  }
  return rows;
}



// ─────────────────────────────────────────────────────────────────────────
// Output CSV schema and resumable I/O
// ─────────────────────────────────────────────────────────────────────────

// One row of --output: a name's current standing plus, when something was
// derived, how. Lives here rather than in cobordismgraph.h because it is
// purely this driver's file format -- the solver itself deals in
// cobordismgraph::Bounds/Verdict and has no opinion about CSV.
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
  long long searchedFaces = 0; // the --max-faces used; 0 means unbounded
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

// Which BoundaryCondition to search a row under; see
// rowsearch::conditionFor() for what each costs.
using rowsearch::BoundaryConditionMode;

constexpr const char *OUTPUT_HEADER =
    "knot,resolved_genus,status,witness_kind,witness_pairsig,via_knot,"
    "via_edge_genus,depends_on,literature_lo,literature_hi,"
    "derived_lo,derived_hi,witness_basis,tubed,searched_faces,search_outcome,"
    "exhausted_depth";

std::string formatOutputRow(const OutputRow &r) {
  std::ostringstream out;
  out << csvField(r.knot) << ',' << r.resolvedGenus << ',' << r.status << ','
      << r.witnessKind << ',' << csvField(r.witnessPairSig) << ','
      << csvField(r.viaKnot) << ',' << r.viaEdgeGenus << ','
      << csvField(r.dependsOn) << ',' << r.literatureLo << ','
      << r.literatureHi << ',' << r.derivedLo << ',' << r.derivedHi << ','
      << r.witnessBasis << ',' << (r.tubed ? "true" : "false") << ','
      << r.searchedFaces << ',' << r.searchOutcome << ',' << r.exhaustedDepth;
  return out.str();
}

// Loads a previously-written --output file, if present, keyed by knot name.
std::unordered_map<std::string, OutputRow>
loadOutputCsv(const std::filesystem::path &path) {
  std::unordered_map<std::string, OutputRow> result;
  std::ifstream in(path);
  if (!in)
    return result;

  std::string line;
  std::getline(in, line); // header
  while (std::getline(in, line)) {
    if (line.empty())
      continue;
    auto f = parseCsvLine(line);
    if (f.size() < 10)
      continue;
    OutputRow r;
    r.knot = f[0];
    try {
      r.resolvedGenus = std::stoi(f[1]);
    } catch (const std::exception &) {
      continue;
    }
    r.status = f[2];
    r.witnessKind = f[3];
    r.witnessPairSig = f[4];
    r.viaKnot = f[5];
    try {
      r.viaEdgeGenus = f[6].empty() ? 0 : std::stoi(f[6]);
    } catch (const std::exception &) {
      r.viaEdgeGenus = 0;
    }
    r.dependsOn = f[7];
    try {
      r.literatureLo = std::stoi(f[8]);
      r.literatureHi = std::stoi(f[9]);
    } catch (const std::exception &) {
      continue;
    }
    // Columns beyond the original ten are optional, so a file written by
    // an older build still loads (its rows simply carry no derived bounds
    // and no search bookkeeping, which is exactly the truth about them).
    if (f.size() > 10)
      r.derivedLo = f[10];
    if (f.size() > 11)
      r.derivedHi = f[11];
    if (f.size() > 12)
      r.witnessBasis = f[12];
    if (f.size() > 13)
      r.tubed = f[13] == "true";
    if (f.size() > 14) {
      try {
        r.searchedFaces = std::stoll(f[14]);
      } catch (const std::exception &) {
        r.searchedFaces = 0;
      }
    }
    if (f.size() > 15)
      r.searchOutcome = f[15];
    if (f.size() > 16) {
      try {
        r.exhaustedDepth = std::stoll(f[16]);
      } catch (const std::exception &) {
        r.exhaustedDepth = -1;
      }
    }
    result[r.knot] = std::move(r);
  }
  return result;
}




// Loads the observed-name -> classical-name table (see --name-aliases).
//
// A far side is named by identify(), from its complement alone, and that
// often lands on something no literature table knows: a bare isomorphism
// signature, a Christy census name ("L108019"), or a SnapPy census manifold
// name ("m129 : #3"). Such an edge bounds nothing. Where we have since
// PROVED what one of those is -- by a Pachner match against a complement
// built from a PD code, or by the peripheral test for a link -- this table
// records it.
//
// Deliberately a separate file rather than a rewrite of cobordisms.csv:
// that file records what the search observed, and must keep doing so. See
// applyNameAliases() for why the distinction has to survive into memory too.
std::unordered_map<std::string, std::string>
loadNameAliases(const std::filesystem::path &path) {
  std::unordered_map<std::string, std::string> aliases;
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("Cannot open name alias table: " + path.string());

  std::string line;
  std::getline(in, line); // header: observed,classical,basis
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    auto f = parseCsvLine(line);
    if (f.size() < 2 || f[0].empty() || f[1].empty())
      continue;
    // An anchor name is an AXIOM to the solver (seedAxioms matches on the
    // string), and identify() only ever emits one from a structural proof.
    // An alias must not be able to manufacture that proof by spelling.
    if (identify::isUnlinkName(f[1]))
      throw std::runtime_error("Name alias table maps '" + f[0] + "' to '" +
                               f[1] +
                               "': an unknot/unlink can only be established "
                               "by identify(), never by alias");
    aliases.emplace(f[0], f[1]);
  }
  return aliases;
}

// `name` as loadNameAliases() keys it: any " : #N" census-hit suffix
// stripped, matching linknames::name()'s own convention.
std::string aliasKey(const std::string &name) {
  return cobordismgraph::stripCensusSuffix(name);
}

// Resolves every witness's far side through the alias table, returning a
// SEPARATE vector for the solver to consume.
//
// Returning a copy rather than mutating in place is the whole point: the
// vector the search holds is what appendWitnesses() writes into
// cobordisms.csv, so aliasing it would quietly bake resolved names into the
// observation record -- exactly what keeping a separate alias table was
// meant to avoid.
//
// `otherCandidates` is re-derived rather than carried across: propagate()
// consumes the stored candidate list, so leaving it keyed to the old name
// would let a witness claim a far side of one name and the variants of
// another.
std::vector<cobordismgraph::Witness>
applyNameAliases(const std::vector<cobordismgraph::Witness> &witnesses,
                 const std::unordered_map<std::string, std::string> &aliases,
                 const cobordismgraph::NameTable &names, size_t &appliedOut) {
  std::vector<cobordismgraph::Witness> resolved = witnesses;
  size_t applied = 0;
  for (cobordismgraph::Witness &w : resolved) {
    if (w.other.empty())
      continue;
    auto it = aliases.find(aliasKey(w.other));
    if (it == aliases.end())
      continue;
    w.other = it->second;
    w.otherCandidates = names.candidates(w.other, w.otherComponents);
    ++applied;
  }
  appliedOut = applied;
  return resolved;
}

// One proved far-side identity, keyed on the witness rather than the name.
struct FarSideResolution {
  std::string boundaryComponent; // "0" or "1", as peripheral_slopes reports it
  std::string name;              // the ORIENTED name we have proved it to be
};

// One row of --far-side-exact: the far side redrawn from the witness's own
// pair signature, oriented by its surface and named with a proof
// (farsidename, exactnaming/).
struct ExactFarSide {
  std::string name;
  bool exact = false; // an identity (may receive a bound), not a description
  int components = 0; // curves drawn: must equal the witness's observed count
};

// Loads --far-side-exact: witness,name,exact,pinned,components,proof. A key
// seen with two different names is a contradiction and is dropped.
std::unordered_map<std::string, ExactFarSide>
loadFarSideExact(const std::filesystem::path &path, size_t &clashes) {
  std::unordered_map<std::string, ExactFarSide> out;
  std::unordered_set<std::string> clash;
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("Cannot open far-side exact names: " + path.string());
  std::string line;
  std::getline(in, line); // header
  while (std::getline(in, line)) {
    auto f = parseCsvLine(line);
    if (f.size() < 5 || f[0].empty() || f[1].empty())
      continue;
    ExactFarSide e{f[1], f[2] == "1", std::stoi(f[4])};
    auto [it, fresh] = out.try_emplace(f[0], e);
    if (!fresh && (it->second.name != e.name || it->second.exact != e.exact))
      clash.insert(f[0]);
  }
  for (const std::string &k : clash)
    out.erase(k);
  clashes = clash.size();
  return out;
}

// Applies --far-side-exact to the solver's copy of the witnesses, last, so it
// outranks aliases and resolutions: it is the far side drawn from this
// witness's own surface. Refused (and counted) when the drawing's curve count
// is not the count the search observed. The candidates are the name alone,
// or its proved alternatives -- never a base's variants.
std::vector<cobordismgraph::Witness>
applyFarSideExact(std::vector<cobordismgraph::Witness> witnesses,
                  const std::unordered_map<std::string, ExactFarSide> &exact,
                  size_t &applied, size_t &refused) {
  applied = refused = 0;
  for (cobordismgraph::Witness &w : witnesses) {
    if (w.kind != cobordismgraph::WitnessKind::cobordism)
      continue;
    if (w.pairSigKey.empty() && !w.pairSig.empty())
      w.pairSigKey = witnesskey::witnessKey(w.pairSig);
    auto it = exact.find(w.pairSigKey);
    if (it == exact.end())
      continue;
    if (it->second.components != w.otherComponents) {
      ++refused;
      continue;
    }
    w.other = it->second.name;
    w.otherCandidates = cobordismgraph::exactCandidates(it->second.name);
    w.farSideProved = true;
    w.farSideExact = it->second.exact;
    ++applied;
  }
  return witnesses;
}

// Loads the per-witness far-side resolution table (see
// --far-side-resolutions).
//
// WHY THIS EXISTS SEPARATELY FROM --name-aliases. An alias is keyed on the
// observed NAME, which is sound only where a name determines the object.
// For a knot it does: Gordon-Luecke makes the complement determine the knot
// up to mirroring, and g_4 is mirror-invariant. For a LINK it does not --
// one complement belongs to infinitely many links (Rolfsen twisting), and
// in our own data one observed census name is a dozen different links
// across different witnesses. A name-keyed row for such a far side would be
// wrong on most of the witnesses it matched.
//
// The pair signature does determine the far side, so link far sides are
// keyed on it (via witnesskey::witnessKey) plus which boundary component of
// that witness is meant.
std::unordered_map<std::string, std::vector<FarSideResolution>>
loadFarSideResolutions(const std::filesystem::path &path) {
  std::unordered_map<std::string, std::vector<FarSideResolution>> resolutions;
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("Cannot open far-side resolution table: " +
                             path.string());

  std::string line;
  std::getline(in, line); // header: witness,boundary_component,resolved_name,...
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    auto f = parseCsvLine(line);
    if (f.size() < 3 || f[0].empty() || f[2].empty())
      continue;
    resolutions[f[0]].push_back({f[1], f[2]});
  }
  return resolutions;
}

// Resolves far sides witness-by-witness, returning a SEPARATE vector for the
// same reason applyNameAliases() does: the search's own vector is what gets
// written to cobordisms.csv, so mutating it in place would bake an
// interpretation into the observation record.
//
// Applied AFTER applyNameAliases(), and strictly more specific than it: a
// resolution names one witness's far side, where an alias can only speak
// about a name. For a KNOT far side the two must agree -- a knot is
// determined by its complement (Gordon-Luecke), so an alias is an identity
// and a disagreement is a bug, and the run stops. For a LINK far side an
// alias can only say which COMPLEMENT was observed, and one complement is
// many links: on 2026-09-24, 21 witnesses whose far side was aliased from a
// census name (e.g. 9^2_55 -> L9n6) were proved per witness, with their own
// meridians, to be another link with the same complement (L9n8). There the
// resolution wins, and the count of such overrides is reported. This mirrors
// cobordism-atlas/tools/frontier.py's load_witnesses(); the two
// implementations are deliberately independent, and `frontier.py --check`
// is only a check while they stay that way.
std::vector<cobordismgraph::Witness> applyFarSideResolutions(
    std::vector<cobordismgraph::Witness> resolved,
    const std::vector<cobordismgraph::Witness> &observed,
    const std::unordered_map<std::string, std::vector<FarSideResolution>>
        &resolutions,
    const cobordismgraph::NameTable &names, size_t &appliedOut) {
  // Taken by value and rewritten in place. (A loaded witness no longer
  // carries its pair signature in memory, only pairSigKey.)
  size_t applied = 0;
  size_t linkAliasesOverridden = 0;
  std::vector<std::string> conflicts;

  for (cobordismgraph::Witness &w : resolved) {
    if (w.other.empty())
      continue;
    if (w.pairSigKey.empty() && !w.pairSig.empty())
      w.pairSigKey = witnesskey::witnessKey(w.pairSig);
    if (w.pairSigKey.empty())
      continue;
    auto it = resolutions.find(w.pairSigKey);
    if (it == resolutions.end())
      continue;

    // A witness has two boundary components and the table names one of them.
    // The component count is what says which: a resolution whose own
    // component count does not match this far side's observed curve count is
    // about the other side, not this one.
    const std::string *match = nullptr;
    for (const FarSideResolution &r : it->second) {
      // The name alone gives the count for a knot, an unlink or a tagged
      // link ("L9a47{0}"). A peripherally proved link arrives as its BASE
      // name -- the meridians pin the link, not its orientation -- and a
      // base name states no count, so ask the table: if it has registered
      // variants with the observed count, the resolution is about this side.
      const bool countFromName =
          cobordismgraph::componentsFromName(r.name) == w.otherComponents;
      // candidates() falls back to {name} itself for an unregistered base,
      // so a real table hit is one whose front is a different (tagged) name.
      // A composite K #_c L has L's components (a knot is summed INTO a
      // component), so it is L's variants that state the count.
      const std::optional<cobordismgraph::CompositeName> cp =
          cobordismgraph::compositeParts(r.name);
      const std::string countable = cp ? cp->link : r.name;
      const std::vector<std::string> variants =
          names.candidates(countable, w.otherComponents);
      const bool countFromTable =
          !variants.empty() && variants.front() != countable &&
          cobordismgraph::componentsFromName(variants.front()) ==
              w.otherComponents;
      if (countFromName || countFromTable) {
        match = &r.name;
        break;
      }
    }
    if (!match)
      continue;

    // The alias layer has already run, so w.other is the aliased name here.
    // Compare base names with any orientation tag stripped: a resolution
    // REFINES `L2a1` to `L2a1{0}`, which is the whole point, but must never
    // turn it into some other link. Only an ALIASED name is worth checking --
    // if no alias fired, w.other is still the raw observed name (an isoSig or
    // a census name), and disagreeing with that is not a contradiction but
    // the entire purpose of the resolution.
    const size_t idx = static_cast<size_t>(&w - resolved.data());
    const bool aliasFired =
        idx < observed.size() && observed[idx].other != w.other;
    const std::string aliasedBase = cobordismgraph::stripOrientationTag(w.other);
    const std::string resolvedBase = cobordismgraph::stripOrientationTag(*match);
    if (aliasFired && !aliasedBase.empty() && aliasedBase != resolvedBase) {
      if (w.otherComponents == 1)
        conflicts.push_back(w.other + " -> " + *match);
      else
        ++linkAliasesOverridden;
    }

    w.other = *match;
    w.otherCandidates = names.candidates(w.other, w.otherComponents);
    w.farSideProved = true;
    ++applied;
  }

  if (!conflicts.empty()) {
    std::ostringstream msg;
    msg << "far-side resolutions contradict name aliases on "
        << conflicts.size() << " KNOT far sides, e.g.";
    for (size_t i = 0; i < conflicts.size() && i < 3; ++i)
      msg << " [" << conflicts[i] << "]";
    msg << ". A knot is determined by its complement, so a resolution may "
           "refine a knot alias but never disagree with it; resolve by hand "
           "before solving.";
    throw std::runtime_error(msg.str());
  }
  if (linkAliasesOverridden)
    std::cerr << "[+] Far-side resolutions: " << linkAliasesOverridden
              << " link far sides named by a complement-level alias were "
                 "proved per witness to be another link with that complement; "
                 "the per-witness proof wins\n";

  appliedOut = applied;
  return resolved;
}


// Reads a witness's pair signature back from its line in the witness file
// (Witness::fileOffset), for the few places that print one: the database's
// memoized reader (the same few bounding witnesses are asked for after
// every row).
witnessstore::PairSigReader &pairSigReader() {
  static witnessstore::PairSigReader reader;
  return reader;
}

// Rewrites the whole --output file from `outputRows` via write-to-temp +
// atomic rename -- so a crash mid-write never corrupts the previous,
// already-durable version. Called once per knot/link processed (resolved,
// attempted-but-unresolved, skipped, or newly propagated).
//
// --output is a single unified table shared across every run ever pointed
// at it, regardless of what any one run's --input covers -- e.g. a
// knots-only run and a links-only run sharing the same --output file both
// resume from and contribute to the same pool of results, so a cobordism
// found this run between (say) a link and a not-yet-resolved knot from an
// earlier knots-only run resolves that knot too (see propagateGraph()),
// right here in the same file. So this writes every row currently in
// `outputRows`, not just ones belonging to this run's own --input: `rows`
// is used only to keep this run's own rows in their familiar
// crossing-count order at the top of the file; every other row (from a
// prior run's --input, sharing this --output) follows after, sorted by
// name for a stable, diffable order.
void writeOutputCsv(const std::filesystem::path &path,
                    const std::vector<InputRow> &rows,
                    const std::unordered_map<std::string, OutputRow>
                        &outputRows) {
  report::atomicWrite(path, [&](std::ostream &out) {
    out << OUTPUT_HEADER << "\n";

    std::unordered_set<std::string> written;
    written.reserve(outputRows.size());
    for (const auto &row : rows) {
      auto it = outputRows.find(row.name);
      if (it == outputRows.end())
        continue; // not yet processed
      out << formatOutputRow(it->second) << "\n";
      written.insert(row.name);
    }

    std::vector<std::string> others;
    others.reserve(outputRows.size());
    for (const auto &[name, unused] : outputRows)
      if (!written.contains(name))
        others.push_back(name);
    std::sort(others.begin(), others.end());
    for (const auto &name : others)
      out << formatOutputRow(outputRows.at(name)) << "\n";
  });
}

// ─────────────────────────────────────────────────────────────────────────
// Turning the solver's verdicts into --output rows
// ─────────────────────────────────────────────────────────────────────────

const char *statusName(cobordismgraph::Status s) {
  using S = cobordismgraph::Status;
  switch (s) {
  case S::verified:
    return "verified";
  case S::verifiedAssisted:
    return "verified-assisted";
  case S::improved:
    return "improved";
  case S::pinned:
    return "pinned";
  case S::bounded:
    return "bounded";
  case S::unresolved:
    return "unresolved";
  case S::contradiction:
    return "contradiction";
  }
  return "unresolved";
}

// Rebuilds `name`'s output row from the solver's current verdict, keeping
// whatever search bookkeeping (searchedFaces/searchOutcome) the row already
// carried -- that records what we DID, which no amount of re-solving
// changes, unlike everything else here.
OutputRow rowFromVerdict(
    const std::string &name, const cobordismgraph::Verdict &v,
    const OutputRow *existing,
    const std::unordered_map<std::string, cobordismgraph::Bounds> &bounds) {
  OutputRow out;
  out.knot = name;
  out.status = statusName(v.status);
  out.literatureLo = v.litLo;
  out.literatureHi = v.litHi;
  out.resolvedGenus = v.value;

  const auto &b = v.bounds;
  if (b.haveUpper()) {
    out.derivedHi = std::to_string(b.hi);
    out.witnessKind =
        b.kind == cobordismgraph::WitnessKind::direct ? "direct" : "cobordism";
    out.witnessPairSig = !b.pairSig.empty()
                             ? b.pairSig
                             : pairSigReader().at(b.pairSigOffset);
    out.viaKnot = b.viaName;
    out.viaEdgeGenus = b.viaGenus;
    out.dependsOn = cobordismgraph::buildDependsOn(b.viaName, bounds);
    out.witnessBasis = b.basis == cobordismgraph::Basis::constructive
                           ? "constructive"
                           : "literature-assisted";
    out.tubed = b.tubed;
  } else {
    out.witnessKind = "none";
  }
  if (b.haveLower())
    out.derivedLo = std::to_string(b.lo);

  if (existing) {
    out.searchedFaces = existing->searchedFaces;
    out.searchOutcome = existing->searchOutcome;
    out.exhaustedDepth = existing->exhaustedDepth;
    // A row that was skipped for crossing count and has still never been
    // searched keeps saying so, rather than being relabelled "unresolved"
    // as though we had tried.
    if (existing->status == "skipped" && !b.haveUpper() && !b.haveLower())
      out.status = "skipped";
  }
  return out;
}

// ─────────────────────────────────────────────────────────────────────────
// Boundary-component identification (--no-cone only)
// ─────────────────────────────────────────────────────────────────────────
// See cobordismgraph.h for BoundarySide/BoundarySplit/splitBoundary().

// ─────────────────────────────────────────────────────────────────────────
// usage()
// ─────────────────────────────────────────────────────────────────────────

void usage(const char *progName, const std::string &error = std::string()) {
  if (!error.empty())
    std::cerr << error << "\n\n";

  std::cerr
      << "Usage:\n    " << progName
      << " --output <csv> [ --input <csv> ] [ --max-crossings N ]\n"
         "    [ --cobordisms <csv> ] [ --solve-only ]\n"
         "    [ --max-faces N ] [ --harvest ] [ --harvest-quiescence S ]\n"
         "    [ --sweep-time-limit S ]\n"
         "    [ --knot-table <csv> ] [ --link-table <csv> ]\n"
         "    [ --name-aliases <csv> ]\n"
         "    [ --far-side-resolutions <csv> ] [ --far-side-exact <csv> ]\n"
         "    [ --link-classes <csv> ] [ --cascade-proofs <csv> ]\n"
         "    [ --knot-symmetry <csv> ] [ --sum-rules ]\n"
         "    [ --threads N ] [ --thicken-layers N ] [ --cone | --no-cone ]\n"
         "    [ --collar-layers N ] [ --iddfs-iterations N --iddfs-step D ]\n"
         "    [ --iddfs-start N ] [ --iddfs-final-threads N ]\n"
         "    [ --per-knot-time-limit S ] [ --surface-target N ]\n"
         "    [ --surface-log <path> ]\n"
         "    [ --surface-stats <path> ]\n"
         "    [ --no-census-updates ] [ --no-retriangulate-on-miss ]\n"
         "    [ --no-diagram-naming ] [ --exact-far-side-names ]\n"
         "    [ --retriangulate-height N ] [ --retriangulate-candidate-budget "
         "N ]\n"
         "    [ --retriangulate-time-budget S ]\n"
         "    [ --census-db <path> ] [ other surfer.cpp-style search-limit "
         "flags ]\n\n"
      << "Verifies the smooth 4D slice genus of every knot in --input by "
         "searching\n"
         "for a certifying surface (or chain of cobordisms) using "
         "surfacesearch.h,\n"
         "the same library surfer.cpp drives -- see that program's --pd "
         "path for\n"
         "the per-knot construction this replicates.\n\n";
  std::cerr
      << "    --output <csv> : Required. Resumable result file -- rewritten "
         "after\n"
         "                     every knot/link processed; rows from a prior "
         "run keep\n"
         "                     their search bookkeeping (every --input row is "
         "searched). A\n"
         "                     single unified table, not scoped to any one "
         "--input: rows\n"
         "                     already here from a different --input (e.g. "
         "a knots run's\n"
         "                     --output reused for a links run) are kept "
         "and can still be\n"
         "                     resolved by cobordism propagation this run "
         "even though\n"
         "                     they're never re-searched -- point multiple "
         "runs with\n"
         "                     different --input tables at the same "
         "--output to build one\n"
         "                     shared knot+link genus table over time.\n";
  std::cerr
      << "    --input <csv>  : The knot/link table (Name,PD Notation,"
         "Genus-4D) to\n"
         "                     search this run. Default:\n"
         "                     4d_smooth_slice_genus_13_crossings_pd_codes."
         "csv\n"
         "                     alongside the source tree. Never written "
         "back to.\n";
  std::cerr
      << "    --cobordisms <csv> : Durable record of every surface found "
         "(default:\n"
         "                     cobordisms.csv). A witness is a FACT (\"a "
         "surface with\n"
         "                     this boundary and this genus exists\"); an "
         "--output row\n"
         "                     is a CONCLUSION. Conclusions improve "
         "whenever the\n"
         "                     solver, the literature tables or the "
         "identification\n"
         "                     improve; facts don't. So searches append "
         "here and are\n"
         "                     never repeated, and --solve-only re-derives "
         "every\n"
         "                     conclusion from them in seconds.\n";
  std::cerr
      << "    --boundary-condition auto|connected|proper : How much boundary "
         "freedom a\n"
         "                     row's search allows (default: proper).\n"
         "                       auto/connected -- a single-component knot "
         "is searched\n"
         "                     under `connected` (at most one surface "
         "boundary curve per\n"
         "                     ambient boundary component), which prunes "
         "hard. Note it\n"
         "                     caps the FAR side's curve count too, so such "
         "a row can\n"
         "                     only ever find single-curve far sides -- it "
         "cannot\n"
         "                     discover a knot-to-LINK cobordism at all. "
         "Multi-component\n"
         "                     rows always fall back to `proper`, which "
         "they must.\n"
         "                       proper -- every row uses `proper`, so knot "
         "rows can\n"
         "                     find link far sides and feed the "
         "knot<->link chains\n"
         "                     this whole approach leans on. Substantially "
         "slower, since\n"
         "                     `connected`'s pruning is what makes knot "
         "rows cheap.\n"
      << "    --skip-drain-on-timeout : When --per-knot-time-limit or "
         "--harvest-\n"
         "                     quiescence stops a row, ALSO abandon the "
         "surfaces still\n"
         "                     queued for boundary identification "
         "(default: off, i.e.\n"
         "                     the drain runs to completion and only the "
         "SEARCH is\n"
         "                     bounded).\n"
         "                       Identification is the product: an "
         "unidentified surface\n"
         "                     says nothing about a slice genus. Cutting "
         "the drain\n"
         "                     short measured 1% identification on a run "
         "that produced\n"
         "                     1,216,027 qualifying surfaces -- 1.2 "
         "million boundaries\n"
         "                     discarded unexamined, which made that "
         "run's \"found\n"
         "                     nothing\" meaningless. Set this only when "
         "throughput\n"
         "                     genuinely matters more than knowing what "
         "was found.\n"
      << "    --solve-only   : Re-derive --output from --cobordisms and "
         "exit. No\n"
         "                     triangulations built, no searching. Run "
         "after any\n"
         "                     solver or literature-table change.\n";
  std::cerr
      << "    --max-faces N  : Hard cap on faces the search may add, "
         "making it\n"
         "                     terminate on its own. WITHOUT this the "
         "final search\n"
         "                     pass is unbounded and the only things that "
         "end a row\n"
         "                     are resolving it, --per-knot-time-limit, or "
         "Ctrl+C.\n"
         "                     With it, a row is finite and EXHAUSTIVE "
         "over exactly\n"
         "                     the surfaces the cap admits -- which is "
         "also what\n"
         "                     makes a result reproducible rather than "
         "\"whatever we\n"
         "                     found in ten minutes\".\n"
         "                       Counts faces ADDED to the collar seed, "
         "not the\n"
         "                     surface's total triangle count (matching "
         "--iddfs-step).\n"
         "                       Recorded per row in --output's "
         "searched_faces, so a\n"
         "                     later run at a LARGER cap re-searches only "
         "the rows\n"
         "                     that could yield something new -- "
         "progressive\n"
         "                     deepening across the whole table (default: "
         "unbounded).\n";
  std::cerr
      << "    --harvest      : Accepted and ignored: every search harvests, "
         "recording\n"
         "                     every distinct cobordism it finds rather than "
         "stopping\n"
         "                     once its row is settled. One expensive search "
         "then yields\n"
         "                     many edges, which bound OTHER rows. So a "
         "search needs a\n"
         "                     stopping rule (--max-faces, --surface-target, "
         "--harvest-\n"
         "                     quiescence or --per-knot-time-limit).\n";
  std::cerr
      << "    --surface-target N : Stop a row's SEARCH once N surfaces "
         "satisfying the\n"
         "                     boundary condition exist, instead of after a "
         "fixed wall\n"
         "                     time. Equalises how much of the root ordering "
         "a row\n"
         "                     covers, which is what the strength of a "
         "negative rests\n"
         "                     on -- wall clock only equalises that if every "
         "host\n"
         "                     explores roots at the same rate, and measured "
         "across\n"
         "                     machines they differ by 3.3x at identical "
         "thread-seconds.\n"
         "                     Pair it with --per-knot-time-limit as a "
         "backstop: a row\n"
         "                     whose space is smaller than N would otherwise "
         "run until\n"
         "                     it exhausts. The drain still runs to "
         "completion either\n"
         "                     way, and the row records which rule stopped "
         "it.\n"
      << "    --harvest-quiescence S : Stop a row once no NEW witness has "
         "been\n"
         "                     recorded for S seconds -- it has stopped "
         "teaching us\n"
         "                     anything, so the rest of its budget is "
         "better spent on\n"
         "                     a row we know nothing about (default: "
         "off).\n";
  std::cerr
      << "    --sweep-time-limit S : Wall-clock cap on the whole run, "
         "checked\n"
         "                     between rows. Everything is resumable via "
         "--output and\n"
         "                     --cobordisms, so a capped sweep loses no "
         "work\n"
         "                     (default: none).\n";
  std::cerr
      << "    --knot-symmetry <csv> : name,symmetry_type (KnotInfo spelling); "
         "makes every\n"
         "                     composite knot whose summands pair off into "
         "concordance\n"
         "                     inverses an anchor (K # m(K^r) is slice).\n"
      << "far-side-resolutions <csv> : per-WITNESS proved far-side "
         "names,\n"
         "                     keyed on sha1(pairsig)[:12] plus boundary "
         "component.\n"
         "                     Use for LINK far sides, where an observed name "
         "is not\n"
         "                     a function of the link and --name-aliases "
         "cannot be\n"
         "                     right on every witness it matches.\n"
      << "    --far-side-exact <csv> : per-WITNESS exact far-side names "
         "(farsidename):\n"
         "                     each far side redrawn from its own pair "
         "signature and\n"
         "                     named with a proof; outranks aliases and "
         "resolutions.\n"
      << "    --link-classes <csv> : the table's link classes (tableclasses): "
         "table\n"
         "                     names that are one oriented link, one graph "
         "node each.\n"
      << "    --cascade-proofs <csv> : certified cascadesearch proofs (the "
         "atlas's\n"
         "                     data/cascade_proofs.csv), each an upper bound "
         "on its\n"
         "                     target resting on the literature values it "
         "names.\n"
      << "    --sum-rules : bound sums along components and splits with link "
         "factors\n"
         "                     from their pieces.\n"
      << "    --name-aliases <csv> : observed -> classical far-side names, "
         "applied\n"
         "        when solving only; cobordisms.csv keeps what identify() "
         "observed.\n"
      << "    --knot-table <csv>, --link-table <csv> : Literature tables "
         "loaded for\n"
         "                     names and bounds ONLY, regardless of which "
         "one --input\n"
         "                     searches. Needed because identify() names a "
         "link by its\n"
         "                     complement, which cannot see component "
         "orientation: a\n"
         "                     far side identified as \"L6a3\" could be "
         "L6a3{0} (genus 2)\n"
         "                     or L6a3{1} (genus 0), so the solver has to "
         "know the full\n"
         "                     set of oriented variants and bound over all "
         "of them.\n";
  std::cerr << "    --max-crossings N : Skip rows whose crossing number "
               "exceeds N\n"
               "                     (default: 13).\n";
  std::cerr
      << "    --per-knot-time-limit S : Wall-clock cap (seconds) per knot's "
         "search;\n"
         "                     stops it (as if by requestStop()) if "
         "exceeded, and\n"
         "                     moves on. Without this, one hard unresolved "
         "knot can\n"
         "                     stall the whole sweep indefinitely (default: "
         "none).\n";
  std::cerr
      << "    --surface-log <path> : Mirrors surfer.cpp's -o: every surface "
         "found\n"
         "                     during the *current* knot's search, "
         "overwritten each\n"
         "                     knot. Debugging only -- lossy on crash "
         "(default: off).\n";
  std::cerr
      << "    --research-settled : Accepted and ignored: every row of --input "
         "is searched,\n"
         "                     settled or not. Choosing which rows to search "
         "is the\n"
         "                     worklist's job, and a row resumed from a "
         "complete frontier,\n"
         "                     or already at its surface target, returns at "
         "once.\n";
  std::cerr
      << "    --audit-linking : Validation only (slow). Compute every petal "
         "linking number\n"
         "                     twice -- by cochains (linkingnumber.h) and by "
         "drilling plus\n"
         "                     homology -- and halt after the row on any "
         "disagreement.\n";
  std::cerr
      << "    --resolve-unlinked : Also accept surfaces that meet themselves "
         "only at\n"
         "                     interior vertices whose trace (all petals "
         "together) is a certified unlink\n"
         "                     (pl_enumeration_draft §4.5): a perturbation "
         "near those\n"
         "                     vertices embeds each one with the same "
         "topology. Such witnesses\n"
         "                     carry resolved_vertices > 0. Changes what "
         "counts toward\n"
         "                     --surface-target, so set it for a whole "
         "campaign or not\n"
         "                     at all.\n"
         "    --no-resolve-unlinked : Accept embedded surfaces only. A search "
         "needs one of\n"
         "                     the two: resolve_unlinked has no default.\n";
  std::cerr
      << "    --self-intersection-census <path> : Measurement only. Appends "
         "one row per\n"
         "                     searched row: how many self-intersecting "
         "candidates of\n"
         "                     each kind the search met, and an audit of "
         "boundary-vertex\n"
         "                     petals on accepted surfaces. Slows the search "
         "(default:\n"
         "                     off).\n";
  std::cerr
      << "    --root-budget-start N : Ration each root's enumeration to N "
         "tryAdd\n"
         "                     attempts per pass, doubling the ration each "
         "pass until\n"
         "                     every root is exhausted. Without this a "
         "time-limited\n"
         "                     search only ever covers a prefix of the "
         "(shallow-first)\n"
         "                     root list, so a longer run is a longer "
         "prefix rather\n"
         "                     than a different sample (default: 0, off).\n";
  std::cerr
      << "    --root-budget-growth N : Factor the ration grows by between "
         "passes;\n"
         "                     >= 2 keeps total work within N/(N-1) of the "
         "final pass\n"
         "                     (default: 2).\n";
  std::cerr
      << "    --surface-stats <path> : Appends, per searched row, the "
         "distribution of\n"
         "                     the surfaces found -- face count against "
         "homeomorphism\n"
         "                     type (genus, punctures, tubed genus, closed "
         "components,\n"
         "                     connectedness). Counts at exactly --max-faces "
         "are the\n"
         "                     ones the cap is truncating, so this is how to "
         "tell\n"
         "                     whether raising it would buy anything. Cheap, "
         "cumulative,\n"
         "                     and safe to leave on (default: off).\n";
  std::cerr << "    --pair-sig-cache <dir> : Keep each row's pair-signature "
               "context\n"
               "                     (the ambient's isoSig and automorphisms: "
               "tens of seconds\n"
               "                     at 10 crossings) in <dir>, keyed by the "
               "exact triangulation\n"
               "                     and verified on loading, so a row searched "
               "again does not\n"
               "                     rebuild it (default: off).\n";
  std::cerr << "    --frontier-dir <dir> : Record each row's search frontier "
               "(searchfrontier.h)\n"
               "                     as <dir>/<row>.frontier: exactly how far "
               "every root got,\n"
               "                     cumulative over any run it resumed. Written "
               "only after the\n"
               "                     row's witnesses are on disk, and only when "
               "its accounting\n"
               "                     balanced and its drain ran to the end "
               "(default: off).\n";
  std::cerr << "    --resume-frontier-dir <dir> : Carry each row on from "
               "<dir>/<row>.frontier,\n"
               "                     if there is one and it is this search's "
               "(same row,\n"
               "                     triangulation, shape and acceptance): the "
               "covered prefix\n"
               "                     is skipped, not re-searched. May be the same "
               "directory as\n"
               "                     --frontier-dir (default: off).\n";
  std::cerr << "    --no-census-updates : Never write the census: neither a "
               "knot row's own\n"
               "                     complement after its search nor a Pachner "
               "search's hit\n"
               "                     (default: both on).\n";
  std::cerr << "    --exact-far-side-names : record a multi-curve far side under "
               "its exact,\n"
               "        ORIENTED name (exactnaming/, fast path), so witnesses are "
               "deduplicated\n"
               "        by the oriented far side; bears nothing until "
               "--far-side-exact.\n";
  std::cerr << "    --no-diagram-naming : Name far sides by drilling their "
               "complement, as before, instead of drawing them as diagrams "
               "(farsidenaming.h; default: draw, falling back to the "
               "complement only for what a diagram cannot name).\n";
  std::cerr
      << "    --no-retriangulate-on-miss : Disable the retriangulate-search "
         "identification\n"
         "                     fallback (default: on for this driver, "
         "opposite of\n"
         "                     surfer.cpp's default).\n";
  std::cerr
      << "    --retriangulate-height N : Pachner-move search depth per "
         "identification\n"
         "                     attempt (default: 2). retriangulate()'s "
         "candidate count\n"
         "                     grows roughly exponentially in this, so it's "
         "the single\n"
         "                     biggest lever on boundary-processing "
         "throughput -- try 1\n"
         "                     if that phase is the bottleneck, at the cost "
         "of catching\n"
         "                     fewer non-canonically-triangulated matches.\n";
  std::cerr
      << "    --retriangulate-candidate-budget N : Max candidate "
         "triangulations tried\n"
         "                     per identification attempt before giving up "
         "(default:\n"
         "                     8000).\n";
  std::cerr
      << "    --retriangulate-time-budget S : Wall-clock cap (seconds) per "
         "identification\n"
         "                     attempt before giving up (default: 20).\n\n";
  std::cerr << "    --threads N, --thicken-layers N, --cone/--no-cone, "
               "--collar-layers N,\n"
               "    --iddfs-iterations/--iddfs-step/--iddfs-start/"
               "--iddfs-final-threads,\n"
               "    --no-simplify, --pending-surface-cap, "
               "--petal-cache-limit,\n"
               "    --recognition-cache-limit, "
               "--boundary-signature-cache-limit,\n"
               "    --boundary-tally-cap, --census-db : as surfer.cpp; the "
               "defaults are\n"
               "    --thicken-layers 2 --collar-layers 2 --no-cone, as the "
               "cascade's.\n";
  std::cerr << "    -h, --help : Display this help\n";
  exit(1);
}

} // namespace

int main(int argc, char *argv[]) {
  std::optional<std::string> outputPath;
  std::string inputPath = "4d_smooth_slice_genus_13_crossings_pd_codes.csv";
  std::string cobordismsPath = "cobordisms.csv";
  // Both literature tables are loaded for metadata regardless of which one
  // --input searches; see the NameTable construction in the main body.
  std::string knotTablePath =
      "4d_smooth_slice_genus_13_crossings_pd_codes.csv";
  std::string linkTablePath =
      "links_4d_smooth_slice_genus_11_crossings_pd_codes.csv";
  // Optional: empty means "no alias table", which is the pre-existing
  // behaviour of taking every far-side name exactly as identify() left it.
  std::string nameAliasPath;
  // Optional: empty means "no per-witness resolution table". Separate from
  // the alias table because it is keyed on the witness, not the name -- see
  // loadFarSideResolutions().
  std::string farSideResolutionPath;
  std::string farSideExactPath;
  std::string linkClassesPath;
  // Optional: the atlas's data/cascade_proofs.csv (certified cascadesearch
  // proofs), each an upper-bound axiom on its target.
  std::string cascadeProofsPath;
  bool sumRules = false;
  // Optional: knot symmetry types (data/knot_symmetry.csv). Without it only
  // the two long-standing slice composites are anchors; with it, every
  // composite whose summands pair off into concordance inverses.
  std::string knotSymmetryPath;
  std::optional<long long> maxFaces;
  std::optional<double> harvestQuiescence;
  std::optional<double> sweepTimeLimit;
  bool solveOnly = false;
  // Identification is the product, so by default a timed-out row still drains
  // its boundary queue to completion; only the SEARCH is bounded.
  bool skipDrainOnTimeout = false;
  BoundaryConditionMode boundaryConditionMode =
      BoundaryConditionMode::proper;
  int maxCrossings = 13;
  std::optional<double> perKnotTimeLimit;
  // Stop the SEARCH once this many boundary-satisfying surfaces exist.
  //
  // --per-knot-time-limit equalises WALL TIME, which only equalises coverage
  // if every host explores roots at the same rate. Measured over 52 rows, they
  // do not: at an identical 3600 thread-seconds one host yielded a median
  // 972,655 qualifying surfaces per row and another 295,178 -- a 3.3x gap,
  // tight on both sides (+/-15%), and a property of the machine rather than of
  // the row. So a "found nothing" from the slower host was a materially weaker
  // claim than the same words from the faster one, and nothing in the record
  // said so. Targeting the surface count instead equalises the thing a
  // negative actually rests on: how much of the root ordering was covered.
  std::optional<long long> surfaceTarget;
  std::optional<std::string> surfaceLogPath;
  std::optional<std::string> surfaceStatsPath;
  std::optional<std::string> frontierDir;       // see --frontier-dir
  std::optional<std::string> pairSigCacheDir;   // see --pair-sig-cache
  std::optional<std::string> resumeFrontierDir; // see --resume-frontier-dir
  long long rootBudgetStart = 0;   // 0 = off, i.e. today's single-pass behaviour
  long long rootBudgetGrowth = 2;
  // Accept surfaces whose only self-intersections are unlinked (paper §4.5,
  // KnottedSurface::isResolvable()). No default (plan divergence 3): it
  // changes which surfaces count toward --surface-target and is part of the
  // frontier fingerprint, so a search run must state it either way.
  std::optional<bool> resolveUnlinked;
  // Measurement only; see SelfIntersectionCensus.
  std::optional<std::string> selfIntersectionCensusPath;
  // Audit trail: see --rejection-sample-log and sampleRejection below.
  std::optional<std::string> rejectionSampleLogPath;
  // Also run the Pachner search on multi-component far sides. Off by
  // default: a link's name bears no bound during a search, and the
  // per-witness far-side pipeline names every far side afterwards anyway.
  bool retriangulateLinksArg = false;
  bool diagramNaming = true;
  bool exactFarSideNames = false;
  // One explicit full rewrite of the witness file (the 12->13 column
  // migration); otherwise the file is only ever appended to.
  bool rewriteWitnesses = false;

  unsigned numThreads = std::thread::hardware_concurrency();
  if (numThreads == 0)
    numThreads = 1;

  // As the cascade's hops (HopShape::layers) and every campaign.
  int thickenLayers = 2;
  bool useCone = false;
  int collarLayers = 2;

  unsigned iddfsIterations = 0;
  long long iddfsStep = 0;
  std::optional<long long> iddfsStart;
  std::optional<unsigned> iddfsFinalThreads;

  SurfaceSearchLimits limits;
  limits.capturePairSig = true;
  size_t recognitionCacheLimitArg = identify::recognitionCacheLimit.load();
  std::string censusPath = SURFER_CENSUS_PATH;
  bool retriangulateOnMissArg = true;
  int retriangulateHeightArg = census::retriangulateHeight.load();
  size_t retriangulateCandidateBudgetArg =
      census::retriangulateCandidateBudget.load();
  long long retriangulateTimeBudgetArg =
      census::retriangulateTimeBudgetSeconds.load();

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      usage(argv[0]);
    } else if (arg == "--output") {
      if (i + 1 >= argc)
        usage(argv[0], "--output requires a value.");
      outputPath = argv[++i];
    } else if (arg == "--input") {
      if (i + 1 >= argc)
        usage(argv[0], "--input requires a value.");
      inputPath = argv[++i];
    } else if (arg == "--cobordisms") {
      if (i + 1 >= argc)
        usage(argv[0], "--cobordisms requires a value.");
      cobordismsPath = argv[++i];
    } else if (arg == "--knot-table") {
      if (i + 1 >= argc)
        usage(argv[0], "--knot-table requires a value.");
      knotTablePath = argv[++i];
    } else if (arg == "--link-table") {
      if (i + 1 >= argc)
        usage(argv[0], "--link-table requires a value.");
      linkTablePath = argv[++i];
    } else if (arg == "--knot-symmetry") {
      if (i + 1 >= argc)
        usage(argv[0], "--knot-symmetry requires a value.");
      knotSymmetryPath = argv[++i];
    } else if (arg == "--far-side-resolutions") {
      if (i + 1 >= argc)
        usage(argv[0], "--far-side-resolutions requires a value.");
      farSideResolutionPath = argv[++i];
    } else if (arg == "--far-side-exact") {
      if (i + 1 >= argc)
        usage(argv[0], "--far-side-exact requires a value.");
      farSideExactPath = argv[++i];
    } else if (arg == "--link-classes") {
      if (i + 1 >= argc)
        usage(argv[0], "--link-classes requires a value.");
      linkClassesPath = argv[++i];
    } else if (arg == "--cascade-proofs") {
      if (i + 1 >= argc)
        usage(argv[0], "--cascade-proofs requires a value.");
      cascadeProofsPath = argv[++i];
    } else if (arg == "--sum-rules") {
      sumRules = true;
    } else if (arg == "--name-aliases") {
      if (i + 1 >= argc)
        usage(argv[0], "--name-aliases requires a value.");
      nameAliasPath = argv[++i];
    } else if (arg == "--max-faces") {
      if (i + 1 >= argc)
        usage(argv[0], "--max-faces requires a value.");
      try {
        maxFaces = std::stoll(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--max-faces requires an integer value.");
      }
    } else if (arg == "--harvest") {
      // Every search harvests (plan divergence 5); still accepted for the
      // tools that pass it, until the config replaces the options.
    } else if (arg == "--harvest-quiescence") {
      if (i + 1 >= argc)
        usage(argv[0], "--harvest-quiescence requires a value.");
      try {
        harvestQuiescence = std::stod(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--harvest-quiescence requires a numeric value.");
      }
    } else if (arg == "--sweep-time-limit") {
      if (i + 1 >= argc)
        usage(argv[0], "--sweep-time-limit requires a value.");
      try {
        sweepTimeLimit = std::stod(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--sweep-time-limit requires a numeric value.");
      }
    } else if (arg == "--boundary-condition") {
      if (i + 1 >= argc)
        usage(argv[0], "--boundary-condition requires a value.");
      std::string mode = argv[++i];
      if (mode == "auto")
        boundaryConditionMode = BoundaryConditionMode::automatic;
      else if (mode == "connected")
        boundaryConditionMode = BoundaryConditionMode::connected;
      else if (mode == "proper")
        boundaryConditionMode = BoundaryConditionMode::proper;
      else
        usage(argv[0],
              "--boundary-condition must be auto, connected or proper.");
    } else if (arg == "--skip-drain-on-timeout") {
      skipDrainOnTimeout = true;
    } else if (arg == "--solve-only") {
      solveOnly = true;
    } else if (arg == "--max-crossings") {
      if (i + 1 >= argc)
        usage(argv[0], "--max-crossings requires a value.");
      try {
        maxCrossings = std::stoi(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--max-crossings requires an integer value.");
      }
    } else if (arg == "--per-knot-time-limit") {
      if (i + 1 >= argc)
        usage(argv[0], "--per-knot-time-limit requires a value.");
      try {
        perKnotTimeLimit = std::stod(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--per-knot-time-limit requires a numeric value.");
      }
    } else if (arg == "--surface-target") {
      if (i + 1 >= argc)
        usage(argv[0], "--surface-target requires a value.");
      surfaceTarget = std::stoll(argv[++i]);
    } else if (arg == "--surface-log") {
      if (i + 1 >= argc)
        usage(argv[0], "--surface-log requires a value.");
      surfaceLogPath = argv[++i];
    } else if (arg == "--pair-sig-cache") {
      if (i + 1 >= argc)
        usage(argv[0], "--pair-sig-cache requires a value.");
      pairSigCacheDir = argv[++i];
    } else if (arg == "--frontier-dir") {
      if (i + 1 >= argc)
        usage(argv[0], "--frontier-dir requires a value.");
      frontierDir = argv[++i];
    } else if (arg == "--resume-frontier-dir") {
      if (i + 1 >= argc)
        usage(argv[0], "--resume-frontier-dir requires a value.");
      resumeFrontierDir = argv[++i];
    } else if (arg == "--research-settled") {
      // Every row is searched (plan divergence 9); still accepted for the
      // tools that pass it, until the config replaces the options.
    } else if (arg == "--resolve-unlinked") {
      resolveUnlinked = true;
    } else if (arg == "--no-resolve-unlinked") {
      resolveUnlinked = false;
    } else if (arg == "--audit-linking") {
      linkingnumber::auditLinkingNumbers.store(true);
    } else if (arg == "--self-intersection-census") {
      if (i + 1 >= argc)
        usage(argv[0], "--self-intersection-census requires a value.");
      selfIntersectionCensusPath = argv[++i];
    } else if (arg == "--root-budget-start") {
      if (i + 1 >= argc)
        usage(argv[0], "--root-budget-start requires a value.");
      try {
        rootBudgetStart = std::stoll(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--root-budget-start requires a numeric value.");
      }
    } else if (arg == "--root-budget-growth") {
      if (i + 1 >= argc)
        usage(argv[0], "--root-budget-growth requires a value.");
      try {
        rootBudgetGrowth = std::stoll(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--root-budget-growth requires a numeric value.");
      }
      if (rootBudgetGrowth < 2)
        usage(argv[0], "--root-budget-growth must be at least 2.");
    } else if (arg == "--surface-stats") {
      if (i + 1 >= argc)
        usage(argv[0], "--surface-stats requires a value.");
      surfaceStatsPath = argv[++i];
    } else if (arg == "--no-census-updates") {
      census::censusUpdates.store(false, std::memory_order_relaxed);
    } else if (arg == "--no-retriangulate-on-miss") {
      retriangulateOnMissArg = false;
    } else if (arg == "--retriangulate-links") {
      retriangulateLinksArg = true;
    } else if (arg == "--no-diagram-naming") {
      diagramNaming = false;
    } else if (arg == "--exact-far-side-names") {
      exactFarSideNames = true;
    } else if (arg == "--rejection-sample-log") {
      if (i + 1 >= argc)
        usage(argv[0], "--rejection-sample-log requires a value.");
      rejectionSampleLogPath = argv[++i];
    } else if (arg == "--rewrite-witnesses") {
      rewriteWitnesses = true;
    } else if (arg == "--retriangulate-height") {
      if (i + 1 >= argc)
        usage(argv[0], "--retriangulate-height requires a value.");
      try {
        retriangulateHeightArg = std::stoi(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--retriangulate-height requires an integer value.");
      }
    } else if (arg == "--retriangulate-candidate-budget") {
      if (i + 1 >= argc)
        usage(argv[0], "--retriangulate-candidate-budget requires a value.");
      try {
        retriangulateCandidateBudgetArg =
            static_cast<size_t>(std::stoul(argv[++i]));
      } catch (const std::exception &) {
        usage(argv[0],
              "--retriangulate-candidate-budget requires an integer value.");
      }
    } else if (arg == "--retriangulate-time-budget") {
      if (i + 1 >= argc)
        usage(argv[0], "--retriangulate-time-budget requires a value.");
      try {
        retriangulateTimeBudgetArg = std::stoll(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0],
              "--retriangulate-time-budget requires an integer (seconds) "
              "value.");
      }
    } else if (arg == "--threads") {
      if (i + 1 >= argc)
        usage(argv[0], "--threads requires a value.");
      try {
        numThreads = static_cast<unsigned>(std::stoul(argv[++i]));
      } catch (const std::exception &) {
        usage(argv[0], "--threads requires an integer value.");
      }
    } else if (arg == "--thicken-layers") {
      if (i + 1 >= argc)
        usage(argv[0], "--thicken-layers requires a value.");
      try {
        thickenLayers = std::stoi(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--thicken-layers requires an integer value.");
      }
    } else if (arg == "--cone") {
      useCone = true;
    } else if (arg == "--no-cone") {
      useCone = false;
    } else if (arg == "--collar-layers") {
      if (i + 1 >= argc)
        usage(argv[0], "--collar-layers requires a value.");
      try {
        collarLayers = std::stoi(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--collar-layers requires an integer value.");
      }
    } else if (arg == "--iddfs-iterations") {
      if (i + 1 >= argc)
        usage(argv[0], "--iddfs-iterations requires a value.");
      try {
        iddfsIterations = static_cast<unsigned>(std::stoul(argv[++i]));
      } catch (const std::exception &) {
        usage(argv[0], "--iddfs-iterations requires an integer value.");
      }
    } else if (arg == "--iddfs-step") {
      if (i + 1 >= argc)
        usage(argv[0], "--iddfs-step requires a value.");
      try {
        iddfsStep = std::stoll(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--iddfs-step requires an integer value.");
      }
    } else if (arg == "--iddfs-start") {
      if (i + 1 >= argc)
        usage(argv[0], "--iddfs-start requires a value.");
      try {
        iddfsStart = std::stoll(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--iddfs-start requires an integer value.");
      }
    } else if (arg == "--iddfs-final-threads") {
      if (i + 1 >= argc)
        usage(argv[0], "--iddfs-final-threads requires a value.");
      try {
        iddfsFinalThreads = static_cast<unsigned>(std::stoul(argv[++i]));
      } catch (const std::exception &) {
        usage(argv[0], "--iddfs-final-threads requires an integer value.");
      }
    } else if (arg == "--no-simplify") {
      simplifyComplements = false;
    } else if (arg == "--pending-surface-cap") {
      if (i + 1 >= argc)
        usage(argv[0], "--pending-surface-cap requires a value.");
      try {
        limits.pendingSurfaceCap = std::stoull(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--pending-surface-cap requires an integer value.");
      }
    } else if (arg == "--petal-cache-limit") {
      if (i + 1 >= argc)
        usage(argv[0], "--petal-cache-limit requires a value.");
      try {
        limits.petalCacheLimit = std::stoull(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--petal-cache-limit requires an integer value.");
      }
    } else if (arg == "--recognition-cache-limit") {
      if (i + 1 >= argc)
        usage(argv[0], "--recognition-cache-limit requires a value.");
      try {
        recognitionCacheLimitArg = std::stoull(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--recognition-cache-limit requires an integer value.");
      }
    } else if (arg == "--boundary-signature-cache-limit") {
      if (i + 1 >= argc)
        usage(argv[0], "--boundary-signature-cache-limit requires a value.");
      try {
        limits.boundarySignatureCacheLimit = std::stoull(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0],
              "--boundary-signature-cache-limit requires an integer value.");
      }
    } else if (arg == "--boundary-tally-cap") {
      if (i + 1 >= argc)
        usage(argv[0], "--boundary-tally-cap requires a value.");
      try {
        limits.boundaryTallyCap = std::stoull(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--boundary-tally-cap requires an integer value.");
      }
    } else if (arg == "--census-db") {
      if (i + 1 >= argc)
        usage(argv[0], "--census-db requires a value.");
      censusPath = argv[++i];
    } else {
      usage(argv[0], "Unknown option: " + arg);
    }
  }

  if (!outputPath)
    usage(argv[0], "--output is required.");
  // The cobordism graph judges every search (divergence 6), reading each
  // find on its row as a goal run reads a stored row's: collared through
  // every layer, not coned.
  if (!solveOnly && (useCone || collarLayers != thickenLayers))
    usage(argv[0],
          "a search runs in a row collared through every layer, not coned (the "
          "cobordism graph reads each find on it): --collar-layers must equal "
          "--thicken-layers, and --cone is retired.");
  if (!solveOnly && !resolveUnlinked)
    usage(argv[0],
          "a search needs --resolve-unlinked or --no-resolve-unlinked: "
          "resolve_unlinked has no default (it changes which surfaces count "
          "toward --surface-target, and every frontier's fingerprint).");
  if (collarLayers < 0)
    usage(argv[0], "--collar-layers requires a value >= 0.");
  if (collarLayers > thickenLayers)
    usage(argv[0], "--collar-layers cannot exceed --thicken-layers.");
  if (iddfsIterations > 0 && iddfsStep <= 0)
    usage(argv[0], "--iddfs-iterations > 0 requires --iddfs-step > 0.");
  if (iddfsStart && *iddfsStart <= 0)
    usage(argv[0], "--iddfs-start requires a value > 0.");
  if (limits.pendingSurfaceCap == 0)
    usage(argv[0], "--pending-surface-cap requires a value > 0.");
  if (limits.petalCacheLimit == 0)
    usage(argv[0], "--petal-cache-limit requires a value > 0.");
  if (recognitionCacheLimitArg == 0)
    usage(argv[0], "--recognition-cache-limit requires a value > 0.");
  if (limits.boundarySignatureCacheLimit == 0)
    usage(argv[0], "--boundary-signature-cache-limit requires a value > 0.");
  if (limits.boundaryTallyCap == 0)
    usage(argv[0], "--boundary-tally-cap requires a value > 0.");
  if (maxFaces && *maxFaces <= 0)
    usage(argv[0], "--max-faces requires a value > 0.");
  if (harvestQuiescence && *harvestQuiescence <= 0)
    usage(argv[0], "--harvest-quiescence requires a value > 0.");
  if (sweepTimeLimit && *sweepTimeLimit <= 0)
    usage(argv[0], "--sweep-time-limit requires a value > 0.");
  if (surfaceTarget && *surfaceTarget <= 0)
    usage(argv[0], "--surface-target requires a value > 0.");
  if (surfaceTarget && harvestQuiescence)
    usage(argv[0],
          "--harvest-quiescence cannot be combined with --surface-target: "
          "quiescence stops a row by how many NEW witnesses it records, so "
          "rows would cover different amounts of the root ordering -- the "
          "very thing --surface-target exists to equalise -- and a row whose "
          "witnesses were being lost would simply look quiescent.");
  if (!solveOnly && !maxFaces && !harvestQuiescence && !perKnotTimeLimit &&
      !surfaceTarget)
    usage(argv[0],
          "a search needs a stopping rule: without --max-faces, "
          "--surface-target, --harvest-quiescence or --per-knot-time-limit, "
          "nothing ends it (every search harvests, so resolving the row does "
          "not) and the sweep will never reach its second row.");
  if (retriangulateHeightArg < 0)
    usage(argv[0], "--retriangulate-height requires a value >= 0.");
  if (retriangulateCandidateBudgetArg == 0)
    usage(argv[0], "--retriangulate-candidate-budget requires a value > 0.");
  if (retriangulateTimeBudgetArg <= 0)
    usage(argv[0], "--retriangulate-time-budget requires a value > 0.");
  identify::recognitionCacheLimit.store(recognitionCacheLimitArg,
                                        std::memory_order_relaxed);
  census::retriangulateOnMiss.store(retriangulateOnMissArg,
                                    std::memory_order_relaxed);
  census::retriangulateHeight.store(retriangulateHeightArg,
                                    std::memory_order_relaxed);
  census::retriangulateCandidateBudget.store(retriangulateCandidateBudgetArg,
                                             std::memory_order_relaxed);
  census::retriangulateTimeBudgetSeconds.store(retriangulateTimeBudgetArg,
                                               std::memory_order_relaxed);
  census::retriangulateLinks.store(retriangulateLinksArg,
                                   std::memory_order_relaxed);
  if (const char *perturb = std::getenv("SURFER_TEST_PERTURB_NAMES");
      perturb && *perturb && std::string(perturb) != "0") {
    identify::perturbNamesForTesting.store(true);
    std::cerr << "[!] SURFER_TEST_PERTURB_NAMES: every identified name is "
                 "perturbed (test mode)\n";
  }
  // Only a multi-curve component's curve COUNT is ever used here, so naming
  // each of its curves separately is work nobody reads.
  limits.nameLinkCurves = false;
  bool censusLoaded = census::setCensusPath(censusPath);

  std::optional<std::ofstream> rejectionSampleLog;
  if (rejectionSampleLogPath) {
    std::error_code ec;
    const bool fresh = !std::filesystem::exists(*rejectionSampleLogPath) ||
                       std::filesystem::file_size(*rejectionSampleLogPath, ec) == 0;
    rejectionSampleLog.emplace(*rejectionSampleLogPath, std::ios::app);
    if (!*rejectionSampleLog)
      usage(argv[0], "--rejection-sample-log: cannot open " +
                         *rejectionSampleLogPath);
    if (fresh)
      *rejectionSampleLog
          << "row,reason,tubed_genus,connected,boundary,pairsig\n";
  }

  std::cout << "------ verifyslicegenus \U0001F30A ------\n\n";
  std::cout << (censusLoaded ? "[+] census: loaded from "
                            : "[+] census: not found at ")
            << censusPath << (censusLoaded ? "\n\n" : ", skipping\n\n");

  std::vector<InputRow> rows;
  try {
    rows = loadInputCsv(inputPath);
  } catch (const std::exception &e) {
    usage(argv[0], e.what());
  }
  std::cout << "[+] Loaded " << rows.size() << " knots from " << inputPath
            << "\n";

  // Lets the incidental-observation path (onSurfaceBoundaryProcessed,
  // below) look up an incidentally-discovered name's own literature
  // bounds in O(1) -- built from the full `rows`, not just `pending`
  // (rows within maxCrossings), so a name outside this run's own
  // maxCrossings cutoff still isn't recorded against, matching that same
  // cutoff's intent.
  std::unordered_map<std::string, const InputRow *> nameToRow;
  nameToRow.reserve(rows.size());
  for (const auto &row : rows)
    nameToRow[row.name] = &row;

  // Literature metadata for every name that could ever appear in the
  // graph, not just this run's --input. A knots run needs the links table
  // to expand an orientation-blind far-side name like "L6a3" into the
  // oriented variants it might be (see NameTable::candidates), and a links
  // run needs the knots table for the same reason in reverse -- so both
  // are always loaded, regardless of which one is being searched.
  cobordismgraph::NameTable names;
  size_t metadataRows = 0;
  for (const auto &row : rows)
    names.addLiterature(row.name, row.lo, row.hi);
  for (const std::string &table : {knotTablePath, linkTablePath}) {
    if (table.empty() || table == inputPath)
      continue;
    try {
      metadataRows += loadNameTable(table, names);
    } catch (const std::exception &e) {
      std::cerr << "[!] could not load name table " << table << ": "
                << e.what() << " (continuing without it)\n";
    }
  }
  std::cout << "[+] Name table: " << names.size() << " names ("
            << metadataRows << " from tables other than --input)\n";

  // Exact diagram signatures of every table knot and link, for naming far
  // sides from their drawings (farsidenaming.h). Only a search needs them.
  std::optional<farside::SignatureTable> signatureTable;
  std::optional<exactnaming::ExactTables> exactTables;
  if (exactFarSideNames && (!diagramNaming || solveOnly))
    std::cerr << "[!] --exact-far-side-names needs diagram naming and a search; "
                 "ignored\n";
  if (exactFarSideNames && diagramNaming && !solveOnly) {
    const auto t0 = std::chrono::steady_clock::now();
    try {
      exactTables = exactnaming::ExactTables::load(knotTablePath, linkTablePath,
                                                   knotSymmetryPath);
      std::cout << "[+] exact far-side names: " << exactTables->size()
                << " table entries ("
                << std::chrono::duration_cast<std::chrono::milliseconds>(
                       std::chrono::steady_clock::now() - t0)
                       .count()
                << " ms)\n";
    } catch (const std::exception &e) {
      std::cerr << "[!] exact far-side names off: " << e.what() << "\n";
    }
  }
  // The tables a search's own cobordism graph names its links by
  // (divergence 6): the exact tables (shared with --exact-far-side-names's,
  // when loaded) and a namer over them, as a goal run names its nodes.
  std::optional<exactnaming::ExactTables> graphTablesOwn;
  const exactnaming::ExactTables *graphTables = exactTables ? &*exactTables : nullptr;
  std::optional<exactnaming::ExactNamer> graphNamer;
  if (!solveOnly) {
    try {
      if (!graphTables)
        graphTables = &graphTablesOwn.emplace(
            exactnaming::ExactTables::load(knotTablePath, linkTablePath, knotSymmetryPath));
      graphNamer.emplace(*graphTables, cascade::NodeAxioms::namerLimits());
    } catch (const std::exception &e) {
      std::cerr << "[!] the cobordism graph cannot load the tables: " << e.what() << "\n";
      return 1;
    }
  }
  if (diagramNaming && !solveOnly) {
    const auto t0 = std::chrono::steady_clock::now();
    try {
      signatureTable = farside::SignatureTable::fromTables(knotTablePath,
                                                           linkTablePath);
      std::cout << "[+] diagram naming: " << signatureTable->knots()
                << " knot and " << signatureTable->links()
                << " link diagram signatures ("
                << std::chrono::duration_cast<std::chrono::milliseconds>(
                       std::chrono::steady_clock::now() - t0)
                       .count()
                << " ms)\n";
    } catch (const std::exception &ex) {
      std::cerr << "[!] diagram naming off: could not load the tables' "
                   "signatures ("
                << ex.what() << ")\n";
    }
  }

  std::unordered_map<std::string, OutputRow> outputRows =
      loadOutputCsv(*outputPath);

  if (rewriteWitnesses && std::filesystem::exists(cobordismsPath)) {
    try {
      rewriteWitnessFile(cobordismsPath);
    } catch (const std::exception &e) {
      std::cerr << "[!] --rewrite-witnesses: " << e.what() << "\n";
      return 2;
    }
  }
  pairSigReader().setPath(cobordismsPath);
  // Both per-witness tables are keyed on the pair signature's key, which is
  // hashed at load only when one of them will be looked up.
  // Every witness the run knows: the file's, and each one its searches
  // record (cobordisms/pending). `witnesses` is what the solver reads.
  cascade::RecordedWitnesses recorded(
      cobordismsPath,
      loadWitnesses(cobordismsPath,
                    !farSideResolutionPath.empty() || !farSideExactPath.empty()));
  const std::vector<cobordismgraph::Witness> &witnesses = recorded.all();
  std::cout << "[+] Resuming with " << witnesses.size()
            << " previously-recorded witnesses from " << cobordismsPath
            << "\n";

  std::unordered_map<std::string, std::string> nameAliases;
  if (!nameAliasPath.empty()) {
    try {
      nameAliases = loadNameAliases(nameAliasPath);
      std::cout << "[+] Name aliases: " << nameAliases.size()
                << " observed names resolved to classical ones from "
                << nameAliasPath << "\n";
    } catch (const std::exception &e) {
      std::cerr << "[!] could not load name aliases " << nameAliasPath << ": "
                << e.what() << " (continuing without them)\n";
    }
  }

  // The witness list the SOLVER sees, which is not the one the search holds:
  // aliases resolve a far side to what we have since proved it to be, while
  // `witnesses` keeps what identify() observed and is what gets written back
  // to cobordisms.csv. Rebuilt on each solve because the search appends to
  // `witnesses` as it goes; the copy costs a few MB against a search measured
  // in minutes.
  if (!knotSymmetryPath.empty()) {
    exactnaming::SymmetryTable types;
    try {
      types = exactnaming::readSymmetryTable(knotSymmetryPath);
    } catch (const std::exception &) {
      std::cerr << "[!] could not open knot symmetry table " << knotSymmetryPath
                << "\n";
      return 1;
    }
    for (const auto &[knot, type] : types)
      names.setSymmetry(knot, type);
    std::cout << "[+] Knot symmetry: " << types.size() << " types from "
              << knotSymmetryPath << "\n";
  }

  std::unordered_map<std::string, std::vector<FarSideResolution>>
      farSideResolutions;
  if (!farSideResolutionPath.empty()) {
    try {
      farSideResolutions = loadFarSideResolutions(farSideResolutionPath);
      std::cout << "[+] Far-side resolutions: " << farSideResolutions.size()
                << " witnesses with a proved far side from "
                << farSideResolutionPath << "\n";
    } catch (const std::exception &e) {
      std::cerr << "[!] could not load far-side resolutions "
                << farSideResolutionPath << ": " << e.what()
                << " (continuing without them)\n";
    }
  }

  std::unordered_map<std::string, ExactFarSide> farSideExact;
  if (!farSideExactPath.empty()) {
    size_t clashes = 0;
    try {
      farSideExact = loadFarSideExact(farSideExactPath, clashes);
    } catch (const std::exception &e) {
      std::cerr << "[!] could not load far-side exact names " << farSideExactPath
                << ": " << e.what() << "\n";
      return 1;
    }
    std::cout << "[+] Far-side exact names: " << farSideExact.size()
              << " witnesses from " << farSideExactPath;
    if (clashes)
      std::cout << " (" << clashes << " with two names dropped)";
    std::cout << "\n";
  }

  // --link-classes: table names that are one oriented link up to mirror and
  // global reversal (tableclasses: a version of one table diagram is the
  // other, or an isometry of the complements carries meridians to meridians
  // with a uniform orientation sign). Every member is read as its class's
  // canonical name, so a class is one node; tableclasses refuses a class
  // whose members' literature values differ.
  std::unordered_map<std::string, std::string> linkClasses;
  if (!linkClassesPath.empty()) {
    std::ifstream in(linkClassesPath);
    if (!in) {
      std::cerr << "[!] could not open link classes " << linkClassesPath << "\n";
      return 1;
    }
    std::string line;
    std::getline(in, line); // name,canonical,proof
    while (std::getline(in, line)) {
      auto f = parseCsvLine(line);
      if (f.size() >= 2 && !f[0].empty() && !f[1].empty())
        linkClasses[f[0]] = f[1];
    }
    std::cout << "[+] Link classes: " << linkClasses.size()
              << " table names read as their class's canonical name, from "
              << linkClassesPath << "\n";
  }
  auto classOf = [&linkClasses](const std::string &name) -> const std::string & {
    auto it = linkClasses.find(name);
    return it == linkClasses.end() ? name : it->second;
  };
  names.setSumRules(sumRules);
  if (sumRules)
    std::cout << "[+] Sum rules: sums along components and splits with link "
                 "factors are bounded from their pieces\n";

  // --cascade-proofs: target,goal,bound,basis,support,verdict,source,...
  // Only CERTIFIED proofs of the connected goal bound g4; names (the target
  // and every literature leaf) read through the link classes, as witnesses'.
  std::vector<cobordismgraph::ExternalProof> externalProofs;
  if (!cascadeProofsPath.empty()) {
    std::ifstream in(cascadeProofsPath);
    if (!in) {
      std::cerr << "[!] could not open cascade proofs " << cascadeProofsPath << "\n";
      return 1;
    }
    std::string line;
    std::getline(in, line);
    const std::vector<std::string> head = parseCsvLine(line);
    auto col = [&head](const std::string &c) {
      return static_cast<size_t>(std::find(head.begin(), head.end(), c) - head.begin());
    };
    const size_t cT = col("target"), cG = col("goal"), cB = col("bound"),
                 cS = col("support"), cV = col("verdict"), cSrc = col("source");
    if (std::max({cT, cG, cB, cS, cV, cSrc}) >= head.size()) {
      std::cerr << "[!] " << cascadeProofsPath
                << ": needs target,goal,bound,support,verdict,source columns\n";
      return 1;
    }
    size_t skipped = 0;
    while (std::getline(in, line)) {
      if (line.empty())
        continue;
      const auto f = parseCsvLine(line);
      if (f.size() < head.size() || f[cV] != "CERTIFIED" || f[cG] != "connected") {
        ++skipped;
        continue;
      }
      cobordismgraph::ExternalProof p;
      p.name = classOf(f[cT]);
      p.genus = std::stoi(f[cB]);
      std::istringstream support(f[cS]);
      for (std::string s; std::getline(support, s, ';');)
        if (!s.empty())
          p.support.push_back(classOf(s));
      p.source = f[cSrc];
      externalProofs.push_back(std::move(p));
    }
    std::cout << "[+] Cascade proofs: " << externalProofs.size()
              << " certified bounds on connected g4 from " << cascadeProofsPath;
    if (skipped)
      std::cout << " (" << skipped << " others skipped)";
    std::cout << "\n";
  }

  size_t aliasesApplied = 0;
  size_t resolutionsApplied = 0;
  size_t exactApplied = 0, exactRefused = 0;
  auto solverWitnesses = [&]() -> std::vector<cobordismgraph::Witness> {
    std::vector<cobordismgraph::Witness> out =
        nameAliases.empty()
            ? witnesses
            : applyNameAliases(witnesses, nameAliases, names, aliasesApplied);
    if (!farSideResolutions.empty())
      out = applyFarSideResolutions(std::move(out), witnesses,
                                    farSideResolutions, names,
                                    resolutionsApplied);
    if (!farSideExact.empty())
      out = applyFarSideExact(std::move(out), farSideExact, exactApplied, exactRefused);
    // Last: whole names only. A name inside a sum or split is a piece,
    // bounded by its literature value, which is the same across a class.
    if (!linkClasses.empty())
      for (cobordismgraph::Witness &w : out) {
        w.subject = classOf(w.subject);
        w.other = classOf(w.other);
        for (std::string &c : w.otherCandidates)
          c = classOf(c);
      }
    return out;
  };

  // Every conclusion is re-derived from the witness set on every run, so a
  // solver fix or a literature-table update takes effect on rows that were
  // searched long ago without re-searching any of them.
  // A resolution/alias contradiction is a data error, not a bug, so it exits
  // rather than aborting; caught here because the table is static -- if it
  // contradicts at all, it does so on this first solve.
  std::vector<cobordismgraph::Witness> initialWitnesses;
  try {
    initialWitnesses = solverWitnesses();
  } catch (const std::exception &e) {
    std::cerr << "[!] " << e.what() << "\n";
    return 1;
  }
  auto bounds = cobordismgraph::propagate(initialWitnesses, names, externalProofs);
  // Release it now: it is a full witness set, pair signatures included, and
  // held for the rest of the run it raised peak memory by ~45% -- enough for
  // a full-master --solve-only to be OOM-killed on yoga (2026-09-24).
  std::vector<cobordismgraph::Witness>().swap(initialWitnesses);
  if (!nameAliases.empty())
    std::cout << "[+] Name aliases: applied to " << aliasesApplied
              << " witness edges\n";
  if (!farSideResolutions.empty())
    std::cout << "[+] Far-side resolutions: applied to " << resolutionsApplied
              << " witness edges\n";
  if (!farSideExact.empty())
    std::cout << "[+] Far-side exact names: applied to " << exactApplied
              << " witness edges"
              << (exactRefused ? " (" + std::to_string(exactRefused) +
                                     " refused: component count differs)"
                               : std::string())
              << "\n";
  std::cout << "[+] Solver: derived bounds for " << bounds.size()
            << " names\n\n";

  std::vector<InputRow> pending;
  pending.reserve(rows.size());
  for (const auto &row : rows) {
    if (row.crossings > maxCrossings) {
      if (!outputRows.contains(row.name)) {
        OutputRow out;
        out.knot = row.name;
        out.status = "skipped";
        out.witnessKind = "none";
        out.literatureLo = row.lo;
        out.literatureHi = row.hi;
        outputRows[row.name] = std::move(out);
      }
      continue;
    }
    pending.push_back(row);
  }
  std::stable_sort(pending.begin(), pending.end(),
                   [](const InputRow &a, const InputRow &b) {
                     return a.crossings < b.crossings;
                   });

  // --solve-only: re-derive every conclusion from the witness file and
  // stop. Emptying `pending` (rather than branching around the loop) means
  // the final solve/write below is the single code path that produces
  // --output, so a solve-only run and a search run can never disagree
  // about how a given witness set is interpreted.
  if (solveOnly) {
    std::cout << "[+] --solve-only: re-deriving from " << witnesses.size()
              << " witnesses, no searching.\n\n";
    pending.clear();
  }

  // Re-solves from the full witness set and rewrites every affected output
  // row. Cheap (pure integer relaxation over the witness list), so it runs
  // after every row rather than only at the end -- which is what lets one
  // row's harvested cobordisms settle a later row before it is ever
  // searched.
  auto resolveAll = [&]() -> std::vector<std::string> {
    bounds = cobordismgraph::propagate(solverWitnesses(), names, externalProofs);
    std::vector<std::string> contradictions;

    // Every name that could need its row rewritten -- crucially including
    // rows that currently HAVE a row but no longer have any derived bound.
    // Iterating `bounds` alone would leave such a row frozen at whatever a
    // previous (possibly buggier) solver wrote, which is exactly the case
    // --solve-only exists to correct.
    std::vector<std::string> toJudge;
    toJudge.reserve(bounds.size() + outputRows.size());
    for (const auto &[name, unused] : bounds)
      toJudge.push_back(name);
    for (const auto &[name, unused] : outputRows)
      if (!bounds.contains(name))
        toJudge.push_back(name);
    // --link-classes: a member is its class's node, so it is judged whenever
    // that node has bounds, row or no row yet.
    for (const auto &[member, canonical] : linkClasses)
      if (bounds.contains(canonical) && !bounds.contains(member) &&
          !outputRows.contains(member))
        toJudge.push_back(member);

    for (const std::string &name : toJudge) {
      cobordismgraph::Bounds b; // default = nothing derived
      // A class member's row reads its class's node (--link-classes).
      if (auto it = bounds.find(classOf(name)); it != bounds.end())
        b = it->second;
      cobordismgraph::Verdict v = cobordismgraph::judge(name, b, names);
      if (v.status == cobordismgraph::Status::contradiction)
        contradictions.push_back(v.reason);
      auto existing = outputRows.find(name);
      // Only names we actually track get a row: the graph is full of
      // incidental nodes (bare isoSigs, unlinks) that are useful for
      // chaining but aren't results in their own right.
      if (existing == outputRows.end() && !names.find(name))
        continue;
      OutputRow updated = rowFromVerdict(
          name, v, existing == outputRows.end() ? nullptr : &existing->second,
          bounds);
      outputRows[name] = std::move(updated);
    }
    return contradictions;
  };

  // How this run's searches differ from the cascade's: one field per
  // divergence the plan removes in phase 4(b) (search/search.h, SearchPolicy).
  cascade::SearchPolicy policy;
  policy.frontierNeedsExamined = true;
  policy.judgeInSearch = true;
  policy.signing = cascade::SearchPolicy::Signing::duringSearch;

  // Every row's search shape, but for its boundary condition (per row,
  // below).
  cascade::SearchShape searchShape;
  searchShape.iddfsIterations = iddfsIterations;
  searchShape.iddfsStep = iddfsStep;
  searchShape.iddfsStart = iddfsStart;
  searchShape.iddfsFinalThreads = iddfsFinalThreads;
  searchShape.maxFaces = maxFaces;
  searchShape.rootBudgetStart = rootBudgetStart;
  searchShape.rootBudgetGrowth = rootBudgetGrowth;
  searchShape.resolveUnlinked = resolveUnlinked;
  searchShape.limits = limits;

  const auto sweepStart = std::chrono::steady_clock::now();
  // Searches whose accounting failed (divergence 2): the run exits 2.
  std::vector<std::string> unaccounted;
  size_t processedThisRun = 0;
  size_t searchedThisRun = 0;
  bool sweepTimedOut = false;

  for (const auto &row : pending) {
    if (sweepTimeLimit &&
        std::chrono::steady_clock::now() - sweepStart >
            std::chrono::duration<double>(*sweepTimeLimit)) {
      sweepTimedOut = true;
      break;
    }

    std::cout << "[+] Searching " << row.name << " (literature [" << row.lo
              << ", " << row.hi << "], " << row.crossings
              << " crossings)...\n";

    // The row's own cobordism graph (divergence 6), which judges its finds,
    // and the row the search runs in: rowsearch::buildRow() of its PD,
    // collared through every layer. buildRow() and the graph's row
    // certification throw for a bad PD, a row map that cannot be built or
    // checked, or a triangulated link that does not redraw as its diagram;
    // letting that escape would abort the whole sweep over one bad row.
    std::unique_ptr<cascade::SearchJudge> judge;
    bool buildFailed = false;
    try {
      judge = std::make_unique<cascade::SearchJudge>(
          row.name, row.pdNotation, thickenLayers, row.lo, *graphTables, *graphNamer,
          names.symmetries(), numThreads);
      if (judge->row().rowBuild().orientation->divergedFromDefaultIsomorphism)
        std::cerr << "[i] " << row.name
                  << ": the diagram's triangulation has a symmetry moving "
                     "L; using the map that takes L onto its own seed "
                     "(isIsomorphicTo() would not have)\n";
    } catch (const std::exception &e) {
      std::cerr << "[!] " << row.name << ": failed to build (" << e.what()
                << "), skipping\n";
      buildFailed = true;
    }
    const rowsearch::RowBuild *rbp = judge ? &judge->row().rowBuild() : nullptr;

    // A row that cannot be built, or that the search refuses: recorded as
    // such, and the run goes on.
    auto recordBuildFailure = [&] {
      OutputRow out;
      out.knot = row.name;
      out.status = "unresolved";
      out.witnessKind = "none";
      out.literatureLo = row.lo;
      out.literatureHi = row.hi;
      out.searchOutcome = "build-failed";
      outputRows[row.name] = std::move(out);
      writeOutputCsv(*outputPath, rows, outputRows);
    };
    if (buildFailed) {
      recordBuildFailure();
      continue;
    }

    // The row's frontier: carried on from, and recorded (see --frontier-dir).
    auto frontierPath = [&](const std::string &dir) {
      return dir + "/" + row.name + ".frontier";
    };
    std::optional<SearchFrontier> resumeFrom;
    if (resumeFrontierDir) {
      try {
        resumeFrom = SearchFrontier::load(frontierPath(*resumeFrontierDir));
      } catch (const std::exception &ex) {
        std::cout << "[!] " << row.name << ": WARNING: frontier not read ("
                  << ex.what() << "); searching from the start\n";
      }
    }

    cascade::SearchRequest request;
    request.name = row.name;
    request.shape = searchShape;
    // A multi-component link's own boundary necessarily puts more than one
    // of the surface's boundary curves on the single search-side ambient
    // component -- impossible under `connected`'s one-curve-per-ambient-
    // component rule, but exactly what `proper` allows. `connected` is the
    // (much tighter-pruning) default for ordinary single-component knots.
    //
    // But note what `connected` also costs, which is easy to miss: it caps
    // the curve count on EVERY ambient boundary component, the far side
    // included. So a knot row searched under `connected` can only ever
    // discover single-curve far sides -- it is structurally incapable of
    // finding a knot-to-LINK cobordism. Measured over a full sweep: all 110
    // witnesses from knot rows had a one-curve far side, and there were
    // zero knot-to-link edges, while link rows produced 62 link-to-knot
    // ones. --boundary-condition proper lifts that, at the cost of much
    // weaker pruning.
    const rowsearch::RowBuild &rb = *rbp;
    request.shape.condition =
        rowsearch::conditionFor(boundaryConditionMode, rb.componentCount);
    // The search stops at the surface target, the per-row or sweep time
    // limit, or quiescence; the boundary drain then finishes, unless
    // --skip-drain-on-timeout.
    request.surfaceTarget = surfaceTarget;
    request.seconds = perKnotTimeLimit;
    request.sweepSeconds = sweepTimeLimit;
    request.sweepStart = sweepStart;
    request.quiescenceSeconds = harvestQuiescence;
    request.skipDrainOnTimeout = skipDrainOnTimeout;
    request.resume = resumeFrom ? &*resumeFrom : nullptr;
    request.recordFrontier = frontierDir.has_value();
    request.pairSigCacheDir = pairSigCacheDir;
    // Every boundary by its complement, unless the row draws its far sides.
    request.diagramNaming = !useCone;
    request.unseeded = true;
    // A knot row's complement goes into the census after its search (when
    // census writes are on).
    request.censusName = cobordismgraph::baseName(row.name);
    request.sweep = {.record = &recorded,
                     .names = &names,
                     .bounds = &bounds,
                     .literatureLo = row.lo,
                     .literatureHi = row.hi,
                     .thickenLayers = thickenLayers};
    // Each find, judged by the row's own cobordism graph as it is kept.
    request.row = &judge->row();
    long long judged = 0;
    request.judge = [&](const cascade::KeptSurface &k) {
      const cascade::SearchJudge::Verdict v =
          judge->add(k.link, k.genus, "find#" + std::to_string(judged++));
      cascade::FindJudgement j;
      if (!v.contradictions.empty())
        j.contradiction = v.contradictions.front();
      return j;
    };
    request.outputs = {.progress = true,
                       .surfaceStats = surfaceStatsPath,
                       .surfaceLog = surfaceLogPath,
                       .rejectionSamples =
                           rejectionSampleLog ? &*rejectionSampleLog : nullptr,
                       .selfIntersectionCensus = selfIntersectionCensusPath};

    // One searcher per row, so each row's exact names start from fresh
    // table caches, as they always have.
    const cascade::HopSearcher searcher(
        signatureTable ? &*signatureTable : nullptr,
        exactTables ? &*exactTables : nullptr, policy, numThreads);
    cascade::HopRun run;
    try {
      run = searcher.run(rb, request);
    } catch (const cascade::SeedInvariantFailure &f) {
      flagFatalBug(row.name + ": " + std::to_string(f.touching) +
                   " searchable non-seed triangles have an edge on the "
                   "search side, so found surfaces could change it.");
      haltIfFatalBugDetected();
    } catch (const cascade::SearchRefused &ex) {
      std::cerr << "[!] " << row.name << ": failed to build (" << ex.what()
                << "), skipping\n";
      recordBuildFailure();
      continue;
    }
    // A contradiction in the row's own cobordism graph: halts once the row's
    // witnesses are written, below.
    if (judge->failures() > 0)
      std::cout << "[!] " << row.name << ": " << judge->failures() << " of "
                << judge->finds() << " finds could not enter the cobordism graph\n";
    if (!run.fatal.empty())
      flagFatalBug(run.fatal);

    // Surfaces were accepted, yet not one reached the witness record. That
    // can be genuine (every one witnesses another oriented variant), but it
    // is also exactly what a broken gate looks like, so it never licenses a
    // negative.
    if (run.nothingExamined)
      std::cout << "[!] " << row.name << ": WARNING: " << run.described
                << " surfaces accepted but none reached the witness record; "
                   "no exhaustion claimed for this search\n";

    if (run.stats.deepestExhaustedCap && run.accountingFailure.empty() &&
        !run.nothingExamined && !run.drainSkipped) {
      OutputRow &out = outputRows[row.name];
      // Never let a shallower run overwrite a deeper exhaustive result.
      out.exhaustedDepth =
          std::max(out.exhaustedDepth, *run.stats.deepestExhaustedCap);
      std::cout << "[+] " << row.name << ": EXHAUSTIVE to "
                << *run.stats.deepestExhaustedCap
                << " added faces -- every root enumerated to completion, so "
                   "no cobordism exists for it at that depth.\n";
    }
    ++searchedThisRun;

    // Record what this search actually cost before anything else, so even
    // a fatal halt below leaves the bookkeeping behind.
    {
      OutputRow &out = outputRows[row.name];
      out.knot = row.name;
      out.literatureLo = row.lo;
      out.literatureHi = row.hi;
      out.searchedFaces = maxFaces.value_or(0);
      out.searchOutcome = run.outcome;
      if (out.status.empty())
        out.status = "unresolved";
    }

    std::vector<std::string> contradictions = resolveAll();
    {
      OutputRow &out = outputRows[row.name];
      out.searchedFaces = maxFaces.value_or(0);
      out.searchOutcome = run.outcome;
    }

    recorded.flush();
    writeOutputCsv(*outputPath, rows, outputRows);

    // The row's breadth, and its frontier: written only now that every
    // witness of the prefix it covers is on disk, and only when the row can
    // vouch for having examined every surface in it -- else a later run
    // would skip surfaces nobody looked at.
    if (resumeFrom || frontierDir) {
      rowsearch::printSweepBreadth(std::cout, row.name, run);
      if (frontierDir && run.recordedFrontier) {
        if (run.frontier) {
          try {
            std::filesystem::create_directories(*frontierDir);
            run.frontier->save(frontierPath(*frontierDir));
          } catch (const std::exception &ex) {
            std::cout << "[!] " << row.name << ": WARNING: frontier not "
                      << "written: " << ex.what() << "\n";
          }
        } else {
          std::cout << "[!] " << row.name << ": frontier not written: the "
                    << "row cannot vouch for every surface in its prefix\n";
        }
      }
    }

    for (const std::string &reason : contradictions)
      flagFatalBug(reason);
    // Divergence 2: a state that cannot occur halts the run, now that this
    // row's witnesses are written; any other accounting failure ends only
    // this search (its outcome is `unaccounted`; no frontier, no exhaustion
    // claim, above) and the run goes on to its other rows; completeness is
    // what a run without a goal is for, so the run then exits 2.
    if (run.impossible > 0) {
      flagFatalBug(row.name + ": surface accounting failed -- " +
                   std::to_string(run.impossible) +
                   " surfaces hit a state that cannot occur (accounting: " +
                   run.accounting + ").");
    } else if (!run.accountingFailure.empty()) {
      unaccounted.push_back(row.name);
      std::cout << "[!!] " << row.name << ": surface accounting failed -- "
                << run.accountingFailure
                << " (its search vouches for nothing: no frontier, no "
                   "exhaustion claim; the run goes on and exits 2)\n";
    }
    if (fatalBugDetected_.load())
      haltIfFatalBugDetected();

    rowsearch::printOutcome(std::cout, row.name, run);
    rowsearch::printIdentification(std::cout, row.name, run);
    rowsearch::printSearchProfile(std::cout, row.name, run);
    // The audit exists to catch exactly this; a wrong linking number prunes
    // (or keeps) surfaces it should not.
    if (run.petals.linkingDisagreements > 0) {
      flagFatalBug(row.name + ": " +
                   std::to_string(run.petals.linkingDisagreements) +
                   " petal linking numbers disagree between the cochain "
                   "and drilling routes (--audit-linking).");
      haltIfFatalBugDetected();
    }
    // Its own line, so the summary line above (parsed by
    // tools/orchestrate/dispatch.py's RE_OUTCOME) is unchanged.
    if (*resolveUnlinked)
      std::cout << "[+] " << row.name << ": " << run.stats.resolvedCount
                << " of " << run.stats.satisfyingCount
                << " accepted surfaces have unlinked self-intersections "
                   "(--resolve-unlinked)\n";

    auto verdict = outputRows.find(row.name);
    if (verdict != outputRows.end()) {
      const OutputRow &out = verdict->second;
      if (out.status == "verified")
        std::cout << "\x1b[1;32m[+] " << row.name << ": VERIFIED at genus "
                  << out.resolvedGenus << " (constructive)\x1b[0m\n";
      else if (out.status == "verified-assisted")
        std::cout << "\x1b[1;36m[+] " << row.name << ": reaches genus "
                  << out.resolvedGenus
                  << ", but via another name's literature value -- NOT an "
                     "independent verification\x1b[0m\n";
      else if (out.status == "improved")
        std::cout << "\x1b[1;33m[!!] " << row.name
                  << ": IMPROVED upper bound to " << out.resolvedGenus
                  << ", beating the literature's " << out.literatureHi
                  << " (" << out.witnessBasis << ")\x1b[0m\n";
      else
        std::cout << "[+] " << row.name << ": " << out.status
                  << (out.derivedHi.empty() ? "" : " (upper bound " +
                                                       out.derivedHi + ")")
                  << "\n";
    }

    ++processedThisRun;
  }

  // One last solve, so a run that searched nothing (or was cut short) still
  // reflects everything its witness file knows.
  for (const std::string &reason : resolveAll())
    flagFatalBug(reason);
  recorded.flush();
  writeOutputCsv(*outputPath, rows, outputRows);
  haltIfFatalBugDetected();

  size_t verified = 0, verifiedAssisted = 0, improved = 0,
         pinnedCount = 0, unresolvedCount = 0;
  for (const auto &[name, out] : outputRows) {
    if (out.status == "verified")
      ++verified;
    else if (out.status == "verified-assisted")
      ++verifiedAssisted;
    else if (out.status == "improved")
      ++improved;
    else if (out.status == "pinned")
      ++pinnedCount;
    else if (out.status == "unresolved")
      ++unresolvedCount;
  }

  std::cout << "\n[+] Done. Searched " << searchedThisRun << " of "
            << processedThisRun << " rows visited this run"
            << (sweepTimedOut ? " (sweep time limit reached)" : "") << ".\n";
  std::cout << "[+] Witness file: " << witnesses.size() << " witnesses in "
            << cobordismsPath << "\n";
  std::cout << "[+] Totals across " << outputRows.size()
            << " tracked names: " << verified << " verified ("
            << verifiedAssisted << " more only with literature help), "
            << improved << " improved, " << pinnedCount << " pinned, "
            << unresolvedCount << " unresolved.\n";
  if (!unaccounted.empty()) {
    std::cerr << "[!] " << unaccounted.size()
              << " search(es) failed their surface accounting (first: "
              << unaccounted.front() << "); exiting 2\n";
    return 2;
  }
  return 0;
}
