//
//  verifyslicegenus.cpp
//
//  Created by John Teague on 07/29/2026.
//

#include <algorithm>
#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
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

#include "cobordismbuilder.h"
#include "cobordismgraph.h"
#include "collar.h"
#include "csvwriter.h"
#include "embeddingsearch.h"
#include "surfacesearch.h"
#include "knotbuilder.h"
#include "linkcomplement.h"
#include "identifycomplement.h"

using namespace cobordismgraph;

namespace {

// Redraws a small status block in place (ANSI cursor-rewind, same erase
// trick as surfer.cpp's RollingReport). Kept file-local and minimal (a
// plain function, not a class) since this driver only ever needs one
// concurrent rolling block -- the live-DFS progress (callbacks.onProgress)
// and the post-search boundary-processing progress
// (callbacks.onBoundaryProcessing*) never run at the same time (the
// latter only starts once the former's DFS phase has fully joined), so
// both can safely share the same redraw region.
size_t progressPrevLines_ = 0;

void redrawProgressBlock(const std::string &text) {
  if (progressPrevLines_ > 0)
    std::cerr << "\x1b[" << progressPrevLines_ << "F\x1b[0J";
  std::cerr << text;
  progressPrevLines_ =
      static_cast<size_t>(std::count(text.begin(), text.end(), '\n'));
}

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

// Fired from callbacks.onProgress once per second while a knot's search runs.
void printProgress(const SearchStats &stats, SurfaceSearch &e) {
  std::ostringstream report;
  report << "[+] elapsed: " << formatElapsed(stats.elapsed)
         << " | candidates examined: " << stats.foundCount
         << " | embedded surfaces found: " << stats.embeddedCount
         << " | satisfying boundary condition: " << stats.satisfyingCount
         << "\n";
  report << "[+] surface homeomorphism types found so far: "
         << e.surfaceTypeTally().summary() << "\n";
  redrawProgressBlock(report.str());
}

// Fired from callbacks.onBoundaryProcessingProgress once per second during
// the post-search boundary-identification phase (processRemainingSurfaceBoundaries).
// `resolvedGenus`, if set, reflects the driver's own resolvedThisRow/target
// for the knot currently being processed (resolution is all-or-nothing --
// see resolvesNow's own comment -- so there is no partial/intermediate
// genus value to show, only "not yet" vs. the exact target once hit), so
// it's visible at a glance how the current genus stands against the
// literature target even while this phase is still churning through
// whatever else was queued.
void printBoundaryProgress(size_t processed, size_t total,
                           std::chrono::steady_clock::duration elapsed,
                           std::optional<int> resolvedGenus, int target) {
  double elapsedSec = std::chrono::duration<double>(elapsed).count();
  double rate = (elapsedSec > 0.0) ? static_cast<double>(processed) / elapsedSec
                                   : 0.0;

  std::ostringstream report;
  report << "[+] boundary processing: elapsed " << formatElapsed(elapsed)
         << " | processed " << processed << "/" << total;
  if (rate > 0.0 && processed < total) {
    auto eta = std::chrono::duration_cast<std::chrono::steady_clock::duration>(
        std::chrono::duration<double>(
            static_cast<double>(total - processed) / rate));
    report << " | ETA " << formatElapsed(eta);
  }
  report << "\n";
  report << "[+] slice genus: "
         << (resolvedGenus ? std::to_string(*resolvedGenus) : "?") << "/"
         << target << " (literature target)"
         << (resolvedGenus ? " -- ACHIEVED" : "") << "\n";
  redrawProgressBlock(report.str());
}

// ─────────────────────────────────────────────────────────────────────────
// CSV parsing helpers
// ─────────────────────────────────────────────────────────────────────────

// One row of the input knot table (Name,PD Notation,Genus-4D). No RFC-4180
// quoting appears in that file (PD Notation uses ';' internally, never a
// literal comma), so a naive two-comma split suffices -- see the input
// table's own format, confirmed during design.
bool splitInputLine(const std::string &line, std::string &name,
                    std::string &pd, std::string &genusField) {
  size_t c1 = line.find(',');
  if (c1 == std::string::npos)
    return false;
  size_t c2 = line.find(',', c1 + 1);
  if (c2 == std::string::npos)
    return false;
  name = line.substr(0, c1);
  pd = line.substr(c1 + 1, c2 - c1 - 1);
  genusField = line.substr(c2 + 1);
  if (!genusField.empty() && genusField.back() == '\r')
    genusField.pop_back();
  return true;
}

// Parses "N" or "[lo;hi]" into lo/hi (lo == hi in the plain-integer case).
void parseGenusField(const std::string &field, int &lo, int &hi) {
  if (!field.empty() && field.front() == '[') {
    size_t semi = field.find(';');
    lo = std::stoi(field.substr(1, semi - 1));
    hi = std::stoi(field.substr(semi + 1, field.size() - semi - 2));
  } else {
    lo = hi = std::stoi(field);
  }
}

std::vector<InputRow> loadInputCsv(const std::filesystem::path &path) {
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("Cannot open input CSV: " + path.string());

  std::vector<InputRow> rows;
  std::string line;
  std::getline(in, line); // header
  while (std::getline(in, line)) {
    if (line.empty())
      continue;
    std::string name, pd, genusField;
    if (!splitInputLine(line, name, pd, genusField))
      continue;
    InputRow row;
    row.name = name;
    row.pdNotation = pd;
    parseGenusField(genusField, row.lo, row.hi);
    // Crossing count is derived from the PD code itself (works uniformly
    // for both knot names like "13n_1109" and link names like "L10a1{0}",
    // which have no leading digit run to parse) rather than from `name`.
    row.crossings =
        static_cast<int>(knotbuilder::parsePDCode(pd).size());
    rows.push_back(std::move(row));
  }
  return rows;
}

// Loads a literature table for its names and bounds only, skipping the PD
// code entirely. Used for tables that aren't this run's --input: we need
// their names (to expand orientation-blind identifications into candidate
// sets) and their bounds, but never build a triangulation from them, so
// there is no reason to pay parsePDCode()'s cost across 12k+ rows.
size_t loadNameTable(const std::filesystem::path &path,
                     cobordismgraph::NameTable &names) {
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("Cannot open name table: " + path.string());

  size_t loaded = 0;
  std::string line;
  std::getline(in, line); // header
  while (std::getline(in, line)) {
    if (line.empty())
      continue;
    std::string name, pd, genusField;
    if (!splitInputLine(line, name, pd, genusField))
      continue;
    int lo = 0, hi = 0;
    try {
      parseGenusField(genusField, lo, hi);
    } catch (const std::exception &) {
      continue;
    }
    names.addLiterature(name, lo, hi);
    ++loaded;
  }
  return loaded;
}

// Minimal RFC-4180 field parser (quotes, doubled-quote escaping) -- needed
// to read our OWN --output file back on resume, since witness_pairsig/
// depends_on may have been written through csvField() and can contain
// commas.
std::vector<std::string> parseCsvLine(const std::string &line) {
  std::vector<std::string> fields;
  size_t i = 0;
  while (i <= line.size()) {
    std::string field;
    if (i < line.size() && line[i] == '"') {
      ++i;
      while (i < line.size()) {
        if (line[i] == '"') {
          if (i + 1 < line.size() && line[i + 1] == '"') {
            field += '"';
            i += 2;
          } else {
            ++i;
            break;
          }
        } else {
          field += line[i++];
        }
      }
    } else {
      while (i < line.size() && line[i] != ',')
        field += line[i++];
    }
    fields.push_back(field);
    if (i < line.size() && line[i] == ',') {
      ++i;
      continue;
    }
    break;
  }
  return fields;
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
};

// Which BoundaryCondition to search a row under; see the switch in the
// main loop for what each costs.
enum class BoundaryConditionMode { automatic, connected, proper };

constexpr const char *OUTPUT_HEADER =
    "knot,resolved_genus,status,witness_kind,witness_pairsig,via_knot,"
    "via_edge_genus,depends_on,literature_lo,literature_hi,"
    "derived_lo,derived_hi,witness_basis,tubed,searched_faces,search_outcome";

std::string formatOutputRow(const OutputRow &r) {
  std::ostringstream out;
  out << csvField(r.knot) << ',' << r.resolvedGenus << ',' << r.status << ','
      << r.witnessKind << ',' << csvField(r.witnessPairSig) << ','
      << csvField(r.viaKnot) << ',' << r.viaEdgeGenus << ','
      << csvField(r.dependsOn) << ',' << r.literatureLo << ','
      << r.literatureHi << ',' << r.derivedLo << ',' << r.derivedHi << ','
      << r.witnessBasis << ',' << (r.tubed ? "true" : "false") << ','
      << r.searchedFaces << ',' << r.searchOutcome;
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
    result[r.knot] = std::move(r);
  }
  return result;
}

// ─────────────────────────────────────────────────────────────────────────
// Witness file (--cobordisms) I/O
// ─────────────────────────────────────────────────────────────────────────
//
// The point of persisting these separately from --output is that a witness
// is a fact ("a surface with this boundary and this genus exists") while an
// --output row is a conclusion. Conclusions get better whenever the solver,
// the literature tables, or the identification improves; facts do not. So
// every fact a search paid for is written here once and never re-searched,
// and `--solve-only` re-derives all the conclusions from them in seconds.

constexpr const char *COBORDISMS_HEADER =
    "kind,subject,subject_components,other,other_candidates,other_components,"
    "genus,tubed,pairsig,source_row,thicken_layers,max_faces";

std::string formatWitness(const cobordismgraph::Witness &w) {
  std::ostringstream candidates;
  for (size_t i = 0; i < w.otherCandidates.size(); ++i) {
    if (i)
      candidates << ';';
    candidates << w.otherCandidates[i];
  }
  std::ostringstream out;
  out << (w.kind == cobordismgraph::WitnessKind::direct ? "direct"
                                                        : "cobordism")
      << ',' << csvField(w.subject) << ',' << w.subjectComponents << ','
      << csvField(w.other) << ',' << csvField(candidates.str()) << ','
      << w.otherComponents << ',' << w.genus << ','
      << (w.tubed ? "true" : "false") << ',' << csvField(w.pairSig) << ','
      << csvField(w.sourceRow) << ',' << w.thickenLayers << ',' << w.maxFaces;
  return out.str();
}

std::vector<cobordismgraph::Witness>
loadWitnesses(const std::filesystem::path &path) {
  std::vector<cobordismgraph::Witness> result;
  std::ifstream in(path);
  if (!in)
    return result;

  std::string line;
  std::getline(in, line); // header
  while (std::getline(in, line)) {
    if (line.empty())
      continue;
    auto f = parseCsvLine(line);
    if (f.size() < 12)
      continue;
    cobordismgraph::Witness w;
    w.kind = f[0] == "direct" ? cobordismgraph::WitnessKind::direct
                              : cobordismgraph::WitnessKind::cobordism;
    w.subject = f[1];
    try {
      w.subjectComponents = std::stoi(f[2]);
      w.otherComponents = std::stoi(f[5]);
      w.genus = std::stoi(f[6]);
      w.thickenLayers = std::stoi(f[10]);
      w.maxFaces = std::stoll(f[11]);
    } catch (const std::exception &) {
      continue;
    }
    w.other = f[3];
    if (!f[4].empty()) {
      std::istringstream candidates(f[4]);
      std::string one;
      while (std::getline(candidates, one, ';'))
        if (!one.empty())
          w.otherCandidates.push_back(one);
    }
    w.tubed = f[7] == "true";
    w.pairSig = f[8];
    w.sourceRow = f[9];
    result.push_back(std::move(w));
  }
  return result;
}

// Rewrites the whole witness file via write-to-temp + atomic rename, the
// same crash-safety pattern writeOutputCsv() uses.
void writeWitnesses(const std::filesystem::path &path,
                    const std::vector<cobordismgraph::Witness> &witnesses) {
  std::filesystem::path tmp = path;
  tmp += ".tmp";
  {
    std::ofstream out(tmp, std::ios::trunc);
    if (!out)
      throw std::runtime_error("Cannot open " + tmp.string() +
                               " for writing");
    out << COBORDISMS_HEADER << "\n";
    for (const auto &w : witnesses)
      out << formatWitness(w) << "\n";
  }
  std::filesystem::rename(tmp, path);
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
  std::filesystem::path tmp = path;
  tmp += ".tmp";
  {
    std::ofstream out(tmp, std::ios::trunc);
    if (!out)
      throw std::runtime_error("Cannot open " + tmp.string() +
                               " for writing");
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
  }
  std::filesystem::rename(tmp, path);
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
    out.witnessPairSig = b.pairSig;
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
         "    [ --threads N ] [ --thicken-layers N ] [ --cone | --no-cone ]\n"
         "    [ --collar-layers N ] [ --iddfs-iterations N --iddfs-step D ]\n"
         "    [ --iddfs-start N ] [ --iddfs-final-threads N ]\n"
         "    [ --per-knot-time-limit S ] [ --surface-log <path> ]\n"
         "    [ --no-census-updates ] [ --no-retriangulate-on-miss ]\n"
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
         "                     every knot/link processed; already-resolved "
         "rows from a\n"
         "                     prior run are trusted directly and not "
         "re-searched. A\n"
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
         "                     row's search allows (default: auto).\n"
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
      << "    --harvest      : Don't stop a row's search once that row is "
         "settled;\n"
         "                     keep recording every distinct cobordism it "
         "finds. One\n"
         "                     expensive search then yields many edges "
         "instead of\n"
         "                     one, and those edges bound OTHER rows via "
         "the solver.\n"
         "                     Requires a stopping rule (--max-faces, "
         "--harvest-\n"
         "                     quiescence, or --per-knot-time-limit), "
         "since resolving\n"
         "                     the row no longer ends it (default: off).\n";
  std::cerr
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
  std::cerr << "    --no-census-updates : Disable live census seeding "
               "(default: on).\n";
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
               "    --boundary-tally-cap, --census-db : as surfer.cpp.\n";
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
  std::optional<long long> maxFaces;
  bool harvest = false;
  std::optional<double> harvestQuiescence;
  std::optional<double> sweepTimeLimit;
  bool solveOnly = false;
  // Identification is the product, so by default a timed-out row still drains
  // its boundary queue to completion; only the SEARCH is bounded.
  bool skipDrainOnTimeout = false;
  BoundaryConditionMode boundaryConditionMode =
      BoundaryConditionMode::automatic;
  int maxCrossings = 13;
  std::optional<double> perKnotTimeLimit;
  std::optional<std::string> surfaceLogPath;
  bool censusUpdates = true;

  unsigned numThreads = std::thread::hardware_concurrency();
  if (numThreads == 0)
    numThreads = 1;

  int thickenLayers = 1;
  bool useCone = true;
  int collarLayers = 1;

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
    } else if (arg == "--max-faces") {
      if (i + 1 >= argc)
        usage(argv[0], "--max-faces requires a value.");
      try {
        maxFaces = std::stoll(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--max-faces requires an integer value.");
      }
    } else if (arg == "--harvest") {
      harvest = true;
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
    } else if (arg == "--surface-log") {
      if (i + 1 >= argc)
        usage(argv[0], "--surface-log requires a value.");
      surfaceLogPath = argv[++i];
    } else if (arg == "--no-census-updates") {
      censusUpdates = false;
    } else if (arg == "--no-retriangulate-on-miss") {
      retriangulateOnMissArg = false;
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
  if (harvest && !maxFaces && !harvestQuiescence && !perKnotTimeLimit)
    usage(argv[0],
          "--harvest needs a stopping rule: without --max-faces, "
          "--harvest-quiescence or --per-knot-time-limit, a search has "
          "nothing to end it (resolving the row no longer does, which is "
          "the point of --harvest) and the sweep will never reach its "
          "second row.");
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
  bool censusLoaded = census::setCensusPath(censusPath);

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

  std::unordered_map<std::string, OutputRow> outputRows =
      loadOutputCsv(*outputPath);

  std::vector<cobordismgraph::Witness> witnesses =
      loadWitnesses(cobordismsPath);
  std::cout << "[+] Resuming with " << witnesses.size()
            << " previously-recorded witnesses from " << cobordismsPath
            << "\n";

  // Every conclusion is re-derived from the witness set on every run, so a
  // solver fix or a literature-table update takes effect on rows that were
  // searched long ago without re-searching any of them.
  auto bounds = cobordismgraph::propagate(witnesses, names);
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
    bounds = cobordismgraph::propagate(witnesses, names);
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

    for (const std::string &name : toJudge) {
      cobordismgraph::Bounds b; // default = nothing derived
      if (auto it = bounds.find(name); it != bounds.end())
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

  // Records one witness if it's genuinely new, returning whether it was.
  // The dedup matters a lot in --harvest mode: a single search reports
  // thousands of near-identical surfaces, and capturing a pair signature
  // for each would dominate the run (see SurfaceFoundInfo::capturePairSig).
  std::mutex witnessMutex;
  std::atomic<long long> lastNewWitnessTick{0};
  auto tickNow = [] {
    return std::chrono::duration_cast<std::chrono::milliseconds>(
               std::chrono::steady_clock::now().time_since_epoch())
        .count();
  };
  auto recordWitness = [&](cobordismgraph::Witness w,
                           const std::function<std::string()> &capturePairSig)
      -> bool {
    std::lock_guard<std::mutex> lock(witnessMutex);
    if (cobordismgraph::haveWitness(witnesses, w))
      return false;
    if (capturePairSig)
      w.pairSig = capturePairSig();
    witnesses.push_back(std::move(w));
    lastNewWitnessTick.store(tickNow(), std::memory_order_relaxed);
    return true;
  };

  const auto sweepStart = std::chrono::steady_clock::now();
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

    // Already settled? `verified` means we constructed a surface meeting
    // the literature's own lower bound, so there is nothing left to find.
    {
      auto it = outputRows.find(row.name);
      if (it != outputRows.end() &&
          (it->second.status == "verified" || it->second.status == "pinned")) {
        continue;
      }
      // T5: don't repeat a search we have already run to exhaustion at
      // this budget -- it would enumerate exactly the same surfaces and
      // learn exactly nothing. A bigger --max-faces does make it worth
      // redoing, which is what turns a re-run into progressive deepening.
      if (it != outputRows.end() && it->second.searchOutcome == "exhausted" &&
          maxFaces && it->second.searchedFaces >= *maxFaces) {
        continue;
      }
    }

    std::cout << "[+] Searching " << row.name << " (literature [" << row.lo
              << ", " << row.hi << "], " << row.crossings
              << " crossings)...\n";

    knotbuilder::PDCode pdcode;
    knotbuilder::TriangulationWithLink link;
    std::optional<CobordismBuilder<3>> cobOpt;
    std::optional<SurfaceSearch> eOpt;
    std::vector<int> seedFaces;
    std::optional<cobordismgraph::RowOrientation> rowOrientation;
    size_t searchSideBC = 0;
    int componentCount = 1;
    std::string rowOwnName;
    regina::Triangulation<4> tri;
    bool buildFailed = false;

    // buildRowOrientation() lives inside this try alongside parsePDCode/
    // buildLink: it throws regina::InvalidArgument too (empty edge set, or
    // a boundary component that isn't isomorphic to the row's own
    // triangulation), and letting that escape would abort the entire
    // sweep over one bad row rather than skipping it.
    try {
      pdcode = knotbuilder::parsePDCode(row.pdNotation);
      link = knotbuilder::buildLink(pdcode);

      auto &[t2, edges2, reversed2] = link;
      Link linkGrouping(t2, edges2);
      componentCount = linkGrouping.countComponents();
      rowOwnName = identify::identify(linkGrouping);

      std::vector<int> edgeIndices;
      edgeIndices.reserve(edges2.size());
      for (const regina::Edge<3> *e : edges2)
        edgeIndices.push_back(static_cast<int>(e->index()));

      cobOpt.emplace(t2);
      CobordismBuilder<3> &cob = *cobOpt;
      CollarBuilder collarBuilder(edgeIndices);
      for (int i = 0; i < thickenLayers; ++i) {
        cob.thicken();
        if (i < collarLayers)
          collarBuilder.addLayer(cob);
      }
      if (useCone)
        cob.cone();

      searchSideBC = cob.baseBoundaryComponent()->index();
      tri = cob.getCobordism();

      if (collarLayers > 0)
        for (regina::Triangle<4> *t : collarBuilder.resolve())
          seedFaces.push_back(static_cast<int>(t->index()));

      if (!seedFaces.empty())
        rowOrientation = cobordismgraph::buildRowOrientation(
            edges2, reversed2, tri.boundaryComponent(searchSideBC)->build());
    } catch (const regina::InvalidArgument &e) {
      std::cerr << "[!] " << row.name << ": failed to build (" << e.what()
                << "), skipping\n";
      buildFailed = true;
    }

    if (buildFailed) {
      OutputRow out;
      out.knot = row.name;
      out.status = "unresolved";
      out.witnessKind = "none";
      out.literatureLo = row.lo;
      out.literatureHi = row.hi;
      out.searchOutcome = "build-failed";
      outputRows[row.name] = std::move(out);
      writeOutputCsv(*outputPath, rows, outputRows);
      continue;
    }

    if (seedFaces.empty())
      eOpt.emplace(tri);
    else
      eOpt.emplace(tri, seedFaces, searchSideBC);
    SurfaceSearch &e = *eOpt;
    e.configureLimits(limits);

    std::optional<CsvWriter> surfaceLog;
    if (surfaceLogPath)
      surfaceLog.emplace(*surfaceLogPath,
                         "orientable,genus,tubed_genus,punctures,triangles,"
                         "pairsig",
                         numThreads);

    // Why the search stopped, for the T5 bookkeeping. Set by whoever
    // requests the stop; "exhausted" means nobody did and the search ran
    // out of candidates on its own, which is only possible with
    // --max-faces (see EmbeddingSearch::search()'s hardFaceCap).
    std::string searchOutcome = "exhausted";
    std::mutex outcomeMutex;
    auto noteStop = [&](const char *why) {
      std::lock_guard<std::mutex> lock(outcomeMutex);
      if (searchOutcome == "exhausted")
        searchOutcome = why;
    };

    std::atomic<long long> orientationRejections{0};
    std::atomic<long long> newWitnessesThisRow{0};
    std::atomic<bool> resolvedThisRow{false};

    lastNewWitnessTick.store(tickNow(), std::memory_order_relaxed);

    SurfaceSearchCallbacks callbacks;
    callbacks.onProgress = [&](const SearchStats &stats) {
      printProgress(stats, e);
    };
    callbacks.onBoundaryProcessingStarted = [&](size_t total,
                                                unsigned threads) {
      progressPrevLines_ = 0;
      std::cerr << "[+] boundary processing: " << total
                << " queued surfaces, " << threads << " threads\n";
    };
    callbacks.onBoundaryProcessingProgress =
        [&](size_t processed, size_t total,
            std::chrono::steady_clock::duration elapsed) {
          printBoundaryProgress(
              processed, total, elapsed,
              resolvedThisRow.load() ? std::optional<int>(row.lo)
                                     : std::nullopt,
              row.hi);
        };
    callbacks.onBoundaryProcessingComplete =
        [&](size_t total, std::chrono::steady_clock::duration elapsed) {
          progressPrevLines_ = 0;
          std::cerr << "[+] boundary processing: done (" << total
                    << " processed in " << formatElapsed(elapsed) << ")\n";
        };

    callbacks.onSurfaceBoundaryProcessed = [&](const SurfaceBoundaryInfo
                                                   &info) {
      if (surfaceLog) {
        std::string pairSig =
            info.capturePairSig ? info.capturePairSig() : std::string{};
        std::ostringstream row2;
        row2 << (info.orientable ? "true" : "false") << ',' << info.genus
             << ',' << info.tubedGenus << ',' << info.punctures << ','
             << info.triangleCount << ',' << csvField(pairSig);
        surfaceLog->writeRow(row2.str());
      }

      if (!info.orientable)
        return; // defensive; orientableOnly=true already prunes these

      // A disconnected find is NOT discarded any more. Its components tube
      // into a single connected surface with the same boundary and genus
      // exactly info.tubedGenus (see KnottedSurface::tubedSurfaceType), so
      // it witnesses precisely what a connected find of that genus would.
      // This is what makes multi-component links tractable at all: their
      // seeded collar starts as one disjoint annulus per component, and
      // nothing forces the DFS to ever bridge them.
      const int witnessGenus = info.tubedGenus;
      const bool tubed = !info.connected;

      auto capturePairSig = [&info] {
        return info.capturePairSig ? info.capturePairSig() : std::string{};
      };

      BoundarySplit split =
          splitBoundary(info.boundaryComponents, searchSideBC, rowOwnName);

      if (split.searchCurveCount != static_cast<size_t>(componentCount))
        return; // this row's own link isn't (wholly) the boundary here, so
                // whatever this surface witnesses, it isn't about this row

      // A mixed orientation match means the surface witnesses a DIFFERENT
      // oriented variant of this same-complement diagram -- exactly the
      // L6a3{0}/L6a3{1} misattribution this check exists to catch.
      if (rowOrientation) {
        std::vector<OrientedCurve> searchSideCurves;
        bool foundSearchSide = false;
        for (auto &[c, curves] : info.captureOrientedBoundaryLinks()) {
          if (c == searchSideBC) {
            searchSideCurves = std::move(curves);
            foundSearchSide = true;
            break;
          }
        }
        if (!foundSearchSide ||
            !cobordismgraph::matchesRowOrientation(*rowOrientation,
                                                   searchSideCurves)) {
          orientationRejections.fetch_add(1, std::memory_order_relaxed);
          return;
        }
      }

      cobordismgraph::Witness w;
      w.subject = row.name;
      w.subjectComponents = componentCount;
      w.genus = witnessGenus;
      w.tubed = tubed;
      w.sourceRow = row.name;
      w.thickenLayers = thickenLayers;
      w.maxFaces = maxFaces.value_or(0);

      if (split.otherSides.empty()) {
        w.kind = cobordismgraph::WitnessKind::direct;
      } else if (split.otherSides.size() == 1) {
        // A genuinely-linked far side is no longer refused. Its
        // orientation is unknowable from a complement, so it is recorded
        // as the SET of oriented variants it could be, and the solver
        // takes the max/min over that set -- a weaker bound than an exact
        // identification would give, but a sound one, where the old
        // behaviour was to discard the cobordism entirely.
        const BoundarySide &far = split.otherSides.front();
        // Normalized, because identify() decorates a translated census hit
        // as "4_1 (m004 : #1)" while the --input tables call it "4_1" --
        // leaving the decoration on would make the far side a different
        // graph node from the row for the very same knot.
        std::string farName = cobordismgraph::normalizeIdentifiedName(far.name);
        w.kind = cobordismgraph::WitnessKind::cobordism;
        w.other = farName;
        w.otherComponents = far.components;
        w.otherCandidates = names.candidates(farName, far.components);
      } else {
        return; // a (k+1)-way cobordism; sound to use, but not yet
                // implemented -- see the plan's out-of-scope note
      }

      if (!recordWitness(w, capturePairSig))
        return;
      newWitnessesThisRow.fetch_add(1, std::memory_order_relaxed);

      // Does this witness alone already settle the row? Checked cheaply
      // against the bounds as of the last full solve rather than by
      // re-running the solver here: the far side's own bound can only have
      // improved since, so this under-reports at worst, and the next
      // resolveAll() picks up anything it missed.
      int implied = cobordismgraph::NO_UPPER_BOUND;
      // Whether that bound leans on a literature value anywhere. It decides
      // whether the row may be considered settled: a conclusion resting on
      // someone else's literature number is not an independent verification
      // (Basis::literatureAssisted), and stopping the search on one would
      // forfeit the chance of finding the surface that upgrades it.
      bool impliedAssisted = false;
      if (w.kind == cobordismgraph::WitnessKind::direct) {
        implied = w.genus;
      } else {
        int worst = 0;
        bool haveAll = !w.otherCandidates.empty();
        for (const std::string &c : w.otherCandidates) {
          int hi = cobordismgraph::NO_UPPER_BOUND;
          bool assisted = false;
          if (auto it = bounds.find(c);
              it != bounds.end() && it->second.haveUpper()) {
            hi = it->second.hi;
            assisted = !it->second.support.empty();
          } else if (const auto *ni = names.find(c); ni && ni->haveLiterature) {
            hi = ni->litHi;
            assisted = true;
          }
          if (hi == cobordismgraph::NO_UPPER_BOUND) {
            haveAll = false;
            break;
          }
          worst = std::max(worst, hi);
          impliedAssisted = impliedAssisted || assisted;
        }
        if (haveAll)
          implied = worst + w.genus + w.otherComponents - 1;
      }

      if (implied != cobordismgraph::NO_UPPER_BOUND && implied < row.lo) {
        std::ostringstream msg;
        msg << row.name << ": witness implies an upper bound of " << implied
            << ", BELOW the literature lower bound " << row.lo << ".";
        flagFatalBug(msg.str());
        noteStop("fatal-bug");
        e.requestStop();
        e.skipRemainingBoundaryProcessing();
        return;
      }

      if (implied != cobordismgraph::NO_UPPER_BOUND && implied <= row.lo &&
          !impliedAssisted) {
        // Only a CONSTRUCTIVE result settles a row. An assisted one is a
        // correct deduction but not an independent verification, and the
        // search might still find the surface that turns it into one -- so
        // it must not end the row. (Moot under --harvest, which never stops
        // early anyway; this matters for a non-harvest run.)
        resolvedThisRow.store(true, std::memory_order_relaxed);
        if (!harvest) {
          // Without --harvest, stop the moment the row is settled: the
          // rest of this cobordism's surfaces would only add edges we
          // aren't going to need. With it, keep going and bank them.
          noteStop("stopped");
          e.requestStop();
          e.skipRemainingBoundaryProcessing();
        }
      }
    };

    std::atomic<bool> searchDone{false};
    std::thread watchdog;
    if (perKnotTimeLimit || harvestQuiescence || sweepTimeLimit) {
      // Stops the DFS, and by default LETS THE BOUNDARY DRAIN FINISH.
      //
      // The time limit bounds the *search*, not the row. Identification is
      // the product: an unidentified surface says nothing whatever about a
      // slice genus, so a surface found and then discarded unexamined is
      // pure waste. Cutting the drain short measured 1% identification on a
      // run that produced 1,216,027 qualifying surfaces -- 1.2 million
      // boundaries thrown away unlooked-at, which is why that run's "found
      // nothing" meant nothing.
      //
      // The cost is that a row overruns its nominal limit by however long
      // the queue takes to drain, which can be minutes. That is the right
      // trade: waiting is cheaper than searching a region and then refusing
      // to look at what it found. --skip-drain-on-timeout restores the old
      // behaviour for when throughput genuinely matters more.
      auto endRowNow = [&](const char *why) {
        noteStop(why);
        e.requestStop();
        if (skipDrainOnTimeout)
          e.skipRemainingBoundaryProcessing();
      };
      watchdog = std::thread([&]() {
        auto knotDeadline =
            std::chrono::steady_clock::now() +
            std::chrono::duration<double>(perKnotTimeLimit.value_or(0));
        while (!searchDone.load(std::memory_order_relaxed)) {
          std::this_thread::sleep_for(std::chrono::milliseconds(200));
          if (searchDone.load(std::memory_order_relaxed))
            break;
          auto now = std::chrono::steady_clock::now();
          if (perKnotTimeLimit && now >= knotDeadline) {
            endRowNow("timeout");
            break;
          }
          if (sweepTimeLimit &&
              now - sweepStart >
                  std::chrono::duration<double>(*sweepTimeLimit)) {
            endRowNow("timeout");
            break;
          }
          // Quiescence: this row has stopped teaching us anything new, so
          // spending the rest of its budget enumerating more of the same
          // is worse than moving on to a row we know nothing about.
          // Meaningful during the boundary drain too, since that is where
          // witnesses are actually identified.
          if (harvestQuiescence) {
            long long idleMs =
                tickNow() - lastNewWitnessTick.load(std::memory_order_relaxed);
            if (idleMs > static_cast<long long>(*harvestQuiescence * 1000)) {
              endRowNow("quiescent");
              break;
            }
          }
        }
      });
    }

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
    BoundaryCondition cond;
    switch (boundaryConditionMode) {
    case BoundaryConditionMode::proper:
      cond = BoundaryCondition::proper;
      break;
    case BoundaryConditionMode::connected:
      // Honoured only where it is actually satisfiable: a multi-component
      // link can never meet `connected` on its own search side, so forcing
      // it there would just guarantee an empty search.
      cond = componentCount == 1 ? BoundaryCondition::connected
                                 : BoundaryCondition::proper;
      break;
    case BoundaryConditionMode::automatic:
    default:
      cond = componentCount == 1 ? BoundaryCondition::connected
                                 : BoundaryCondition::proper;
      break;
    }
    e.search(numThreads, cond, callbacks, iddfsIterations, iddfsStep,
             iddfsStart, iddfsFinalThreads, /*orientableOnly=*/true, maxFaces);

    searchDone.store(true, std::memory_order_relaxed);
    if (watchdog.joinable())
      watchdog.join();

    if (surfaceLog)
      surfaceLog->finalize();
    progressPrevLines_ = 0;
    ++searchedThisRun;

    // Record what this search actually cost before anything else, so even
    // a fatal halt below leaves the bookkeeping behind.
    {
      OutputRow &out = outputRows[row.name];
      out.knot = row.name;
      out.literatureLo = row.lo;
      out.literatureHi = row.hi;
      out.searchedFaces = maxFaces.value_or(0);
      out.searchOutcome = searchOutcome;
      if (out.status.empty())
        out.status = "unresolved";
    }

    std::vector<std::string> contradictions = resolveAll();
    {
      OutputRow &out = outputRows[row.name];
      out.searchedFaces = maxFaces.value_or(0);
      out.searchOutcome = searchOutcome;
    }

    writeWitnesses(cobordismsPath, witnesses);
    writeOutputCsv(*outputPath, rows, outputRows);

    for (const std::string &reason : contradictions)
      flagFatalBug(reason);
    if (fatalBugDetected_.load())
      haltIfFatalBugDetected();

    long long rejected = orientationRejections.load();
    std::cout << "[+] " << row.name << ": "
              << newWitnessesThisRow.load() << " new witnesses, outcome "
              << searchOutcome;
    if (rejected > 0)
      std::cout << ", " << rejected
                << " surfaces rejected on orientation mismatch";
    std::cout << "\n";

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

    if (censusUpdates) {
      // The complement cannot distinguish oriented variants, so inserting
      // an orientation-TAGGED name against a complement isoSig would make
      // identify() confidently return (say) L6a3{0} for L6a3{1} forever
      // after -- and insertCensusEntry is INSERT OR IGNORE, so the first
      // variant processed would own that isoSig permanently. Store the
      // base name; the candidate-set machinery restores the variants.
      auto &[t2, edges2, reversed2] = link;
      Link linkGrouping(t2, edges2);
      regina::Triangulation<3> complement = linkGrouping.buildComplement();
      census::insertCensusEntry(complement.isoSig(),
                                cobordismgraph::baseName(row.name));
    }

    ++processedThisRun;
  }

  // One last solve, so a run that searched nothing (or was cut short) still
  // reflects everything its witness file knows.
  for (const std::string &reason : resolveAll())
    flagFatalBug(reason);
  writeWitnesses(cobordismsPath, witnesses);
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
  return 0;
}
