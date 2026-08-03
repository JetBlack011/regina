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

constexpr const char *OUTPUT_HEADER =
    "knot,resolved_genus,status,witness_kind,witness_pairsig,via_knot,"
    "via_edge_genus,depends_on,literature_lo,literature_hi";

std::string formatOutputRow(const OutputRow &r) {
  std::ostringstream out;
  out << csvField(r.knot) << ',' << r.resolvedGenus << ',' << r.status << ','
      << r.witnessKind << ',' << csvField(r.witnessPairSig) << ','
      << csvField(r.viaKnot) << ',' << r.viaEdgeGenus << ','
      << csvField(r.dependsOn) << ',' << r.literatureLo << ','
      << r.literatureHi;
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
    result[r.knot] = std::move(r);
  }
  return result;
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

  std::unordered_map<std::string, OutputRow> outputRows =
      loadOutputCsv(*outputPath);
  std::unordered_map<std::string, int> knownGenus;
  knownGenus["Unknot"] = 0;
  for (const auto &[name, out] : outputRows)
    if (out.status == "resolved" || out.status == "range")
      knownGenus[name] = out.resolvedGenus;
  std::cout << "[+] Resuming with " << knownGenus.size() - 1
            << " already-resolved knots from " << *outputPath << "\n\n";

  Graph graph;

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

  size_t processedThisRun = 0;
  for (const auto &row : pending) {
    if (knownGenus.contains(row.name))
      continue;

    int target = row.hi;
    bool resolvedThisRow = false;

    if (auto info = resolvable(row.name, target, graph, knownGenus)) {
      knownGenus[row.name] = target;
      OutputRow out;
      out.knot = row.name;
      out.resolvedGenus = target;
      out.status = (row.lo == row.hi) ? "resolved" : "range";
      out.witnessKind = "propagated";
      out.viaKnot = info->viaKnot;
      out.viaEdgeGenus = info->viaEdgeGenus;
      out.dependsOn = buildDependsOn(info->viaKnot, outputRows);
      out.literatureLo = row.lo;
      out.literatureHi = row.hi;
      outputRows[row.name] = std::move(out);
      resolvedThisRow = true;
      std::cout << "[+] " << row.name
                << ": resolved via propagation (genus " << target << ")\n";
    } else {
      std::cout << "[+] Searching " << row.name << " (target genus "
                << target << ", " << row.crossings << " crossings)...\n";

      knotbuilder::PDCode pdcode;
      knotbuilder::TriangulationWithLink link;
      bool buildFailed = false;
      try {
        pdcode = knotbuilder::parsePDCode(row.pdNotation);
        link = knotbuilder::buildLink(pdcode);
      } catch (const regina::InvalidArgument &e) {
        std::cerr << "[!] " << row.name << ": failed to build from PD code ("
                  << e.what() << "), skipping\n";
        buildFailed = true;
      }

      if (!buildFailed) {
        auto &[t2, edges2] = link;

        // Splits edges2 into connected components -- 1 for an ordinary
        // knot, >1 for a genuine multi-component link. Reused below both
        // to pick the search's BoundaryCondition and to build the census
        // complement.
        Link linkGrouping(t2, edges2);
        int componentCount = linkGrouping.countComponents();

        std::vector<int> edgeIndices;
        edgeIndices.reserve(edges2.size());
        for (const regina::Edge<3> *e : edges2)
          edgeIndices.push_back(static_cast<int>(e->index()));

        CobordismBuilder<3> cob(t2);
        CollarBuilder collarBuilder(edgeIndices);
        for (int i = 0; i < thickenLayers; ++i) {
          cob.thicken();
          if (i < collarLayers)
            collarBuilder.addLayer(cob);
        }
        if (useCone)
          cob.cone();

        // Must be read before getCobordism() is copied into `tri` below:
        // baseBoundaryComponent() resolves against cob's own internal
        // triangulation, whose boundary-component indices are preserved
        // across the copy (same "indices are preserved" guarantee
        // CobordismBuilder::baseTriangulation()'s own doc comment relies
        // on elsewhere in this file).
        size_t knotSideBC = cob.baseBoundaryComponent()->index();

        regina::Triangulation<4> tri = cob.getCobordism();

        std::vector<int> seedFaces;
        if (collarLayers > 0)
          for (regina::Triangle<4> *t : collarBuilder.resolve())
            seedFaces.push_back(static_cast<int>(t->index()));

        std::optional<SurfaceSearch> eOpt;
        if (seedFaces.empty())
          eOpt.emplace(tri);
        else
          eOpt.emplace(tri, seedFaces);
        SurfaceSearch &e = *eOpt;
        e.configureLimits(limits);

        std::optional<CsvWriter> surfaceLog;
        if (surfaceLogPath)
          surfaceLog.emplace(*surfaceLogPath,
                             "orientable,genus,punctures,triangles,pairsig",
                             numThreads);

        SurfaceSearchCallbacks callbacks;
        callbacks.onProgress = [&](const SearchStats &stats) {
          printProgress(stats, e);
        };
        callbacks.onBoundaryProcessingStarted = [&](size_t total,
                                                    unsigned threads) {
          // Commits (rather than redraws over) whatever progress block the
          // DFS phase left on screen, so the transition into this phase is
          // visible rather than silently overwritten by the first tick.
          progressPrevLines_ = 0;
          std::cerr << "[+] boundary processing: " << total
                    << " queued surfaces, " << threads << " threads\n";
        };
        callbacks.onBoundaryProcessingProgress =
            [&](size_t processed, size_t total,
                std::chrono::steady_clock::duration elapsed) {
              printBoundaryProgress(
                  processed, total, elapsed,
                  resolvedThisRow ? std::optional<int>(target) : std::nullopt,
                  target);
            };
        callbacks.onBoundaryProcessingComplete =
            [&](size_t total, std::chrono::steady_clock::duration elapsed) {
              progressPrevLines_ = 0;
              std::cerr << "[+] boundary processing: done (" << total
                        << " processed in " << formatElapsed(elapsed)
                        << ")\n";
            };
        callbacks.onSurfaceBoundaryProcessed =
            [&](const SurfaceBoundaryInfo &info) {
              if (surfaceLog) {
                std::string pairSig =
                    info.capturePairSig ? info.capturePairSig() : std::string{};
                std::ostringstream row2;
                row2 << (info.orientable ? "true" : "false") << ','
                    << info.genus << ',' << info.punctures << ','
                    << info.triangleCount << ',' << csvField(pairSig);
                surfaceLog->writeRow(row2.str());
              }

              if (!info.orientable)
                return; // defensive; orientableOnly=true already prunes these
              if (!info.connected)
                return; // a disconnected find isn't a valid single witness
                        // surface -- info.genus has no real meaning for it
                        // (see SurfaceFoundInfo::genus's own doc comment).
                        // Most likely with componentCount > 1, since a
                        // multi-component link's seeded collar starts out
                        // as one disjoint piece per component and nothing
                        // requires the search to ever bridge them.

              auto capturePairSig = [&info] {
                return info.capturePairSig ? info.capturePairSig()
                                           : std::string{};
              };
              auto handleOutcome = [&](const WitnessOutcome &outcome) {
                if (outcome.fatalBug) {
                  flagFatalBug(outcome.fatalBugMessage);
                  e.requestStop();
                  e.skipRemainingBoundaryProcessing();
                }
              };
              // Like handleOutcome, but also announces a genuine
              // resolution -- unlike row.name's own resolution (announced
              // separately, after e.search() returns, once this row's
              // search has actually stopped), an incidental resolution has
              // no other announcement point: this callback is the only
              // place that ever learns about it.
              auto handleIncidentalOutcome =
                  [&](const std::string &name, const WitnessOutcome &outcome) {
                    handleOutcome(outcome);
                    if (outcome.resolved)
                      std::cout << "\x1b[1;36m[+] " << name
                                << ": resolved incidentally (genus "
                                << knownGenus[name] << ") while searching "
                                << row.name << "\x1b[0m\n";
                  };

              BoundarySplit split =
                  splitBoundary(info.boundaryComponents, knotSideBC);

              if (split.mineCurveCount == static_cast<size_t>(componentCount)) {
                // This row's own diagram is fully witnessed -- the
                // original reason this search is running.
                if (!resolvedThisRow) {
                  WitnessOutcome outcome;
                  if (split.otherSides.empty()) {
                    // Direct witness: generalizes the old punctures==1
                    // knot-only case to any component count.
                    outcome = recordWitness(row.name, row.lo, row.hi, target,
                                            "", info.genus, capturePairSig,
                                            "direct", graph, knownGenus,
                                            outputRows);
                  } else if (split.otherSides.size() == 1 &&
                            split.otherSides.front().safe) {
                    // Cobordism witness: generalizes the old punctures==2
                    // knot-vs-knot case to fire for any component count on
                    // either side. An unsafe (genuinely linked multi-
                    // component) far side is skipped -- it's still shown
                    // in boundaryDescription for a human to read, just
                    // never used for a genus deduction (see
                    // identify::isOrientationSafeName()).
                    outcome = recordWitness(
                        row.name, row.lo, row.hi, target,
                        split.otherSides.front().name, info.genus,
                        capturePairSig, "cobordism", graph, knownGenus,
                        outputRows);
                  }
                  // More than one other side: an unhandled multi-way
                  // cobordism -- skipped, not guessed at.
                  handleOutcome(outcome);
                  if (outcome.resolved) {
                    resolvedThisRow = true;
                    e.requestStop();
                    e.skipRemainingBoundaryProcessing();
                  }
                }
                return;
              }

              if (split.mineCurveCount != 0)
                return; // partial "mine" presence isn't cleanly nameable

              // Incidental observation: this row's own diagram is
              // COMPLETELY absent from this surface's boundary -- a side
              // effect of searching row.name's own cobordism, unrelated to
              // it. describeBoundary_() already paid for identifying every
              // side regardless of whether "mine" is present, so this is
              // free information that would otherwise just be discarded;
              // recording it here can let a LATER row's own search be
              // skipped entirely (see the `if (knownGenus.contains(...))`
              // check above the main loop, and propagateGraph(), which
              // runs after every row and will pick up anything recorded
              // here). Never calls e.requestStop() here except on a fatal
              // bug -- this row's own search keeps running toward its own
              // goal regardless of what gets harvested along the way.
              if (split.otherSides.size() == 1 &&
                  split.otherSides.front().safe) {
                const std::string &y = split.otherSides.front().name;
                auto it = nameToRow.find(y);
                if (it != nameToRow.end()) {
                  const InputRow &yRow = *it->second;
                  handleIncidentalOutcome(
                      y, recordWitness(y, yRow.lo, yRow.hi, yRow.hi, "",
                                       info.genus, capturePairSig,
                                       "incidental-direct", graph, knownGenus,
                                       outputRows));
                }
              } else if (split.otherSides.size() == 2 &&
                        split.otherSides[0].safe &&
                        split.otherSides[1].safe) {
                const std::string &y = split.otherSides[0].name;
                const std::string &z = split.otherSides[1].name;

                // Tries to resolve `subject` via a cobordism to `via`,
                // only if `subject` is actually a tracked --input row.
                // addEdge()'s own dedup (hasEdgeWithGenus(), inside
                // recordWitness()) makes calling this for both directions
                // a cheap no-op on the second call once the first has
                // already recorded the edge.
                auto tryDirection = [&](const std::string &subject,
                                        const std::string &via) {
                  auto it = nameToRow.find(subject);
                  if (it == nameToRow.end())
                    return false;
                  const InputRow &subjectRow = *it->second;
                  WitnessOutcome outcome = recordWitness(
                      subject, subjectRow.lo, subjectRow.hi, subjectRow.hi,
                      via, info.genus, capturePairSig, "incidental-cobordism",
                      graph, knownGenus, outputRows);
                  handleIncidentalOutcome(subject, outcome);
                  return outcome.fatalBug;
                };
                if (!tryDirection(y, z))
                  tryDirection(z, y);
              }
            };

        std::atomic<bool> searchDone{false};
        std::thread watchdog;
        if (perKnotTimeLimit) {
          watchdog = std::thread([&]() {
            auto deadline = std::chrono::steady_clock::now() +
                            std::chrono::duration<double>(*perKnotTimeLimit);
            while (!searchDone.load(std::memory_order_relaxed) &&
                  std::chrono::steady_clock::now() < deadline)
              std::this_thread::sleep_for(std::chrono::milliseconds(200));
            if (!searchDone.load(std::memory_order_relaxed))
              e.requestStop();
          });
        }

        // A multi-component link's own boundary necessarily puts more than
        // one of the surface's own boundary curves on the single knot-side
        // ambient boundary component (see knotSideBC above) -- impossible
        // under `connected`'s distinct-ambient-component-per-curve
        // requirement, but exactly what `proper` allows. `connected` stays
        // the (tighter-pruning) default for ordinary single-component
        // knots.
        BoundaryCondition cond = componentCount == 1
                                     ? BoundaryCondition::connected
                                     : BoundaryCondition::proper;
        e.search(numThreads, cond, callbacks,
                iddfsIterations, iddfsStep, iddfsStart, iddfsFinalThreads,
                /*orientableOnly=*/true);

        searchDone.store(true, std::memory_order_relaxed);
        if (watchdog.joinable())
          watchdog.join();
        if (fatalBugDetected_.load()) {
          // Durably record whatever legitimate results already exist
          // before halting -- the outer loop's own writeOutputCsv() call
          // at the end of this row never runs, since
          // haltIfFatalBugDetected() exits the process below.
          writeOutputCsv(*outputPath, rows, outputRows);
          haltIfFatalBugDetected();
        }
        // Leave this knot's final progress block in place (rather than
        // erasing it on the next tick, which won't come until the *next*
        // knot's search starts) so the subsequent resolved/unresolved
        // line prints cleanly underneath it.
        progressPrevLines_ = 0;

        if (surfaceLog)
          surfaceLog->finalize();

        if (!resolvedThisRow) {
          OutputRow out;
          out.knot = row.name;
          out.resolvedGenus = target;
          out.status = "unresolved";
          out.witnessKind = "none";
          out.literatureLo = row.lo;
          out.literatureHi = row.hi;
          outputRows[row.name] = std::move(out);
          std::cout << "[+] " << row.name << ": unresolved this run\n";
        } else {
          std::cout << "[+] " << row.name << ": resolved (genus " << target
                    << ")\n";
        }

        if (censusUpdates) {
          regina::Triangulation<3> complement = linkGrouping.buildComplement();
          census::insertCensusEntry(complement.isoSig(), row.name);
        }
      } else {
        OutputRow out;
        out.knot = row.name;
        out.status = "unresolved";
        out.witnessKind = "none";
        out.literatureLo = row.lo;
        out.literatureHi = row.hi;
        outputRows[row.name] = std::move(out);
      }
    }

    for (const std::string &name :
        propagateGraph(pending, graph, knownGenus, outputRows))
      std::cout << "\x1b[1;36m[+] " << name
                << ": resolved incidentally (genus " << knownGenus[name]
                << ") via graph propagation, no search needed\x1b[0m\n";
    writeOutputCsv(*outputPath, rows, outputRows);
    ++processedThisRun;
  }

  std::cout << "\n[+] Done. Processed " << processedThisRun
            << " knots this run.\n";
  return 0;
}
