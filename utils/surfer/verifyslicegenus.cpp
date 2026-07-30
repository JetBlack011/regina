//
//  verifyslicegenus.cpp
//
//  Created by John Teague on 07/29/2026.
//

#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
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
#include "collar.h"
#include "csvwriter.h"
#include "embeddingsearch.h"
#include "surfacesearch.h"
#include "knotbuilder.h"
#include "linkcomplement.h"
#include "identifycomplement.h"

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
// `resolved` reflects the driver's own resolvedThisRow for the knot
// currently being processed, so it's visible whether the goal (matching
// the literature slice genus) has already been met even while this phase
// is still churning through whatever else was queued.
void printBoundaryProgress(size_t processed, size_t total,
                           std::chrono::steady_clock::duration elapsed,
                           bool resolved) {
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
  report << "[+] goal (literature slice genus) "
         << (resolved ? "ACHIEVED" : "not yet achieved") << "\n";
  redrawProgressBlock(report.str());
}

// ─────────────────────────────────────────────────────────────────────────
// CSV parsing helpers
// ─────────────────────────────────────────────────────────────────────────

// One row of the input knot table (Name,PD Notation,Genus-4D). No RFC-4180
// quoting appears in that file (PD Notation uses ';' internally, never a
// literal comma), so a naive two-comma split suffices -- see the input
// table's own format, confirmed during design.
struct InputRow {
  std::string name;
  std::string pdNotation;
  int lo = 0, hi = 0; // literature genus bounds; lo == hi except for a
                      // handful of [lo;hi] range rows
  int crossings = 0;
};

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

// The leading digit run of `name` (e.g. "13n_1109" -> 13, "8_6" -> 8).
int crossingsFromName(const std::string &name) {
  size_t i = 0;
  while (i < name.size() &&
        std::isdigit(static_cast<unsigned char>(name[i])))
    ++i;
  return std::stoi(name.substr(0, i));
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
    row.crossings = crossingsFromName(name);
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

struct OutputRow {
  std::string knot;
  int resolvedGenus = 0;
  std::string status;       // resolved | range | unresolved | skipped
  std::string witnessKind;  // direct | cobordism | propagated | none
  std::string witnessPairSig;
  std::string viaKnot;
  int viaEdgeGenus = 0;
  std::string dependsOn;
  int literatureLo = 0;
  int literatureHi = 0;
};

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

// Rewrites the whole --output file from `outputRows`, in `rows`' original
// order (skipping any row not yet processed), via write-to-temp +
// atomic rename -- so a crash mid-write never corrupts the previous,
// already-durable version. Called once per knot processed (resolved,
// attempted-but-unresolved, skipped, or newly propagated).
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
    for (const auto &row : rows) {
      auto it = outputRows.find(row.name);
      if (it == outputRows.end())
        continue; // not yet processed
      out << formatOutputRow(it->second) << "\n";
    }
  }
  std::filesystem::rename(tmp, path);
}

// ─────────────────────────────────────────────────────────────────────────
// The genus-implication graph
// ─────────────────────────────────────────────────────────────────────────

struct GraphEdge {
  std::string other;
  int genus;
  std::string pairSig;
};

using Graph = std::unordered_map<std::string, std::vector<GraphEdge>>;

void addEdge(Graph &graph, const std::string &a, const std::string &b,
            int genus, const std::string &pairSig) {
  graph[a].push_back({b, genus, pairSig});
  graph[b].push_back({a, genus, pairSig});
}

// Whether `graph` already has an (a, b) edge at exactly `genus` -- see the
// caller for why exact-genus (not "any genus <= this one") is the right
// dedup key: resolvable() checks an exact equality (target == g + h or
// target == h - g), not an inequality, so two edges to the same neighbor
// with different genus are NOT interchangeable, but two with the SAME
// genus (however many redundant surfaces happened to witness it) are.
bool hasEdgeWithGenus(const Graph &graph, const std::string &a,
                      const std::string &b, int genus) {
  auto it = graph.find(a);
  if (it == graph.end())
    return false;
  for (const auto &e : it->second)
    if (e.other == b && e.genus == genus)
      return true;
  return false;
}

struct ResolveInfo {
  std::string viaKnot;
  int viaEdgeGenus;
  std::string pairSig; // the resolving edge's own pairSig, if any
};

// Returns the edge (if any) that pins `name`'s genus to exactly `target`,
// given `knownGenus`'s currently-established values: an edge (name, other,
// g) with other's genus known as h pins target only when target == g + h
// (constructive direction) or target == h - g (contrapositive direction) --
// merely falling inside [h-g, h+g] is consistent, not conclusive. See
// identifycomplement.h-adjacent design notes (the plan document) for the
// derivation.
std::optional<ResolveInfo>
resolvable(const std::string &name, int target, const Graph &graph,
          const std::unordered_map<std::string, int> &knownGenus) {
  auto it = graph.find(name);
  if (it == graph.end())
    return std::nullopt;
  for (const auto &e : it->second) {
    auto hIt = knownGenus.find(e.other);
    if (hIt == knownGenus.end())
      continue;
    int h = hIt->second;
    if (target == e.genus + h || target == h - e.genus)
      return ResolveInfo{e.other, e.genus, e.pairSig};
  }
  return std::nullopt;
}

// Walks from `viaKnot` back through outputRows' own via_knot chain to
// whatever ultimately bottomed out (a direct witness, or the Unknot, which
// has no outputRows entry of its own), joining the path with ';'.
std::string buildDependsOn(
    const std::string &viaKnot,
    const std::unordered_map<std::string, OutputRow> &outputRows) {
  std::vector<std::string> chain;
  std::unordered_set<std::string> seen;
  std::string cur = viaKnot;
  while (!cur.empty() && seen.insert(cur).second) {
    chain.push_back(cur);
    if (cur == "Unknot")
      break;
    auto it = outputRows.find(cur);
    if (it == outputRows.end() || it->second.viaKnot.empty())
      break;
    cur = it->second.viaKnot;
  }
  std::ostringstream out;
  for (size_t i = 0; i < chain.size(); ++i) {
    if (i)
      out << ';';
    out << chain[i];
  }
  return out.str();
}

// Sweeps every still-pending row for new resolutions unlockable purely from
// the current knownGenus/graph state, repeating until a full pass makes no
// further progress. Each newly-resolved row is recorded with witness_kind
// "propagated" (no search run for it this session).
void propagateGraph(const std::vector<InputRow> &rows, const Graph &graph,
                    std::unordered_map<std::string, int> &knownGenus,
                    std::unordered_map<std::string, OutputRow> &outputRows) {
  bool changed = true;
  while (changed) {
    changed = false;
    for (const auto &row : rows) {
      if (knownGenus.contains(row.name))
        continue;
      int target = row.hi;
      auto info = resolvable(row.name, target, graph, knownGenus);
      if (!info)
        continue;

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
      changed = true;
    }
  }
}

// ─────────────────────────────────────────────────────────────────────────
// Boundary-component identification (--no-cone only)
// ─────────────────────────────────────────────────────────────────────────

// Splits a SurfaceBoundaryInfo::boundaryDescription string into
// {ambient-component-number -> curve name}. Only valid when the search ran
// under BoundaryCondition::connected AND the surface has exactly 2 of its
// own boundary curves across exactly 2 ambient boundary components (our
// --no-cone, punctures==2 case) -- connected's "each surface boundary
// component maps to a distinct ambient one" guarantees each ambient
// component's segment is then exactly one curve name, with no comma inside
// it and no "(<link>)" parenthetical (that only appears when a component
// holds more than one curve), so a top-level split on ", " is safe.
std::unordered_map<int, std::string>
parseTwoComponentBoundary(const std::string &desc) {
  std::unordered_map<int, std::string> result;
  size_t start = 0;
  while (start <= desc.size()) {
    size_t comma = desc.find(", ", start);
    std::string piece = desc.substr(
        start, comma == std::string::npos ? std::string::npos
                                          : comma - start);
    size_t colon = piece.find(": ");
    if (colon != std::string::npos) {
      int comp = std::stoi(piece.substr(0, colon));
      result[comp] = piece.substr(colon + 2);
    }
    if (comma == std::string::npos)
      break;
    start = comma + 2;
  }
  return result;
}

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
         "                     every knot processed; already-resolved rows "
         "from a\n"
         "                     prior run are trusted directly and not "
         "re-searched.\n";
  std::cerr
      << "    --input <csv>  : The knot table (Name,PD Notation,Genus-4D). "
         "Default:\n"
         "                     4d_smooth_slice_genus_13_crossings_pd_codes."
         "csv\n"
         "                     alongside the source tree.\n";
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
              printBoundaryProgress(processed, total, elapsed,
                                   resolvedThisRow);
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
              if (resolvedThisRow)
                return;

              if (info.punctures == 1) {
                if (info.genus == target) {
                  knownGenus[row.name] = target;
                  OutputRow out;
                  out.knot = row.name;
                  out.resolvedGenus = target;
                  out.status = (row.lo == row.hi) ? "resolved" : "range";
                  out.witnessKind = "direct";
                  out.witnessPairSig =
                      info.capturePairSig ? info.capturePairSig() : std::string{};
                  out.literatureLo = row.lo;
                  out.literatureHi = row.hi;
                  outputRows[row.name] = std::move(out);
                  resolvedThisRow = true;
                  e.requestStop();
                  e.skipRemainingBoundaryProcessing();
                } else if (info.genus < row.lo) {
                  std::cerr
                      << "[!] " << row.name
                      << ": found a direct genus-" << info.genus
                      << " surface, BELOW the literature lower bound "
                      << row.lo
                      << " -- possible pipeline bug or literature error\n";
                }
              } else if (info.punctures == 2) {
                auto parts = parseTwoComponentBoundary(info.boundaryDescription);
                std::string farName;
                for (const auto &[comp, nm] : parts)
                  if (comp != static_cast<int>(knotSideBC) + 1)
                    farName = nm;
                if (farName.empty())
                  return;

                // Skip the (expensive -- see capturePairSig's own doc
                // comment) pairSig capture entirely when it wouldn't teach
                // the graph anything new: resolvable()'s check is an exact
                // equality on genus, so two edges to the same neighbor at
                // the SAME genus are fully interchangeable (whichever
                // witness got recorded first is just as good), and with
                // thousands of near-duplicate surfaces per knot pair, this
                // is the overwhelmingly common case.
                bool haveThisGenusEdge =
                    hasEdgeWithGenus(graph, row.name, farName, info.genus);

                auto hIt = knownGenus.find(farName);
                bool resolvesNow =
                    hIt != knownGenus.end() &&
                    (target == info.genus + hIt->second ||
                     target == hIt->second - info.genus);

                if (haveThisGenusEdge && !resolvesNow)
                  return;

                std::string pairSig =
                    info.capturePairSig ? info.capturePairSig() : std::string{};
                if (!haveThisGenusEdge)
                  addEdge(graph, row.name, farName, info.genus, pairSig);

                if (resolvesNow) {
                  knownGenus[row.name] = target;
                  OutputRow out;
                  out.knot = row.name;
                  out.resolvedGenus = target;
                  out.status = (row.lo == row.hi) ? "resolved" : "range";
                  out.witnessKind = "cobordism";
                  out.witnessPairSig = pairSig;
                  out.viaKnot = farName;
                  out.viaEdgeGenus = info.genus;
                  out.dependsOn = buildDependsOn(farName, outputRows);
                  out.literatureLo = row.lo;
                  out.literatureHi = row.hi;
                  outputRows[row.name] = std::move(out);
                  resolvedThisRow = true;
                  e.requestStop();
                  e.skipRemainingBoundaryProcessing();
                }
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

        e.search(numThreads, BoundaryCondition::connected, callbacks,
                iddfsIterations, iddfsStep, iddfsStart, iddfsFinalThreads,
                /*orientableOnly=*/true);

        searchDone.store(true, std::memory_order_relaxed);
        if (watchdog.joinable())
          watchdog.join();
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
          EdgeComplement complementBuilder(t2, edges2);
          regina::Triangulation<3> complement =
              complementBuilder.buildComplement();
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

    propagateGraph(pending, graph, knownGenus, outputRows);
    writeOutputCsv(*outputPath, rows, outputRows);
    ++processedThisRun;
  }

  std::cout << "\n[+] Done. Processed " << processedThisRun
            << " knots this run.\n";
  return 0;
}
