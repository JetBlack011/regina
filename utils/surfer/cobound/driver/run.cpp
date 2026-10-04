//
//  run.cpp
//
//  cobound run: a search from each target, scheduled by the goal. A run
//  without a goal (verifyslicegenus's sweep, created by John Teague on
//  07/29/2026) searches each target once.
//

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "cobound/bounds/axioms.h"
#include "cobound/bounds/searchjudge.h"
#include "cobound/cobordisms/database.h"
#include "cobound/cobordisms/pending.h"
#include "cobound/driver/commands.h"
#include "cobound/driver/config.h"
#include "cobound/driver/fatal.h"
#include "cobound/driver/scheduler.h"
#include "cobound/driver/setup.h"
#include "cobound/driver/signals.h"
#include "cobound/driver/targets.h"
#include "cobound/search/search.h"
#include "cobound/search/searchreport.h"
#include "cobound/solver/literature.h"
#include "cobound/solver/verdicts.h"
#include "cobound/frozen.h"
#include "linknaming/names.h"
#include "linknaming/tables.h"
#include "surfer/report/atomicwrite.h"

namespace {

using solver::InputRow;
using search::BoundaryConditionMode;
using verdicts::OutputRow;

BoundaryConditionMode boundaryCondition(const std::string &mode) {
  if (mode == "auto") return BoundaryConditionMode::automatic;
  if (mode == "connected") return BoundaryConditionMode::connected;
  return BoundaryConditionMode::proper;
}

/**
 * A run without a goal (context run, depth 0): each target -- the rows of
 * `targets` within max_crossings, in crossing order, or the one diagram
 * `target_pd` -- searched once, each by the one search (search/search.h),
 * judged by its own cobordism graph (SearchJudge), its finds written to its
 * pending file (<work>/hop_<k>_n0/kept.csv) and all of them signed into
 * `cobordisms` at the run's end. Each target's search record goes to its
 * row of the verdicts file; its status and bounds are `solve`'s (plan divergence
 * 10). The log lines are verifyslicegenus's, every one (frozen).
 *
 * Returns the exit code: 0 (a first SIGINT or SIGTERM included: the run
 * searches no further target, records `interrupted` and ends as usual), 1 when
 * its tables or symmetry types cannot be read, or 2 when any search failed
 * its accounting (after the run's other targets) or the run halted (fatal.h).
 * \throws config::Error for a configuration a search cannot run with.
 */
int runWithoutGoal(const config::Config &cfg) {
  const std::string outputPath = cfg.text("verdicts");
  const std::string cobordismsPath = cfg.text("cobordisms");
  const std::string knotTablePath = cfg.text("knot_table");
  const std::string linkTablePath = cfg.text("link_table");
  const std::string knotSymmetryPath = cfg.text("knot_symmetry");
  const std::optional<std::string> targetPD = cfg.optionalText("target_pd");
  // The rows to search: the table's, or the one diagram.
  const std::string inputPath = targetPD ? std::string() : cfg.text("targets");
  const std::optional<long long> maxFaces = cfg.optionalInteger("max_faces");
  const int maxCrossings = static_cast<int>(cfg.integer("max_crossings"));
  const std::optional<double> perKnotTimeLimit = cfg.optionalReal("search_seconds");
  // Stop the SEARCH once this many boundary-satisfying surfaces exist.
  //
  // A wall-clock limit equalises WALL TIME, which only equalises coverage if
  // every host explores roots at the same rate. Measured over 52 searches,
  // they do not: at an identical 3600 thread-seconds one host yielded a median
  // 972,655 qualifying surfaces per search and another 295,178 -- a 3.3x gap,
  // tight on both sides (+/-15%), and a property of the machine rather than of
  // the link searched. Targeting the surface count instead equalises the thing a
  // negative actually rests on: how much of the root ordering was covered.
  const std::optional<long long> surfaceTarget = cfg.optionalInteger("surface_target");
  const std::optional<std::string> surfaceLogPath = cfg.optionalText("surface_log");
  const std::optional<std::string> surfaceStatsPath = cfg.optionalText("surface_stats");
  const std::optional<std::string> frontierDir = cfg.optionalText("frontier_dir");
  const std::optional<std::string> pairSigCacheDir = cfg.optionalText("pair_sig_cache");
  const std::optional<std::string> resumeFrontierDir = cfg.optionalText("resume_frontier_dir");
  const std::optional<std::string> selfIntersectionCensusPath =
      cfg.optionalText("self_intersection_census");
  const std::optional<std::string> rejectionSampleLogPath =
      cfg.optionalText("rejection_sample_log");
  // Accept surfaces whose only self-intersections are unlinked (paper §4.5,
  // KnottedSurface::isResolvable()). No default (plan divergence 3).
  const bool resolveUnlinked = cfg.flag("resolve_unlinked");
  const bool outgoingNames = cfg.flag("outgoing_names");
  const std::string workDir = cfg.text("work");
  const unsigned numThreads = cfg.threads();
  // Thickened and collared through every layer (divergence 6: the cobordism
  // graph reads each find on its thickening as a goal run reads a stored
  // cobordism's).
  const int thickenLayers = static_cast<int>(cfg.integer("layers"));
  const BoundaryConditionMode boundaryConditionMode =
      boundaryCondition(cfg.text("boundary_condition"));
  const unsigned iddfsIterations = static_cast<unsigned>(cfg.integer("iddfs_iterations"));
  const long long iddfsStep = cfg.integer("iddfs_step");
  const std::optional<long long> iddfsStart = cfg.optionalInteger("iddfs_start");
  const long long rootBudgetStart = cfg.integer("root_budget_start");
  const long long rootBudgetGrowth = cfg.integer("root_budget_growth");

  SurfaceSearchLimits limits;
  limits.capturePairSig = true;
  limits.pendingSurfaceCap = static_cast<size_t>(cfg.integer("pending_surface_cap"));
  limits.petalCacheLimit = static_cast<size_t>(cfg.integer("petal_cache_limit"));
  limits.boundarySignatureCacheLimit =
      static_cast<size_t>(cfg.integer("boundary_signature_cache_limit"));
  // Only a multi-curve component's curve COUNT is ever used here, so naming
  // each of its curves separately is work nobody reads.
  limits.nameLinkCurves = false;

  auto refuse = [](const std::string &why) { throw config::Error(why); };
  if (thickenLayers < 1) refuse("layers must be at least 1");
  if (iddfsIterations > 0 && iddfsStep <= 0)
    refuse("iddfs_iterations > 0 requires iddfs_step > 0");
  if (iddfsStart && *iddfsStart <= 0) refuse("iddfs_start must be > 0");
  if (limits.pendingSurfaceCap == 0) refuse("pending_surface_cap must be > 0");
  if (limits.petalCacheLimit == 0) refuse("petal_cache_limit must be > 0");
  if (cfg.integer("complement_cache_limit") <= 0) refuse("complement_cache_limit must be > 0");
  if (limits.boundarySignatureCacheLimit == 0)
    refuse("boundary_signature_cache_limit must be > 0");
  if (maxFaces && *maxFaces <= 0) refuse("max_faces must be > 0");
  if (surfaceTarget && *surfaceTarget <= 0) refuse("surface_target must be > 0");
  if (perKnotTimeLimit && *perKnotTimeLimit <= 0) refuse("search_seconds must be > 0");
  if (rootBudgetGrowth < 2) refuse("root_budget_growth must be at least 2");
  if (cfg.integer("retriangulate_time_budget") <= 0)
    refuse("retriangulate_time_budget must be > 0");
  if (!maxFaces && !perKnotTimeLimit && !surfaceTarget)
    refuse("a search needs a stopping rule: without max_faces, surface_target or "
           "search_seconds nothing ends it (every search harvests, so resolving the row "
           "does not) and the run will never reach its second row");

  const bool censusLoaded = setup::applyRunSettings(cfg);

  std::optional<std::ofstream> rejectionSampleLog;
  if (rejectionSampleLogPath) {
    std::error_code ec;
    const bool fresh = !std::filesystem::exists(*rejectionSampleLogPath) ||
                       std::filesystem::file_size(*rejectionSampleLogPath, ec) == 0;
    rejectionSampleLog.emplace(*rejectionSampleLogPath, std::ios::app);
    if (!*rejectionSampleLog) refuse("rejection_sample_log: cannot open " + *rejectionSampleLogPath);
    if (fresh) *rejectionSampleLog << kFrozenRejectionSampleHeader;
  }

  // One policy for SIGINT and SIGTERM (divergence 8): the first ends the
  // current target's search cleanly and the run searches no further target; a
  // second ends the process.
  runsignals::install();

  std::cout << "------ verifyslicegenus \U0001F30A ------\n\n";
  std::cout << (censusLoaded ? "[+] census: loaded from " : "[+] census: not found at ")
            << cfg.text("census") << (censusLoaded ? "\n\n" : ", skipping\n\n");

  std::vector<InputRow> rows;
  if (targetPD)
    rows = {targets::oneDiagram(cfg.text("target_name"), *targetPD,
                                {knotTablePath, linkTablePath})};
  else
    rows = targets::loadInputCsv(inputPath);
  std::cout << "[+] Loaded " << rows.size() << " knots from "
            << (targetPD ? "target_pd" : inputPath) << "\n";

  // Literature metadata for every name that could ever appear in the
  // graph, not just this run's targets. A knots run needs the links table
  // to expand an orientation-blind outgoing name like "L6a3" into the
  // oriented variants it might be (see NameTable::candidates), and a links
  // run needs the knots table for the same reason in reverse -- so both
  // are always loaded, regardless of which one is being searched.
  solver::NameTable names;
  size_t metadataRows = 0;
  for (const auto &row : rows)
    names.addLiterature(row.name, row.lo, row.hi);
  for (const std::string &table : {knotTablePath, linkTablePath}) {
    if (table.empty() || table == inputPath)
      continue;
    try {
      metadataRows += solver::loadNameTable(table, names);
    } catch (const std::exception &e) {
      std::cerr << "[!] could not load name table " << table << ": " << e.what()
                << " (continuing without it)\n";
    }
  }
  std::cout << "[+] Name table: " << names.size() << " names (" << metadataRows
            << " from tables other than --input)\n";

  // The tables, loaded once: the tables, which a search's own
  // cobordism graph names its links by (divergence 6) and its outgoing links
  // are named by (outgoing_names), and diagram naming's
  // signature table, drawn from them.
  std::optional<linknaming::Tables> outgoingTables;
  if (outgoingNames) {
    const auto t0 = std::chrono::steady_clock::now();
    try {
      outgoingTables =
          linknaming::Tables::load(knotTablePath, linkTablePath, knotSymmetryPath);
      std::cout << kFrozenExactFarSideNamesLine << outgoingTables->size() << " table entries ("
                << std::chrono::duration_cast<std::chrono::milliseconds>(
                       std::chrono::steady_clock::now() - t0)
                       .count()
                << " ms)\n";
    } catch (const std::exception &e) {
      std::cerr << kFrozenExactFarSideNamesOff << e.what() << "\n";
    }
  }
  std::optional<linknaming::Tables> graphTablesOwn;
  const linknaming::Tables *graphTables = outgoingTables ? &*outgoingTables : nullptr;
  std::optional<linknaming::LinkNamer> graphNamer;
  try {
    if (!graphTables)
      graphTables = &graphTablesOwn.emplace(
          linknaming::Tables::load(knotTablePath, linkTablePath, knotSymmetryPath));
    graphNamer.emplace(*graphTables, bounds::LinkAxioms::namerLimits());
  } catch (const std::exception &e) {
    std::cerr << "[!] the cobordism graph cannot load the tables: " << e.what() << "\n";
    return 1;
  }
  std::optional<linknaming::SignatureTable> signatureTable;
  {
    const auto t0 = std::chrono::steady_clock::now();
    try {
      // From the tables the graph already loaded: one table load (phase 5).
      signatureTable = linknaming::SignatureTable::fromTables(*graphTables);
      std::cout << "[+] diagram naming: " << signatureTable->knots() << " knot and "
                << signatureTable->links() << " link diagram signatures ("
                << std::chrono::duration_cast<std::chrono::milliseconds>(
                       std::chrono::steady_clock::now() - t0)
                       .count()
                << " ms)\n";
    } catch (const std::exception &ex) {
      std::cerr << "[!] diagram naming off: could not load the tables' signatures ("
                << ex.what() << ")\n";
    }
  }

  std::unordered_map<std::string, OutputRow> outputRows = verdicts::loadOutputCsv(outputPath);

  // The database the run signs into, loaded: a search keeps no cobordism
  // whose identity it holds (one per identity across the database).
  cobordisms::LoadedDatabase recorded(cobordismsPath,
                                      cobordisms::loadCobordisms(cobordismsPath, false));
  const std::vector<cobordisms::Cobordism> &cobordisms = recorded.all();
  std::cout << "[+] Resuming with " << cobordisms.size()
            << " previously-recorded witnesses from " << cobordismsPath << "\n";

  // The run directory (plan divergence 7): each search's pending cobordisms
  // go to <work>/hop_<k>_n0/kept.csv, as a goal run's searches' do, and the run's
  // end signs every pending file there into the database, on all the run's
  // threads (`cobound sign` does the same for a run that was killed).
  int nextSearch = 0;
  if (std::filesystem::is_directory(workDir))
    for (const auto &d : std::filesystem::directory_iterator(workDir)) {
      const std::string n = d.path().filename().string();
      if (n.rfind(kFrozenHopDirPrefix, 0) != 0) continue;
      try {
        nextSearch = std::max(nextSearch, std::stoi(n.substr(sizeof kFrozenHopDirPrefix - 1)) + 1);
      } catch (const std::exception &) {
      }
    }
  size_t signedAppended = 0;
  bool signedPending = false;
  // Whether the work directory holds a search's pending directory (it also
  // holds the run's cobound.conf).
  auto anyPending = [&] {
    if (std::filesystem::is_directory(workDir))
      for (const auto &d : std::filesystem::directory_iterator(workDir))
        if (d.is_directory() && d.path().filename().string().rfind(kFrozenHopDirPrefix, 0) == 0) return true;
    return false;
  };
  auto signPending = [&] {
    if (signedPending || !anyPending()) return;
    signedPending = true;
    // A cobordism searched on its subject's table PD as written needs no
    // `.rows.csv` line (plan, "The .rows.csv sidecar"): every target of a
    // table-driven run.
    std::unordered_map<std::string, std::string> tablePD;
    for (const std::string &table : {knotTablePath, linkTablePath})
      if (!table.empty() && std::filesystem::exists(table))
        for (const linknaming::TableRow &r : linknaming::readTableRows(table))
          tablePD.emplace(r.name, r.pd);
    const auto sidecarLine = [&](const cobordisms::PendingCobordism &p) {
      auto it = tablePD.find(p.cobordism.subject);
      return it == tablePD.end() || it->second != p.incomingPD;
    };
    const cobordisms::SignResult s = cobordisms::signPending(
        workDir, cobordismsPath, {}, names, numThreads, pairSigCacheDir.value_or(""),
        recorded.loaded(), sidecarLine);
    signedAppended = s.appended;
    std::cout << kFrozenWitnessStoreLine << s.kept << " kept, " << s.fresh << " new, "
              << s.appended << " appended to " << cobordismsPath << " (signed in "
              << std::fixed << std::setprecision(0) << s.signSeconds << " s)\n"
              << std::defaultfloat;
  };
  // A halt writes what was found first: the pending cobordisms are signed.
  fatal::beforeHalt([&] {
    try {
      signPending();
    } catch (const std::exception &e) {
      std::cerr << "[!] the pending cobordisms were not signed (" << e.what()
                << "); `cobound sign` with work = " << workDir << " signs them\n";
    }
  });

  if (!knotSymmetryPath.empty()) {
    linknaming::SymmetryTable types;
    try {
      types = linknaming::readSymmetryTable(knotSymmetryPath);
    } catch (const std::exception &) {
      std::cerr << "[!] could not open knot symmetry table " << knotSymmetryPath << "\n";
      return 1;
    }
    for (const auto &[knot, type] : types)
      names.setSymmetry(knot, type);
    std::cout << "[+] Knot symmetry: " << types.size() << " types from " << knotSymmetryPath
              << "\n";
  }

  std::vector<InputRow> pending = targets::searchOrder(rows, maxCrossings, outputRows);

  // Every target's search shape, but for its boundary condition (per
  // target, below).
  search::SearchShape searchShape;
  searchShape.iddfsIterations = iddfsIterations;
  searchShape.iddfsStep = iddfsStep;
  searchShape.iddfsStart = iddfsStart;
  searchShape.maxFaces = maxFaces;
  searchShape.rootBudgetStart = rootBudgetStart;
  searchShape.rootBudgetGrowth = rootBudgetGrowth;
  searchShape.resolveUnlinked = resolveUnlinked;
  searchShape.limits = limits;

  // Searches whose accounting failed (divergence 2), and searches whose
  // outputs could not be written (an I/O error): the run exits 2.
  std::vector<std::string> unaccounted;
  std::vector<std::string> ioFailed;
  size_t processedThisRun = 0;
  size_t searchedThisRun = 0;

  for (const auto &row : pending) {
    if (runsignals::interrupted()) {
      std::cout << "[!] " << runsignals::name() << ": the run searches no further row\n";
      break;
    }

    std::cout << "[+] Searching " << row.name << " (literature [" << row.lo << ", " << row.hi
              << "], " << row.crossings << " crossings)...\n";

    // The target's own cobordism graph (divergence 6), which judges its finds,
    // and the thickening the search runs in: search::buildIncoming() of its
    // PD, collared through every layer. buildIncoming() and the graph's
    // certification of the incoming link throw for a bad PD, an incoming map
    // that cannot be built or checked, or a triangulated link that does not
    // redraw as its diagram; letting that escape would abort the whole run
    // over one bad target.
    std::unique_ptr<bounds::SearchJudge> judge;
    bool buildFailed = false;
    try {
      judge = std::make_unique<bounds::SearchJudge>(row.name, row.pdNotation, thickenLayers,
                                                     row.lo, *graphTables, *graphNamer,
                                                     names.symmetries(), numThreads);
      if (judge->reader().thickened().orientation->divergedFromDefaultIsomorphism)
        std::cerr << "[i] " << row.name
                  << ": the diagram's triangulation has a symmetry moving "
                     "L; using the map that takes L onto its own seed "
                     "(isIsomorphicTo() would not have)\n";
    } catch (const std::exception &e) {
      std::cerr << "[!] " << row.name << ": failed to build (" << e.what() << "), skipping\n";
      buildFailed = true;
    }
    const search::IncomingThickening *thickenedOrNull = judge ? &judge->reader().thickened() : nullptr;

    // A target that cannot be built, or that the search refuses: recorded as
    // such, and the run goes on.
    auto recordBuildFailure = [&] {
      OutputRow out;
      out.knot = row.name;
      out.status = "unresolved";
      out.cobordismKind = "none";
      out.literatureLo = row.lo;
      out.literatureHi = row.hi;
      out.searchOutcome = "build-failed";
      outputRows[row.name] = std::move(out);
      verdicts::writeOutputCsv(outputPath, rows, outputRows);
    };
    if (buildFailed) {
      recordBuildFailure();
      continue;
    }

    // The search's frontier: carried on from, and recorded (see frontier_dir).
    //
    // INVARIANT: a torn, partial or unparsable frontier loads as ABSENT. The
    // search starts fresh and the log names the reason; it never resumes from
    // such a file. Frontiers are written atomically but not fsynced (phase 5),
    // so a crash can leave one empty or cut short; read() refuses any file
    // without its 'end' line or with a malformed line (submanifoldsearch_test
    // cuts one at every byte; frontier_pending_test resumes from a truncated
    // and a garbage one). A whole frontier of another search is refused
    // later, by its fingerprint ("resumed no").
    auto frontierPath = [&](const std::string &dir) { return dir + "/" + row.name + ".frontier"; };
    std::optional<SearchFrontier> resumeFrom;
    if (resumeFrontierDir) {
      try {
        resumeFrom = SearchFrontier::load(frontierPath(*resumeFrontierDir));
      } catch (const std::exception &ex) {
        std::cout << "[!] " << row.name << ": WARNING: frontier not read (" << ex.what()
                  << "); searching from the start\n";
      }
    }

    search::SearchRequest request;
    request.name = row.name;
    request.shape = searchShape;
    // A multi-component link's own boundary necessarily puts more than one
    // of the surface's boundary curves on the single incoming ambient
    // component -- impossible under `connected`'s one-curve-per-ambient-
    // component rule, but exactly what `proper` allows. `connected` caps the
    // curve count on EVERY ambient boundary component, the outgoing link
    // included, so a knot searched under it can only ever discover
    // single-curve outgoing links: proper, the default, lifts that.
    const search::IncomingThickening &thickened = *thickenedOrNull;
    request.shape.condition = search::conditionFor(boundaryConditionMode, thickened.componentCount);
    // The search stops at the surface target or the per-search time limit; the
    // boundary drain then finishes.
    request.surfaceTarget = surfaceTarget;
    request.seconds = perKnotTimeLimit;
    request.resume = resumeFrom ? &*resumeFrom : nullptr;
    request.recordFrontier = frontierDir.has_value();
    request.pairSigCacheDir = pairSigCacheDir;
    // A knot target's complement goes into the census after its search (when
    // census writes are on).
    request.censusName = linknaming::baseName(row.name);
    request.literature = {.literatureLo = row.lo, .literatureHi = row.hi};
    // Its finds: kept one per cobordism (none the database holds), written to
    // its pending file as it runs, and signed at the run's end (divergence 7).
    const std::string searchDir = workDir + "/" + kFrozenHopDirPrefix +
                               std::to_string(nextSearch++) + kFrozenHopDirNodeMark + "0";
    std::filesystem::create_directories(searchDir);
    request.pending = searchDir + "/kept.csv";
    request.incomingPD = row.pdNotation;
    request.layers = thickenLayers;
    request.knownIdentities = &recorded.identities();
    request.runDirectory = workDir;
    // Each find, judged by the target's own cobordism graph as it is kept.
    request.reader = &judge->reader();
    long long judged = 0;
    request.judge = [&](const search::KeptSurface &k) {
      const bounds::SearchJudge::Verdict v =
          judge->add(k.link, k.genus, "find#" + std::to_string(judged++));
      search::FindJudgement j;
      if (!v.contradictions.empty()) j.contradiction = v.contradictions.front();
      j.constructive = v.constructive;
      return j;
    };
    request.outputs = {.progress = true,
                       .surfaceStats = surfaceStatsPath,
                       .surfaceLog = surfaceLogPath,
                       .rejectionSamples = rejectionSampleLog ? &*rejectionSampleLog : nullptr,
                       .selfIntersectionCensus = selfIntersectionCensusPath};

    // One searcher per target, so each search's names start from fresh
    // table caches, as they always have.
    const search::Searcher searcher(signatureTable ? &*signatureTable : nullptr,
                                        outgoingTables ? &*outgoingTables : nullptr, numThreads);
    search::SearchResult run;
    try {
      run = searcher.run(thickened, request);
    } catch (const search::SeedInvariantFailure &f) {
      fatal::flag(row.name + ": " + std::to_string(f.touching) +
                  " searchable non-seed triangles have an edge on the "
                  "search side, so found surfaces could change it.");
      fatal::haltIfFlagged();
    } catch (const search::SearchRefused &ex) {
      std::cerr << "[!] " << row.name << ": failed to build (" << ex.what() << "), skipping\n";
      recordBuildFailure();
      continue;
    }
    // A contradiction in the target's own cobordism graph: halts once the
    // search's cobordisms are written, below.
    if (judge->failures() > 0)
      std::cout << "[!] " << row.name << ": " << judge->failures() << " of " << judge->finds()
                << " finds could not enter the cobordism graph\n";
    if (!run.fatal.empty())
      fatal::flag(run.fatal);

    // The frontier, written before anything is claimed from the search: one
    // that cannot be written makes the search an I/O error (no frontier, no
    // exhaustion claim), as a failed write during the search does. Written
    // only when the search can vouch for every surface in its prefix: every
    // cobordism of that prefix is in its pending file, fsynced -- else a later
    // run would skip surfaces nobody looked at.
    std::string frontierFailure;
    if (frontierDir && run.frontier) {
      try {
        std::filesystem::create_directories(*frontierDir);
        run.frontier->save(frontierPath(*frontierDir));
      } catch (const std::exception &ex) {
        frontierFailure = ex.what();
        run.frontier.reset();
        if (run.ioFailure.empty()) run.ioFailure = "frontier: " + frontierFailure;
        if (run.outcome != "fatal-bug") run.outcome = "io-error";
      }
    }

    // Surfaces were accepted, yet not one reached the cobordism record. That
    // can be genuine (every one bounds another oriented variant), but it
    // is also exactly what a broken gate looks like, so it never licenses a
    // negative.
    if (run.nothingExamined)
      std::cout << "[!] " << row.name << ": WARNING: " << run.described
                << " surfaces accepted but none reached the witness record; "
                   "no exhaustion claimed for this search\n";

    if (run.stats.deepestExhaustedCap && run.accountingFailure.empty() &&
        !run.nothingExamined && !run.drainSkipped && run.ioFailure.empty()) {
      OutputRow &out = outputRows[row.name];
      // Never let a shallower run overwrite a deeper exhaustive result.
      out.exhaustedDepth = std::max(out.exhaustedDepth, *run.stats.deepestExhaustedCap);
      std::cout << "[+] " << row.name << ": EXHAUSTIVE to " << *run.stats.deepestExhaustedCap
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

    verdicts::writeOutputCsv(outputPath, rows, outputRows);

    // The search's breadth, and what became of its frontier (written above).
    if (resumeFrom || frontierDir) {
      search::printBreadth(std::cout, row.name, run);
      if (frontierDir && run.recordedFrontier) {
        if (!frontierFailure.empty()) {
          std::cout << "[!] " << row.name << ": frontier not written: " << frontierFailure
                    << "\n";
        } else if (!run.frontier) {
          std::cout << "[!] " << row.name << ": frontier not written: the "
                    << "row cannot vouch for every surface in its prefix\n";
        }
      }
    }

    // Divergence 2: a state that cannot occur halts the run, now that this
    // search's cobordisms are written; any other accounting failure ends only
    // this search (its outcome is `unaccounted`; no frontier, no exhaustion
    // claim, above) and the run goes on to its other targets; completeness is
    // what a run without a goal is for, so the run then exits 2.
    if (run.impossible > 0) {
      fatal::flag(row.name + ": surface accounting failed -- " + std::to_string(run.impossible) +
                  " surfaces hit a state that cannot occur (accounting: " + run.accounting +
                  ").");
    } else if (!run.accountingFailure.empty()) {
      unaccounted.push_back(row.name);
      std::cout << "[!!] " << row.name << ": surface accounting failed -- "
                << run.accountingFailure
                << " (its search vouches for nothing: no frontier, no "
                   "exhaustion claim; the run goes on and exits 2)\n";
    }
    if (!run.ioFailure.empty()) {
      ioFailed.push_back(row.name);
      std::cout << "[!!] " << row.name << ": an output write failed -- " << run.ioFailure
                << " (its search vouches for nothing: no frontier, no exhaustion claim; "
                   "the run goes on and exits 2)\n";
    }
    if (fatal::flagged())
      fatal::haltIfFlagged();

    search::printOutcome(std::cout, row.name, run);
    search::printComplementNaming(std::cout, row.name, run);
    search::printSearchProfile(std::cout, row.name, run);
    // The audit exists to catch exactly this; a wrong linking number prunes
    // (or keeps) surfaces it should not.
    if (run.petals.linkingDisagreements > 0) {
      fatal::flag(row.name + ": " + std::to_string(run.petals.linkingDisagreements) +
                  " petal linking numbers disagree between the cochain "
                  "and drilling routes (audit_linking).");
      fatal::haltIfFlagged();
    }
    // Its own line, so the summary line above (parsed by
    // tools/orchestrate/dispatch.py's RE_OUTCOME) is unchanged.
    if (resolveUnlinked)
      std::cout << "[+] " << row.name << ": " << run.stats.resolvedCount << " of "
                << run.stats.satisfyingCount
                << " accepted surfaces have unlinked self-intersections "
                   "(--resolve-unlinked)\n";

    ++processedThisRun;
  }

  verdicts::writeOutputCsv(outputPath, rows, outputRows);
  // The run's pending cobordisms, signed into the database (divergence 7).
  signPending();
  fatal::haltIfFlagged();

  std::cout << "\n[+] Done. Searched " << searchedThisRun << " of " << processedThisRun
            << " rows visited this run.\n";
  std::cout << "[+] Witness file: " << cobordisms.size() + signedAppended << " witnesses in "
            << cobordismsPath << "\n";
  if (!unaccounted.empty())
    std::cerr << "[!] " << unaccounted.size()
              << " search(es) failed their surface accounting (first: " << unaccounted.front()
              << "); exiting 2\n";
  if (!ioFailed.empty())
    std::cerr << "[!] " << ioFailed.size()
              << " search(es) failed an output write (first: " << ioFailed.front()
              << "); exiting 2\n";
  if (!unaccounted.empty() || !ioFailed.empty()) return 2;
  return 0;
}

} // namespace

// `run` (was verifyslicegenus's search runs and cascadesearch): the config
// (driver/config.h) names the targets, the work directory and, optionally, a
// goal. Without one each target is searched once (runWithoutGoal() above);
// with one the scheduler searches outwards from the target until the goal
// has a proof or the limits are spent (scheduler.h). The configuration it
// ran with is written to <work>/cobound.conf first: atomically, but not
// fsynced. It is a record of the run, rewritten by every run, and nothing
// reads it back to resume; an fsync before the search can stall for a
// second or more while the system writes back other data (phase 5's C2
// measurement).
int commands::run(const std::vector<std::string> &args) {
  try {
    const config::Config cfg = config::forCommand("run", std::nullopt, args);
    const std::string work = cfg.text("work");
    std::filesystem::create_directories(work);
    report::atomicWrite(
        work + "/cobound.conf", [&](std::ostream &out) { cfg.writeEffective(out, "run"); },
        report::Durability::cache);
    if (cfg.context() == config::Context::run)
      return runWithoutGoal(cfg);
    scheduler::GoalOptions options = scheduler::goalOptions(cfg);
    // The process-wide settings, from the config, before the first naming.
    options.censusLoaded = setup::applyRunSettings(cfg);
    // One policy for SIGINT and SIGTERM (divergence 8): the first ends the
    // running search cleanly and the run searches no further; a second ends
    // the process.
    runsignals::install();
    return scheduler::runToGoal(options);
  } catch (const std::exception &e) {
    std::cerr << "cobound run: " << e.what() << "\n";
    return 2;
  }
}
