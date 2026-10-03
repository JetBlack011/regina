//
//  sweep.cpp
//
//  A run without a goal (verifyslicegenus's sweep, created by John Teague
//  on 07/29/2026): each target searched once.
//

#include "cobound/driver/sweep.h"

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
#include "cobound/driver/fatal.h"
#include "cobound/driver/setup.h"
#include "cobound/driver/signals.h"
#include "cobound/driver/targets.h"
#include "cobound/search/search.h"
#include "cobound/search/searchreport.h"
#include "cobound/solver/literature.h"
#include "cobound/solver/verdicts.h"
#include "linknaming/names.h"
#include "linknaming/tables.h"

namespace sweep {

using cobordismgraph::InputRow;
using rowsearch::BoundaryConditionMode;
using verdicts::OutputRow;

namespace {

BoundaryConditionMode boundaryCondition(const std::string &mode) {
  if (mode == "auto") return BoundaryConditionMode::automatic;
  if (mode == "connected") return BoundaryConditionMode::connected;
  return BoundaryConditionMode::proper;
}

} // namespace

int run(const config::Config &cfg) {
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
  // every host explores roots at the same rate. Measured over 52 rows, they
  // do not: at an identical 3600 thread-seconds one host yielded a median
  // 972,655 qualifying surfaces per row and another 295,178 -- a 3.3x gap,
  // tight on both sides (+/-15%), and a property of the machine rather than of
  // the row. Targeting the surface count instead equalises the thing a
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
  const bool exactFarSideNames = cfg.flag("exact_far_side_names");
  const std::string workDir = cfg.text("work");
  const unsigned numThreads = cfg.threads();
  // Thickened and collared through every layer (divergence 6: the cobordism
  // graph reads each find on its row as a goal run reads a stored row's).
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
    if (fresh) *rejectionSampleLog << "row,reason,tubed_genus,connected,boundary,pairsig\n";
  }

  // One policy for SIGINT and SIGTERM (divergence 8): the first ends the
  // current row's search cleanly and the run searches no further row; a
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
      metadataRows += witnessstore::loadNameTable(table, names);
    } catch (const std::exception &e) {
      std::cerr << "[!] could not load name table " << table << ": " << e.what()
                << " (continuing without it)\n";
    }
  }
  std::cout << "[+] Name table: " << names.size() << " names (" << metadataRows
            << " from tables other than --input)\n";

  // The tables, loaded once: the exact tables, which a search's own
  // cobordism graph names its links by (divergence 6) and its outgoing links
  // are named exactly by (exact_far_side_names), and diagram naming's
  // signature table, drawn from them.
  std::optional<exactnaming::ExactTables> exactTables;
  if (exactFarSideNames) {
    const auto t0 = std::chrono::steady_clock::now();
    try {
      exactTables =
          exactnaming::ExactTables::load(knotTablePath, linkTablePath, knotSymmetryPath);
      std::cout << "[+] exact far-side names: " << exactTables->size() << " table entries ("
                << std::chrono::duration_cast<std::chrono::milliseconds>(
                       std::chrono::steady_clock::now() - t0)
                       .count()
                << " ms)\n";
    } catch (const std::exception &e) {
      std::cerr << "[!] exact far-side names off: " << e.what() << "\n";
    }
  }
  std::optional<exactnaming::ExactTables> graphTablesOwn;
  const exactnaming::ExactTables *graphTables = exactTables ? &*exactTables : nullptr;
  std::optional<exactnaming::ExactNamer> graphNamer;
  try {
    if (!graphTables)
      graphTables = &graphTablesOwn.emplace(
          exactnaming::ExactTables::load(knotTablePath, linkTablePath, knotSymmetryPath));
    graphNamer.emplace(*graphTables, cascade::NodeAxioms::namerLimits());
  } catch (const std::exception &e) {
    std::cerr << "[!] the cobordism graph cannot load the tables: " << e.what() << "\n";
    return 1;
  }
  std::optional<farside::SignatureTable> signatureTable;
  {
    const auto t0 = std::chrono::steady_clock::now();
    try {
      // From the tables the graph already loaded: one table load (phase 5).
      signatureTable = farside::SignatureTable::fromTables(*graphTables);
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
  cascade::RecordedWitnesses recorded(cobordismsPath,
                                      witnessstore::loadWitnesses(cobordismsPath, false));
  const std::vector<cobordismgraph::Witness> &witnesses = recorded.all();
  std::cout << "[+] Resuming with " << witnesses.size()
            << " previously-recorded witnesses from " << cobordismsPath << "\n";

  // The run directory (plan divergence 7): each search's pending cobordisms
  // go to <work>/hop_<k>_n0/kept.csv, as a goal run's hops' do, and the run's
  // end signs every pending file there into the database, on all the run's
  // threads (`cobound sign` does the same for a run that was killed).
  int nextHop = 0;
  if (std::filesystem::is_directory(workDir))
    for (const auto &d : std::filesystem::directory_iterator(workDir)) {
      const std::string n = d.path().filename().string();
      if (n.rfind("hop_", 0) != 0) continue;
      try {
        nextHop = std::max(nextHop, std::stoi(n.substr(4)) + 1);
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
        if (d.is_directory() && d.path().filename().string().rfind("hop_", 0) == 0) return true;
    return false;
  };
  auto signPending = [&] {
    if (signedPending || !anyPending()) return;
    signedPending = true;
    // A cobordism searched on its subject's table PD as written needs no
    // `.rows.csv` line (plan, "The .rows.csv sidecar"): every row of a
    // table-driven run.
    std::unordered_map<std::string, std::string> tablePD;
    for (const std::string &table : {knotTablePath, linkTablePath})
      if (!table.empty() && std::filesystem::exists(table))
        for (const exactnaming::TableRow &r : exactnaming::readTableRows(table))
          tablePD.emplace(r.name, r.pd);
    const auto sidecarLine = [&](const cascade::PendingWitness &p) {
      auto it = tablePD.find(p.witness.subject);
      return it == tablePD.end() || it->second != p.rowPD;
    };
    const cascade::StoreResult s = cascade::signPending(
        workDir, cobordismsPath, {}, names, numThreads, pairSigCacheDir.value_or(""),
        recorded.loaded(), sidecarLine);
    signedAppended = s.appended;
    std::cout << "[+] witness store: " << s.kept << " kept, " << s.fresh << " new, "
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
    exactnaming::SymmetryTable types;
    try {
      types = exactnaming::readSymmetryTable(knotSymmetryPath);
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

  // Every row's search shape, but for its boundary condition (per row,
  // below).
  cascade::SearchShape searchShape;
  searchShape.iddfsIterations = iddfsIterations;
  searchShape.iddfsStep = iddfsStep;
  searchShape.iddfsStart = iddfsStart;
  searchShape.maxFaces = maxFaces;
  searchShape.rootBudgetStart = rootBudgetStart;
  searchShape.rootBudgetGrowth = rootBudgetGrowth;
  searchShape.resolveUnlinked = resolveUnlinked;
  searchShape.limits = limits;

  // Searches whose accounting failed (divergence 2): the run exits 2.
  std::vector<std::string> unaccounted;
  size_t processedThisRun = 0;
  size_t searchedThisRun = 0;

  for (const auto &row : pending) {
    if (runsignals::interrupted()) {
      std::cout << "[!] " << runsignals::name() << ": the run searches no further row\n";
      break;
    }

    std::cout << "[+] Searching " << row.name << " (literature [" << row.lo << ", " << row.hi
              << "], " << row.crossings << " crossings)...\n";

    // The row's own cobordism graph (divergence 6), which judges its finds,
    // and the row the search runs in: rowsearch::buildRow() of its PD,
    // collared through every layer. buildRow() and the graph's row
    // certification throw for a bad PD, a row map that cannot be built or
    // checked, or a triangulated link that does not redraw as its diagram;
    // letting that escape would abort the whole sweep over one bad row.
    std::unique_ptr<cascade::SearchJudge> judge;
    bool buildFailed = false;
    try {
      judge = std::make_unique<cascade::SearchJudge>(row.name, row.pdNotation, thickenLayers,
                                                     row.lo, *graphTables, *graphNamer,
                                                     names.symmetries(), numThreads);
      if (judge->row().rowBuild().orientation->divergedFromDefaultIsomorphism)
        std::cerr << "[i] " << row.name
                  << ": the diagram's triangulation has a symmetry moving "
                     "L; using the map that takes L onto its own seed "
                     "(isIsomorphicTo() would not have)\n";
    } catch (const std::exception &e) {
      std::cerr << "[!] " << row.name << ": failed to build (" << e.what() << "), skipping\n";
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
      verdicts::writeOutputCsv(outputPath, rows, outputRows);
    };
    if (buildFailed) {
      recordBuildFailure();
      continue;
    }

    // The row's frontier: carried on from, and recorded (see frontier_dir).
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

    cascade::SearchRequest request;
    request.name = row.name;
    request.shape = searchShape;
    // A multi-component link's own boundary necessarily puts more than one
    // of the surface's boundary curves on the single search-side ambient
    // component -- impossible under `connected`'s one-curve-per-ambient-
    // component rule, but exactly what `proper` allows. `connected` caps the
    // curve count on EVERY ambient boundary component, the far side
    // included, so a knot row searched under it can only ever discover
    // single-curve far sides: proper, the default, lifts that.
    const rowsearch::RowBuild &rb = *rbp;
    request.shape.condition = rowsearch::conditionFor(boundaryConditionMode, rb.componentCount);
    // The search stops at the surface target or the per-row time limit; the
    // boundary drain then finishes.
    request.surfaceTarget = surfaceTarget;
    request.seconds = perKnotTimeLimit;
    request.resume = resumeFrom ? &*resumeFrom : nullptr;
    request.recordFrontier = frontierDir.has_value();
    request.pairSigCacheDir = pairSigCacheDir;
    // A knot row's complement goes into the census after its search (when
    // census writes are on).
    request.censusName = cobordismgraph::baseName(row.name);
    request.sweep = {.literatureLo = row.lo, .literatureHi = row.hi};
    // Its finds: kept one per cobordism (none the database holds), written to
    // its pending file as it runs, and signed at the run's end (divergence 7).
    const std::string hopDir = workDir + "/hop_" + std::to_string(nextHop++) + "_n0";
    std::filesystem::create_directories(hopDir);
    request.pending = hopDir + "/kept.csv";
    request.rowPD = row.pdNotation;
    request.layers = thickenLayers;
    request.knownIdentities = &recorded.identities();
    request.runDirectory = workDir;
    // Each find, judged by the row's own cobordism graph as it is kept.
    request.row = &judge->row();
    long long judged = 0;
    request.judge = [&](const cascade::KeptSurface &k) {
      const cascade::SearchJudge::Verdict v =
          judge->add(k.link, k.genus, "find#" + std::to_string(judged++));
      cascade::FindJudgement j;
      if (!v.contradictions.empty()) j.contradiction = v.contradictions.front();
      j.constructive = v.constructive;
      return j;
    };
    request.outputs = {.progress = true,
                       .surfaceStats = surfaceStatsPath,
                       .surfaceLog = surfaceLogPath,
                       .rejectionSamples = rejectionSampleLog ? &*rejectionSampleLog : nullptr,
                       .selfIntersectionCensus = selfIntersectionCensusPath};

    // One searcher per row, so each row's exact names start from fresh
    // table caches, as they always have.
    const cascade::HopSearcher searcher(signatureTable ? &*signatureTable : nullptr,
                                        exactTables ? &*exactTables : nullptr, numThreads);
    cascade::HopRun run;
    try {
      run = searcher.run(rb, request);
    } catch (const cascade::SeedInvariantFailure &f) {
      fatal::flag(row.name + ": " + std::to_string(f.touching) +
                  " searchable non-seed triangles have an edge on the "
                  "search side, so found surfaces could change it.");
      fatal::haltIfFlagged();
    } catch (const cascade::SearchRefused &ex) {
      std::cerr << "[!] " << row.name << ": failed to build (" << ex.what() << "), skipping\n";
      recordBuildFailure();
      continue;
    }
    // A contradiction in the row's own cobordism graph: halts once the row's
    // witnesses are written, below.
    if (judge->failures() > 0)
      std::cout << "[!] " << row.name << ": " << judge->failures() << " of " << judge->finds()
                << " finds could not enter the cobordism graph\n";
    if (!run.fatal.empty())
      fatal::flag(run.fatal);

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

    // The row's breadth, and its frontier: written only now that every
    // cobordism of the prefix it covers is in its pending file, fsynced, and
    // only when the row can vouch for having examined every surface in it --
    // else a later run would skip surfaces nobody looked at.
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

    // Divergence 2: a state that cannot occur halts the run, now that this
    // row's witnesses are written; any other accounting failure ends only
    // this search (its outcome is `unaccounted`; no frontier, no exhaustion
    // claim, above) and the run goes on to its other rows; completeness is
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
    if (fatal::flagged())
      fatal::haltIfFlagged();

    rowsearch::printOutcome(std::cout, row.name, run);
    rowsearch::printIdentification(std::cout, row.name, run);
    rowsearch::printSearchProfile(std::cout, row.name, run);
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
  std::cout << "[+] Witness file: " << witnesses.size() + signedAppended << " witnesses in "
            << cobordismsPath << "\n";
  if (!unaccounted.empty()) {
    std::cerr << "[!] " << unaccounted.size()
              << " search(es) failed their surface accounting (first: " << unaccounted.front()
              << "); exiting 2\n";
    return 2;
  }
  return 0;
}

} // namespace sweep
