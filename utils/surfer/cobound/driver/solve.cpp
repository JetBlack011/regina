//
//  solve.cpp
//

#include "cobound/driver/commands.h"

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "cobound/cobordisms/database.h"
#include "cobound/driver/config.h"
#include "cobound/driver/fatal.h"
#include "cobound/driver/targets.h"
#include "cobound/solver/literature.h"
#include "cobound/solver/solver.h"
#include "cobound/solver/solverinputs.h"
#include "cobound/solver/verdicts.h"
#include "linknaming/census/censusnaming.h"
#include "linknaming/tables.h"

namespace {

using cobordismgraph::InputRow;
using verdicts::OutputRow;

// `solve` (was verifyslicegenus --solve-only): every conclusion re-derived
// from the database, the tables and the solver's other inputs
// (solver/solverinputs.h), into the verdicts file (solver/verdicts.h). It
// never writes the database, and never searches. Its log lines are
// verifyslicegenus's (verify.sh greps `Totals across`). Returns 0, or 1 when
// an input cannot be read; halts with 2 on a contradiction (fatal.h), after
// writing the verdicts.
int solveWith(const config::Config &cfg) {
  const std::string outputPath = cfg.text("verdicts");
  const std::string inputPath = cfg.text("targets");
  const std::string cobordismsPath = cfg.text("cobordisms");
  const std::string knotTablePath = cfg.text("knot_table");
  const std::string linkTablePath = cfg.text("link_table");
  const std::string knotSymmetryPath = cfg.text("knot_symmetry");
  const std::string nameAliasPath = cfg.text("name_aliases");
  const std::string farSideResolutionPath = cfg.text("far_side_resolutions");
  const std::string farSideExactPath = cfg.text("far_side_exact");
  const std::string linkClassesPath = cfg.text("link_classes");
  const std::string cascadeProofsPath = cfg.text("cascade_proofs");
  const bool sumRules = cfg.flag("sum_rules");
  const int maxCrossings = static_cast<int>(cfg.integer("max_crossings"));
  const std::string censusPath = cfg.text("census");
  const bool censusLoaded = census::setCensusPath(censusPath);

  std::cout << "------ verifyslicegenus \U0001F30A ------\n\n";
  std::cout << (censusLoaded ? "[+] census: loaded from " : "[+] census: not found at ")
            << censusPath << (censusLoaded ? "\n\n" : ", skipping\n\n");

  std::vector<InputRow> rows;
  try {
    rows = targets::loadInputCsv(inputPath);
  } catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    return 1;
  }
  std::cout << "[+] Loaded " << rows.size() << " knots from " << inputPath << "\n";

  // Literature metadata for every name that could ever appear in the
  // graph, not just this run's targets: both tables, always.
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

  std::unordered_map<std::string, OutputRow> outputRows = verdicts::loadOutputCsv(outputPath);

  // A bound's witness pair signature, read back from its line when the
  // verdicts print it.
  witnessstore::PairSigReader pairSigReader;
  pairSigReader.setPath(cobordismsPath);
  // Both per-witness tables are keyed on the pair signature's key, which is
  // hashed at load only when one of them will be looked up.
  const std::vector<cobordismgraph::Witness> witnesses = witnessstore::loadWitnesses(
      cobordismsPath, !farSideResolutionPath.empty() || !farSideExactPath.empty());
  std::cout << "[+] Resuming with " << witnesses.size()
            << " previously-recorded witnesses from " << cobordismsPath << "\n";

  std::unordered_map<std::string, std::string> nameAliases;
  if (!nameAliasPath.empty()) {
    try {
      nameAliases = solverinputs::loadNameAliases(nameAliasPath);
      std::cout << "[+] Name aliases: " << nameAliases.size()
                << " observed names resolved to classical ones from " << nameAliasPath << "\n";
    } catch (const std::exception &e) {
      std::cerr << "[!] could not load name aliases " << nameAliasPath << ": " << e.what()
                << " (continuing without them)\n";
    }
  }

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

  std::unordered_map<std::string, std::vector<solverinputs::FarSideResolution>>
      farSideResolutions;
  if (!farSideResolutionPath.empty()) {
    try {
      farSideResolutions = solverinputs::loadFarSideResolutions(farSideResolutionPath);
      std::cout << "[+] Far-side resolutions: " << farSideResolutions.size()
                << " witnesses with a proved far side from " << farSideResolutionPath << "\n";
    } catch (const std::exception &e) {
      std::cerr << "[!] could not load far-side resolutions " << farSideResolutionPath << ": "
                << e.what() << " (continuing without them)\n";
    }
  }

  std::unordered_map<std::string, solverinputs::ExactFarSide> farSideExact;
  if (!farSideExactPath.empty()) {
    size_t clashes = 0;
    try {
      farSideExact = solverinputs::loadFarSideExact(farSideExactPath, clashes);
    } catch (const std::exception &e) {
      std::cerr << "[!] could not load far-side exact names " << farSideExactPath << ": "
                << e.what() << "\n";
      return 1;
    }
    std::cout << "[+] Far-side exact names: " << farSideExact.size() << " witnesses from "
              << farSideExactPath;
    if (clashes)
      std::cout << " (" << clashes << " with two names dropped)";
    std::cout << "\n";
  }

  // link_classes: table names that are one oriented link up to mirror and
  // global reversal (tableclasses: a version of one table diagram is the
  // other, or an isometry of the complements carries meridians to meridians
  // with a uniform orientation sign). Every member is read as its class's
  // canonical name, so a class is one node; tableclasses refuses a class
  // whose members' literature values differ.
  std::unordered_map<std::string, std::string> linkClasses;
  if (!linkClassesPath.empty()) {
    try {
      linkClasses = solverinputs::loadLinkClasses(linkClassesPath);
    } catch (const std::exception &e) {
      std::cerr << "[!] " << e.what() << "\n";
      return 1;
    }
    std::cout << "[+] Link classes: " << linkClasses.size()
              << " table names read as their class's canonical name, from " << linkClassesPath
              << "\n";
  }
  auto classOf = [&linkClasses](const std::string &name) -> const std::string & {
    auto it = linkClasses.find(name);
    return it == linkClasses.end() ? name : it->second;
  };
  names.setSumRules(sumRules);
  if (sumRules)
    std::cout << "[+] Sum rules: sums along components and splits with link "
                 "factors are bounded from their pieces\n";

  // cascade_proofs: only CERTIFIED proofs of the connected goal bound g4;
  // names (the target and every literature leaf) read through the link
  // classes, as witnesses'.
  std::vector<cobordismgraph::ExternalProof> externalProofs;
  if (!cascadeProofsPath.empty()) {
    solverinputs::CertifiedBounds certified;
    try {
      certified = solverinputs::loadCascadeProofs(cascadeProofsPath, classOf);
    } catch (const std::exception &e) {
      std::cerr << "[!] " << e.what() << "\n";
      return 1;
    }
    externalProofs = std::move(certified.proofs);
    std::cout << "[+] Cascade proofs: " << externalProofs.size()
              << " certified bounds on connected g4 from " << cascadeProofsPath;
    if (certified.skipped)
      std::cout << " (" << certified.skipped << " others skipped)";
    std::cout << "\n";
  }

  // The witness list the SOLVER sees, which is not the database's: aliases,
  // resolutions and exact names resolve a far side to what we have since
  // proved it to be, while `witnesses` keeps what the search observed.
  size_t aliasesApplied = 0;
  size_t resolutionsApplied = 0;
  size_t exactApplied = 0, exactRefused = 0;
  auto solverWitnesses = [&]() -> std::vector<cobordismgraph::Witness> {
    std::vector<cobordismgraph::Witness> out =
        nameAliases.empty()
            ? witnesses
            : solverinputs::applyNameAliases(witnesses, nameAliases, names, aliasesApplied);
    if (!farSideResolutions.empty())
      out = solverinputs::applyFarSideResolutions(std::move(out), witnesses, farSideResolutions,
                                                  names, resolutionsApplied);
    if (!farSideExact.empty())
      out = solverinputs::applyFarSideExact(std::move(out), farSideExact, exactApplied,
                                            exactRefused);
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
  // searched long ago without re-searching any of them. A resolution/alias
  // contradiction is a data error, not a bug, so it exits rather than
  // aborting; caught here because the table is static -- if it contradicts
  // at all, it does so on this first solve.
  std::unordered_map<std::string, cobordismgraph::Bounds> bounds;
  {
    std::vector<cobordismgraph::Witness> initialWitnesses;
    try {
      initialWitnesses = solverWitnesses();
    } catch (const std::exception &e) {
      std::cerr << "[!] " << e.what() << "\n";
      return 1;
    }
    bounds = cobordismgraph::propagate(initialWitnesses, names, externalProofs);
    // Release it now: it is a full witness set, pair signatures included, and
    // held for the rest of the run it raised peak memory by ~45% -- enough for
    // a full-master solve to be OOM-killed on yoga (2026-09-24).
    std::vector<cobordismgraph::Witness>().swap(initialWitnesses);
  }
  if (!nameAliases.empty())
    std::cout << "[+] Name aliases: applied to " << aliasesApplied << " witness edges\n";
  if (!farSideResolutions.empty())
    std::cout << "[+] Far-side resolutions: applied to " << resolutionsApplied
              << " witness edges\n";
  if (!farSideExact.empty())
    std::cout << "[+] Far-side exact names: applied to " << exactApplied << " witness edges"
              << (exactRefused ? " (" + std::to_string(exactRefused) +
                                     " refused: component count differs)"
                               : std::string())
              << "\n";
  std::cout << "[+] Solver: derived bounds for " << bounds.size() << " names\n\n";

  // Rows above max_crossings that the verdicts do not hold yet are recorded
  // as skipped; nothing is searched.
  targets::searchOrder(rows, maxCrossings, outputRows);
  std::cout << "[+] --solve-only: re-deriving from " << witnesses.size()
            << " witnesses, no searching.\n\n";

  // The solve, from everything the database knows: every affected verdict
  // rewritten.
  bounds = cobordismgraph::propagate(solverWitnesses(), names, externalProofs);
  // Every name that could need its row rewritten -- crucially including
  // rows that currently HAVE a row but no longer have any derived bound.
  // Iterating `bounds` alone would leave such a row frozen at whatever a
  // previous (possibly buggier) solver wrote, which is exactly the case
  // solve exists to correct.
  std::vector<std::string> toJudge;
  toJudge.reserve(bounds.size() + outputRows.size());
  for (const auto &[name, unused] : bounds)
    toJudge.push_back(name);
  for (const auto &[name, unused] : outputRows)
    if (!bounds.contains(name))
      toJudge.push_back(name);
  // link_classes: a member is its class's node, so it is judged whenever
  // that node has bounds, row or no row yet.
  for (const auto &[member, canonical] : linkClasses)
    if (bounds.contains(canonical) && !bounds.contains(member) && !outputRows.contains(member))
      toJudge.push_back(member);
  for (const std::string &name : toJudge) {
    cobordismgraph::Bounds b; // default = nothing derived
    // A class member's row reads its class's node (link_classes).
    if (auto it = bounds.find(classOf(name)); it != bounds.end())
      b = it->second;
    cobordismgraph::Verdict v = cobordismgraph::judge(name, b, names);
    if (v.status == cobordismgraph::Status::contradiction)
      fatal::flag(v.reason);
    auto existing = outputRows.find(name);
    // Only names we actually track get a row: the graph is full of
    // incidental nodes (bare isoSigs, unlinks) that are useful for
    // chaining but aren't results in their own right.
    if (existing == outputRows.end() && !names.find(name))
      continue;
    OutputRow updated = verdicts::rowFromVerdict(
        name, v, existing == outputRows.end() ? nullptr : &existing->second, bounds,
        pairSigReader);
    outputRows[name] = std::move(updated);
  }
  verdicts::writeOutputCsv(outputPath, rows, outputRows);
  fatal::haltIfFlagged();

  size_t verified = 0, verifiedAssisted = 0, improved = 0, pinnedCount = 0, unresolvedCount = 0;
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
  std::cout << "\n[+] Done. Searched 0 of 0 rows visited this run.\n";
  std::cout << "[+] Witness file: " << witnesses.size() << " witnesses in " << cobordismsPath
            << "\n";
  std::cout << "[+] Totals across " << outputRows.size() << " tracked names: " << verified
            << " verified (" << verifiedAssisted << " more only with literature help), "
            << improved << " improved, " << pinnedCount << " pinned, " << unresolvedCount
            << " unresolved.\n";
  return 0;
}

} // namespace

int commands::solve(const std::vector<std::string> &args) {
  try {
    return solveWith(config::forCommand("solve", config::Context::solve, args));
  } catch (const config::Error &e) {
    std::cerr << "cobound solve: " << e.what() << "\n";
    return 2;
  }
}
