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

using solver::InputRow;
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
  const std::string outgoingResolutionsPath = cfg.text("outgoing_resolutions");
  const std::string outgoingNamesPath = cfg.text("outgoing_names_file");
  const std::string linkClassesPath = cfg.text("link_classes");
  const std::string certifiedBoundsPath = cfg.text("certified_bounds");
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

  std::unordered_map<std::string, OutputRow> outputRows = verdicts::loadOutputCsv(outputPath);

  // A bound's cobordism pair signature, read back from its line when the
  // verdicts print it.
  cobordisms::PairSigReader pairSigReader;
  pairSigReader.setPath(cobordismsPath);
  // Both per-cobordism tables are keyed on the pair signature's key, which is
  // hashed at load only when one of them will be looked up.
  const std::vector<cobordisms::Cobordism> cobordisms = cobordisms::loadCobordisms(
      cobordismsPath, !outgoingResolutionsPath.empty() || !outgoingNamesPath.empty());
  std::cout << "[+] Resuming with " << cobordisms.size()
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

  std::unordered_map<std::string, std::vector<solverinputs::OutgoingResolution>>
      outgoingResolutions;
  if (!outgoingResolutionsPath.empty()) {
    try {
      outgoingResolutions = solverinputs::loadOutgoingResolutions(outgoingResolutionsPath);
      std::cout << "[+] Far-side resolutions: " << outgoingResolutions.size()
                << " witnesses with a proved far side from " << outgoingResolutionsPath << "\n";
    } catch (const std::exception &e) {
      std::cerr << "[!] could not load far-side resolutions " << outgoingResolutionsPath << ": "
                << e.what() << " (continuing without them)\n";
    }
  }

  std::unordered_map<std::string, solverinputs::OutgoingName> outgoingNamed;
  if (!outgoingNamesPath.empty()) {
    size_t clashes = 0;
    try {
      outgoingNamed = solverinputs::loadOutgoingNames(outgoingNamesPath, clashes);
    } catch (const std::exception &e) {
      std::cerr << "[!] could not load far-side exact names " << outgoingNamesPath << ": "
                << e.what() << "\n";
      return 1;
    }
    std::cout << "[+] Far-side exact names: " << outgoingNamed.size() << " witnesses from "
              << outgoingNamesPath;
    if (clashes)
      std::cout << " (" << clashes << " with two names dropped)";
    std::cout << "\n";
  }

  // link_classes: table names that are one oriented link up to mirror and
  // global reversal (tableclasses: a version of one table diagram is the
  // other, or an isometry of the complements carries meridians to meridians
  // with a uniform orientation sign). Every member is read as its class's
  // canonical name, so a class is one name to the solver; tableclasses refuses a class
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

  // certified_bounds (the atlas's cascade_proofs.csv): only CERTIFIED proofs
  // of the connected goal bound g4; names (the target and every literature
  // leaf) read through the link classes, as cobordisms'.
  std::vector<solver::ExternalProof> externalProofs;
  if (!certifiedBoundsPath.empty()) {
    solverinputs::CertifiedBounds certified;
    try {
      certified = solverinputs::loadCertifiedBounds(certifiedBoundsPath, classOf);
    } catch (const std::exception &e) {
      std::cerr << "[!] " << e.what() << "\n";
      return 1;
    }
    externalProofs = std::move(certified.proofs);
    std::cout << "[+] Cascade proofs: " << externalProofs.size()
              << " certified bounds on connected g4 from " << certifiedBoundsPath;
    if (certified.skipped)
      std::cout << " (" << certified.skipped << " others skipped)";
    std::cout << "\n";
  }

  // The cobordism list the SOLVER sees, which is not the database's: aliases,
  // resolutions and names resolve an outgoing link to what we have since
  // proved it to be, while `cobordisms` keeps what the search observed.
  size_t aliasesApplied = 0;
  size_t resolutionsApplied = 0;
  size_t namesApplied = 0, namesRefused = 0;
  size_t alternativesSplit = 0;
  auto solverCobordisms = [&]() -> std::vector<cobordisms::Cobordism> {
    std::vector<cobordisms::Cobordism> out =
        nameAliases.empty()
            ? cobordisms
            : solverinputs::applyNameAliases(cobordisms, nameAliases, names, aliasesApplied);
    if (!outgoingResolutions.empty())
      out = solverinputs::applyOutgoingResolutions(std::move(out), cobordisms, outgoingResolutions,
                                                  names, resolutionsApplied);
    if (!outgoingNamed.empty())
      out = solverinputs::applyOutgoingNames(std::move(out), outgoingNamed, namesApplied,
                                            namesRefused);
    // A stored candidate of alternatives "A|B" is its alternatives, whichever
    // build signed it.
    alternativesSplit = solverinputs::splitStoredAlternatives(out);
    // Last: whole names only. A name inside a sum or split is a piece,
    // bounded by its literature value, which is the same across a class.
    if (!linkClasses.empty())
      for (cobordisms::Cobordism &w : out) {
        w.subject = classOf(w.subject);
        w.other = classOf(w.other);
        for (std::string &c : w.otherCandidates)
          c = classOf(c);
      }
    return out;
  };

  // Every conclusion is re-derived from the cobordism set on every run, so a
  // solver fix or a literature-table update takes effect on rows that were
  // searched long ago without re-searching any of them. A resolution/alias
  // contradiction is a data error, not a bug, so it exits rather than
  // aborting; caught here because the table is static -- if it contradicts
  // at all, it does so on this first solve.
  std::unordered_map<std::string, solver::Bounds> bounds;
  {
    std::vector<cobordisms::Cobordism> initialCobordisms;
    try {
      initialCobordisms = solverCobordisms();
    } catch (const std::exception &e) {
      std::cerr << "[!] " << e.what() << "\n";
      return 1;
    }
    bounds = solver::propagate(initialCobordisms, names, externalProofs);
    // Release it now: it is a full cobordism set, pair signatures included, and
    // held for the rest of the run it raised peak memory by ~45% -- enough for
    // a full-master solve to be OOM-killed on yoga (2026-09-24).
    std::vector<cobordisms::Cobordism>().swap(initialCobordisms);
  }
  if (!nameAliases.empty())
    std::cout << "[+] Name aliases: applied to " << aliasesApplied << " witness edges\n";
  if (!outgoingResolutions.empty())
    std::cout << "[+] Far-side resolutions: applied to " << resolutionsApplied
              << " witness edges\n";
  if (!outgoingNamed.empty())
    std::cout << "[+] Far-side exact names: applied to " << namesApplied << " witness edges"
              << (namesRefused ? " (" + std::to_string(namesRefused) +
                                     " refused: component count differs)"
                               : std::string())
              << "\n";
  if (alternativesSplit)
    std::cout << "[+] Stored candidates: " << alternativesSplit
              << " witness edges recorded alternatives \"A|B\" as one candidate, read as "
                 "their alternatives\n";
  std::cout << "[+] Solver: derived bounds for " << bounds.size() << " names\n\n";

  // Rows above max_crossings that the verdicts do not hold yet are recorded
  // as skipped; nothing is searched.
  targets::searchOrder(rows, maxCrossings, outputRows);
  std::cout << "[+] --solve-only: re-deriving from " << cobordisms.size()
            << " witnesses, no searching.\n\n";

  // The solve, from everything the database knows: every affected verdict
  // rewritten.
  bounds = solver::propagate(solverCobordisms(), names, externalProofs);
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
  // link_classes: a member reads as its class's canonical name, so it is
  // judged whenever that name has bounds, row or no row yet.
  for (const auto &[member, canonical] : linkClasses)
    if (bounds.contains(canonical) && !bounds.contains(member) && !outputRows.contains(member))
      toJudge.push_back(member);
  for (const std::string &name : toJudge) {
    solver::Bounds b; // default = nothing derived
    // A class member's row reads its class's canonical name (link_classes).
    if (auto it = bounds.find(classOf(name)); it != bounds.end())
      b = it->second;
    solver::Verdict v = solver::judge(name, b, names);
    if (v.status == solver::Status::contradiction)
      fatal::flag(v.reason);
    auto existing = outputRows.find(name);
    // Only names we actually track get a row: the graph is full of
    // incidental names (bare isoSigs, unlinks) that are useful for
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
  std::cout << "[+] Witness file: " << cobordisms.size() << " witnesses in " << cobordismsPath
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
