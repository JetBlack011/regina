//
//  verdicts.cpp
//

#include "cobound/solver/verdicts.h"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <unordered_set>

#include "surfer/report/atomicwrite.h"
#include "surfer/report/csvwriter.h"

namespace verdicts {

const char *const OUTPUT_HEADER =
    "knot,resolved_genus,status,witness_kind,witness_pairsig,via_knot,"
    "via_edge_genus,depends_on,literature_lo,literature_hi,"
    "derived_lo,derived_hi,witness_basis,tubed,searched_faces,search_outcome,"
    "exhausted_depth";

std::string formatOutputRow(const OutputRow &r) {
  std::ostringstream out;
  out << csvField(r.knot) << ',' << r.resolvedGenus << ',' << r.status << ','
      << r.cobordismKind << ',' << csvField(r.cobordismPairSig) << ','
      << csvField(r.viaKnot) << ',' << r.viaCobordismGenus << ','
      << csvField(r.dependsOn) << ',' << r.literatureLo << ','
      << r.literatureHi << ',' << r.derivedLo << ',' << r.derivedHi << ','
      << r.cobordismBasis << ',' << (r.tubed ? "true" : "false") << ','
      << r.searchedFaces << ',' << r.searchOutcome << ',' << r.exhaustedDepth;
  return out.str();
}

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
    r.cobordismKind = f[3];
    r.cobordismPairSig = f[4];
    r.viaKnot = f[5];
    try {
      r.viaCobordismGenus = f[6].empty() ? 0 : std::stoi(f[6]);
    } catch (const std::exception &) {
      r.viaCobordismGenus = 0;
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
      r.cobordismBasis = f[12];
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

void writeOutputCsv(const std::filesystem::path &path,
                    const std::vector<solver::InputRow> &rows,
                    const std::unordered_map<std::string, OutputRow> &outputRows) {
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

const char *statusName(solver::Status s) {
  using S = solver::Status;
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

OutputRow rowFromVerdict(
    const std::string &name, const solver::Verdict &v, const OutputRow *existing,
    const std::unordered_map<std::string, solver::Bounds> &bounds,
    cobordisms::PairSigReader &reader) {
  OutputRow out;
  out.knot = name;
  out.status = statusName(v.status);
  out.literatureLo = v.litLo;
  out.literatureHi = v.litHi;
  out.resolvedGenus = v.value;

  const auto &b = v.bounds;
  if (b.haveUpper()) {
    out.derivedHi = std::to_string(b.hi);
    out.cobordismKind =
        b.kind == cobordisms::CobordismKind::direct ? "direct" : "cobordism";
    out.cobordismPairSig = !b.pairSig.empty() ? b.pairSig : reader.at(b.pairSigOffset);
    out.viaKnot = b.viaName;
    out.viaCobordismGenus = b.viaGenus;
    out.dependsOn = solver::buildDependsOn(b.viaName, bounds);
    out.cobordismBasis = b.basis == solver::Basis::constructive
                           ? "constructive"
                           : "literature-assisted";
    out.tubed = b.tubed;
  } else {
    out.cobordismKind = "none";
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

} // namespace verdicts
