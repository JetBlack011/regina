//
//  solverinputs.cpp
//

#include "cobound/solver/solverinputs.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

#include "cobound/cobordisms/cobordismkey.h"
#include "linknaming/complement/unlinknaming.h"
#include "linknaming/names.h"
#include "surfer/report/csvwriter.h"

namespace solverinputs {

std::unordered_map<std::string, std::string> loadNameAliases(const std::filesystem::path &path) {
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
    if (complement::isUnlinkName(f[1]))
      throw std::runtime_error("Name alias table maps '" + f[0] + "' to '" + f[1] +
                               "': an unknot/unlink can only be established "
                               "by identify(), never by alias");
    aliases.emplace(f[0], f[1]);
  }
  return aliases;
}

std::string aliasKey(const std::string &name) {
  return linknaming::stripCensusSuffix(name);
}

std::vector<cobordisms::Witness>
applyNameAliases(const std::vector<cobordisms::Witness> &witnesses,
                 const std::unordered_map<std::string, std::string> &aliases,
                 const solver::NameTable &names, size_t &appliedOut) {
  std::vector<cobordisms::Witness> resolved = witnesses;
  size_t applied = 0;
  for (cobordisms::Witness &w : resolved) {
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

std::vector<cobordisms::Witness>
applyFarSideExact(std::vector<cobordisms::Witness> witnesses,
                  const std::unordered_map<std::string, ExactFarSide> &exact, size_t &applied,
                  size_t &refused) {
  applied = refused = 0;
  for (cobordisms::Witness &w : witnesses) {
    if (w.kind != cobordisms::WitnessKind::cobordism)
      continue;
    if (w.pairSigKey.empty() && !w.pairSig.empty())
      w.pairSigKey = cobordisms::witnessKey(w.pairSig);
    auto it = exact.find(w.pairSigKey);
    if (it == exact.end())
      continue;
    if (it->second.components != w.otherComponents) {
      ++refused;
      continue;
    }
    w.other = it->second.name;
    w.otherCandidates = linknaming::exactCandidates(it->second.name);
    w.farSideProved = true;
    w.farSideExact = it->second.exact;
    ++applied;
  }
  return witnesses;
}

std::unordered_map<std::string, std::vector<FarSideResolution>>
loadFarSideResolutions(const std::filesystem::path &path) {
  std::unordered_map<std::string, std::vector<FarSideResolution>> resolutions;
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("Cannot open far-side resolution table: " + path.string());

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

std::vector<cobordisms::Witness> applyFarSideResolutions(
    std::vector<cobordisms::Witness> resolved,
    const std::vector<cobordisms::Witness> &observed,
    const std::unordered_map<std::string, std::vector<FarSideResolution>> &resolutions,
    const solver::NameTable &names, size_t &appliedOut) {
  // Taken by value and rewritten in place. (A loaded witness no longer
  // carries its pair signature in memory, only pairSigKey.)
  size_t applied = 0;
  size_t linkAliasesOverridden = 0;
  std::vector<std::string> conflicts;

  for (cobordisms::Witness &w : resolved) {
    if (w.other.empty())
      continue;
    if (w.pairSigKey.empty() && !w.pairSig.empty())
      w.pairSigKey = cobordisms::witnessKey(w.pairSig);
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
          linknaming::componentsFromName(r.name) == w.otherComponents;
      // candidates() falls back to {name} itself for an unregistered base,
      // so a real table hit is one whose front is a different (tagged) name.
      // A composite K #_c L has L's components (a knot is summed INTO a
      // component), so it is L's variants that state the count.
      const std::optional<linknaming::CompositeName> cp =
          linknaming::compositeParts(r.name);
      const std::string countable = cp ? cp->link : r.name;
      const std::vector<std::string> variants =
          names.candidates(countable, w.otherComponents);
      const bool countFromTable =
          !variants.empty() && variants.front() != countable &&
          linknaming::componentsFromName(variants.front()) == w.otherComponents;
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
    const bool aliasFired = idx < observed.size() && observed[idx].other != w.other;
    const std::string aliasedBase = linknaming::stripOrientationTag(w.other);
    const std::string resolvedBase = linknaming::stripOrientationTag(*match);
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
    msg << "far-side resolutions contradict name aliases on " << conflicts.size()
        << " KNOT far sides, e.g.";
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

std::unordered_map<std::string, std::string> loadLinkClasses(const std::filesystem::path &path) {
  std::unordered_map<std::string, std::string> linkClasses;
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("could not open link classes " + path.string());
  std::string line;
  std::getline(in, line); // name,canonical,proof
  while (std::getline(in, line)) {
    auto f = parseCsvLine(line);
    if (f.size() >= 2 && !f[0].empty() && !f[1].empty())
      linkClasses[f[0]] = f[1];
  }
  return linkClasses;
}

CertifiedBounds loadCascadeProofs(const std::filesystem::path &path,
                                  const std::function<const std::string &(const std::string &)> &classOf) {
  CertifiedBounds out;
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("could not open cascade proofs " + path.string());
  std::string line;
  std::getline(in, line);
  const std::vector<std::string> head = parseCsvLine(line);
  auto col = [&head](const std::string &c) {
    return static_cast<size_t>(std::find(head.begin(), head.end(), c) - head.begin());
  };
  const size_t cT = col("target"), cG = col("goal"), cB = col("bound"), cS = col("support"),
               cV = col("verdict"), cSrc = col("source");
  if (std::max({cT, cG, cB, cS, cV, cSrc}) >= head.size())
    throw std::runtime_error(path.string() +
                             ": needs target,goal,bound,support,verdict,source columns");
  while (std::getline(in, line)) {
    if (line.empty())
      continue;
    const auto f = parseCsvLine(line);
    if (f.size() < head.size() || f[cV] != "CERTIFIED" || f[cG] != "connected") {
      ++out.skipped;
      continue;
    }
    solver::ExternalProof p;
    p.name = classOf(f[cT]);
    p.genus = std::stoi(f[cB]);
    std::istringstream support(f[cS]);
    for (std::string s; std::getline(support, s, ';');)
      if (!s.empty())
        p.support.push_back(classOf(s));
    p.source = f[cSrc];
    out.proofs.push_back(std::move(p));
  }
  return out;
}

} // namespace solverinputs
