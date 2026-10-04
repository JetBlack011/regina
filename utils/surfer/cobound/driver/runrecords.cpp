//
//  runrecords.cpp
//

#include "cobound/driver/runrecords.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

#include "cobound/bounds/searchcobordisms.h"
#include "cobound/json.h"
#include "cobound/parallelfor.h"
#include "cobound/frozen.h"
#include "linknaming/complement/unlinknaming.h"
#include "linknaming/diagrams/gaussdiagram.h"
#include "surfer/report/csvwriter.h"

namespace runrecords {

using namespace bounds;

void append(const std::string &work, const std::string &line) {
  std::ofstream(work + "/" + kFrozenCascadeJsonl, std::ios::app) << line << "\n";
}

void writePartitionGenera(const std::string &work, const GraphView &v,
                   const std::function<std::string(LinkId)> &subjectName,
                   const std::function<bool(LinkId)> &searched) {
  const CobordismGraph &g = v.g;
  const LinkRegistry &reg = v.reg;
  std::ofstream out(work + "/" + kFrozenProfilesJsonl);
  for (LinkId n = 0; n < static_cast<LinkId>(g.linkCount()); ++n) {
    // A crossingless link is named as the atlas's recorder names it (cascade_record.py):
    // it is never a search's subject, so subjectName() has no better name.
    std::string name = subjectName(n);
    if (reg.known(n) && reg.info(n).diagram.signs.empty() && n != v.target) {
      const int k = g.link(n).components;
      name = complement::unlinkName(static_cast<size_t>(k));
    }
    out << "{\"node\":" << n << ",\"name\":\"" << json::escape(name) << '"';
    if (auto it = v.tableName.find(n); it != v.tableName.end())
      out << ",\"table\":\"" << json::escape(it->second) << '"';
    out << ",\"label\":\"" << json::escape(g.link(n).label) << '"';
    if (auto it = v.depth.find(n); it != v.depth.end())
      out << ",\"depth\":" << it->second;
    // A split outgoing link's whole is added to the graph, not the registry,
    // and has no diagram of its own.
    if (reg.known(n))
      out << ",\"crossings\":" << reg.info(n).diagram.signs.size();
    out << ",\"searched\":" << (searched(n) ? "true" : "false") << ','
        << g.partitionGeneraFields(n) << "}\n";
  }
}

void writeLinkBounds(const std::string &work, const GraphView &v) {
  using Kind = CobordismGraph::LowerReason::Kind;
  const CobordismGraph &g = v.g;
  const LinkRegistry &reg = v.reg;
  std::ofstream o(work + "/" + kFrozenNodeBoundsJsonl);
  auto kindName = [](Kind k) {
    switch (k) {
    case Kind::literature: return "literature";
    case Kind::linking: return "linking";
    case Kind::cobordism: return kFrozenLowerKindWitness;
    case Kind::splitWhole: return "split-whole";
    case Kind::splitPiece: return "split-piece";
    case Kind::seed: return "seed";
    case Kind::sumPiece: return "sum-piece";
    default: return "none";
    }
  };
  for (size_t i = 0; i < g.linkCount(); ++i) {
    const LinkId n = static_cast<LinkId>(i);
    const GraphLink &link = g.link(n);
    o << "{\"node\":" << n << ",\"label\":\"" << json::escape(link.label) << "\",\"components\":"
      << link.components << ",\"target\":" << (n == v.target ? "true" : "false");
    if (auto t = v.tableName.find(n); t != v.tableName.end())
      o << ",\"table\":\"" << json::escape(t->second) << "\"";
    if (reg.known(n)) {
      const LinkInfo &ni = reg.info(n);
      const linknaming::GaussDiagram &d = ni.diagram;
      o << ",\"crossings\":" << d.crossings() << ",\"hyperbolic\":"
        << (ni.hyperbolic ? "true" : "false");
      if (ni.hyperbolic) o << ",\"volume\":" << std::setprecision(12) << ni.volume;
      o << ",\"pd\":\"" << json::escape(diagramPD(d)) << "\",\"signs\":" << json::array(d.signs) << ",\"gauss\":[";
      for (size_t c = 0; c < d.comps.size(); ++c) o << (c ? "," : "") << json::array(d.comps[c]);
      o << "],\"linking\":[";
      for (size_t a = 0; a < ni.linking.size(); ++a) o << (a ? "," : "") << json::array(ni.linking[a]);
      o << "]";
    }
    o << ",\"upper\":[";
    bool first = true;
    for (const PartitionGenus &e : link.partitionGenera.entries()) {
      bool constructive = true;
      for (DerivationId r : g.proof(e.derivation))
        if (g.derivation(r).kind == DerivationKind::leaf &&
            g.derivation(r).source.rfind("literature", 0) == 0)
          constructive = false;
      o << (first ? "" : ",") << "{\"partition\":\"" << e.partition.str() << "\",\"genus\":"
        << e.genus << ",\"record\":" << e.derivation << ",\"constructive\":"
        << (constructive ? "true" : "false") << "}";
      first = false;
    }
    o << "],\"lower\":[";
    first = true;
    if (link.components <= CobordismGraph::kMaxLowerComponents)
      for (const Partition &p : allPartitions(link.components)) {
        const auto f = g.lowerWhy(n, p);
        if (f.value <= 0 || !(f.storedFor == p)) continue; // stored bounds only, once each
        o << (first ? "" : ",") << "{\"partition\":\"" << p.str() << "\",\"value\":"
          << (f.value >= CobordismGraph::kNoSurface ? std::string("\"inf\"") : std::to_string(f.value))
          << ",\"kind\":\"" << kindName(f.reason.kind) << "\"";
        if (f.reason.kind == Kind::literature)
          o << ",\"source\":\"" << json::escape(link.lowerBoundSource) << "\"";
        o << "}";
        first = false;
      }
    o << "]}\n";
  }
}

void writeLowerReport(const std::string &work, const GraphView &v,
                      const linknaming::Tables &tables, const std::string &targetName,
                      const std::map<std::string, bool> &special, unsigned threads) {
  // For every tabulated link Y: the least charge of carrying a lower bound
  // from Y to the target, over every path the graph holds. Measured by
  // seeding Y alone at a large M in a copy and reading what reaches the
  // target (propagateLower() takes the maximum over sources, and M dwarfs
  // every real bound), so charge = M - lower(target). A literature bound
  // lo(Y) carries lo(Y) - charge; hi(Y) - charge is the most Y could ever
  // carry. README.md, "Lower bounds": only a charge-0 path to a source whose
  // bound is not Lipschitz (lower_bound_sources.csv, `special`) can beat the
  // target's own literature bound.
  const CobordismGraph &g = v.g;
  const LinkId target = v.target;
  const Partition goal = v.goal;
  const int targetLower = g.lower(target, goal);
  int litLo = -1;
  if (const linknaming::TableEntry *e = tables.entry(targetName))
    if (auto g4 = linknaming::parseTableG4(e->g4)) litLo = g4->first;
  std::ofstream out(work + "/lower_report.jsonl");
  out << "{\"target\":\"" << json::escape(targetName) << "\",\"target_lower\":" << targetLower
      << ",\"lit_lo\":" << litLo << ",\"nodes\":" << g.linkCount() << "}\n";
  // What link n's lower bound `seed` alone carries to the target: every other
  // lower bound forgotten (clearLowerBounds()), n seeded, relaxed. Only
  // consistent facts are ever seeded -- a value the link could really have,
  // at most its best proved genus -- or the split rules, which read proved
  // surfaces, would pump bounds without limit. -1 when the what-if itself
  // meets a contradiction (then nothing it says is used).
  auto carried = [&](LinkId n, int seed) {
    CobordismGraph what = g;
    what.clearLowerBounds();
    const size_t before = what.contradictions().size();
    what.setGenusLowerBound(n, seed, "lower-report what-if");
    what.propagateLower();
    if (what.contradictions().size() != before) return -1;
    return what.lower(target, goal);
  };
  // Each what-if copies the graph and relaxes it alone, so they run on the
  // run's threads (serially they were a third of a run's driver time,
  // 2026-09-30); the lines are then written in the same order as before.
  struct Job {
    LinkId n;
    std::string name;
    int lo, hi, could;
    int carries = 0, couldCarry = 0;
  };
  std::vector<Job> jobs;
  for (const auto &[n, name] : v.tableName) {
    if (n == target) continue;
    const linknaming::TableEntry *e = tables.entry(name);
    if (!e) continue;
    auto g4 = linknaming::parseTableG4(e->g4);
    if (!g4) continue;
    // The most n could be: its literature upper end, or less if a surface
    // for it is already proved.
    int could = g4->second;
    if (auto b = g.bestConnected(n)) could = std::min(could, b->genus);
    jobs.push_back({n, name, g4->first, g4->second, could});
  }
  parallelFor(jobs.size(), std::max(threads, 1u), [&](size_t i) {
    Job &j = jobs[i];
    j.carries = j.lo > 0 ? carried(j.n, j.lo) : 0;
    j.couldCarry = j.could > 0 ? (j.could == j.lo ? j.carries : carried(j.n, j.could)) : 0;
  });
  int bestCarry = 0, bestCould = 0;
  std::string bestName, bestCouldName;
  for (const Job &j : jobs) {
    const std::string &name = j.name;
    const int carries = j.carries, couldCarry = j.couldCarry;
    if (carries <= 0 && couldCarry <= 0) continue; // n reaches the target with nothing
    const auto sp = special.find(name);
    out << "{\"node\":" << j.n << ",\"name\":\"" << json::escape(name) << "\",\"lit_lo\":"
        << j.lo << ",\"lit_hi\":" << j.hi << ",\"could\":" << j.could
        << ",\"special\":" << (sp == special.end() ? "null" : sp->second ? "true" : "false")
        << ",\"carries\":" << carries << ",\"could_carry\":" << couldCarry << "}\n";
    if (carries > bestCarry) {
      bestCarry = carries;
      bestName = name;
    }
    if (couldCarry > bestCould) {
      bestCould = couldCarry;
      bestCouldName = name;
    }
  }
  std::cout << "[+] lower report: target lower " << targetLower << " (literature " << litLo
            << "); best carried " << (bestName.empty() ? std::string("none")
                                                       : std::to_string(bestCarry) + " from " + bestName)
            << "; most any tabulated node could carry "
            << (bestCouldName.empty() ? std::string("none")
                                      : std::to_string(bestCould) + " from " + bestCouldName)
            << "\n";
}

void writeLinksCsv(const std::string &work, const std::map<LinkId, std::string> &subjects,
                   const LinkRegistry &reg) {
  // The `cascade:` subjects, as the atlas's results/cascade/nodes.csv lists
  // them (cascade_record.py), so a later identity can be attached to each.
  std::ofstream links(work + "/" + kFrozenNodesCsv);
  links << "name,components,crossings,pd,signs,gauss,label\n";
  for (const auto &[n, name] : subjects) {
    if (name.rfind(kFrozenCascadeSubjectPrefix, 0) != 0) continue;
    const linknaming::GaussDiagram &d = reg.info(n).diagram;
    std::ostringstream signs, gauss;
    signs << '[';
    for (size_t i = 0; i < d.signs.size(); ++i) signs << (i ? ", " : "") << d.signs[i];
    signs << ']';
    gauss << '[';
    for (size_t c = 0; c < d.comps.size(); ++c) {
      gauss << (c ? ", [" : "[");
      for (size_t j = 0; j < d.comps[c].size(); ++j) gauss << (j ? ", " : "") << d.comps[c][j];
      gauss << ']';
    }
    gauss << ']';
    links << csvField(name) << ',' << d.components() << ',' << d.crossings() << ','
          << csvField(diagramPD(d)) << ',' << csvField(signs.str()) << ','
          << csvField(gauss.str()) << ',' << kFrozenNodesCsvLabel << n << '\n';
  }
}

} // namespace runrecords
