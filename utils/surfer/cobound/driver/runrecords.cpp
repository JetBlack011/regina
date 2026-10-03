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
#include "linknaming/complement/unlinknaming.h"
#include "linknaming/diagrams/gaussdiagram.h"
#include "surfer/report/csvwriter.h"

namespace runrecords {

using namespace cascade;

void append(const std::string &work, const std::string &line) {
  std::ofstream(work + "/cascade.jsonl", std::ios::app) << line << "\n";
}

void writeProfiles(const std::string &work, const GraphView &v,
                   const std::function<std::string(NodeId)> &subjectName,
                   const std::function<bool(NodeId)> &searched) {
  const ProofGraph &g = v.g;
  const NodeRegistry &reg = v.reg;
  std::ofstream out(work + "/profiles.jsonl");
  for (NodeId n = 0; n < static_cast<NodeId>(g.nodeCount()); ++n) {
    // A crossingless node is named as the store names it (cascade_record.py):
    // it is never a hop's subject, so subjectName() has no better name.
    std::string name = subjectName(n);
    if (reg.known(n) && reg.info(n).diagram.signs.empty() && n != v.target) {
      const int k = g.node(n).components;
      name = identify::unlinkName(static_cast<size_t>(k));
    }
    out << "{\"node\":" << n << ",\"name\":\"" << json::escape(name) << '"';
    if (auto it = v.tableName.find(n); it != v.tableName.end())
      out << ",\"table\":\"" << json::escape(it->second) << '"';
    out << ",\"label\":\"" << json::escape(g.node(n).label) << '"';
    if (auto it = v.depth.find(n); it != v.depth.end())
      out << ",\"depth\":" << it->second;
    // A split far side's whole is added to the graph, not the registry,
    // and has no diagram of its own.
    if (reg.known(n))
      out << ",\"crossings\":" << reg.info(n).diagram.signs.size();
    out << ",\"searched\":" << (searched(n) ? "true" : "false") << ','
        << g.profileFields(n) << "}\n";
  }
}

void writeNodeBounds(const std::string &work, const GraphView &v) {
  using Kind = ProofGraph::LowerReason::Kind;
  const ProofGraph &g = v.g;
  const NodeRegistry &reg = v.reg;
  std::ofstream o(work + "/node_bounds.jsonl");
  auto kindName = [](Kind k) {
    switch (k) {
    case Kind::literature: return "literature";
    case Kind::linking: return "linking";
    case Kind::witness: return "witness";
    case Kind::splitWhole: return "split-whole";
    case Kind::splitPiece: return "split-piece";
    case Kind::seed: return "seed";
    case Kind::sumPiece: return "sum-piece";
    default: return "none";
    }
  };
  for (size_t i = 0; i < g.nodeCount(); ++i) {
    const NodeId n = static_cast<NodeId>(i);
    const Node &node = g.node(n);
    o << "{\"node\":" << n << ",\"label\":\"" << json::escape(node.label) << "\",\"components\":"
      << node.components << ",\"target\":" << (n == v.target ? "true" : "false");
    if (auto t = v.tableName.find(n); t != v.tableName.end())
      o << ",\"table\":\"" << json::escape(t->second) << "\"";
    if (reg.known(n)) {
      const NodeInfo &ni = reg.info(n);
      const exactnaming::GaussDiagram &d = ni.diagram;
      o << ",\"crossings\":" << d.crossings() << ",\"hyperbolic\":"
        << (ni.hyperbolic ? "true" : "false");
      if (ni.hyperbolic) o << ",\"volume\":" << std::setprecision(12) << ni.volume;
      o << ",\"pd\":\"" << json::escape(rowPD(d)) << "\",\"signs\":" << json::array(d.signs) << ",\"gauss\":[";
      for (size_t c = 0; c < d.comps.size(); ++c) o << (c ? "," : "") << json::array(d.comps[c]);
      o << "],\"linking\":[";
      for (size_t a = 0; a < ni.linking.size(); ++a) o << (a ? "," : "") << json::array(ni.linking[a]);
      o << "]";
    }
    o << ",\"upper\":[";
    bool first = true;
    for (const ProfileEntry &e : node.profile.entries()) {
      bool constructive = true;
      for (RecordId r : g.proof(e.record))
        if (g.record(r).kind == RecordKind::leaf &&
            g.record(r).source.rfind("literature", 0) == 0)
          constructive = false;
      o << (first ? "" : ",") << "{\"partition\":\"" << e.partition.str() << "\",\"genus\":"
        << e.genus << ",\"record\":" << e.record << ",\"constructive\":"
        << (constructive ? "true" : "false") << "}";
      first = false;
    }
    o << "],\"lower\":[";
    first = true;
    if (node.components <= ProofGraph::kMaxLowerComponents)
      for (const Partition &p : allPartitions(node.components)) {
        const auto f = g.lowerWhy(n, p);
        if (f.value <= 0 || !(f.storedFor == p)) continue; // stored bounds only, once each
        o << (first ? "" : ",") << "{\"partition\":\"" << p.str() << "\",\"value\":"
          << (f.value >= ProofGraph::kNoSurface ? std::string("\"inf\"") : std::to_string(f.value))
          << ",\"kind\":\"" << kindName(f.reason.kind) << "\"";
        if (f.reason.kind == Kind::literature)
          o << ",\"source\":\"" << json::escape(node.lowerBoundSource) << "\"";
        o << "}";
        first = false;
      }
    o << "]}\n";
  }
}

void writeLowerReport(const std::string &work, const GraphView &v,
                      const exactnaming::ExactTables &tables, const std::string &targetName,
                      const std::map<std::string, bool> &special, unsigned threads) {
  // For every tabulated node Y: the least charge of carrying a lower bound
  // from Y to the target, over every path the graph holds. Measured by
  // seeding Y alone at a large M in a copy and reading what reaches the
  // target (propagateLower() takes the maximum over sources, and M dwarfs
  // every real bound), so charge = M - lower(target). A literature bound
  // lo(Y) carries lo(Y) - charge; hi(Y) - charge is the most Y could ever
  // carry. README.md, "Lower bounds": only a charge-0 path to a source whose
  // bound is not Lipschitz (lower_bound_sources.csv, `special`) can beat the
  // target's own literature bound.
  const ProofGraph &g = v.g;
  const NodeId target = v.target;
  const Partition goal = v.goal;
  const int targetLower = g.lower(target, goal);
  int litLo = -1;
  if (const exactnaming::TableEntry *e = tables.entry(targetName))
    if (auto g4 = exactnaming::parseTableG4(e->g4)) litLo = g4->first;
  std::ofstream out(work + "/lower_report.jsonl");
  out << "{\"target\":\"" << json::escape(targetName) << "\",\"target_lower\":" << targetLower
      << ",\"lit_lo\":" << litLo << ",\"nodes\":" << g.nodeCount() << "}\n";
  // What node n's lower bound `seed` alone carries to the target: every other
  // lower bound forgotten (clearLowerBounds()), n seeded, relaxed. Only
  // consistent facts are ever seeded -- a value the node could really have,
  // at most its best proved genus -- or the split rules, which read proved
  // surfaces, would pump bounds without limit. -1 when the what-if itself
  // meets a contradiction (then nothing it says is used).
  auto carried = [&](NodeId n, int seed) {
    ProofGraph what = g;
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
    NodeId n;
    std::string name;
    int lo, hi, could;
    int carries = 0, couldCarry = 0;
  };
  std::vector<Job> jobs;
  for (const auto &[n, name] : v.tableName) {
    if (n == target) continue;
    const exactnaming::TableEntry *e = tables.entry(name);
    if (!e) continue;
    auto g4 = exactnaming::parseTableG4(e->g4);
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

void writeNodesCsv(const std::string &work, const std::map<NodeId, std::string> &subjects,
                   const NodeRegistry &reg) {
  // The cascade: subjects, as the atlas's results/cascade/nodes.csv lists
  // them (cascade_record.py), so a later identity can be attached to each.
  std::ofstream nodes(work + "/nodes.csv");
  nodes << "name,components,crossings,pd,signs,gauss,label\n";
  for (const auto &[n, name] : subjects) {
    if (name.rfind("cascade:", 0) != 0) continue;
    const exactnaming::GaussDiagram &d = reg.info(n).diagram;
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
    nodes << csvField(name) << ',' << d.components() << ',' << d.crossings() << ','
          << csvField(rowPD(d)) << ',' << csvField(signs.str()) << ','
          << csvField(gauss.str()) << ",node " << n << '\n';
  }
}

} // namespace runrecords
