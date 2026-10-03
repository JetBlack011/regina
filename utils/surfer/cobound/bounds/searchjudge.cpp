// searchjudge.cpp

#include "cobound/bounds/searchjudge.h"

#include <numeric>
#include <stdexcept>

#include "linknaming/diagrams/simplification.h"

using exactnaming::GaussDiagram;

namespace cascade {

namespace {

// Names compared as names: a depth-0 graph computes no link classes (plan,
// "Startup per process"). Only which links may take a literature UPPER bound
// as a leaf depends on it, and a duplicate of the searched link taking its
// own class-mate's would give an assisted bound at least its literature
// lower bound: never a contradiction, never constructive.
NodeAxioms::Options judgeOptions(unsigned threads) {
  NodeAxioms::Options o;
  o.literature = true;
  o.classes = false;
  o.threads = threads;
  o.log = nullptr;
  return o;
}

} // namespace

SearchJudge::SearchJudge(const std::string &name, const std::string &pd, int layers,
                         int literatureLo, const exactnaming::ExactTables &tables,
                         const exactnaming::ExactNamer &namer,
                         const exactnaming::SymmetryTable &symmetries, unsigned threads)
    : reg_(g_), axioms_(g_, reg_, tables, namer, symmetries, judgeOptions(threads)),
      literatureLo_(literatureLo) {
  // The searched link as a node, as a goal run's graph takes a table row's
  // stored cobordisms: interned from its simplified diagram, assembled
  // against the row's own (HopAssembler certifies that the triangulated
  // link redraws as it).
  const regina::Link link = exactnaming::linkFromTablePD(pd);
  std::vector<size_t> origin(link.countComponents());
  std::iota(origin.begin(), origin.end(), 0);
  const GaussDiagram diagram = GaussDiagram::of(link, origin);
  const NodeMatch m = reg_.intern(simplifyKeepingComponents(diagram), "target " + name);
  target_ = m.node;
  components_ = g_.node(target_).components;
  axioms_.target = target_;
  axioms_.targetClass = name;
  axioms_.depth.emplace(target_, 0);
  // Its literature lower bound, always: what the gates hold every proof to.
  g_.setGenusLowerBound(target_, literatureLo, "literature " + name);
  HopRow row;
  row.node = target_;
  row.diagram = diagram;
  row.nodeMap = m.componentMap;
  row.pd = pd;
  row.layers = layers;
  hop_ = std::make_unique<HopAssembler>(g_, reg_, row);
}

SearchJudge::Verdict SearchJudge::add(const farside::OutgoingLink &link, int genus,
                                      const std::string &key) {
  ++finds_;
  const size_t before = g_.nodeCount();
  HopEdge e;
  try {
    e = hop_->addRead(link, genus, key);
  } catch (const std::exception &ex) {
    e.ok = false;
    e.why = ex.what();
  }
  if (!e.ok) ++failures_;
  std::vector<NodeId> fresh;
  for (size_t n = before; n < g_.nodeCount(); ++n) fresh.push_back(static_cast<NodeId>(n));
  axioms_.name(fresh, 1);
  g_.propagate();
  g_.propagateLower();

  Verdict v;
  v.contradictions = g_.contradictions();
  if (auto best = g_.best(target_, Partition::coarsest(components_));
      best && best->genus <= literatureLo_) {
    bool constructive = true;
    for (RecordId r : g_.proof(best->record))
      if (g_.record(r).kind == RecordKind::leaf &&
          g_.record(r).source.rfind("literature", 0) == 0)
        constructive = false;
    if (constructive) v.constructive = best->genus;
  }
  return v;
}

} // namespace cascade
