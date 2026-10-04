// searchjudge.cpp

#include "cobound/bounds/searchjudge.h"

#include <numeric>
#include <stdexcept>

#include "linknaming/diagrams/simplification.h"

using linknaming::GaussDiagram;

namespace bounds {

namespace {

// Names compared as names: a depth-0 graph computes no link classes (plan,
// "Startup per process"). Only which links may take a literature UPPER bound
// as a leaf depends on it, and a duplicate of the searched link taking its
// own class-mate's would give an assisted bound at least its literature
// lower bound: never a contradiction, never constructive.
LinkAxioms::Options judgeOptions(unsigned threads) {
  LinkAxioms::Options o;
  o.literature = true;
  o.classes = false;
  o.threads = threads;
  o.log = nullptr;
  return o;
}

} // namespace

SearchJudge::SearchJudge(const std::string &name, const std::string &pd, int layers,
                         int literatureLo, const linknaming::Tables &tables,
                         const linknaming::LinkNamer &namer,
                         const linknaming::SymmetryTable &symmetries, unsigned threads)
    : reg_(g_), axioms_(g_, reg_, tables, namer, symmetries, judgeOptions(threads)),
      literatureLo_(literatureLo) {
  // The searched link as a graph link, as a goal run's graph takes a table link's
  // stored cobordisms: interned from its simplified diagram, assembled
  // against the searched diagram (CobordismAssembler certifies that the triangulated
  // link redraws as it).
  const regina::Link link = linknaming::linkFromTablePD(pd);
  std::vector<size_t> origin(link.countComponents());
  std::iota(origin.begin(), origin.end(), 0);
  const GaussDiagram diagram = GaussDiagram::of(link, origin);
  const LinkMatch m = reg_.intern(linknaming::simplifyKeepingComponents(diagram), "target " + name);
  target_ = m.link;
  components_ = g_.link(target_).components;
  axioms_.target = target_;
  axioms_.targetClass = name;
  axioms_.depth.emplace(target_, 0);
  // Its literature lower bound, always: what the gates hold every proof to.
  g_.setGenusLowerBound(target_, literatureLo, "literature " + name);
  SearchedLink searched;
  searched.link = target_;
  searched.diagram = diagram;
  searched.linkMap = m.componentMap;
  searched.pd = pd;
  searched.layers = layers;
  assembler_ = std::make_unique<CobordismAssembler>(g_, reg_, searched);
}

SearchJudge::Verdict SearchJudge::add(const outgoing::OutgoingLink &link, int genus,
                                      const std::string &key) {
  ++finds_;
  const size_t before = g_.linkCount();
  AddedCobordism e;
  try {
    e = assembler_->addRead(link, genus, key);
  } catch (const std::exception &ex) {
    e.ok = false;
    e.why = ex.what();
  }
  if (!e.ok) ++failures_;
  std::vector<LinkId> fresh;
  for (size_t n = before; n < g_.linkCount(); ++n) fresh.push_back(static_cast<LinkId>(n));
  axioms_.name(fresh, 1);
  g_.propagate();
  g_.propagateLower();

  Verdict v;
  v.contradictions = g_.contradictions();
  if (auto best = g_.best(target_, Partition::coarsest(components_));
      best && best->genus <= literatureLo_) {
    bool constructive = true;
    for (DerivationId r : g_.proof(best->derivation))
      if (g_.derivation(r).kind == DerivationKind::leaf &&
          g_.derivation(r).source.rfind("literature", 0) == 0)
        constructive = false;
    if (constructive) v.constructive = best->genus;
  }
  return v;
}

} // namespace bounds
