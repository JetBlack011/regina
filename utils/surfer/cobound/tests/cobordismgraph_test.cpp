// cobordismgraph_test.cpp
//
// Tests for bounds/cobordismgraph.h: the fixed point over a graph with cycles,
// both directions of every cobordism, splits, component maps, derivations,
// and the contradiction gates. Each block names the assumption it pins
// (README.md, "Assumptions and their tests").

#include <algorithm>
#include <map>
#include <random>
#include <set>
#include <vector>

#include "cobound/bounds/cobordismgraph.h"
#include "linknaming/tests/check.h"

using namespace bounds;

namespace {

CobordismShape pants() { // one component: 1 incoming curve, 2 outgoing
  CobordismShape c;
  c.components = 1;
  c.genus = 0;
  c.inComponent = {0};
  c.outComponent = {0, 0};
  return c;
}

CobordismShape annulusShape(int genus) { // knot to knot, genus g
  CobordismShape c;
  c.components = 1;
  c.genus = genus;
  c.inComponent = {0};
  c.outComponent = {0};
  return c;
}

void checkAllDerivations(const CobordismGraph &g, const char *what) {
  for (DerivationId r = 0; r < static_cast<DerivationId>(g.derivationCount()); ++r) {
    std::string why = g.recheck(r);
    CHECK(why.empty(), std::string(what) + ": record " + std::to_string(r) +
                           " rechecks (" + why + ")");
    for (DerivationId c : g.derivation(r).children)
      CHECK(c < r, std::string(what) + ": children are older");
  }
}

// ---------------------------------------------------------------------------

void testBandToDisjointDiscs() {
  // B1: knot K -> (band) -> 2-component link L; L gets "disjoint discs" as a
  // leaf; K becomes slice. With only an annulus for L, K gets genus 1.
  CobordismGraph g;
  LinkId K = g.addLink(1, "K");
  LinkId L = g.addLink(2, "L", std::vector<std::vector<int>>{{0, 0}, {0, 0}});
  g.addCobordism(K, L, pants(), {0}, {0, 1}, "w");
  g.addLeaf(L, Partition::coarsest(2), 0, "annulus");
  g.propagate();
  CHECK_EQ(g.bestConnected(K)->genus, 1, "annulus for L gives K genus 1");
  g.addLeaf(L, Partition::singletons(2), 0, "disjoint discs");
  g.propagate();
  CHECK_EQ(g.bestConnected(K)->genus, 0, "disjoint discs for L make K slice");
  // And the reverse direction gave L an annulus from nothing new: L's own
  // leaf already implied it, so no record was created for it.
  checkAllDerivations(g, "band");
}

void testCycleImprovesAncestor() {
  // B2 (John's example): a genus-1 cobordism L0 -> L1 is found first, L1 has a
  // leaf (slice), so L0 <= 1. Then a search from L1 finds a genus-0 cobordism
  // L1 -> L0. Read backwards, it improves L0 to 0 -- without re-expanding
  // L0, and with a well-founded proof.
  CobordismGraph g;
  LinkId L0 = g.addLink(1, "L0"), L1 = g.addLink(1, "L1");
  g.addLeaf(L1, Partition::coarsest(1), 0, "L1 is slice");
  g.addCobordism(L0, L1, annulusShape(1), {0}, {0}, "g1 edge");
  g.propagate();
  CHECK_EQ(g.bestConnected(L0)->genus, 1, "first L0 <= 1");
  g.addCobordism(L1, L0, annulusShape(0), {0}, {0}, "g0 edge back");
  g.propagate();
  auto best = g.bestConnected(L0);
  CHECK_EQ(best->genus, 0, "the back edge improves L0 to 0");
  CHECK(g.derivation(best->derivation).kind == DerivationKind::cobordismReverse,
        "via the back edge, read backwards");
  auto pf = g.proof(best->derivation);
  CHECK(std::is_sorted(pf.begin(), pf.end()), "proof is in creation order");
  CHECK_EQ(g.derivation(pf.front()).kind == DerivationKind::leaf, true,
           "proof starts at a leaf");
  checkAllDerivations(g, "cycle");
}

void testCycleWithoutLeafGivesNothing() {
  // B3: a cycle of cobordisms with no leaf anywhere derives nothing.
  CobordismGraph g;
  LinkId a = g.addLink(1, "a"), b = g.addLink(1, "b"), c = g.addLink(1, "c");
  g.addCobordism(a, b, annulusShape(0), {0}, {0}, "ab");
  g.addCobordism(b, c, annulusShape(0), {0}, {0}, "bc");
  g.addCobordism(c, a, annulusShape(0), {0}, {0}, "ca");
  g.propagate();
  CHECK_EQ(g.derivationCount(), static_cast<size_t>(0), "no leaves, no records");
}

void testCycleCannotSelfImprove() {
  // B4: a leaf at a, then a cycle a -> b -> a of genus-0 cobordisms: a gets
  // nothing better than its leaf, and nothing is created for a at all.
  CobordismGraph g;
  LinkId a = g.addLink(1, "a"), b = g.addLink(1, "b");
  g.addLeaf(a, Partition::coarsest(1), 2, "a <= 2");
  g.addCobordism(a, b, annulusShape(0), {0}, {0}, "ab");
  g.addCobordism(b, a, annulusShape(0), {0}, {0}, "ba");
  g.propagate();
  CHECK_EQ(g.bestConnected(a)->genus, 2, "a keeps its leaf");
  CHECK(g.derivation(g.bestConnected(a)->derivation).kind == DerivationKind::leaf,
        "and its best record is still the leaf");
  CHECK_EQ(g.bestConnected(b)->genus, 2, "b inherits 2");
  CHECK_EQ(g.saturate(), 0L, "fixed point: nothing more to derive");
}

void testComponentMapsMatter() {
  // B5: a cobordism whose outgoing curves are listed in a different order from
  // the outgoing link's components. L (3 components) bounds pieces {0,1} and {2}
  // with genus 0 (a leaf). The cobordism's outgoing curve j is L's component
  // outMap[j] = {2, 0, 1}: curves 1 and 2 share a piece, curve 0 is alone.
  // Its shape: component A carries the incoming curve and outgoing curves
  // 1, 2; component B is an annulus... not possible from one incoming curve,
  // so use a 2-component incoming link: A: in 0, out 1 and 2; B: in 1, out 0.
  CobordismGraph g;
  LinkId N = g.addLink(2, "N");
  LinkId L = g.addLink(3, "L");
  CobordismShape c;
  c.components = 2;
  c.genus = 0;
  c.inComponent = {0, 1};
  c.outComponent = {1, 0, 0};
  g.addLeaf(L, Partition::fromLabels({0, 0, 1}), 0, "pieces {0,1},{2}");
  g.addCobordism(N, L, c, {0, 1}, {2, 0, 1}, "permuted");
  g.propagate();
  // Curves 1,2 (L's 0,1) are one piece of G: joined to C's component A by
  // two edges -> one cycle, genus 1, and B is alone with curve 0 (L's 2).
  auto sep = g.best(N, Partition::singletons(2));
  CHECK(sep.has_value(), "the permuted map keeps N's components apart");
  CHECK_EQ(sep->genus, 1, "with the cycle through curves 1 and 2");
  // With the identity map instead, curves 0 and 1 (L's 0, 1) share a piece,
  // joining A and B: connected, genus 0. Different answer: the map matters.
  CobordismGraph h;
  LinkId N2 = h.addLink(2, "N"), L2 = h.addLink(3, "L");
  h.addLeaf(L2, Partition::fromLabels({0, 0, 1}), 0, "pieces {0,1},{2}");
  h.addCobordism(N2, L2, c, {0, 1}, {0, 1, 2}, "identity map");
  h.propagate();
  CHECK(!h.best(N2, Partition::singletons(2)).has_value(),
        "identity map: no separated surface");
  CHECK_EQ(h.bestConnected(N2)->genus, 0, "identity map: connected, genus 0");
  checkAllDerivations(g, "maps");
  checkAllDerivations(h, "maps identity");
}

void testSplits() {
  // B6: a split union A u B of two slice knots bounds two discs; a band
  // from K to A u B then makes K slice. And the restriction direction: a
  // surface for A u B that keeps A and B apart gives each piece a bound.
  CobordismGraph g;
  LinkId K = g.addLink(1, "K"), W = g.addLink(2, "A u B"),
         A = g.addLink(1, "A"), B = g.addLink(1, "B");
  g.addSplit(W, {A, B}, {{0}, {1}});
  g.addCobordism(K, W, pants(), {0}, {0, 1}, "band");
  g.addLeaf(A, Partition::coarsest(1), 0, "A slice");
  g.propagate();
  CHECK(!g.bestConnected(K).has_value(), "one piece alone gives nothing");
  g.addLeaf(B, Partition::coarsest(1), 0, "B slice");
  g.propagate();
  CHECK_EQ(g.best(W, Partition::singletons(2))->genus, 0, "A u B: two discs");
  CHECK_EQ(g.bestConnected(K)->genus, 0, "K slice through the split");
  // Restriction: C u D with a leaf keeping pieces apart at genus 3.
  CobordismGraph h;
  LinkId W2 = h.addLink(2, "C u D"), C = h.addLink(1, "C"),
         D = h.addLink(1, "D");
  h.addSplit(W2, {C, D}, {{0}, {1}});
  h.addLeaf(W2, Partition::singletons(2), 3, "separated, genus 3");
  h.addLeaf(W2, Partition::coarsest(2), 0, "annulus");
  h.propagate();
  CHECK_EQ(h.bestConnected(C)->genus, 3, "restriction gives C <= total");
  CHECK_EQ(h.bestConnected(D)->genus, 3, "and D <= total");
  // The annulus (K u -K style) mixes pieces, so it restricts to nothing:
  // no false "each piece is slice".
  CHECK(h.bestConnected(C)->genus != 0, "an annulus across pieces restricts to nothing");
  checkAllDerivations(g, "split");
  checkAllDerivations(h, "restrict");
}

void testContradictionGates() {
  // B7: a derived bound below a proved lower bound is reported; so is a
  // partition the linking numbers forbid (a wrong component map would do
  // this: pretend the Hopf link bounds two discs).
  CobordismGraph g;
  LinkId K = g.addLink(1, "K"), L = g.addLink(2, "L");
  g.setGenusLowerBound(K, 1, "literature");
  g.addCobordism(K, L, pants(), {0}, {0, 1}, "band");
  g.addLeaf(L, Partition::singletons(2), 0, "discs");
  g.propagate();
  CHECK_EQ(static_cast<int>(g.contradictions().size()), 1,
           "below-literature bound reported");
  CobordismGraph h;
  LinkId H = h.addLink(2, "Hopf", std::vector<std::vector<int>>{{0, 1}, {1, 0}});
  h.addLeaf(H, Partition::singletons(2), 0, "wrong");
  CHECK_EQ(static_cast<int>(h.contradictions().size()), 1,
           "linking-forbidden partition reported");
  CobordismGraph ok;
  LinkId H2 = ok.addLink(2, "Hopf", std::vector<std::vector<int>>{{0, 1}, {1, 0}});
  ok.addLeaf(H2, Partition::coarsest(2), 0, "annulus");
  CHECK(ok.contradictions().empty(), "an annulus for the Hopf link is fine");
}

void testSaturationAndMaps() {
  // B8: addCobordism/addSplit reject maps that are not bijections.
  CobordismGraph g;
  LinkId K = g.addLink(1, "K"), L = g.addLink(2, "L");
  bool threw = false;
  try {
    g.addCobordism(K, L, pants(), {0}, {0, 0}, "bad");
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  CHECK(threw, "a non-bijective map is refused");
  threw = false;
  try {
    CobordismShape bad = pants();
    bad.components = 2; // component 1 carries no curve
    g.addCobordism(K, L, bad, {0}, {0, 1}, "closed component");
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  CHECK(threw, "a shape with a boundaryless component is refused");
}

// ---------------------------------------------------------------------------
// B9: random graphs. The fixed point must not depend on the order in which
// edges and leaves arrive, must equal a naive closure computed from
// scratch, must be saturated, and every record must recheck.

struct RandomGraph {
  int links;
  std::vector<int> comps;
  struct W { int in, out; CobordismShape shape; std::vector<int> inMap, outMap; };
  struct Lf { int link; Partition p; int genus; };
  std::vector<W> cobordisms;
  std::vector<Lf> leaves;
};

RandomGraph randomGraph(std::mt19937 &rng) {
  auto u = [&](int lo, int hi) {
    return std::uniform_int_distribution<int>(lo, hi)(rng);
  };
  RandomGraph g;
  g.links = u(2, 6);
  for (int i = 0; i < g.links; ++i)
    g.comps.push_back(u(1, 3));
  int nw = u(1, 9);
  for (int k = 0; k < nw; ++k) {
    RandomGraph::W w;
    w.in = u(0, g.links - 1);
    w.out = u(0, g.links - 1);
    int nin = g.comps[w.in], nout = g.comps[w.out];
    w.shape.components = u(1, std::min(nin, 3));
    w.shape.genus = u(0, 2);
    // Every component gets an incoming curve (as a cobordism's must).
    w.shape.inComponent.resize(nin);
    for (int i = 0; i < nin; ++i)
      w.shape.inComponent[i] = i < w.shape.components ? i : u(0, w.shape.components - 1);
    w.shape.outComponent.resize(nout);
    for (int &x : w.shape.outComponent)
      x = u(0, w.shape.components - 1);
    w.inMap.resize(nin);
    w.outMap.resize(nout);
    std::iota(w.inMap.begin(), w.inMap.end(), 0);
    std::iota(w.outMap.begin(), w.outMap.end(), 0);
    std::shuffle(w.inMap.begin(), w.inMap.end(), rng);
    std::shuffle(w.outMap.begin(), w.outMap.end(), rng);
    g.cobordisms.push_back(w);
  }
  int nl = u(1, 4);
  for (int k = 0; k < nl; ++k) {
    int n = u(0, g.links - 1);
    auto parts = allPartitions(g.comps[n]);
    g.leaves.push_back({n, parts[u(0, static_cast<int>(parts.size()) - 1)], u(0, 3)});
  }
  return g;
}

// The partition genera of every link after building `rg` in the given order.
std::vector<std::set<std::pair<std::vector<int>, int>>>
build(const RandomGraph &rg, const std::vector<int> &order, bool stepwise,
      CobordismGraph *out = nullptr) {
  CobordismGraph g;
  for (int i = 0; i < rg.links; ++i)
    g.addLink(rg.comps[i], "n" + std::to_string(i));
  // order: indices into cobordisms (0..W-1) then leaves (W..W+L-1), mixed.
  const int W = static_cast<int>(rg.cobordisms.size());
  for (int idx : order) {
    if (idx < W) {
      const auto &w = rg.cobordisms[idx];
      g.addCobordism(w.in, w.out, w.shape, w.inMap, w.outMap, "w");
    } else {
      const auto &l = rg.leaves[idx - W];
      g.addLeaf(l.link, l.p, l.genus, "leaf");
    }
    if (stepwise)
      g.propagate();
  }
  g.propagate();
  std::vector<std::set<std::pair<std::vector<int>, int>>> ans(rg.links);
  for (int i = 0; i < rg.links; ++i)
    for (const PartitionGenus &e : g.link(i).partitionGenera.entries())
      ans[i].insert({e.partition.labels(), e.genus});
  if (out)
    *out = std::move(g);
  return ans;
}

// A naive closure, independent of propagate(): repeatedly apply every
// cobordism in both directions to every fact until nothing new appears, then
// take Pareto minima. Uses glue() (tested against its own model above) but
// none of CobordismGraph's bookkeeping.
std::vector<std::set<std::pair<std::vector<int>, int>>>
naiveClosure(const RandomGraph &rg) {
  std::vector<std::set<std::pair<std::vector<int>, int>>> facts(rg.links);
  for (const auto &l : rg.leaves)
    facts[l.link].insert({l.p.labels(), l.genus});
  const int maxGenus = 40; // genus only grows along derivations; cap for termination
  bool changed = true;
  while (changed) {
    changed = false;
    for (const auto &w : rg.cobordisms)
      for (int dir = 0; dir < 2; ++dir) {
        const int from = dir == 0 ? w.out : w.in, to = dir == 0 ? w.in : w.out;
        const auto &gm = dir == 0 ? w.outMap : w.inMap;
        const auto &fm = dir == 0 ? w.inMap : w.outMap;
        auto snapshot = facts[from];
        for (const auto &[labels, genus] : snapshot) {
          Partition p = Partition::fromLabels(labels);
          std::vector<int> cl(gm.size());
          for (size_t j = 0; j < gm.size(); ++j) cl[j] = p.blockOf(gm[j]);
          auto r = glue(w.shape, dir == 0 ? Side::outgoing : Side::incoming,
                        Partition::fromLabels(cl), genus);
          if (!r || r->genus > maxGenus) continue;
          std::vector<int> nl(fm.size());
          for (size_t i = 0; i < fm.size(); ++i) nl[fm[i]] = r->partition.blockOf(i);
          auto np = Partition::fromLabels(nl);
          // Keep the fact if not implied by an existing one.
          bool implied = false;
          for (const auto &[l2, g2] : facts[to])
            if (g2 <= r->genus && Partition::fromLabels(l2).refines(np))
              implied = true;
          if (!implied && facts[to].insert({np.labels(), r->genus}).second)
            changed = true;
        }
      }
  }
  // Pareto minima.
  std::vector<std::set<std::pair<std::vector<int>, int>>> ans(rg.links);
  for (int i = 0; i < rg.links; ++i)
    for (const auto &[l, gnum] : facts[i]) {
      bool dominated = false;
      for (const auto &[l2, g2] : facts[i])
        if (std::make_pair(l2, g2) != std::make_pair(l, gnum) && g2 <= gnum &&
            Partition::fromLabels(l2).refines(Partition::fromLabels(l)))
          dominated = true;
      if (!dominated)
        ans[i].insert({l, gnum});
    }
  return ans;
}

void testRandomFixedPoints() {
  std::mt19937 rng(424242);
  int graphs = 0, withCycles = 0;
  for (int t = 0; t < 1500; ++t) {
    RandomGraph rg = randomGraph(rng);
    const int total = static_cast<int>(rg.cobordisms.size() + rg.leaves.size());
    std::vector<int> order(total);
    std::iota(order.begin(), order.end(), 0);
    CobordismGraph g0;
    auto base = build(rg, order, /*stepwise=*/false, &g0);
    for (int k = 0; k < 3; ++k) {
      std::shuffle(order.begin(), order.end(), rng);
      auto other = build(rg, order, /*stepwise=*/k % 2 == 0);
      CHECK(other == base, "fixed point is independent of arrival order");
    }
    CHECK(naiveClosure(rg) == base, "fixed point equals the naive closure");
    CHECK_EQ(g0.saturate(), 0L, "propagate() reached a fixed point");
    checkAllDerivations(g0, "random");
    for (DerivationId r = 0; r < static_cast<DerivationId>(g0.derivationCount()); ++r) {
      auto pf = g0.proof(r);
      CHECK(pf.back() == r, "proof ends at its record");
      for (DerivationId x : pf)
        for (DerivationId c : g0.derivation(x).children)
          CHECK(std::binary_search(pf.begin(), pf.end(), c),
                "proof is closed under children");
    }
    ++graphs;
    // Count graphs with a directed cycle through cobordisms, to be sure
    // cycles are exercised.
    std::vector<std::vector<int>> adj(rg.links);
    for (const auto &w : rg.cobordisms) adj[w.in].push_back(w.out);
    std::vector<int> state(rg.links, 0);
    bool cyc = false;
    std::function<void(int)> dfs = [&](int v) {
      state[v] = 1;
      for (int x : adj[v]) {
        if (state[x] == 1) cyc = true;
        else if (state[x] == 0) dfs(x);
      }
      state[v] = 2;
    };
    for (int v = 0; v < rg.links; ++v)
      if (!state[v]) dfs(v);
    if (cyc) ++withCycles;
  }
  CHECK(graphs == 1500, "all random graphs ran");
  CHECK(withCycles > 300, "random graphs include cycles");
}

// ---------------------------------------------------------------------------
// Lower bounds (README.md, "Lower bounds").

void testLowerConcordance() {
  // L1: K --annulus--> K', literature g4(K') >= 1: then g4(K) >= 1 (genus 0
  // edge), and nothing from a genus-1 edge.
  for (int g : {0, 1}) {
    CobordismGraph pg;
    LinkId K = pg.addLink(1, "K"), K2 = pg.addLink(1, "K'");
    pg.addCobordism(K, K2, annulusShape(g), {0}, {0}, "w");
    pg.setGenusLowerBound(K2, 1, "literature");
    pg.propagateLower();
    CHECK_EQ(pg.lower(K, Partition::coarsest(1)), 1 - g, "lower bound across an annulus");
    CHECK(pg.contradictions().empty(), "no contradiction");
  }
}

void testLowerPaperCases() {
  // L2: the paper's reverse inequality g4(L0) >= g4(L1) - g - n0 + 1 as
  // special cases. A band K -> L (n0 = 1): g4(K) >= g4(L).
  {
    CobordismGraph pg;
    LinkId K = pg.addLink(1, "K"), L = pg.addLink(2, "L");
    pg.addCobordism(K, L, pants(), {0}, {0, 1}, "band");
    pg.setGenusLowerBound(L, 1, "literature");
    pg.propagateLower();
    CHECK_EQ(pg.lower(K, Partition::coarsest(1)), 1, "band: g4(K) >= g4(L)");
  }
  // A band merging the two components of L0 into K (n0 = 2): the real
  // penalty 1 applies, g4(L0) >= g4(K) - 1.
  {
    CobordismGraph pg;
    LinkId L0 = pg.addLink(2, "L0"), K = pg.addLink(1, "K");
    CobordismShape merge;
    merge.components = 1;
    merge.genus = 0;
    merge.inComponent = {0, 0};
    merge.outComponent = {0};
    pg.addCobordism(L0, K, merge, {0, 1}, {0}, "merge");
    pg.setGenusLowerBound(K, 2, "literature");
    pg.propagateLower();
    CHECK_EQ(pg.lower(L0, Partition::coarsest(2)), 1, "merging band: penalty 1");
  }
}

void testLowerNoPenaltyForAnnuli() {
  // L3: two annuli L0 -> L1 (each piece carries one component of L0): the
  // lower bound transports with NO penalty, where the solvers' n0 - 1 = 1
  // would lose it.
  CobordismGraph pg;
  LinkId L0 = pg.addLink(2, "L0"), L1 = pg.addLink(2, "L1");
  pg.addCobordism(L0, L1, CobordismShape::product(2), {0, 1}, {0, 1}, "annuli");
  pg.setGenusLowerBound(L1, 1, "literature");
  pg.propagateLower();
  CHECK_EQ(pg.lower(L0, Partition::coarsest(2)), 1, "annuli: no penalty");
}

void testLowerSplit() {
  // L4: W = A u B. Whole from pieces (non-mixing partitions only), and a
  // piece from the whole minus a proved surface for the other piece.
  CobordismGraph pg;
  LinkId W = pg.addLink(2, "A u B"), A = pg.addLink(1, "A"), B = pg.addLink(1, "B");
  pg.addSplit(W, {A, B}, {{0}, {1}});
  pg.setGenusLowerBound(A, 2, "lit A");
  pg.setGenusLowerBound(B, 1, "lit B");
  pg.propagateLower();
  CHECK_EQ(pg.lower(W, Partition::singletons(2)), 3, "separated pieces: bounds add");
  CHECK_EQ(pg.lower(W, Partition::coarsest(2)), 0,
           "a connected surface may mix pieces: no additive bound (K u -K)");
  CobordismGraph ph;
  LinkId W2 = ph.addLink(2, "C u D"), C = ph.addLink(1, "C"), D = ph.addLink(1, "D");
  ph.addSplit(W2, {C, D}, {{0}, {1}});
  ph.addLeaf(D, Partition::coarsest(1), 1, "D bounds genus 1");
  ph.setGenusLowerBound(W2, 3, "lit whole");
  ph.propagate();
  ph.propagateLower();
  CHECK_EQ(ph.lower(C, Partition::coarsest(1)), 2,
           "piece from whole: g4(C) >= g4(C u D) - g(D-surface)");
  CHECK(ph.contradictions().empty(), "no contradiction");
}

void testLowerWhatIf() {
  // L6: the lower report's what-if (runrecords::writeLowerReport()):
  // forget every lower bound, seed one link, and read what reaches the
  // target. A concordance carries the seed whole, a merging band loses 1,
  // and the literature bound elsewhere no longer contributes.
  CobordismGraph pg;
  LinkId T = pg.addLink(2, "T"), C = pg.addLink(2, "C"), K = pg.addLink(1, "K"),
         O = pg.addLink(1, "O");
  pg.addCobordism(T, C, CobordismShape::product(2), {0, 1}, {0, 1}, "annuli");
  CobordismShape merge;
  merge.components = 1;
  merge.genus = 0;
  merge.inComponent = {0, 0};
  merge.outComponent = {0};
  pg.addCobordism(T, K, merge, {0, 1}, {0}, "merge");
  pg.setGenusLowerBound(O, 5, "literature, unconnected");
  pg.setGenusLowerBound(T, 1, "the target's own literature");
  pg.propagateLower();
  auto reach = [&](LinkId n, int seed) {
    CobordismGraph what = pg;
    what.clearLowerBounds();
    what.setGenusLowerBound(n, seed, "what-if");
    what.propagateLower();
    return what.lower(T, Partition::coarsest(2));
  };
  CHECK_EQ(reach(C, 1), 1, "what-if: a concordance carries the seed whole");
  CHECK_EQ(reach(K, 2), 1, "what-if: a merging band loses 1");
  CHECK_EQ(reach(O, 5), 0, "what-if: an unconnected node carries nothing");
  CobordismGraph cleared = pg;
  cleared.clearLowerBounds();
  CHECK_EQ(cleared.lower(T, Partition::coarsest(2)), 0,
           "cleared: the target's own literature bound is forgotten too");
  CHECK_EQ(pg.lower(T, Partition::coarsest(2)), 1, "the original graph is untouched");
}

void testLowerTransportMonotone() {
  // L7 (lem:transport-monotone): what a cobordism transports for surfaces
  // with EXACTLY partition p never falls when p is refined, so the bound for
  // surfaces refining q is transportedLower() at q itself. lowerAcross()
  // relies on this instead of minimising over refinements; if the claim
  // ever broke (a shape or a map for which refining the cap raised the
  // addition), this test fails. Checked at the fixed point of random worlds,
  // both directions of every edge, every q, with the lower bounds relaxed
  // from random literature values the world allows, not only true minima.
  std::mt19937 rng(4242);
  long checked = 0, strict = 0;
  for (int t = 0; t < 400; ++t) {
    RandomGraph rg = randomGraph(rng);
    const int total = static_cast<int>(rg.cobordisms.size() + rg.leaves.size());
    std::vector<int> order(total);
    std::iota(order.begin(), order.end(), 0);
    CobordismGraph g;
    build(rg, order, false, &g);
    std::uniform_int_distribution<int> lit(0, 3);
    for (int n = 0; n < rg.links; ++n) {
      int cap = 3;
      if (auto b = g.bestConnected(n)) cap = std::min(cap, b->genus);
      const int v = std::min(lit(rng), cap);
      if (v > 0) g.setGenusLowerBound(n, v, "literature");
    }
    g.propagateLower();
    if (!g.contradictions().empty()) continue;
    for (RelationId e = 0; e < static_cast<RelationId>(rg.cobordisms.size()); ++e)
      for (bool toIsIn : {true, false}) {
        const LinkCobordism &w = g.cobordism(e);
        const LinkId to = toIsIn ? w.in : w.out;
        const int k = g.link(to).components;
        if (k > CobordismGraph::kMaxLowerComponents) continue;
        for (const Partition &q : allPartitions(k)) {
          const int atQ = g.transportedLower(w, toIsIn, q);
          int minOverRefinements = CobordismGraph::kNoSurface;
          for (const Partition &p : allPartitions(k))
            if (p.refines(q)) {
              const int v = g.transportedLower(w, toIsIn, p);
              minOverRefinements = std::min(minOverRefinements, v);
              if (v > atQ) ++strict;
            }
          CHECK_EQ(minOverRefinements, atQ,
                   "transport is monotone: the min over refinements is at q");
          ++checked;
        }
      }
  }
  CHECK(checked > 10000, "monotonicity was checked on many (edge, q) pairs");
  CHECK(strict > 100, "refinements often transport strictly more, so the check bites");
}

void testLowerIf() {
  // L8: lowerIf() seeds a copy that KEEPS the graph's own bounds (unlike the
  // report's cleared what-if in L6): a seed on one route and a literature
  // bound on another can complete a bound together; a seed above a proved
  // surface is refused; the graph is untouched.
  CobordismGraph pg;
  LinkId T = pg.addLink(2, "T"), C = pg.addLink(2, "C"), K = pg.addLink(1, "K");
  pg.addCobordism(T, C, CobordismShape::product(2), {0, 1}, {0, 1}, "annuli");
  CobordismShape merge;
  merge.components = 1;
  merge.genus = 0;
  merge.inComponent = {0, 0};
  merge.outComponent = {0};
  pg.addCobordism(T, K, merge, {0, 1}, {0}, "merge");
  pg.setGenusLowerBound(K, 2, "literature");
  pg.propagate();
  pg.propagateLower();
  const Partition goal = Partition::coarsest(2);
  CHECK_EQ(pg.lower(T, goal), 1, "K's bound loses 1 across the band");
  auto r = pg.lowerIf({{C, goal, 3}}, T, goal);
  CHECK(r && *r == 3, "a seed on C carries whole across the annuli, K's bound kept");
  r = pg.lowerIf({{C, goal, 0}}, T, goal);
  CHECK(r && *r == 1, "a seed below what is known changes nothing");
  pg.addLeaf(C, goal, 1, "proved");
  pg.propagate();
  r = pg.lowerIf({{C, goal, 3}}, T, goal);
  CHECK(!r, "a seed above a proved surface is refused");
  r = pg.lowerIf({{C, goal, 1}}, T, goal);
  CHECK(r && *r == 1, "a consistent seed is read");
  CHECK_EQ(pg.lower(T, goal), 1, "the graph is untouched");
  CHECK(pg.contradictions().empty(), "no contradiction leaks into the graph");
}

void testLowerWhy() {
  // L9: every raised lower bound remembers the fact that raised it, and
  // following those facts from the target reaches a literature leaf: T's
  // bound came across the merging band from K's literature value, with the
  // band's addition of 1 recorded; lowerVersion() changes on every raise.
  CobordismGraph pg;
  LinkId T = pg.addLink(2, "T"), K = pg.addLink(1, "K");
  CobordismShape merge;
  merge.components = 1;
  merge.genus = 0;
  merge.inComponent = {0, 0};
  merge.outComponent = {0};
  const RelationId e = pg.addCobordism(T, K, merge, {0, 1}, {0}, "merge");
  const long v0 = pg.lowerVersion();
  pg.setGenusLowerBound(K, 2, "literature");
  pg.propagateLower();
  CHECK(pg.lowerVersion() > v0, "raising a bound changes the version");
  auto why = pg.lowerWhy(T, Partition::coarsest(2));
  CHECK_EQ(why.value, 1, "T's bound is 1");
  CHECK(why.reason.kind == CobordismGraph::LowerReason::Kind::cobordism, "it came across a witness");
  CHECK(why.reason.relation == e && why.reason.toIsIn, "the merging band, read at its incoming end");
  CHECK_EQ(why.reason.from, 2, "from K's bound 2");
  CHECK_EQ(why.reason.addition, 1, "the cap added 1");
  auto leaf = pg.lowerWhy(K, Partition::fromLabels(why.reason.fromPartition));
  CHECK(leaf.reason.kind == CobordismGraph::LowerReason::Kind::literature, "K's is a literature leaf");
  CHECK_EQ(leaf.value, 2, "with value 2");
  auto none = pg.lowerWhy(T, Partition::singletons(2));
  CHECK_EQ(none.value, 2, "the singleton partition transports 2 (no addition)");
  CobordismGraph empty;
  LinkId X = empty.addLink(1, "X");
  CHECK(empty.lowerWhy(X, Partition::coarsest(1)).reason.kind ==
            CobordismGraph::LowerReason::Kind::none,
        "nothing known: kind none");
}

void testSums() {
  // S1 (paper lem:sum-partitions): a sum along components. Two Hopf links
  // H1, H2 (each realizes {0,1} at genus 0: an annulus; singletons are
  // forbidden by lk = 1) summed along a component form a 3-component chain
  // W: W's components 0,1 from H1, 1,2 from H2 (component 1 shared). The
  // annuli boundary-sum to a pair of pants: W realizes {0,1,2} at genus 0.
  // A knot K of genus 1 summed into a component of H1 (K #_c H1) realizes
  // {0,1} at genus 1: the paper's lem:sum-along-components(i). The lower
  // rule (cor:sum-pieces): lower(K) = 1 and H1's annulus (genus 0, 2
  // components) give lower(K #_c H1) >= 1 - (0 + 2 - 1) = 0, nothing; with
  // lower(K) = 3 it gives 2. Every record rechecks.
  CobordismGraph pg;
  LinkId H1 = pg.addLink(2, "H1", std::vector<std::vector<int>>{{0, 1}, {1, 0}});
  LinkId H2 = pg.addLink(2, "H2", std::vector<std::vector<int>>{{0, 1}, {1, 0}});
  LinkId W = pg.addLink(3, "W");
  pg.addLeaf(H1, Partition::coarsest(2), 0, "annulus");
  pg.addLeaf(H2, Partition::coarsest(2), 0, "annulus");
  pg.addSum(W, {H1, H2}, {{0, 1}, {1, 2}});
  pg.propagate();
  auto b = pg.best(W, Partition::coarsest(3));
  CHECK(b && b->genus == 0, "a chain of two Hopf links bounds a planar surface");
  CHECK(pg.derivation(b->derivation).kind == DerivationKind::sumCombine, "by the sum rule");
  LinkId K = pg.addLink(1, "K"), KH = pg.addLink(2, "K#H1");
  pg.addLeaf(K, Partition::coarsest(1), 1, "genus 1");
  pg.addSum(KH, {K, H1}, {{0}, {0, 1}});
  pg.propagate();
  auto c = pg.best(KH, Partition::coarsest(2));
  CHECK(c && c->genus == 1, "K #_c H1 realizes genus 1 with one piece");
  pg.setGenusLowerBound(K, 1, "literature");
  pg.propagateLower();
  CHECK_EQ(pg.lower(KH, Partition::coarsest(2)), 0, "lower 1 - (0 + 2 - 1) = 0");
  // A knot K2 with literature lower bound 3 and no proved surface, summed
  // into H1: lower(K2 #_c H1) >= 3 - (0 + 2 - 1) = 2.
  LinkId K2 = pg.addLink(1, "K2"), K2H = pg.addLink(2, "K2#H1");
  pg.addSum(K2H, {K2, H1}, {{0}, {0, 1}});
  pg.setGenusLowerBound(K2, 3, "literature");
  pg.propagate();
  pg.propagateLower();
  CHECK_EQ(pg.lower(K2H, Partition::coarsest(2)), 2, "lower 3 - 1 = 2");
  auto why = pg.lowerWhy(K2H, Partition::coarsest(2));
  CHECK(why.reason.kind == CobordismGraph::LowerReason::Kind::sumPiece && why.reason.piece == 0,
        "reason: the sum rule from piece 0");
  checkAllDerivations(pg, "sum edges");
  CHECK(pg.contradictions().empty(), "no contradiction");
  // Pieces that are not summed at all (a map missing a whole component) are refused.
  bool refused = false;
  try {
    pg.addSum(W, {H1, H2}, {{0, 1}, {0, 1}});
  } catch (const std::invalid_argument &) {
    refused = true;
  }
  CHECK(refused, "a sum whose pieces miss a component of the whole is refused");
  // Two Hopf links summed along BOTH pairs of components would be a cycle
  // of sum sites, which no sphere decomposition gives: refused.
  LinkId W2 = pg.addLink(2, "W2");
  refused = false;
  try {
    pg.addSum(W2, {H1, H2}, {{0, 1}, {0, 1}});
  } catch (const std::invalid_argument &) {
    refused = true;
  }
  CHECK(refused, "sum sites forming a cycle are refused");
}

void testLowerSoundOnRandomWorlds() {
  // L5: take a random graph's upper-bound closure as the whole world, give
  // every link its true connected minimum as a literature lower bound, and
  // propagate: no lower bound may exceed a surface the world contains.
  std::mt19937 rng(777);
  int worlds = 0, raised = 0;
  for (int t = 0; t < 600; ++t) {
    RandomGraph rg = randomGraph(rng);
    const int total = static_cast<int>(rg.cobordisms.size() + rg.leaves.size());
    std::vector<int> order(total);
    std::iota(order.begin(), order.end(), 0);
    CobordismGraph g;
    build(rg, order, false, &g);
    for (int n = 0; n < rg.links; ++n)
      if (auto b = g.bestConnected(n))
        g.setGenusLowerBound(n, b->genus, "true minimum");
    raised += g.propagateLower();
    CHECK(g.contradictions().empty(), "sound: no contradiction in a consistent world");
    for (int n = 0; n < rg.links; ++n)
      for (const PartitionGenus &e : g.link(n).partitionGenera.entries())
        CHECK(g.lower(n, e.partition) <= e.genus,
              "no lower bound exceeds an achieved surface");
    ++worlds;
  }
  CHECK(worlds == 600, "all worlds ran");
  CHECK(raised > 100, "lower bounds were actually transported");
}

void testPartitionGeneraFields() {
  // P1: profiles.jsonl's fields. B1's graph: the band's Pareto sets (L's
  // disjoint discs dominate its annulus, K's disc its genus-1 surface), a
  // Hopf link whose split partition the linking numbers forbid, and a
  // literature bound read back per partition.
  CobordismGraph g;
  LinkId K = g.addLink(1, "K");
  LinkId L = g.addLink(2, "L", std::vector<std::vector<int>>{{0, 0}, {0, 0}});
  g.addCobordism(K, L, pants(), {0}, {0, 1}, "w");
  g.addLeaf(L, Partition::coarsest(2), 0, "annulus");
  g.propagate();
  g.addLeaf(L, Partition::singletons(2), 0, "disjoint discs");
  g.propagate();
  LinkId H = g.addLink(2, "Hopf", std::vector<std::vector<int>>{{0, 1}, {1, 0}});
  g.setGenusLowerBound(H, 0, "table");
  LinkId X = g.addLink(1, "X");
  g.setGenusLowerBound(X, 2, "table");
  g.propagateLower();
  const long rK = g.bestConnected(K)->derivation;
  const long rL = g.best(L, Partition::singletons(2))->derivation;
  CHECK_EQ(g.partitionGeneraFields(K),
           std::string("\"components\":1,\"entries\":[{\"p\":\"{0}\",\"g\":0,\"r\":") +
               std::to_string(rK) + "}],\"lower\":[]",
           "a slice knot: its disc alone, no lower bound");
  CHECK_EQ(g.partitionGeneraFields(L),
           std::string("\"components\":2,\"linking\":[[0,0],[0,0]],\"entries\":[{\"p\":"
                       "\"{0}{1}\",\"g\":0,\"r\":") +
               std::to_string(rL) + "}],\"lower\":[]",
           "disjoint discs dominate the annulus");
  CHECK_EQ(g.partitionGeneraFields(H),
           std::string("\"components\":2,\"linking\":[[0,1],[1,0]],\"genus_lower\":0,"
                       "\"entries\":[],\"lower\":[{\"p\":\"{0}{1}\",\"forbidden\":true}]"),
           "the Hopf link's components cannot bound disjoint surfaces");
  CHECK_EQ(g.partitionGeneraFields(X),
           std::string("\"components\":1,\"genus_lower\":2,\"entries\":[],"
                       "\"lower\":[{\"p\":\"{0}\",\"lo\":2}]"),
           "a literature bound, per partition");
}

} // namespace

int main() {
  testPartitionGeneraFields();
  testLowerConcordance();
  testLowerPaperCases();
  testLowerNoPenaltyForAnnuli();
  testLowerSplit();
  testLowerWhatIf();
  testLowerTransportMonotone();
  testLowerIf();
  testLowerWhy();
  testSums();
  testLowerSoundOnRandomWorlds();
  testBandToDisjointDiscs();
  testCycleImprovesAncestor();
  testCycleWithoutLeafGivesNothing();
  testCycleCannotSelfImprove();
  testComponentMapsMatter();
  testSplits();
  testContradictionGates();
  testSaturationAndMaps();
  testRandomFixedPoints();
  return checks::finish("proofgraph_test");
}
