// nodes_test.cpp
//
// Tests for cascade/nodes.h: node identity is exact, and the component map
// it returns is right (README.md, "Component maps").

#include <algorithm>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include <link/examplelink.h>
#include <link/link.h>

#include "../diagramiso.h"
#include "../hopedges.h"
#include "../nodes.h"
#include "check.h"
#include "exactnaming/exacttables.h"

using exactnaming::GaussDiagram;
using namespace cascade;

namespace {

GaussDiagram of(const regina::Link &l) {
  std::vector<size_t> origin(l.countComponents());
  std::iota(origin.begin(), origin.end(), 0);
  return GaussDiagram::of(l, origin);
}

// Real table rows (the same strings exactnaming_test uses).
const char *L4a1_0 = "PD[X[6; 1; 7; 2]; X[8; 3; 5; 4]; X[2; 5; 3; 6]; X[4; 7; 1; 8]]";
const char *L4a1_1 = "PD[X[6; 2; 7; 1]; X[8; 4; 5; 3]; X[2; 8; 3; 7]; X[4; 6; 1; 5]]";
const char *L7n1_0 = "PD[X[6; 1; 7; 2]; X[12; 7; 13; 8]; X[4; 13; 1; 14]; X[5; 10; 6; 11]; "
                     "X[3; 8; 4; 9]; X[9; 14; 10; 5]; X[11; 2; 12; 3]]";
const char *L7n1_1 = "PD[X[12; 2; 13; 1]; X[6; 11; 7; 12]; X[4; 6; 1; 5]; X[13; 8; 14; 9]; "
                     "X[3; 11; 4; 10]; X[9; 14; 10; 5]; X[7; 3; 8; 2]]";

// Random classical Reidemeister II/III moves, then simplify(): usually a
// different diagram of the same link, with `origin` carried along.
GaussDiagram scramble(const GaussDiagram &d, std::mt19937 &rng, int moves) {
  regina::Link l = d.link();
  auto u = [&](size_t n) { return std::uniform_int_distribution<size_t>(0, n - 1)(rng); };
  for (int k = 0; k < moves && l.size() > 0; ++k) {
    regina::Crossing *c = l.crossing(u(l.size()));
    regina::StrandRef a = c->strand(u(2));
    if (u(2)) {
      regina::Crossing *c2 = l.crossing(u(l.size()));
      regina::StrandRef b = c2->strand(u(2));
      l.r2(a, static_cast<int>(u(2)), b, static_cast<int>(u(2)));
    } else {
      l.r3(a, static_cast<int>(u(2)));
    }
  }
  GaussDiagram s = GaussDiagram::of(l, d.origin);
  return simplifyKeepingComponents(s);
}

void testLinking() {
  using regina::ExampleLink;
  auto hopf = linkingMatrix(of(ExampleLink::hopf()));
  CHECK(std::abs(hopf[0][1]) == 1 && hopf[0][1] == hopf[1][0], "Hopf: lk = +-1");
  auto t24 = linkingMatrix(of(ExampleLink::torus(2, 4)));
  CHECK(std::abs(t24[0][1]) == 2, "T(2,4): lk = +-2");
  auto wh = linkingMatrix(of(ExampleLink::whitehead()));
  CHECK_EQ(wh[0][1], 0, "Whitehead: lk = 0");
}

void testSimplifyKeepsComponents() {
  // simplify() never permutes or reverses components: pairwise linking
  // numbers (checked inside simplifyKeepingComponents) and origins survive
  // hundreds of random scrambles.
  using regina::ExampleLink;
  std::mt19937 rng(21);
  std::vector<regina::Link> ls = {ExampleLink::torus(3, 6), ExampleLink::borromean(),
                                  ExampleLink::whitehead(), ExampleLink::torus(2, 6)};
  int ok = 0;
  for (const auto &l : ls)
    for (int t = 0; t < 60; ++t) {
      GaussDiagram d = of(l);
      bool threw = false;
      try {
        GaussDiagram s = scramble(d, rng, 12);
        CHECK(s.origin == d.origin, "origins carried through simplify");
        ++ok;
      } catch (const std::logic_error &e) {
        threw = true;
      }
      CHECK(!threw, "simplify kept every pairwise linking number");
    }
  CHECK(ok > 200, "scrambles ran");
}

void testDiagramHits() {
  // The same diagram relabelled: a "diagram" hit whose map is the relabelling.
  using regina::ExampleLink;
  ProofGraph g;
  NodeRegistry reg(g);
  regina::Link u = ExampleLink::hopf();
  u.insertLink(ExampleLink::torus(2, 4));
  // Two split pieces here; intern one connected piece: T(2,4) alone and a
  // three-component connected link, T(3,6).
  GaussDiagram t36 = of(ExampleLink::torus(3, 6));
  NodeMatch first = reg.intern(simplifyKeepingComponents(t36), "T(3,6)");
  CHECK(first.created, "first sight creates a node");
  // Reverse the component order: the map must follow.
  GaussDiagram perm = t36;
  std::reverse(perm.comps.begin(), perm.comps.end());
  std::reverse(perm.origin.begin(), perm.origin.end());
  NodeMatch again = reg.intern(perm, "T(3,6) permuted");
  CHECK(!again.created && again.node == first.node, "same diagram, same node");
  CHECK_EQ(again.method, std::string("diagram"), "found as a diagram");
  // perm's component i is t36's component 2-i; the node's components are
  // t36's. T(3,6) has symmetries permuting components, so check the map is
  // an isomorphism rather than a specific permutation:
  GaussDiagram mapped;
  mapped.signs = perm.signs;
  mapped.comps.assign(3, {});
  for (int i = 0; i < 3; ++i) mapped.comps[again.componentMap[i]] = perm.comps[i];
  CHECK(findDiagramIsomorphism(reg.info(first.node).diagram, mapped, false, true)
            .has_value(),
        "the returned map realises an isomorphism onto the node");
}

void testDifferentDiagramsSameLink() {
  // Scrambled diagrams of one hyperbolic link: the same node, by diagram or
  // by isometry, with a map consistent with the tracked origins up to a
  // symmetry of the link.
  using regina::ExampleLink;
  std::mt19937 rng(22);
  int isometryHits = 0, trials = 0;
  for (const regina::Link &l : {ExampleLink::whitehead(), ExampleLink::borromean(),
                                ExampleLink::conway()}) {
    ProofGraph g;
    NodeRegistry reg(g);
    GaussDiagram base = simplifyKeepingComponents(of(l));
    NodeMatch n0 = reg.intern(base, "base");
    for (int t = 0; t < 25; ++t) {
      GaussDiagram s = scramble(base, rng, 20);
      for (const GaussDiagram &q : exactnaming::splitPieces(s)) {
        if (q.components() != base.components()) continue; // came apart: skip
        NodeMatch m = reg.intern(q, "scrambled");
        ++trials;
        CHECK_EQ(m.node, n0.node, "scrambled diagram is the same node");
        if (m.method == "isometry") ++isometryHits;
        // Linking numbers must transport through the map (up to the global
        // sign a mirror gives).
        auto lq = linkingMatrix(q), ln = reg.info(n0.node).linking;
        bool okMap = true;
        const int sgn = m.mirrored ? -1 : 1;
        for (size_t i = 0; i < q.components(); ++i)
          for (size_t j = 0; j < q.components(); ++j)
            if (i != j && lq[i][j] * sgn != ln[m.componentMap[i]][m.componentMap[j]])
              okMap = false;
        CHECK(okMap, "linking numbers transport through the component map");
      }
    }
  }
  CHECK(trials > 30, "enough scrambled diagrams");
  CHECK(isometryHits > 0, "some matches needed the isometry test");
}

void testOrientationVariantsAreDifferentNodes() {
  // L7n1{0} (g4 2) and L7n1{1} (g4 0), L4a1{0} (g4 0) and L4a1{1} (g4 1):
  // one link with different component orientations. Merging them would be
  // unsound, so they must be different nodes.
  ProofGraph g;
  NodeRegistry reg(g);
  for (auto [a, b] : {std::pair{L7n1_0, L7n1_1}, std::pair{L4a1_0, L4a1_1}}) {
    GaussDiagram da = simplifyKeepingComponents(of(exactnaming::linkFromTablePD(a)));
    GaussDiagram db = simplifyKeepingComponents(of(exactnaming::linkFromTablePD(b)));
    NodeMatch ma = reg.intern(da, a), mb = reg.intern(db, b);
    CHECK(ma.node != mb.node, "orientation variants are different nodes");
    // And each variant, re-interned, finds itself.
    CHECK_EQ(reg.intern(da, a).node, ma.node, "variant 0 finds itself");
    CHECK_EQ(reg.intern(db, b).node, mb.node, "variant 1 finds itself");
  }
}

void testUnknot() {
  ProofGraph g;
  NodeRegistry reg(g);
  GaussDiagram u;
  u.comps = {{}};
  u.origin = {3};
  NodeMatch m = reg.intern(u, "crossingless");
  CHECK_EQ(m.node, reg.unknot(), "a crossingless piece is the unknot");
  CHECK_EQ(g.bestConnected(m.node)->genus, 0, "the unknot bounds a disc");
  NodeMatch m2 = reg.intern(u, "again");
  CHECK_EQ(m2.node, m.node, "one unknot node");
}

// `a` summed along its component `onComponent` with the knot `k` through a
// twist: that component's word becomes `c(over) a... c(under) k...`, where
// crossing 0 is the twist. Nugatory by construction, with all of `a` on the
// side it cuts off (a's other components included, being linked only with
// that side).
GaussDiagram sumThroughTwist(const GaussDiagram &a, size_t onComponent,
                             const GaussDiagram &k, int sign) {
  GaussDiagram d;
  d.signs.push_back(sign);
  d.signs.insert(d.signs.end(), a.signs.begin(), a.signs.end());
  d.signs.insert(d.signs.end(), k.signs.begin(), k.signs.end());
  auto shift = [](long x, long by) { return x > 0 ? x + by : x - by; };
  const long na = static_cast<long>(a.crossings());
  for (size_t c = 0; c < a.components(); ++c) {
    std::vector<long> w;
    if (c == onComponent) w.push_back(+1);
    for (long x : a.comps[c]) w.push_back(shift(x, 1));
    if (c == onComponent) {
      w.push_back(-1);
      for (long x : k.comps[0]) w.push_back(shift(x, 1 + na));
    }
    d.comps.push_back(w);
  }
  d.origin = a.origin;
  return d;
}

bool sameDiagram(const GaussDiagram &a, const GaussDiagram &b) {
  return a.signs == b.signs && a.comps == b.comps && a.origin == b.origin;
}

// removeNugatoryCrossings(): every nugatory crossing goes, nothing else
// changes -- the link (Jones polynomial), every linking number, the
// component order -- and a reduced diagram comes back untouched.
void testNugatoryCrossings() {
  using regina::ExampleLink;
  const GaussDiagram trefoil = of(ExampleLink::trefoilLeft());
  const GaussDiagram fig8 = of(ExampleLink::figureEight());
  const GaussDiagram whitehead = of(ExampleLink::whitehead());
  const GaussDiagram borromean = of(ExampleLink::borromean());

  for (const GaussDiagram *r : {&trefoil, &fig8, &whitehead, &borromean}) {
    CHECK(!nugatoryCrossing(*r), "a reduced diagram has no nugatory crossing");
    CHECK(sameDiagram(removeNugatoryCrossings(*r), *r),
          "a reduced diagram comes back untouched");
  }

  // Connected sums through a twist, knots and links, both twist signs.
  struct Case { const char *name; const GaussDiagram *a; size_t on; const GaussDiagram *k; };
  const Case cases[] = {{"3_1 # 4_1", &trefoil, 0, &fig8},
                        {"4_1 # 3_1", &fig8, 0, &trefoil},
                        {"whitehead(0) # 3_1", &whitehead, 0, &trefoil},
                        {"whitehead(1) # 4_1", &whitehead, 1, &fig8},
                        {"borromean(2) # 3_1", &borromean, 2, &trefoil}};
  for (const Case &cs : cases)
    for (int sign : {+1, -1}) {
      const std::string name = std::string(cs.name) + (sign > 0 ? " (+)" : " (-)");
      const GaussDiagram d = sumThroughTwist(*cs.a, cs.on, *cs.k, sign);
      auto c = nugatoryCrossing(d);
      CHECK(c && *c == 0, name + ": the twist is found nugatory");
      const GaussDiagram r = removeNugatoryCrossings(d);
      CHECK_EQ(r.crossings(), cs.a->crossings() + cs.k->crossings(),
               name + ": exactly the twist goes");
      CHECK(!nugatoryCrossing(r), name + ": the result is reduced");
      CHECK_EQ(r.components(), d.components(), name + ": the components stay");
      CHECK(r.origin == d.origin, name + ": origins stay");
      CHECK(linkingMatrix(r) == linkingMatrix(d), name + ": every linking number stays");
      CHECK(r.link().jones() == d.link().jones(), name + ": the same link (Jones)");
      CHECK(r.link().jones() == cs.a->link().jones() * cs.k->link().jones(),
            name + ": the connected sum (Jones is multiplicative)");
    }

  // Kinks, from Regina's own type I moves, anywhere, either way round.
  std::mt19937 rng(7);
  for (const GaussDiagram *base : {&fig8, &whitehead}) {
    for (int trial = 0; trial < 20; ++trial) {
      regina::Link l = base->link();
      for (int t = 0; t < 3; ++t) {
        regina::Crossing *x = l.crossing(rng() % l.size());
        l.r1(rng() % 2 ? x->upper() : x->lower(), static_cast<int>(rng() % 2),
             rng() % 2 ? 1 : -1);
      }
      GaussDiagram d = GaussDiagram::of(l, base->origin);
      CHECK(nugatoryCrossing(d).has_value(), "kinks are nugatory");
      const GaussDiagram r = removeNugatoryCrossings(d);
      CHECK_EQ(r.crossings(), base->crossings(), "every kink goes, nothing else");
      CHECK(linkingMatrix(r) == linkingMatrix(*base), "kinks: linking numbers stay");
      CHECK(r.link().jones() == base->link().jones(), "kinks: the same link (Jones)");
    }
  }
}

// Why reduction matters: a hop's row is certified by drawing knotbuilder's
// link back (HopAssembler), and knotbuilder's drawer cannot draw a diagram
// with a nugatory crossing. The reduced diagram's row certifies.
void testReducedRowsCertify() {
  using regina::ExampleLink;
  const GaussDiagram d =
      sumThroughTwist(of(ExampleLink::trefoilLeft()), 0, of(ExampleLink::figureEight()), 1);
  auto certifies = [](const GaussDiagram &diagram) {
    ProofGraph g;
    NodeRegistry reg(g);
    HopRow row;
    row.node = reg.intern(diagram, "row").node;
    row.diagram = diagram;
    row.nodeMap = {0};
    row.pd = rowPD(diagram);
    row.layers = 2;
    try {
      HopAssembler hop(g, reg, row);
      return true;
    } catch (const std::exception &) {
      return false;
    }
  };
  CHECK(!certifies(d), "a row with a nugatory crossing cannot be certified");
  CHECK(certifies(removeNugatoryCrossings(d)), "its reduced diagram's row certifies");
}

// A far side a hop refused (11a_239's run, 2026-09-28, node 46): component
// 1 passes over both crossings it meets, so its PD code cannot carry its
// orientation and the row did not certify. Lifted off, it is a split unknot:
// the link is unchanged (Jones polynomial, linking numbers), the diagram
// comes apart, and every piece's row certifies.
void testLiftedRowsCertify() {
  GaussDiagram d;
  d.signs = {1, 1, 1, -1, -1, -1};
  d.comps = {{5, -6}, {4, 2}, {1, -4, -5, 6, -2, -3}, {-1, 3}};
  d.origin = {0, 1, 2, 3};
  auto certifies = [](const GaussDiagram &diagram) {
    ProofGraph g;
    NodeRegistry reg(g);
    HopRow row;
    row.node = reg.intern(diagram, "row").node;
    row.diagram = diagram;
    row.nodeMap.resize(diagram.components());
    std::iota(row.nodeMap.begin(), row.nodeMap.end(), 0);
    row.pd = rowPD(diagram);
    row.layers = 2;
    try {
      HopAssembler hop(g, reg, row);
      return true;
    } catch (const std::exception &) {
      return false;
    }
  };
  CHECK(d.link().pdAmbiguous(), "the refused row's PD code is ambiguous");
  CHECK(!certifies(d), "the refused row does not certify");

  const GaussDiagram lifted = liftSplitComponents(d);
  CHECK(lifted.components() == 4 && lifted.comps[1].empty(),
        "the over-everywhere component is lifted off, crossingless");
  CHECK(lifted.crossings() == 4, "its two crossings go, and only those");
  CHECK(lifted.link().jones() == d.link().jones(), "lifting keeps the link (Jones)");
  CHECK(linkingMatrix(lifted) == linkingMatrix(d), "lifting keeps every linking number");
  CHECK(liftSplitComponents(lifted).comps == lifted.comps, "nothing further to lift");

  const GaussDiagram s = simplifyKeepingComponents(d);
  const auto pieces = exactnaming::splitPieces(s);
  CHECK(pieces.size() >= 2, "simplified, the far side comes apart");
  bool all = true;
  for (const GaussDiagram &p : pieces)
    if (p.crossings() > 0 && !certifies(p)) all = false;
  CHECK(all, "every piece's row certifies");

  // A component over everywhere but linked cannot exist in a planar
  // diagram; one mixing over and under is never lifted.
  GaussDiagram hopf;
  hopf.signs = {1, 1};
  hopf.comps = {{1, -2}, {-1, 2}};
  hopf.origin = {0, 1};
  CHECK(liftSplitComponents(hopf).comps == hopf.comps, "the Hopf link is left alone");
}

} // namespace

int main() {
  testNugatoryCrossings();
  testReducedRowsCertify();
  testLiftedRowsCertify();
  testLinking();
  testSimplifyKeepsComponents();
  testDiagramHits();
  testDifferentDiagramsSameLink();
  testOrientationVariantsAreDifferentNodes();
  testUnknot();
  return cascadetest::finish("nodes_test");
}
