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

} // namespace

int main() {
  testLinking();
  testSimplifyKeepsComponents();
  testDiagramHits();
  testDifferentDiagramsSameLink();
  testOrientationVariantsAreDifferentNodes();
  testUnknot();
  return cascadetest::finish("nodes_test");
}
