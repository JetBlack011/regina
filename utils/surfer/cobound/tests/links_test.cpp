// links_test.cpp
//
// Tests for bounds/links.h: graph link identity is exact, and the component map
// it returns is right (README.md, "Component maps").

#include <algorithm>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include <link/examplelink.h>
#include <link/link.h>

#include "linknaming/diagrams/diagramiso.h"
#include "cobound/bounds/searchcobordisms.h"
#include "cobound/bounds/links.h"
#include "linknaming/tests/check.h"
#include "linknaming/tests/gaussfixtures.h"
#include "linknaming/tables.h"

using linknaming::GaussDiagram;
using namespace bounds;
using namespace checks;

namespace {

// Real table rows (the same strings linknamer_test uses).
const char *L4a1_0 = "PD[X[6; 1; 7; 2]; X[8; 3; 5; 4]; X[2; 5; 3; 6]; X[4; 7; 1; 8]]";
const char *L4a1_1 = "PD[X[6; 2; 7; 1]; X[8; 4; 5; 3]; X[2; 8; 3; 7]; X[4; 6; 1; 5]]";
const char *L7n1_0 = "PD[X[6; 1; 7; 2]; X[12; 7; 13; 8]; X[4; 13; 1; 14]; X[5; 10; 6; 11]; "
                     "X[3; 8; 4; 9]; X[9; 14; 10; 5]; X[11; 2; 12; 3]]";
const char *L7n1_1 = "PD[X[12; 2; 13; 1]; X[6; 11; 7; 12]; X[4; 6; 1; 5]; X[13; 8; 14; 9]; "
                     "X[3; 11; 4; 10]; X[9; 14; 10; 5]; X[7; 3; 8; 2]]";

void testDiagramHits() {
  // The same diagram relabelled: a "diagram" hit whose map is the relabelling.
  using regina::ExampleLink;
  CobordismGraph g;
  LinkRegistry reg(g);
  regina::Link u = ExampleLink::hopf();
  u.insertLink(ExampleLink::torus(2, 4));
  // Two split pieces here; intern one connected piece: T(2,4) alone and a
  // three-component connected link, T(3,6).
  GaussDiagram t36 = of(ExampleLink::torus(3, 6));
  LinkMatch first = reg.intern(simplifyKeepingComponents(t36), "T(3,6)");
  CHECK(first.created, "first sight creates a node");
  // Reverse the component order: the map must follow.
  GaussDiagram perm = t36;
  std::reverse(perm.comps.begin(), perm.comps.end());
  std::reverse(perm.origin.begin(), perm.origin.end());
  LinkMatch again = reg.intern(perm, "T(3,6) permuted");
  CHECK(!again.created && again.link == first.link, "same diagram, same node");
  CHECK_EQ(again.method, std::string("diagram"), "found as a diagram");
  // perm's component i is t36's component 2-i; the link's components are
  // t36's. T(3,6) has symmetries permuting components, so check the map is
  // an isomorphism rather than a specific permutation:
  GaussDiagram mapped;
  mapped.signs = perm.signs;
  mapped.comps.assign(3, {});
  for (int i = 0; i < 3; ++i) mapped.comps[again.componentMap[i]] = perm.comps[i];
  CHECK(findDiagramIsomorphism(reg.info(first.link).diagram, mapped, false, true)
            .has_value(),
        "the returned map realises an isomorphism onto the node");
}

void testDifferentDiagramsSameLink() {
  // Scrambled diagrams of one hyperbolic link: the same link, by diagram or
  // by isometry, with a map consistent with the tracked origins up to a
  // symmetry of the link.
  using regina::ExampleLink;
  std::mt19937 rng(22);
  int isometryHits = 0, trials = 0;
  for (const regina::Link &l : {ExampleLink::whitehead(), ExampleLink::borromean(),
                                ExampleLink::conway()}) {
    CobordismGraph g;
    LinkRegistry reg(g);
    GaussDiagram base = simplifyKeepingComponents(of(l));
    LinkMatch n0 = reg.intern(base, "base");
    for (int t = 0; t < 25; ++t) {
      GaussDiagram s = scramble(base, rng, 20);
      for (const GaussDiagram &q : linknaming::splitPieces(s)) {
        if (q.components() != base.components()) continue; // came apart: skip
        LinkMatch m = reg.intern(q, "scrambled");
        ++trials;
        CHECK_EQ(m.link, n0.link, "scrambled diagram is the same node");
        if (m.method == "isometry") ++isometryHits;
        // Linking numbers must transport through the map (up to the global
        // sign a mirror gives).
        auto lq = linkingMatrix(q), ln = reg.info(n0.link).linking;
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

void testOrientationVariantsAreDifferentLinks() {
  // L7n1{0} (g4 2) and L7n1{1} (g4 0), L4a1{0} (g4 0) and L4a1{1} (g4 1):
  // one link with different component orientations. Merging them would be
  // unsound, so they must be different links.
  CobordismGraph g;
  LinkRegistry reg(g);
  for (auto [a, b] : {std::pair{L7n1_0, L7n1_1}, std::pair{L4a1_0, L4a1_1}}) {
    GaussDiagram da = simplifyKeepingComponents(of(linknaming::linkFromTablePD(a)));
    GaussDiagram db = simplifyKeepingComponents(of(linknaming::linkFromTablePD(b)));
    LinkMatch ma = reg.intern(da, a), mb = reg.intern(db, b);
    CHECK(ma.link != mb.link, "orientation variants are different nodes");
    // And each variant, re-interned, finds itself.
    CHECK_EQ(reg.intern(da, a).link, ma.link, "variant 0 finds itself");
    CHECK_EQ(reg.intern(db, b).link, mb.link, "variant 1 finds itself");
  }
}

void testUnknot() {
  CobordismGraph g;
  LinkRegistry reg(g);
  GaussDiagram u;
  u.comps = {{}};
  u.origin = {3};
  LinkMatch m = reg.intern(u, "crossingless");
  CHECK_EQ(m.link, reg.unknot(), "a crossingless piece is the unknot");
  CHECK_EQ(g.bestConnected(m.link)->genus, 0, "the unknot bounds a disc");
  LinkMatch m2 = reg.intern(u, "again");
  CHECK_EQ(m2.link, m.link, "one unknot node");
}

// Why reduction matters: a search's incoming link is certified by drawing
// T's link back (CobordismAssembler), and the drawer (todiagram.h)
// cannot draw a diagram with a nugatory crossing. The reduced diagram
// certifies.
void testReducedDiagramsCertify() {
  using regina::ExampleLink;
  const GaussDiagram d =
      sumThroughTwist(of(ExampleLink::trefoilLeft()), 0, of(ExampleLink::figureEight()), 1);
  auto certifies = [](const GaussDiagram &diagram) {
    CobordismGraph g;
    LinkRegistry reg(g);
    SearchedLink searched;
    searched.link = reg.intern(diagram, "row").link;
    searched.diagram = diagram;
    searched.linkMap = {0};
    searched.pd = diagramPD(diagram);
    searched.layers = 2;
    try {
      CobordismAssembler assembler(g, reg, searched);
      return true;
    } catch (const std::exception &) {
      return false;
    }
  };
  CHECK(!certifies(d), "a row with a nugatory crossing cannot be certified");
  CHECK(certifies(removeNugatoryCrossings(d)), "its reduced diagram's row certifies");
}

// An outgoing link a search refused (11a_239's run, 2026-09-28, graph link 46):
// component
// 1 passes over both crossings it meets, so its PD code cannot carry its
// orientation and it did not certify. Lifted off, it is a split unknot:
// the link is unchanged (Jones polynomial, linking numbers), the diagram
// comes apart, and every piece certifies.
void testLiftedDiagramsCertify() {
  GaussDiagram d;
  d.signs = {1, 1, 1, -1, -1, -1};
  d.comps = {{5, -6}, {4, 2}, {1, -4, -5, 6, -2, -3}, {-1, 3}};
  d.origin = {0, 1, 2, 3};
  auto certifies = [](const GaussDiagram &diagram) {
    CobordismGraph g;
    LinkRegistry reg(g);
    SearchedLink searched;
    searched.link = reg.intern(diagram, "row").link;
    searched.diagram = diagram;
    searched.linkMap.resize(diagram.components());
    std::iota(searched.linkMap.begin(), searched.linkMap.end(), 0);
    searched.pd = diagramPD(diagram);
    searched.layers = 2;
    try {
      CobordismAssembler assembler(g, reg, searched);
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
  const auto pieces = linknaming::splitPieces(s);
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
  testReducedDiagramsCertify();
  testLiftedDiagramsCertify();
  testDiagramHits();
  testDifferentDiagramsSameLink();
  testOrientationVariantsAreDifferentLinks();
  testUnknot();
  return checks::finish("links_test");
}
