// simplification_test.cpp
//
// Tests for linknaming/diagrams/simplification.h: linking numbers,
// simplify() keeping components, and removing nugatory crossings (moved
// from cobound's links_test, refactor phase 3).

#include <algorithm>
#include <random>
#include <string>
#include <vector>

#include <link/examplelink.h>
#include <link/link.h>

#include "linknaming/diagrams/simplification.h"
#include "linknaming/tests/check.h"
#include "linknaming/tests/gaussfixtures.h"

using exactnaming::GaussDiagram;
using namespace cascade;
using namespace cascadetest;

namespace {

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

} // namespace

int main() {
  testNugatoryCrossings();
  testLinking();
  testSimplifyKeepsComponents();
  return cascadetest::finish("simplification_test");
}
