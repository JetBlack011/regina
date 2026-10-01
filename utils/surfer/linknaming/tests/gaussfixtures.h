// gaussfixtures.h: Gauss-diagram helpers the diagram tests share -- linknaming's
// simplification_test and cobound's links_test (moved from links_test,
// refactor phase 3).

#pragma once

#include <numeric>
#include <random>
#include <vector>

#include <link/link.h>

#include "linknaming/diagrams/gaussdiagram.h"
#include "linknaming/diagrams/simplification.h"

namespace cascadetest {

using exactnaming::GaussDiagram;

inline GaussDiagram of(const regina::Link &l) {
  std::vector<size_t> origin(l.countComponents());
  std::iota(origin.begin(), origin.end(), 0);
  return GaussDiagram::of(l, origin);
}

// Random classical Reidemeister II/III moves, then simplify(): usually a
// different diagram of the same link, with `origin` carried along.
inline GaussDiagram scramble(const GaussDiagram &d, std::mt19937 &rng, int moves) {
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
  return cascade::simplifyKeepingComponents(s);
}

// `a` summed along its component `onComponent` with the knot `k` through a
// twist: that component's word becomes `c(over) a... c(under) k...`, where
// crossing 0 is the twist. Nugatory by construction, with all of `a` on the
// side it cuts off (a's other components included, being linked only with
// that side).
inline GaussDiagram sumThroughTwist(const GaussDiagram &a, size_t onComponent,
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

inline bool sameDiagram(const GaussDiagram &a, const GaussDiagram &b) {
  return a.signs == b.signs && a.comps == b.comps && a.origin == b.origin;
}

} // namespace cascadetest
