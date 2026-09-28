// diagramiso.cpp

#include "diagramiso.h"

#include <algorithm>
#include <cstdlib>
#include <functional>

using exactnaming::GaussDiagram;

namespace cascade {

GaussDiagram reverseAll(const GaussDiagram &d) {
  // Reversing every strand keeps each crossing's sign (both strands turn
  // round) and its over/under.
  GaussDiagram r = d;
  for (auto &w : r.comps)
    std::reverse(w.begin(), w.end());
  return r;
}

GaussDiagram mirrorImage(const GaussDiagram &d) {
  GaussDiagram r = d;
  for (int &s : r.signs)
    s = -s;
  for (auto &w : r.comps)
    for (long &x : w)
      x = -x;
  return r;
}

GaussDiagram reverseComponent(const GaussDiagram &d, size_t c) {
  GaussDiagram r = d;
  std::reverse(r.comps[c].begin(), r.comps[c].end());
  // A crossing between component c and another changes sign; a self-crossing
  // of c does not (both of its strands turn round).
  std::vector<int> onC(d.crossings(), 0);
  for (long x : d.comps[c])
    ++onC[std::labs(x) - 1];
  for (size_t k = 0; k < d.crossings(); ++k)
    if (onC[k] == 1)
      r.signs[k] = -r.signs[k];
  return r;
}

namespace {

// Orientation-preserving, non-mirroring isomorphism a -> b, or nullopt; one
// realising `required` (a's component i -> b's required[i]) when given.
std::optional<std::vector<int>> strictIsomorphism(const GaussDiagram &a,
                                                  const GaussDiagram &b,
                                                  const std::vector<int> *required) {
  const size_t m = a.components(), n = a.crossings();
  if (m != b.components() || n != b.crossings())
    return std::nullopt;
  {
    std::vector<size_t> la, lb;
    for (const auto &w : a.comps) la.push_back(w.size());
    for (const auto &w : b.comps) lb.push_back(w.size());
    std::sort(la.begin(), la.end());
    std::sort(lb.begin(), lb.end());
    if (la != lb)
      return std::nullopt;
    long pa = 0, pb = 0;
    for (int s : a.signs) pa += s;
    for (int s : b.signs) pb += s;
    if (pa != pb)
      return std::nullopt; // writhe
  }
  // Order a's components longest first: long words pin crossings early.
  std::vector<size_t> order(m);
  for (size_t i = 0; i < m; ++i) order[i] = i;
  std::stable_sort(order.begin(), order.end(), [&](size_t x, size_t y) {
    return a.comps[x].size() > a.comps[y].size();
  });

  std::vector<long> crossMap(n, -1), crossUsed(n, 0); // a crossing -> b crossing
  std::vector<int> compMap(m, -1);
  std::vector<char> compUsed(m, 0);

  std::function<bool(size_t)> rec = [&](size_t idx) -> bool {
    if (idx == m)
      return true;
    const size_t ca = order[idx];
    const auto &wa = a.comps[ca];
    for (size_t cb = 0; cb < m; ++cb) {
      if (compUsed[cb] || b.comps[cb].size() != wa.size())
        continue;
      if (required && (*required)[ca] != static_cast<int>(cb))
        continue;
      const auto &wb = b.comps[cb];
      const size_t len = wa.size();
      for (size_t rot = 0; rot < std::max<size_t>(len, 1); ++rot) {
        // Try mapping wa[i] -> wb[(i + rot) % len].
        std::vector<size_t> newly;
        bool ok = true;
        for (size_t i = 0; i < len && ok; ++i) {
          const long xa = wa[i], xb = wb[(i + rot) % len];
          if ((xa > 0) != (xb > 0)) { ok = false; break; } // over vs under
          const size_t ka = std::labs(xa) - 1, kb = std::labs(xb) - 1;
          if (a.signs[ka] != b.signs[kb]) { ok = false; break; }
          if (crossMap[ka] < 0) {
            if (crossUsed[kb]) { ok = false; break; }
            crossMap[ka] = static_cast<long>(kb);
            crossUsed[kb] = 1;
            newly.push_back(ka);
          } else if (crossMap[ka] != static_cast<long>(kb)) {
            ok = false;
          }
        }
        if (ok) {
          compMap[ca] = static_cast<int>(cb);
          compUsed[cb] = 1;
          if (rec(idx + 1))
            return true;
          compUsed[cb] = 0;
          compMap[ca] = -1;
        }
        for (size_t ka : newly) {
          crossUsed[crossMap[ka]] = 0;
          crossMap[ka] = -1;
        }
        if (len == 0)
          break; // a crossingless component: one "rotation"
      }
    }
    return false;
  };
  if (!rec(0))
    return std::nullopt;
  return compMap;
}

} // namespace

std::optional<DiagramIsomorphism>
findDiagramIsomorphism(const GaussDiagram &a, const GaussDiagram &b,
                       bool allowMirror, bool allowReverse,
                       const std::vector<int> *componentMap) {
  if (componentMap && componentMap->size() != a.components())
    return std::nullopt;
  for (int mir = 0; mir <= (allowMirror ? 1 : 0); ++mir)
    for (int rev = 0; rev <= (allowReverse ? 1 : 0); ++rev) {
      GaussDiagram t = a;
      if (mir) t = mirrorImage(t);
      if (rev) t = reverseAll(t);
      if (auto map = strictIsomorphism(t, b, componentMap))
        return DiagramIsomorphism{*map, mir == 1, rev == 1};
    }
  return std::nullopt;
}

} // namespace cascade
