// partitiongenera.cpp

#include "cobound/bounds/partitiongenera.h"

#include <numeric>
#include <stdexcept>

namespace bounds {

namespace {

// Plain union-find; the graphs here have a handful of vertices.
struct UnionFind {
  std::vector<int> parent;
  explicit UnionFind(int n) : parent(n) {
    std::iota(parent.begin(), parent.end(), 0);
  }
  int find(int x) {
    while (parent[x] != x)
      x = parent[x] = parent[parent[x]];
    return x;
  }
  void unite(int a, int b) { parent[find(a)] = find(b); }
};

} // namespace

void CobordismShape::validate() const {
  if (components < 0 || genus < 0)
    throw std::invalid_argument("CobordismShape: negative count");
  std::vector<int> curves(components, 0);
  for (const auto *side : {&inComponent, &outComponent})
    for (int c : *side) {
      if (c < 0 || c >= components)
        throw std::invalid_argument("CobordismShape: component out of range");
      ++curves[c];
    }
  for (int k : curves)
    if (k == 0)
      throw std::invalid_argument(
          "CobordismShape: a component carries no boundary curve");
}

CobordismShape CobordismShape::reversed() const {
  CobordismShape r = *this;
  std::swap(r.inComponent, r.outComponent);
  return r;
}

CobordismShape CobordismShape::product(int n) {
  CobordismShape c;
  c.components = n;
  c.genus = 0;
  c.inComponent.resize(n);
  c.outComponent.resize(n);
  for (int i = 0; i < n; ++i)
    c.inComponent[i] = c.outComponent[i] = i;
  return c;
}

std::optional<Glued> glue(const CobordismShape &c, Side glued,
                          const Partition &gluedPartition, int gluedGenus) {
  c.validate();
  const std::vector<int> &gluedCurves =
      glued == Side::outgoing ? c.outComponent : c.inComponent;
  const std::vector<int> &freeCurves =
      glued == Side::outgoing ? c.inComponent : c.outComponent;
  if (gluedPartition.size() != static_cast<int>(gluedCurves.size()))
    throw std::invalid_argument("glue: partition size != glued curve count");
  if (gluedGenus < 0)
    throw std::invalid_argument("glue: negative genus");
  if (freeCurves.empty())
    return std::nullopt;

  // Vertices: C's components 0..c.components-1, then G's pieces.
  const int nv = c.components + gluedPartition.blocks();
  UnionFind uf(nv);
  for (size_t j = 0; j < gluedCurves.size(); ++j)
    uf.unite(gluedCurves[j], c.components + gluedPartition.blockOf(j));

  // Sum over graph components of b1 = E - V + 1, i.e. E - V + (#components).
  // Every vertex has an edge or a free curve (validate() for C's; every
  // block is nonempty for G's), so isolated vertices are genuine pieces.
  std::vector<char> isRoot(nv, 0);
  for (int v = 0; v < nv; ++v)
    isRoot[uf.find(v)] = 1;
  const int graphComponents = std::accumulate(isRoot.begin(), isRoot.end(), 0);
  const int b1 = static_cast<int>(gluedCurves.size()) - nv + graphComponents;

  std::vector<int> labels(freeCurves.size());
  for (size_t i = 0; i < freeCurves.size(); ++i)
    labels[i] = uf.find(freeCurves[i]);

  return Glued{Partition::fromLabels(labels), c.genus + gluedGenus + b1};
}

bool linkingAllows(const Partition &p,
                   const std::vector<std::vector<int>> &lk) {
  const int n = p.size();
  if (static_cast<int>(lk.size()) != n)
    throw std::invalid_argument("linkingAllows: matrix size != partition size");
  const int b = p.blocks();
  // total[a][c]: the linking number between blocks a and c.
  std::vector<std::vector<long>> total(b, std::vector<long>(b, 0));
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      if (i != j)
        total[p.blockOf(i)][p.blockOf(j)] += lk[i][j];
  for (int a = 0; a < b; ++a)
    for (int c = 0; c < b; ++c)
      if (a != c && total[a][c] != 0)
        return false;
  return true;
}

bool PartitionGenera::implies(const Partition &p, int g) const {
  for (const PartitionGenus &e : entries_)
    if (e.genus <= g && e.partition.refines(p))
      return true;
  return false;
}

bool PartitionGenera::insert(const Partition &p, int g, long derivation) {
  if (p.size() != n_)
    throw std::invalid_argument("PartitionGenera::insert: partition size != link");
  if (implies(p, g))
    return false;
  std::erase_if(entries_, [&](const PartitionGenus &e) {
    return g <= e.genus && p.refines(e.partition);
  });
  entries_.push_back({p, g, derivation});
  return true;
}

std::optional<PartitionGenus> PartitionGenera::best(const Partition &target) const {
  std::optional<PartitionGenus> ans;
  for (const PartitionGenus &e : entries_)
    if (e.partition.refines(target) && (!ans || e.genus < ans->genus))
      ans = e;
  return ans;
}

} // namespace bounds
