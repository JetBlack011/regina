// profile_test.cpp
//
// Tests for cascade/partition.h and cascade/profile.h. Each block names the
// assumption it pins (README.md, "Assumptions and their tests").

#include <map>
#include <random>
#include <set>
#include <vector>

#include "cobound/bounds/partition.h"
#include "cobound/bounds/partitiongenera.h"
#include "cobound/tests/check.h"

using namespace cascade;

namespace {

// ---------------------------------------------------------------------------
// An independent model of glue(): explicit pieces with per-piece genera,
// connected components by breadth-first search, and the genus of each glued
// surface read off its Euler characteristic. glue() instead uses union-find
// and the cycle-rank formula; the two must agree on every input.
struct ModelResult {
  std::map<int, std::vector<int>> freeCurvesByPiece; // piece id -> free curves
  int genus = 0;
  bool integral = true;
  int closedPieces = 0; // glued pieces with no free curve
  int cycles = 0;       // total cycle rank, to be sure b1 > 0 is exercised
};

ModelResult modelGlue(const CobordismShape &c, Side gluedSide,
                      const Partition &gp, const std::vector<int> &cGenus,
                      const std::vector<int> &gGenus) {
  const auto &glued = gluedSide == Side::outgoing ? c.outComponent
                                                  : c.inComponent;
  const auto &freeSide = gluedSide == Side::outgoing ? c.inComponent
                                                     : c.outComponent;
  const int nc = c.components, ng = gp.blocks(), nv = nc + ng;
  std::vector<std::vector<int>> adj(nv);
  for (size_t j = 0; j < glued.size(); ++j) {
    int a = glued[j], b = nc + gp.blockOf(j);
    adj[a].push_back(b);
    adj[b].push_back(a);
  }
  // Boundary circles each vertex's piece has before gluing.
  std::vector<int> circles(nv, 0);
  for (int x : c.inComponent) ++circles[x];
  for (int x : c.outComponent) ++circles[x];
  for (int j = 0; j < gp.size(); ++j) ++circles[nc + gp.blockOf(j)];
  std::vector<int> comp(nv, -1);
  int ncomp = 0;
  for (int s = 0; s < nv; ++s) {
    if (comp[s] >= 0) continue;
    std::vector<int> queue{s};
    comp[s] = ncomp;
    for (size_t q = 0; q < queue.size(); ++q)
      for (int y : adj[queue[q]])
        if (comp[y] < 0) {
          comp[y] = ncomp;
          queue.push_back(y);
        }
    ++ncomp;
  }
  std::vector<int> chi(ncomp, 0), freeCount(ncomp, 0);
  for (int v = 0; v < nv; ++v) {
    int g = v < nc ? cGenus[v] : gGenus[v - nc];
    chi[comp[v]] += 2 - 2 * g - circles[v];
  }
  ModelResult r;
  for (size_t i = 0; i < freeSide.size(); ++i) {
    int k = comp[freeSide[i]];
    ++freeCount[k];
    r.freeCurvesByPiece[k].push_back(static_cast<int>(i));
  }
  for (int k = 0; k < ncomp; ++k) {
    // A closed piece (no free curve) has chi = 2 - 2g; an open one
    // chi = 2 - 2g - (#free curves). glue() counts closed pieces too.
    int twiceG = 2 - freeCount[k] - chi[k];
    if (twiceG % 2 != 0 || twiceG < 0)
      r.integral = false;
    r.genus += twiceG / 2;
    if (freeCount[k] == 0)
      ++r.closedPieces;
  }
  r.cycles = static_cast<int>(glued.size()) - nv + ncomp;
  return r;
}

Partition partitionFromGroups(const std::map<int, std::vector<int>> &g, int n) {
  std::vector<int> labels(n, -1);
  for (const auto &[k, members] : g)
    for (int m : members)
      labels[m] = k;
  return Partition::fromLabels(labels);
}

// Splits `total` into `parts` nonnegative integers at random.
std::vector<int> randomSplit(int total, int parts, std::mt19937 &rng) {
  std::vector<int> v(parts, 0);
  for (int i = 0; i < total; ++i)
    ++v[std::uniform_int_distribution<int>(0, parts - 1)(rng)];
  return v;
}

CobordismShape randomShape(std::mt19937 &rng) {
  auto u = [&](int lo, int hi) {
    return std::uniform_int_distribution<int>(lo, hi)(rng);
  };
  CobordismShape c;
  c.components = u(1, 4);
  int nin = u(0, 5), nout = u(0, 5);
  c.inComponent.resize(nin);
  c.outComponent.resize(nout);
  for (int &x : c.inComponent) x = u(0, c.components - 1);
  for (int &x : c.outComponent) x = u(0, c.components - 1);
  // Make sure every component carries a curve.
  std::vector<int> count(c.components, 0);
  for (int x : c.inComponent) ++count[x];
  for (int x : c.outComponent) ++count[x];
  for (int k = 0; k < c.components; ++k)
    if (!count[k]) {
      if (u(0, 1)) c.inComponent.push_back(k);
      else c.outComponent.push_back(k);
    }
  c.genus = u(0, 3);
  return c;
}

Partition randomPartition(int n, std::mt19937 &rng) {
  std::vector<int> labels(n);
  for (int &l : labels) l = std::uniform_int_distribution<int>(0, n)(rng);
  return Partition::fromLabels(labels);
}

// ---------------------------------------------------------------------------

void testPartition() {
  // A1: normal form makes equal partitions compare equal.
  CHECK(Partition::fromLabels({5, 5, 2}) == Partition::fromLabels({0, 0, 7}),
        "normal form: relabelled partitions are equal");
  CHECK_EQ(Partition::fromLabels({3, 1, 3, 2}).str(), std::string("{0,2}{1}{3}"),
           "str()");
  // A2: refines is reflexive, and singletons refine everything, which
  // refines the one-block partition.
  for (int n = 0; n <= 5; ++n)
    for (const Partition &p : allPartitions(n)) {
      CHECK(p.refines(p), "refines is reflexive");
      CHECK(Partition::singletons(n).refines(p), "singletons refine all");
      CHECK(p.refines(Partition::coarsest(n)), "all refine the one block");
    }
  // Bell numbers.
  const int bell[] = {1, 1, 2, 5, 15, 52, 203};
  for (int n = 0; n <= 6; ++n)
    CHECK_EQ(static_cast<int>(allPartitions(n).size()), bell[n],
             "allPartitions gives Bell(n)");
  // refines is antisymmetric and transitive on all of Bell(4).
  auto all4 = allPartitions(4);
  for (const auto &a : all4)
    for (const auto &b : all4) {
      if (a.refines(b) && b.refines(a))
        CHECK(a == b, "refines is antisymmetric");
      for (const auto &c : all4)
        if (a.refines(b) && b.refines(c))
          CHECK(a.refines(c), "refines is transitive");
    }
}

void testGlueAgainstModel() {
  // A3: glue()'s genus (sum of genera + cycle rank of the gluing graph) and
  // partition agree with the Euler-characteristic model on random inputs,
  // for any split of the genera among pieces, on both sides.
  std::mt19937 rng(20260928);
  int trials = 0, closedCases = 0, cycleCases = 0;
  for (int t = 0; t < 20000; ++t) {
    CobordismShape c = randomShape(rng);
    for (Side side : {Side::outgoing, Side::incoming}) {
      const auto &glued = side == Side::outgoing ? c.outComponent
                                                 : c.inComponent;
      const auto &freeSide = side == Side::outgoing ? c.inComponent
                                                    : c.outComponent;
      Partition gp = randomPartition(static_cast<int>(glued.size()), rng);
      int gg = std::uniform_int_distribution<int>(0, 3)(rng);
      auto result = glue(c, side, gp, gg);
      if (freeSide.empty()) {
        CHECK(!result.has_value(), "no free curves: nullopt");
        continue;
      }
      CHECK(result.has_value(), "free curves: a result");
      auto cg = randomSplit(c.genus, c.components, rng);
      auto ggs = randomSplit(gg, std::max(1, gp.blocks()), rng);
      ggs.resize(gp.blocks());
      if (gp.blocks() == 0 && gg > 0)
        continue; // no pieces to carry genus: not a meaningful input
      ModelResult m = modelGlue(c, side, gp, cg, ggs);
      CHECK(m.integral, "model: every glued piece has integral genus");
      CHECK_EQ(result->genus, m.genus, "glue genus == Euler-characteristic model");
      CHECK(result->partition ==
                partitionFromGroups(m.freeCurvesByPiece,
                                    static_cast<int>(freeSide.size())),
            "glue partition == model's pieces");
      ++trials;
      if (m.closedPieces > 0)
        ++closedCases;
      if (m.cycles > 0)
        ++cycleCases;
    }
  }
  CHECK(trials > 20000, "enough random glue trials");
  CHECK(closedCases > 500, "the trials include pieces that close up");
  CHECK(cycleCases > 500, "the trials include gluing graphs with cycles");
}

void testPaperSpecialCases() {
  // A4: the paper's lem:cobordism-inequality is the special case "C and the
  // capping surface both connected": g4(L0) <= g4(L1) + g + n1 - 1.
  for (int n1 = 1; n1 <= 4; ++n1)
    for (int g = 0; g <= 2; ++g)
      for (int h = 0; h <= 2; ++h) {
        CobordismShape c;
        c.components = 1;
        c.genus = g;
        c.inComponent = {0};
        c.outComponent.assign(n1, 0);
        auto r = glue(c, Side::outgoing, Partition::coarsest(n1), h);
        CHECK_EQ(r->genus, h + g + n1 - 1, "cobordism inequality");
        CHECK_EQ(r->partition.blocks(), 1, "a knot's surface is connected");
      }
  // A5: the unlink rule: an n-component unlink bounds disjoint discs, so a
  // connected genus-g cobordism to it gives genus g, not g + n - 1.
  for (int n1 = 1; n1 <= 4; ++n1) {
    CobordismShape c;
    c.components = 1;
    c.genus = 1;
    c.inComponent = {0};
    c.outComponent.assign(n1, 0);
    auto r = glue(c, Side::outgoing, Partition::singletons(n1), 0);
    CHECK_EQ(r->genus, 1, "unlink rule: no n-1 penalty with disjoint discs");
  }
  // A6: a band move K -> L (pair of pants): disjoint discs for L make K
  // slice; an annulus for L (g4(L) = 0 in LinkInfo's sense) gives only 1.
  CobordismShape pants;
  pants.components = 1;
  pants.genus = 0;
  pants.inComponent = {0};
  pants.outComponent = {0, 0};
  CHECK_EQ(glue(pants, Side::outgoing, Partition::singletons(2), 0)->genus, 0,
           "band to disjoint discs: slice");
  CHECK_EQ(glue(pants, Side::outgoing, Partition::coarsest(2), 0)->genus, 1,
           "band to an annulus: genus 1 only");
  // A7: read backwards, the same pants bounds L by a connected surface of
  // K's genus: g4(L) <= g4(K) (the reverse inequality, n0 = 1).
  auto back = glue(pants, Side::incoming, Partition::coarsest(1), 0);
  CHECK_EQ(back->genus, 0, "reverse band: annulus for L from a disc for K");
  CHECK(back->partition == Partition::coarsest(2), "reverse band: connected");
  // A8: the product cobordism transports every (partition, genus) exactly.
  for (int n = 1; n <= 4; ++n)
    for (const Partition &p : allPartitions(n))
      for (Side s : {Side::outgoing, Side::incoming}) {
        auto r = glue(CobordismShape::product(n), s, p, 2);
        CHECK(r->partition == p && r->genus == 2, "product is the identity");
      }
  // A9: disconnected witnesses (the ~70% tubed ones): two annuli from a
  // 2-component link to a 2-component link keep blocks apart.
  CobordismShape annuli = CobordismShape::product(2);
  auto sep = glue(annuli, Side::outgoing, Partition::singletons(2), 0);
  CHECK(sep->partition == Partition::singletons(2) && sep->genus == 0,
        "two annuli preserve disjoint discs");
  // A10: a piece that closes up is discarded but its genus still counted
  // (so the bound stays an upper bound): C = {disc on in-curve 0} u
  // {annulus in 1 -> out 0}; read backwards with G = two discs, the disc
  // piece closes into a sphere.
  CobordismShape cap;
  cap.components = 2;
  cap.genus = 0;
  cap.inComponent = {0, 1};
  cap.outComponent = {1};
  auto capBack = glue(cap, Side::incoming, Partition::singletons(2), 0);
  CHECK_EQ(capBack->genus, 0, "closed sphere piece adds no genus");
  CHECK_EQ(capBack->partition.size(), 1, "only the open piece has free curves");
  auto capBackG = glue(cap, Side::incoming, Partition::singletons(2), 1);
  CHECK_EQ(capBackG->genus, 1, "closed piece's genus is still counted");
}

void testLinking() {
  // A11: the linking condition is on block TOTALS, not on pairs.
  std::vector<std::vector<int>> hopf = {{0, -1}, {-1, 0}};
  CHECK(!linkingAllows(Partition::singletons(2), hopf),
        "Hopf link bounds no disjoint discs");
  CHECK(linkingAllows(Partition::coarsest(2), hopf), "Hopf: an annulus is fine");
  std::vector<std::vector<int>> cancel = {{0, 1, -1}, {1, 0, 0}, {-1, 0, 0}};
  CHECK(linkingAllows(Partition::fromLabels({0, 1, 1}), cancel),
        "blocks {0},{1,2} with lk(0,1)+lk(0,2) = 0 are allowed");
  CHECK(!linkingAllows(Partition::singletons(3), cancel),
        "but singletons are not");
  std::vector<std::vector<int>> unlink(3, std::vector<int>(3, 0));
  for (const Partition &p : allPartitions(3))
    CHECK(linkingAllows(p, unlink), "unlink: every partition allowed");
}

void testProfile() {
  // A12: the Pareto set keeps exactly the non-implied entries.
  Profile pr(3);
  CHECK(pr.insert(Partition::coarsest(3), 2, 0), "first entry");
  CHECK(!pr.insert(Partition::coarsest(3), 3, 1), "worse genus is implied");
  CHECK(pr.insert(Partition::fromLabels({0, 0, 1}), 2, 2),
        "finer partition, same genus: new");
  CHECK_EQ(static_cast<int>(pr.entries().size()), 1,
           "and it removes the coarser entry it implies");
  CHECK(pr.insert(Partition::coarsest(3), 1, 3), "coarser but better: new");
  CHECK_EQ(static_cast<int>(pr.entries().size()), 2, "two incomparable entries");
  CHECK_EQ(pr.best(Partition::coarsest(3))->genus, 1, "connected best");
  CHECK_EQ(pr.best(Partition::fromLabels({0, 0, 1}))->genus, 2, "block best");
  CHECK(!pr.best(Partition::singletons(3)).has_value(), "no disc entry");
  // Random: after any insertion sequence, no entry implies another, and
  // every inserted pair is implied by the set.
  std::mt19937 rng(7);
  for (int t = 0; t < 300; ++t) {
    Profile q(4);
    std::vector<std::pair<Partition, int>> inserted;
    for (int k = 0; k < 12; ++k) {
      Partition p = randomPartition(4, rng);
      int g = std::uniform_int_distribution<int>(0, 3)(rng);
      q.insert(p, g, k);
      inserted.push_back({p, g});
    }
    for (const auto &[p, g] : inserted)
      CHECK(q.implies(p, g), "every inserted fact stays implied");
    for (const auto &a : q.entries())
      for (const auto &b : q.entries())
        if (&a != &b)
          CHECK(!(a.genus <= b.genus && a.partition.refines(b.partition)),
                "no entry implies another");
  }
}

void testGlueMonotone() {
  // A13: glue() is monotone: a finer glued partition and a smaller glued
  // genus never give a worse result. This is what lets propagate() skip
  // records that were later superseded.
  std::mt19937 rng(99);
  for (int t = 0; t < 5000; ++t) {
    CobordismShape c = randomShape(rng);
    if (c.inComponent.empty())
      continue;
    int n = static_cast<int>(c.outComponent.size());
    Partition coarse = randomPartition(n, rng);
    // A random refinement of `coarse`: split each block at random.
    std::vector<int> labels(n);
    for (int j = 0; j < n; ++j)
      labels[j] = coarse.blockOf(j) * 10 +
                  std::uniform_int_distribution<int>(0, 1)(rng);
    Partition fine = Partition::fromLabels(labels);
    auto a = glue(c, Side::outgoing, coarse, 2);
    auto b = glue(c, Side::outgoing, fine, 1);
    CHECK(b->genus <= a->genus && b->partition.refines(a->partition),
          "glue is monotone in the glued entry");
    // And never below its input genus (so cycles cannot self-improve).
    CHECK(a->genus >= 2, "glue never lowers genus below its input");
  }
}

} // namespace

int main() {
  testPartition();
  testGlueAgainstModel();
  testPaperSpecialCases();
  testLinking();
  testProfile();
  testGlueMonotone();
  return cascadetest::finish("profile_test");
}
