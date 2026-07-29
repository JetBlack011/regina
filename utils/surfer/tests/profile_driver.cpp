//
//  profile_driver.cpp
//
//  Standalone driver for profiling EmbeddingSearch<dim,2>/SurfaceSearch
//  (embeddingsearch.h) against a range of target triangulations. Not part
//  of CTest -- build and run manually, typically under `perf record`:
//
//    ./profile_driver bary3 <subdivisions> <cond> <threads>
//    ./profile_driver bary4 <subdivisions> <cond> <threads>
//    ./profile_driver isosig4 <isosig> <cond> <threads>
//    ./profile_driver pipeline <pdcode> <thickenLayers> <cond> <threads>
//
//  cond in {all, closed, proper, connected}. search() no longer prints
//  anything on its own -- passing no callbacks (as every mode but
//  pipeline-loud does) means steady-state work runs with no output at all,
//  so there's nothing left to suppress for perf's benefit.

#include <algorithm>
#include <atomic>
#include <chrono>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <triangulation/example4.h>

#include "cobordismbuilder.h"
#include "embeddingsearch.h"
#include "surfacesearch.h"
#include "knotbuilder.h"

namespace {

BoundaryCondition parseCond(const std::string &s) {
  if (s == "all")
    return BoundaryCondition::all;
  if (s == "closed")
    return BoundaryCondition::closed;
  if (s == "proper")
    return BoundaryCondition::proper;
  if (s == "connected")
    return BoundaryCondition::connected;
  throw std::runtime_error("Unknown condition: " + s);
}

// Times search(), passing no callbacks -- so it runs with no output at all,
// keeping perf's samples on the actual search rather than any reporting.
template <typename Search>
void runTimed(Search &search, unsigned threads, BoundaryCondition cond) {
  const auto start = std::chrono::steady_clock::now();
  search.search(threads, cond);
  const auto elapsed = std::chrono::steady_clock::now() - start;
  std::cerr << "[driver] wall time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed)
                   .count()
            << " ms\n";
}

// DIAGNOSTIC: exposes EmbeddingSearch<4,2>'s protected skeleton_/graph_ and
// EmbeddednessPredicate so we can drive the DFS single-threaded ourselves and
// count, at fine grain, how many visited nodes actually satisfy a given
// BoundaryCondition (i.e. how often surfaceTypeKey()/boundaryComponentsMap-
// Injectively() actually fire) -- independent of embeddingsearch.cpp's own
// FLUSH_EVERY-batched reporting, which is too coarse-grained when the hit
// rate is low.
class DiagSearch : public EmbeddingSearch<4, 2> {
public:
  using EmbeddingSearch<4, 2>::EmbeddingSearch;

  void diagRun(BoundaryCondition cond, long long reportEvery) {
    ConnectedInducedSubgraphEnumerator enumerator(graph_.adjList.first,
                                                   graph_.adjList.second);
    EmbeddedSubmanifold<4, 2> embedding(skeleton_);
    EmbeddednessPredicate predicate(embedding, graph_.graphToSkel);

    long long visited = 0, satisfying = 0;
    const auto start = std::chrono::steady_clock::now();

    for (int s = 1; s <= graph_.adjList.first; ++s) {
      enumerator.enumerateFromRootFiltered(
          s,
          [&](const std::vector<int> &) {
            ++visited;
            if (embedding.satisfies(cond))
              ++satisfying;
            if (visited % reportEvery == 0) {
              double pct = 100.0 * static_cast<double>(satisfying) /
                           static_cast<double>(visited);
              std::cerr << "[diag] elapsed="
                        << formatElapsed(std::chrono::steady_clock::now() -
                                         start)
                        << " visited=" << visited
                        << " satisfying=" << satisfying << " (" << pct
                        << "%)\n";
            }
          },
          predicate);
    }
  }
};

// Returns this process's current peak resident set size, in bytes (Linux
// only, via /proc/self/status' VmHWM -- more direct than getrusage(), which
// on Linux reports the same VmHWM value but in KB via a less obviously-named
// field).
long long peakRssBytes() {
  std::ifstream status("/proc/self/status");
  std::string line;
  while (std::getline(status, line)) {
    if (!line.starts_with("VmHWM:"))
      continue;
    std::istringstream iss(line.substr(6));
    long long kb = 0;
    iss >> kb;
    return kb * 1024;
  }
  return -1;
}

// DIAGNOSTIC: single-threaded analogue of DiagSearch, but driving
// KnottedSurface (not the base EmbeddedSubmanifold<4,2>) so the
// local-flatness/transverse-self-intersection checks -- and PetalCache, the
// memoization added to keep them off the hot path -- actually run. Reports,
// at each interval: how many addFace() calls each check has rejected (i.e.
// how much of the search tree these checks are really pruning) and
// PetalCache's cumulative hit/miss counts (i.e. how much of that pruning
// cost memoization is actually avoiding recomputing).
class SurfaceDiagSearch : public EmbeddingSearch<4, 2> {
public:
  using EmbeddingSearch<4, 2>::EmbeddingSearch;

  void diagRun(BoundaryCondition cond, long long reportEvery) {
    ConnectedInducedSubgraphEnumerator enumerator(graph_.adjList.first,
                                                   graph_.adjList.second);
    KnottedSurface embedding(skeleton_);
    EmbeddednessPredicate<KnottedSurface> predicate(embedding,
                                                     graph_.graphToSkel);

    long long visited = 0, satisfying = 0;
    const auto start = std::chrono::steady_clock::now();

    auto report = [&] {
      const auto &stats = embedding.petalCacheStats();
      std::cerr << "[surfacediag] elapsed="
                << formatElapsed(std::chrono::steady_clock::now() - start)
                << " visited=" << visited << " satisfying=" << satisfying
                << " localFlatnessRejections=" << stats.localFlatnessRejections
                << " transverseRejections=" << stats.transverseRejections
                << " unknotChecks=" << stats.unknotChecks
                << " unknotCacheHits=" << stats.unknotCacheHits
                << " linkingChecks=" << stats.linkingChecks
                << " linkingCacheHits=" << stats.linkingCacheHits
                << " peakRssMB=" << (peakRssBytes() / (1024 * 1024)) << "\n";
    };

    for (int s = 1; s <= graph_.adjList.first; ++s) {
      enumerator.enumerateFromRootFiltered(
          s,
          [&](const std::vector<int> &) {
            ++visited;
            if (embedding.satisfies(cond))
              ++satisfying;
            if (visited % reportEvery == 0)
              report();
          },
          predicate);
    }
    report();
  }
};

// PHASE 0 DIAGNOSTIC (see the VertexLinkData plan): estimates, per ambient
// vertex, how many candidate "closed petals" would need to be precomputed --
// before writing any of the real precomputation machinery.
//
// Reduction (see the plan for the derivation): ambient spokes at v <-> link
// vertices of Lk(v); ambient triangle-corners at v <-> link edges of Lk(v);
// a closed petal <-> a simple cycle in Lk(v)'s 1-skeleton. So this reuses
// ConnectedInducedSubgraphEnumerator directly over Lk(v)'s
// edges-share-a-vertex graph, with a small predicate mirroring the
// monotonic-count idiom EmbeddedSubmanifold already uses for
// badProperCount_ (each link vertex may be used by at most 2 currently-added
// link edges). Counts closed results only -- no curve/isUnknot/linking work
// yet, that's Phase 2.
namespace vertexstats {

// Hereditary (removing an edge only lowers counts, never raises them) --
// matches ConnectedInducedSubgraphEnumerator::enumerateFiltered()'s
// precondition. "Closed" itself is NOT hereditary and is checked separately
// via isClosed(), not baked into tryAdd()/undo().
class SpokeCountPredicate : public ConditionalPredicate {
  const std::vector<std::pair<int, int>> &edgeEndpoints_;
      /**< 0-indexed link-vertex ids per link-edge, indexed by (graph node id - 1). */
  std::vector<int> count_; /**< Per link-vertex: how many added link-edges currently use it (0/1/2). */

public:
  bool sawLoopEdge = false;
      /**< Whether any link-edge with both endpoints equal was encountered
           (would mean a triangle whose two spokes are the same ambient
           edge) -- excluded from closure for now; Phase 1 decides whether
           this can genuinely occur here and how to handle it if so. */

  SpokeCountPredicate(const std::vector<std::pair<int, int>> &endpoints,
                       size_t numLinkVertices)
      : edgeEndpoints_(endpoints), count_(numLinkVertices, 0) {}

  bool tryAdd(int v) override {
    auto [a, b] = edgeEndpoints_[static_cast<size_t>(v - 1)];
    if (a == b) {
      sawLoopEdge = true;
      return false;
    }
    if (count_[static_cast<size_t>(a)] >= 2 || count_[static_cast<size_t>(b)] >= 2)
      return false;
    ++count_[static_cast<size_t>(a)];
    ++count_[static_cast<size_t>(b)];
    return true;
  }

  void undo(int v) override {
    auto [a, b] = edgeEndpoints_[static_cast<size_t>(v - 1)];
    --count_[static_cast<size_t>(a)];
    --count_[static_cast<size_t>(b)];
  }

  // Valid only immediately after U was reported (i.e. count_ reflects
  // exactly U's members) -- true iff every link-vertex U touches is at
  // count 2 (not left dangling at 1).
  bool isClosed(const std::vector<int> &U) const {
    for (int edgeId : U) {
      auto [a, b] = edgeEndpoints_[static_cast<size_t>(edgeId - 1)];
      if (count_[static_cast<size_t>(a)] != 2 || count_[static_cast<size_t>(b)] != 2)
        return false;
    }
    return true;
  }
};

struct VertexStats {
  size_t index = 0;
  size_t linkTetrahedra = 0; /**< = vertex degree (number of incident pentachora). */
  size_t linkVertices = 0; /**< = number of ambient spokes at this vertex. */
  size_t linkEdges = 0;    /**< = number of ambient triangle-corners at this vertex. */
  int maxLinkVertexDegree = 0;
  long long closedCount = 0;
  long long visitedCount = 0;
  bool budgetExceeded = false;
  bool sawLoopEdge = false;
  std::chrono::steady_clock::duration elapsed{};
};

VertexStats analyzeVertex(const regina::Vertex<4> *v, long long budget) {
  VertexStats stats;
  stats.index = v->index();

  const auto &link = v->buildLink();
  stats.linkTetrahedra = link.size();
  stats.linkVertices = link.countVertices();
  stats.linkEdges = link.countEdges();
  if (stats.linkEdges == 0)
    return stats;

  std::vector<std::pair<int, int>> endpoints(stats.linkEdges);
  std::vector<std::vector<int>> incident(stats.linkVertices);
  for (size_t e = 0; e < stats.linkEdges; ++e) {
    const auto *edge = link.edge(e);
    int a = static_cast<int>(edge->vertex(0)->index());
    int b = static_cast<int>(edge->vertex(1)->index());
    endpoints[e] = {a, b};
    incident[static_cast<size_t>(a)].push_back(static_cast<int>(e));
    if (b != a)
      incident[static_cast<size_t>(b)].push_back(static_cast<int>(e));
  }

  for (const auto &lst : incident)
    stats.maxLinkVertexDegree =
        std::max(stats.maxLinkVertexDegree, static_cast<int>(lst.size()));

  // 1-indexed graph over link-edges (nodes), adjacent iff they share a link
  // vertex -- what ConnectedInducedSubgraphEnumerator expects.
  std::vector<std::vector<int>> graphAdj(stats.linkEdges + 1);
  for (const auto &lst : incident) {
    for (size_t i = 0; i < lst.size(); ++i) {
      for (size_t j = i + 1; j < lst.size(); ++j) {
        int e1 = lst[i] + 1, e2 = lst[j] + 1;
        graphAdj[static_cast<size_t>(e1)].push_back(e2);
        graphAdj[static_cast<size_t>(e2)].push_back(e1);
      }
    }
  }
  for (auto &lst : graphAdj) {
    std::ranges::sort(lst);
    lst.erase(std::ranges::unique(lst).begin(), lst.end());
  }

  SpokeCountPredicate predicate(endpoints, stats.linkVertices);
  std::atomic<bool> stop{false};
  std::atomic<bool> pause{false}; // never set here; this driver has no pause/drain concept
  InterruptiblePredicate interruptible(predicate, stop, pause);
  ConnectedInducedSubgraphEnumerator enumerator(static_cast<int>(stats.linkEdges),
                                                 graphAdj);

  const auto start = std::chrono::steady_clock::now();
  enumerator.enumerateFiltered(
      [&](const std::vector<int> &U) {
        ++stats.visitedCount;
        if (predicate.isClosed(U))
          ++stats.closedCount;
        if (stats.visitedCount >= budget) {
          stats.budgetExceeded = true;
          stop.store(true, std::memory_order_relaxed);
        }
      },
      interruptible);
  stats.elapsed = std::chrono::steady_clock::now() - start;
  stats.sawLoopEdge = predicate.sawLoopEdge;

  return stats;
}

// Direct simple-cycle counter over Lk(v)'s own 1-skeleton, for comparison
// against analyzeVertex() above. analyzeVertex() drives
// ConnectedInducedSubgraphEnumerator over the "link-edges share a
// link-vertex" graph with only a degree<=2 predicate -- but that predicate
// accepts every simple PATH as well as every simple CYCLE (a degree-<=2
// subgraph is just a disjoint union of paths and cycles), and simple paths
// vastly outnumber simple cycles in any graph with real branching. This
// walks simple paths directly via backtracking DFS and only counts closures
// back to the anchor, so path-only dead ends are never separately reported
// -- much closer to what an output-sensitive cycle enumerator (e.g.
// Johnson's algorithm) would cost, without yet implementing the "blocked
// set" optimization that makes those algorithms truly output-sensitive.
//
// Anchor convention: only vertices with id > the cycle's minimum-id vertex
// are ever extended into, so every cycle is discovered exactly twice (once
// per traversal direction) rooted at its own minimum vertex -- halve the
// raw count for the true number of distinct cycles.
struct CycleCountResult {
  long long rawCycleCount = 0; // before halving for direction double-counting
  long long stepsTaken = 0;
  bool budgetExceeded = false;
  std::chrono::steady_clock::duration elapsed{};
  std::vector<long long> lengthHistogram; // raw (pre-halving) count by cycle edge-length
};

CycleCountResult countSimpleCyclesDirect(const regina::Triangulation<3> &link,
                                          long long budget,
                                          int maxLen = -1) {
  CycleCountResult result;
  const size_t n = link.countVertices();
  if (n == 0)
    return result;

  std::vector<std::vector<std::pair<int, int>>> adj(n); // vertex -> (neighbor, edgeId)
  for (size_t e = 0; e < link.countEdges(); ++e) {
    const auto *edge = link.edge(e);
    int a = static_cast<int>(edge->vertex(0)->index());
    int b = static_cast<int>(edge->vertex(1)->index());
    if (a == b)
      continue; // loop edge -- can never close a valid simple cycle
    adj[static_cast<size_t>(a)].emplace_back(b, static_cast<int>(e));
    adj[static_cast<size_t>(b)].emplace_back(a, static_cast<int>(e));
  }

  std::vector<uint8_t> visited(n, 0);
  std::vector<uint8_t> usedEdge(link.countEdges(), 0);
  const auto start = std::chrono::steady_clock::now();
  bool stop = false;

  // path length counted in edges; anchor is the fixed root, cur the current
  // frontier vertex.
  std::function<void(int anchor, int cur, int pathEdges)> dfs =
      [&](int anchor, int cur, int pathEdges) {
        if (stop)
          return;
        ++result.stepsTaken;
        if (result.stepsTaken >= budget) {
          result.budgetExceeded = true;
          stop = true;
          return;
        }
        for (auto [nbr, edgeId] : adj[static_cast<size_t>(cur)]) {
          if (stop)
            return;
          if (usedEdge[static_cast<size_t>(edgeId)])
            continue;
          if (nbr == anchor) {
            // pathEdges counts edges from anchor to cur; closing needs a
            // *distinct* edge back to anchor (already guaranteed by the
            // usedEdge check above), so pathEdges >= 1 is enough -- this
            // also correctly counts a digon (2 parallel link-edges between
            // the same pair of vertices) as a valid 2-edge cycle.
            if (pathEdges >= 1) {
              ++result.rawCycleCount;
              auto len = static_cast<size_t>(pathEdges) + 1;
              if (len >= result.lengthHistogram.size())
                result.lengthHistogram.resize(len + 1, 0);
              ++result.lengthHistogram[len];
            }
            continue; // never step onto the anchor as a live frontier
          }
          if (nbr < anchor || visited[static_cast<size_t>(nbr)])
            continue;
          if (maxLen >= 0 && pathEdges + 1 >= maxLen)
            continue; // one more edge would already exceed the cap
          visited[static_cast<size_t>(nbr)] = 1;
          usedEdge[static_cast<size_t>(edgeId)] = 1;
          dfs(anchor, nbr, pathEdges + 1);
          usedEdge[static_cast<size_t>(edgeId)] = 0;
          visited[static_cast<size_t>(nbr)] = 0;
        }
      };

  for (size_t s = 0; s < n && !stop; ++s) {
    visited[s] = 1;
    dfs(static_cast<int>(s), static_cast<int>(s), 0);
    visited[s] = 0;
  }

  result.elapsed = std::chrono::steady_clock::now() - start;
  return result;
}

// Johnson's algorithm (1975), generalized from simple directed graphs to our
// undirected multigraph: CIRCUIT(v) explores v's neighbors (restricted to
// id >= s, so each cycle is found anchored at its own minimum-id vertex,
// same convention as countSimpleCyclesDirect), blocking v for the duration
// of its (possibly fruitless) exploration and only re-enabling it -- via
// UNBLOCK's cascade -- once something reachable through it actually closes a
// circuit back to s. That's the source of the output-sensitivity
// countSimpleCyclesDirect's plain backtracking lacks: a vertex whose subtree
// is proven fruitless under the current blocked set is never re-explored
// until an ancestor's success proves it might not be fruitless anymore.
//
// The only adaptation multi-edges need beyond the classic vertex-blocking
// scheme: track the single edge id used to arrive at the current vertex
// (`arrivalEdge`) and never immediately step back out over that same edge --
// this is what stops a lone undirected edge from being misread as a 2-cycle,
// while still correctly counting a genuine digon formed by two *distinct*
// parallel link-edges between the same pair of link-vertices. Every other
// repeat (revisiting a vertex further back on the path) is already excluded
// by `blocked[]`, so no full used-edge set is needed, unlike
// countSimpleCyclesDirect.
//
// Undirected cycles are still found twice (once per traversal direction),
// exactly as in countSimpleCyclesDirect -- halve rawCycleCount for the true
// count.
CycleCountResult johnsonCountCycles(const regina::Triangulation<3> &link,
                                     long long budget) {
  CycleCountResult result;
  const size_t n = link.countVertices();
  if (n == 0)
    return result;

  std::vector<std::vector<std::pair<int, int>>> adj(n); // vertex -> (neighbor, edgeId)
  for (size_t e = 0; e < link.countEdges(); ++e) {
    const auto *edge = link.edge(e);
    int a = static_cast<int>(edge->vertex(0)->index());
    int b = static_cast<int>(edge->vertex(1)->index());
    if (a == b)
      continue; // loop edge -- can never close a valid simple cycle
    adj[static_cast<size_t>(a)].emplace_back(b, static_cast<int>(e));
    adj[static_cast<size_t>(b)].emplace_back(a, static_cast<int>(e));
  }

  std::vector<uint8_t> blocked(n, 0);
  std::vector<std::vector<int>> B(n);
  const auto start = std::chrono::steady_clock::now();
  bool stop = false;
  int anchorS = 0;
  std::vector<int> pathLen(n, 0); // depth (in edges from s) at which each vertex currently sits, for the histogram

  std::function<void(int)> unblock = [&](int v) {
    blocked[static_cast<size_t>(v)] = 0;
    auto &bv = B[static_cast<size_t>(v)];
    for (int w : bv)
      if (blocked[static_cast<size_t>(w)])
        unblock(w);
    bv.clear();
  };

  std::function<bool(int, int, int)> circuit =
      [&](int v, int arrivalEdge, int depth) -> bool {
        if (stop)
          return false;
        ++result.stepsTaken;
        if (result.stepsTaken >= budget) {
          result.budgetExceeded = true;
          stop = true;
          return false;
        }
        bool f = false;
        blocked[static_cast<size_t>(v)] = 1;
        for (auto [w, edgeId] : adj[static_cast<size_t>(v)]) {
          if (stop)
            break;
          if (w < anchorS || edgeId == arrivalEdge)
            continue;
          if (w == anchorS) {
            ++result.rawCycleCount;
            auto len = static_cast<size_t>(depth) + 1;
            if (len >= result.lengthHistogram.size())
              result.lengthHistogram.resize(len + 1, 0);
            ++result.lengthHistogram[len];
            f = true;
          } else if (!blocked[static_cast<size_t>(w)]) {
            if (circuit(w, edgeId, depth + 1))
              f = true;
          }
        }
        if (f) {
          unblock(v);
        } else {
          for (auto [w, edgeId] : adj[static_cast<size_t>(v)]) {
            (void)edgeId;
            if (w < anchorS)
              continue;
            auto &bw = B[static_cast<size_t>(w)];
            if (std::ranges::find(bw, v) == bw.end())
              bw.push_back(v);
          }
        }
        return f;
      };

  for (size_t s = 0; s < n && !stop; ++s) {
    anchorS = static_cast<int>(s);
    std::ranges::fill(blocked, 0);
    for (auto &bv : B)
      bv.clear();
    circuit(anchorS, -1, 0);
  }

  result.elapsed = std::chrono::steady_clock::now() - start;
  return result;
}

void run(const regina::Triangulation<4> &tri, long long budgetPerVertex) {
  long long totalClosed = 0;
  size_t verticesAnalyzed = 0, verticesSkippedBoundary = 0,
         verticesExceededBudget = 0, verticesWithLoopEdges = 0;
  size_t maxLinkEdges = 0;
  int maxLinkVertexDegree = 0;
  const auto start = std::chrono::steady_clock::now();

  for (const auto &v : tri.vertices()) {
    if (v->isBoundary()) {
      ++verticesSkippedBoundary;
      continue;
    }
    VertexStats s = analyzeVertex(v, budgetPerVertex);
    ++verticesAnalyzed;
    totalClosed += s.closedCount;
    maxLinkEdges = std::max(maxLinkEdges, s.linkEdges);
    maxLinkVertexDegree = std::max(maxLinkVertexDegree, s.maxLinkVertexDegree);
    if (s.budgetExceeded)
      ++verticesExceededBudget;
    if (s.sawLoopEdge)
      ++verticesWithLoopEdges;

    CycleCountResult direct = countSimpleCyclesDirect(v->buildLink(), budgetPerVertex);

    std::cerr << "[vertexstats] v=" << s.index
               << " linkTetrahedra=" << s.linkTetrahedra
               << " linkVertices=" << s.linkVertices
               << " linkEdges=" << s.linkEdges
               << " maxLinkVertexDegree=" << s.maxLinkVertexDegree
               << " visited=" << s.visitedCount
               << " closedPetals(CIS)=" << s.closedCount
               << (s.budgetExceeded ? " [CIS BUDGET EXCEEDED]" : "")
               << (s.sawLoopEdge ? " [SAW LOOP EDGE]" : "")
               << " cycles(direct)=" << (direct.rawCycleCount / 2)
               << " steps(direct)=" << direct.stepsTaken
               << (direct.budgetExceeded ? " [DIRECT BUDGET EXCEEDED]" : "")
               << " elapsed(CIS)=" << formatElapsed(s.elapsed)
               << " elapsed(direct)=" << formatElapsed(direct.elapsed) << "\n";

    std::cerr << "[vertexstats]   length histogram (raw, pre-halving):";
    for (size_t len = 0; len < direct.lengthHistogram.size(); ++len)
      if (direct.lengthHistogram[len] > 0)
        std::cerr << " " << len << ":" << direct.lengthHistogram[len];
    std::cerr << "\n";

    // TEMPORARY (Phase 0 exploration): exhaustive counts under a cycle-length
    // cap, to see whether bounding petal length tames the blowup while still
    // covering the common case.
    for (int cap : {6, 8, 10, 12, 14, 16, 20}) {
      CycleCountResult capped =
          countSimpleCyclesDirect(v->buildLink(), 200000000, cap);
      std::cerr << "[vertexstats]   maxLen=" << cap
                 << " cycles=" << (capped.rawCycleCount / 2)
                 << " steps=" << capped.stepsTaken
                 << (capped.budgetExceeded ? " [EXCEEDED]" : " [EXHAUSTIVE]")
                 << " elapsed=" << formatElapsed(capped.elapsed) << "\n";
    }

    // Johnson's algorithm, uncapped: cross-check against countSimpleCyclesDirect
    // (both total count and per-length histogram should agree wherever
    // neither hit budget), and see how much further the blocking
    // optimization gets within the same step budget.
    CycleCountResult johnson = johnsonCountCycles(v->buildLink(), 200000000);
    std::cerr << "[vertexstats]   johnson: cycles=" << (johnson.rawCycleCount / 2)
               << " steps=" << johnson.stepsTaken
               << (johnson.budgetExceeded ? " [EXCEEDED]" : " [EXHAUSTIVE]")
               << " elapsed=" << formatElapsed(johnson.elapsed) << "\n";
    std::cerr << "[vertexstats]   johnson length histogram (raw, pre-halving):";
    for (size_t len = 0; len < johnson.lengthHistogram.size(); ++len)
      if (johnson.lengthHistogram[len] > 0)
        std::cerr << " " << len << ":" << johnson.lengthHistogram[len];
    std::cerr << "\n";
  }

  const auto elapsed = std::chrono::steady_clock::now() - start;
  std::cerr << "[vertexstats] SUMMARY: vertices=" << tri.countVertices()
             << " analyzed=" << verticesAnalyzed
             << " skippedBoundary=" << verticesSkippedBoundary
             << " exceededBudget=" << verticesExceededBudget
             << " withLoopEdges=" << verticesWithLoopEdges
             << " maxLinkEdges=" << maxLinkEdges
             << " maxLinkVertexDegree=" << maxLinkVertexDegree
             << " totalClosedPetals=" << totalClosed
             << " totalElapsed=" << formatElapsed(elapsed) << "\n";
}

} // namespace vertexstats

} // namespace

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage:\n"
                 "  profile_driver bary3 <subdivisions> <cond> <threads>\n"
                 "  profile_driver bary4 <subdivisions> <cond> <threads>\n"
                 "  profile_driver isosig4 <isosig> <cond> <threads>\n"
                 "  profile_driver pipeline <pdcode> <thickenLayers> <cond> "
                 "<threads>\n"
                 "  profile_driver diag <pdcode> <thickenLayers> <cond> "
                 "<reportEvery>\n"
                 "  profile_driver surfacediag <pdcode> <thickenLayers> "
                 "<cond> <reportEvery>\n"
                 "  profile_driver vertexstats bary4 <subdivisions> "
                 "<budgetPerVertex>\n"
                 "  profile_driver vertexstats isosig4 <isosig> "
                 "<budgetPerVertex>\n"
                 "  profile_driver vertexstats pipeline <pdcode> "
                 "<thickenLayers> <budgetPerVertex>\n";
    return 1;
  }

  std::string mode = argv[1];

  if (mode == "bary3") {
    int subdivisions = std::stoi(argv[2]);
    BoundaryCondition cond = parseCond(argv[3]);
    unsigned threads = static_cast<unsigned>(std::stoi(argv[4]));

    regina::Triangulation<3> tri;
    tri.newSimplex();
    for (int i = 0; i < subdivisions; ++i)
      tri.subdivide();

    std::cerr << "[driver] dim3 tetrahedra = " << tri.size()
              << ", triangles = " << tri.countTriangles() << "\n";

    EmbeddingSearch<3, 2> search(tri);
    std::cerr << "[driver] embeddable faces = " << search.numEmbeddableFaces()
              << "\n";
    runTimed(search, threads, cond);
  } else if (mode == "bary4") {
    int subdivisions = std::stoi(argv[2]);
    BoundaryCondition cond = parseCond(argv[3]);
    unsigned threads = static_cast<unsigned>(std::stoi(argv[4]));

    regina::Triangulation<4> tri;
    tri.newSimplex();
    for (int i = 0; i < subdivisions; ++i)
      tri.subdivide();

    std::cerr << "[driver] dim4 pentachora = " << tri.size()
              << ", triangles = " << tri.countTriangles() << "\n";

    SurfaceSearch search(tri);
    std::cerr << "[driver] embeddable faces = " << search.numEmbeddableFaces()
              << "\n";
    runTimed(search, threads, cond);
  } else if (mode == "isosig4") {
    std::string isosig = argv[2];
    BoundaryCondition cond = parseCond(argv[3]);
    unsigned threads = static_cast<unsigned>(std::stoi(argv[4]));

    regina::Triangulation<4> tri(isosig);
    std::cerr << "[driver] dim4 pentachora = " << tri.size()
              << ", triangles = " << tri.countTriangles() << "\n";

    SurfaceSearch search(tri);
    std::cerr << "[driver] embeddable faces = " << search.numEmbeddableFaces()
              << "\n";
    runTimed(search, threads, cond);
  } else if (mode == "pipeline" || mode == "pipeline-loud") {
    std::string pdcodeStr = argv[2];
    int layers = std::stoi(argv[3]);
    BoundaryCondition cond = parseCond(argv[4]);
    unsigned threads = static_cast<unsigned>(std::stoi(argv[5]));

    knotbuilder::PDCode pdcode = knotbuilder::parsePDCode(pdcodeStr);
    auto [t2, edges0] = knotbuilder::buildLink(pdcode);

    CobordismBuilder<3> cob(t2);
    if (layers > 0)
      cob.thicken(layers);
    regina::Triangulation<4> tri = cob.cone();

    std::cerr << "[driver] dim4 pentachora = " << tri.size()
              << ", triangles = " << tri.countTriangles() << "\n";

    SurfaceSearch search(tri);
    std::cerr << "[driver] embeddable faces = " << search.numEmbeddableFaces()
              << "\n";
    if (mode == "pipeline-loud") {
      SurfaceSearchCallbacks callbacks;
      callbacks.onProgress = [](const SearchStats &stats) {
        std::cerr << "[driver] elapsed=" << formatElapsed(stats.elapsed)
                  << " roots=" << stats.rootsCompleted << "/"
                  << stats.totalRoots << " found=" << stats.foundCount
                  << " satisfying=" << stats.satisfyingCount << "\n";
      };
      search.search(threads, cond, callbacks); // live progress reports
    } else {
      runTimed(search, threads, cond);
    }
  } else if (mode == "diag") {
    std::string pdcodeStr = argv[2];
    int layers = std::stoi(argv[3]);
    BoundaryCondition cond = parseCond(argv[4]);
    long long reportEvery = std::stoll(argv[5]);

    knotbuilder::PDCode pdcode = knotbuilder::parsePDCode(pdcodeStr);
    auto [t2, edges0] = knotbuilder::buildLink(pdcode);

    CobordismBuilder<3> cob(t2);
    if (layers > 0)
      cob.thicken(layers);
    regina::Triangulation<4> tri = cob.cone();

    std::cerr << "[driver] dim4 pentachora = " << tri.size()
              << ", triangles = " << tri.countTriangles() << "\n";

    DiagSearch search(tri);
    std::cerr << "[driver] embeddable faces = " << search.numEmbeddableFaces()
              << "\n";
    search.diagRun(cond, reportEvery);
  } else if (mode == "surfacediag") {
    std::string pdcodeStr = argv[2];
    int layers = std::stoi(argv[3]);
    BoundaryCondition cond = parseCond(argv[4]);
    long long reportEvery = std::stoll(argv[5]);

    knotbuilder::PDCode pdcode = knotbuilder::parsePDCode(pdcodeStr);
    auto [t2, edges0] = knotbuilder::buildLink(pdcode);

    CobordismBuilder<3> cob(t2);
    if (layers > 0)
      cob.thicken(layers);
    regina::Triangulation<4> tri = cob.cone();

    std::cerr << "[driver] dim4 pentachora = " << tri.size()
              << ", triangles = " << tri.countTriangles() << "\n";

    SurfaceDiagSearch search(tri);
    std::cerr << "[driver] embeddable faces = " << search.numEmbeddableFaces()
              << "\n";
    search.diagRun(cond, reportEvery);
  } else if (mode == "vertexstats") {
    std::string submode = argv[2];

    if (submode == "bary4") {
      int subdivisions = std::stoi(argv[3]);
      long long budget = std::stoll(argv[4]);

      regina::Triangulation<4> tri;
      tri.newSimplex();
      for (int i = 0; i < subdivisions; ++i)
        tri.subdivide();

      std::cerr << "[driver] dim4 pentachora = " << tri.size()
                << ", triangles = " << tri.countTriangles() << "\n";
      vertexstats::run(tri, budget);
    } else if (submode == "isosig4") {
      std::string isosig = argv[3];
      long long budget = std::stoll(argv[4]);

      regina::Triangulation<4> tri(isosig);
      std::cerr << "[driver] dim4 pentachora = " << tri.size()
                << ", triangles = " << tri.countTriangles() << "\n";
      vertexstats::run(tri, budget);
    } else if (submode == "pipeline") {
      std::string pdcodeStr = argv[3];
      int layers = std::stoi(argv[4]);
      long long budget = std::stoll(argv[5]);

      knotbuilder::PDCode pdcode = knotbuilder::parsePDCode(pdcodeStr);
      auto [t2, edges0] = knotbuilder::buildLink(pdcode);

      CobordismBuilder<3> cob(t2);
      if (layers > 0)
        cob.thicken(layers);
      regina::Triangulation<4> tri = cob.cone();

      std::cerr << "[driver] dim4 pentachora = " << tri.size()
                << ", triangles = " << tri.countTriangles() << "\n";
      vertexstats::run(tri, budget);
    } else if (submode == "tiny") {
      // Deliberately small closed 4-manifold (2-pentachoron 4-sphere), so
      // both countSimpleCyclesDirect and johnsonCountCycles can run to
      // completion (no budget/cap needed) for a direct correctness
      // cross-check between the two algorithms.
      long long budget = std::stoll(argv[3]);
      regina::Triangulation<4> tri = regina::Example<4>::fourSphere();
      std::cerr << "[driver] dim4 pentachora = " << tri.size()
                << ", triangles = " << tri.countTriangles() << "\n";
      vertexstats::run(tri, budget);
    } else {
      std::cerr << "Unknown vertexstats submode: " << submode << "\n";
      return 1;
    }
  } else {
    std::cerr << "Unknown mode: " << mode << "\n";
    return 1;
  }

  return 0;
}
