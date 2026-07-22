#ifndef SEEDED_SEARCH_H
#define SEEDED_SEARCH_H

// Enumerate all connected induced subgraphs of G that are supersets of a
// given connected seed S, by contracting S into a single vertex with the
// globally smallest id, then delegating entirely to the existing,
// unmodified ConnectedInducedSubgraphEnumerator on the contracted graph.
//
// Why this works: enumerateFromRoot(s, ...) finds every connected set whose
// anchor (smallest-id member) is s. If s is the smallest id in the WHOLE
// graph, "anchor is s" and "contains s" are the same statement -- so
// enumerateFromRoot(1, ...) on the contracted graph finds every connected
// set containing the contracted seed, which (since S is connected) is
// exactly every connected superset of S in the original graph.

#include "enumerate_cis.h"
#include <unordered_set>
#include <set>
#include <algorithm>
#include <cassert>
#include <queue>

struct SeededGraph {
    int n;                            // vertex count of the contracted graph
    std::vector<std::vector<int>> adj; // 1-indexed; vertex 1 = contracted seed
    std::vector<int> originalOf;       // originalOf[i], i>=2: original id of contracted vertex i
};

// Precondition: `seed` must itself induce a connected subgraph of
// (origN, origAdj). Checked with an assert (BFS within the seed).
inline SeededGraph buildSeededGraph(int origN,
                                     const std::vector<std::vector<int>>& origAdj,
                                     const std::vector<int>& seed) {
    std::unordered_set<int> S(seed.begin(), seed.end());

    // Sanity check: seed must be connected as an induced subgraph.
    {
        std::vector<char> seen(origN + 1, 0);
        std::queue<int> q;
        seen[seed[0]] = 1; q.push(seed[0]);
        int count = 1;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int v : origAdj[u])
                if (S.count(v) && !seen[v]) { seen[v] = 1; ++count; q.push(v); }
        }
        assert(count == (int)seed.size() && "seed must be a connected induced subgraph");
    }

    // External vertices, in a fixed (ascending) order -> new ids 2..n'
    std::vector<int> externalVerts;
    for (int v = 1; v <= origN; ++v) if (!S.count(v)) externalVerts.push_back(v);

    std::vector<int> newIdOf(origN + 1, 0);   // original id -> new id (external only)
    for (size_t i = 0; i < externalVerts.size(); ++i) newIdOf[externalVerts[i]] = (int)i + 2;

    int nPrime = 1 + (int)externalVerts.size();
    SeededGraph g;
    g.n = nPrime;
    g.adj.assign(nPrime + 1, {});
    g.originalOf.assign(nPrime + 1, -1);   // originalOf[1] left as -1: "this is the contracted seed"
    for (size_t i = 0; i < externalVerts.size(); ++i) g.originalOf[i + 2] = externalVerts[i];

    // Seed (vertex 1) <-> external neighbors: union over all seed vertices,
    // deduplicated (a vertex adjacent to several seed members still only
    // gets ONE edge to the contracted vertex).
    std::set<int> seedNeighbors;
    for (int v : seed)
        for (int u : origAdj[v])
            if (!S.count(u)) seedNeighbors.insert(newIdOf[u]);
    for (int u : seedNeighbors) { g.adj[1].push_back(u); g.adj[u].push_back(1); }

    // External <-> external edges: unchanged, just relabeled.
    for (int v : externalVerts) {
        int vNew = newIdOf[v];
        for (int u : origAdj[v]) {
            if (S.count(u)) continue;             // handled above
            int uNew = newIdOf[u];
            if (uNew > vNew) { g.adj[vNew].push_back(uNew); g.adj[uNew].push_back(vNew); }
        }
    }
    for (auto& lst : g.adj) {
        std::sort(lst.begin(), lst.end());
        lst.erase(std::unique(lst.begin(), lst.end()), lst.end());
    }
    return g;
}

// Enumerates every connected induced subgraph of (origN, origAdj) that is a
// superset of `seed`, calling visit(U) once per result (U given as original
// vertex ids, unsorted: seed vertices first, then any added external ones).
inline void enumerateSeeded(int origN, const std::vector<std::vector<int>>& origAdj,
                             const std::vector<int>& seed,
                             const std::function<void(const std::vector<int>&)>& visit) {
    SeededGraph g = buildSeededGraph(origN, origAdj, seed);
    ConnectedInducedSubgraphEnumerator enumerator(g.n, g.adj);
    enumerator.enumerateFromRoot(1, [&](const std::vector<int>& U) {
        std::vector<int> full = seed;
        for (int v : U) if (v != 1) full.push_back(g.originalOf[v]);
        visit(full);
    });
}

// Same, but filtered by a hereditary predicate (stateless or
// ConditionalPredicate), exactly as with the unmodified enumerator.
template <class Predicate>
inline void enumerateSeededFiltered(int origN, const std::vector<std::vector<int>>& origAdj,
                                     const std::vector<int>& seed,
                                     const std::function<void(const std::vector<int>&)>& visit,
                                     Predicate& predicate) {
    SeededGraph g = buildSeededGraph(origN, origAdj, seed);
    ConnectedInducedSubgraphEnumerator enumerator(g.n, g.adj);
    enumerator.enumerateFromRootFiltered(1, [&](const std::vector<int>& U) {
        std::vector<int> full = seed;
        for (int v : U) if (v != 1) full.push_back(g.originalOf[v]);
        visit(full);
    }, predicate);
}

#endif // SEEDED_SEARCH_H
