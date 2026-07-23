//
//  enumerate_cis.cpp
//

#include "enumerate_cis.h"

#include <algorithm>
#include <cassert>
#include <queue>
#include <set>
#include <unordered_set>

ConnectedInducedSubgraphEnumerator::SeededGraph
ConnectedInducedSubgraphEnumerator::contractSeed(
    int origN, const std::vector<std::vector<int>> &origAdj,
    const std::vector<int> &seed) {
    std::unordered_set<int> S(seed.begin(), seed.end());

    // Sanity check: seed must be connected as an induced subgraph.
    {
        std::vector<char> seen(origN + 1, 0);
        std::queue<int> q;
        seen[seed[0]] = 1;
        q.push(seed[0]);
        int count = 1;
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            for (int v : origAdj[u])
                if (S.count(v) && !seen[v]) {
                    seen[v] = 1;
                    ++count;
                    q.push(v);
                }
        }
        assert(count == (int)seed.size() &&
              "seed must be a connected induced subgraph");
    }

    // External vertices, in a fixed (ascending) order -> new ids 2..n'
    std::vector<int> externalVerts;
    for (int v = 1; v <= origN; ++v)
        if (!S.count(v))
            externalVerts.push_back(v);

    std::vector<int> newIdOf(origN + 1, 0); // original id -> new id (external only)
    for (size_t i = 0; i < externalVerts.size(); ++i)
        newIdOf[externalVerts[i]] = (int)i + 2;

    int nPrime = 1 + (int)externalVerts.size();
    SeededGraph g;
    g.n = nPrime;
    g.adj.assign(nPrime + 1, {});
    g.originalOf.assign(nPrime + 1, -1); // originalOf[1] left as -1: "this is the seed"
    for (size_t i = 0; i < externalVerts.size(); ++i)
        g.originalOf[i + 2] = externalVerts[i];

    // Seed (vertex 1) <-> external neighbors: union over all seed vertices,
    // deduplicated (a vertex adjacent to several seed members still only
    // gets ONE edge to the contracted vertex).
    std::set<int> seedNeighbors;
    for (int v : seed)
        for (int u : origAdj[v])
            if (!S.count(u))
                seedNeighbors.insert(newIdOf[u]);
    for (int u : seedNeighbors) {
        g.adj[1].push_back(u);
        g.adj[u].push_back(1);
    }

    // External <-> external edges: unchanged, just relabeled.
    for (int v : externalVerts) {
        int vNew = newIdOf[v];
        for (int u : origAdj[v]) {
            if (S.count(u))
                continue; // handled above
            int uNew = newIdOf[u];
            if (uNew > vNew) {
                g.adj[vNew].push_back(uNew);
                g.adj[uNew].push_back(vNew);
            }
        }
    }
    for (auto &lst : g.adj) {
        std::sort(lst.begin(), lst.end());
        lst.erase(std::unique(lst.begin(), lst.end()), lst.end());
    }
    return g;
}

void ConnectedInducedSubgraphEnumerator::enumerateFromRoot(
    int s, const std::function<void(const std::vector<int> &)> &visit) {
    if (isSeeded_) {
        // s is a sibling of the seed, already sitting in the candidate list
        // from seedFastForward_(). Anchor stays 1 throughout -- see
        // seedFastForward_'s comment and the class-level isSeeded_ comment.
        report = &visit;
        const int w = s;

        detachToU(w);
        U.push_back(w);
        inU[w] = true;

        std::vector<int> &introduced = introducedBuf[U.size()];
        introduced.clear();
        for (int x : adj[w])
            if (x > 1 && !inU[x] && !inC[x]) {
                addCandidate(x, w, dist[w] + 1);
                introduced.push_back(x);
            }

        (*report)(U);
        extend(1);

        for (int x : introduced)
            removeCandidate(x);
        U.pop_back();
        inU[w] = false;
        addCandidate(w, 1, 1); // restore as a candidate for the next root

        return;
    }

    report = &visit;

    U.push_back(s);
    inU[s] = true;
    dist[s] = 0;
    (*report)(U); // output {s}

    for (int w : adj[s])
        if (w > s && !inU[w] && !inC[w])
            addCandidate(w, s, 1);

    extend(s);

    while (!listEmpty())
        removeCandidate(listFront()); // reset before next root
    U.pop_back();
    inU[s] = false;
    dist[s] = -1;
}

void ConnectedInducedSubgraphEnumerator::enumerateFromRootFiltered(
    int s, const std::function<void(const std::vector<int> &)> &visit,
    ConditionalPredicate &predicate) {
    if (isSeeded_) {
        report = &visit;
        const int w = s;

        detachToU(w);
        U.push_back(w);
        inU[w] = true;

        std::vector<int> &introduced = introducedBuf[U.size()];
        introduced.clear();
        for (int x : adj[w])
            if (x > 1 && !inU[x] && !inC[x]) {
                addCandidate(x, w, dist[w] + 1);
                introduced.push_back(x);
            }

        if (predicate.tryAdd(w)) {
            (*report)(U);
            extendFiltered(1, predicate);
            predicate.undo(w);
        }

        for (int x : introduced)
            removeCandidate(x);
        U.pop_back();
        inU[w] = false;
        addCandidate(w, 1, 1); // restore as a candidate for the next root

        return;
    }

    report = &visit;

    U.push_back(s);
    inU[s] = true;
    dist[s] = 0;

    if (predicate.tryAdd(s)) {
        (*report)(U); // output {s}

        for (int w : adj[s])
            if (w > s && !inU[w] && !inC[w])
                addCandidate(w, s, 1);

        extendFiltered(s, predicate);

        while (!listEmpty())
            removeCandidate(listFront()); // reset before next root
        predicate.undo(s);
    }

    U.pop_back();
    inU[s] = false;
    dist[s] = -1;
}

void ConnectedInducedSubgraphEnumerator::seedFastForward_(
    ConditionalPredicate *predicate) {
    U.push_back(1);
    inU[1] = true;
    dist[1] = 0;

    if (predicate) {
        // The caller (e.g. EmbeddingSearch) is expected to have already
        // validated that the seed can jointly be committed before ever
        // constructing a seeded enumerator -- see EmbeddingSearch's seeded
        // constructor, which does exactly this validation eagerly. This is
        // a defensive check for standalone/direct use of this class, not
        // the primary place that validation happens.
        [[maybe_unused]] bool committed = predicate->tryAdd(1);
        assert(committed && "seed must be jointly committable via tryAdd(1)");
    }

    for (int w : adj[1])
        if (w > 1 && !inU[w] && !inC[w])
            addCandidate(w, 1, 1);

    roots_.clear();
    for (int w = candNext[0]; w != 0; w = candNext[w]) {
        if (!predicate) {
            roots_.push_back(w);
            continue;
        }
        if (predicate->tryAdd(w)) {
            predicate->undo(w);
            roots_.push_back(w);
        }
    }
}

void ConnectedInducedSubgraphEnumerator::addCandidate(int v, int par, int d) {
    inC[v] = true;
    dist[v] = d;
    parentOf[v] = par;
    listPushBack(v);
}

void ConnectedInducedSubgraphEnumerator::removeCandidate(int v) {
    listErase(v);
    inC[v] = false;
    dist[v] = -1;
    parentOf[v] = 0;
}

void ConnectedInducedSubgraphEnumerator::detachToU(int v) {
    listErase(v);
    inC[v] = false;
}

void ConnectedInducedSubgraphEnumerator::extend(int s) {
    const int u = U.back(); // utmost(U)
    const int du = dist[u];

    // Freeze the sibling list: vertices discovered while branching into
    // one candidate belong to a deeper recursion level and must not be
    // tried as siblings at this level. Reuse this depth's scratch buffer
    // (see siblingBuf's declaration) instead of heap-allocating a fresh
    // vector on every call.
    std::vector<int> &siblings = siblingBuf[U.size()];
    siblings.clear();
    for (int v = candNext[0]; v != 0; v = candNext[v])
        siblings.push_back(v);

    for (int w : siblings) {
        const int dw = dist[w];
        const bool validChild = (w > s) && (dw > du || (dw == du && w > u));
        if (!validChild)
            continue;

        const int wDist = dw;
        const int wParent = parentOf[w];

        // --- descend: move w from C into U ---
        detachToU(w);
        U.push_back(w);
        inU[w] = true;

        std::vector<int> &introduced = introducedBuf[U.size()];
        introduced.clear();
        for (int x : adj[w])
            if (x > s && !inU[x] && !inC[x]) {
                addCandidate(x, w, wDist + 1);
                introduced.push_back(x);
            }

        (*report)(U); // output U (parent set plus w)
        extend(s);

        // --- backtrack ---
        for (int x : introduced)
            removeCandidate(x);
        U.pop_back();
        inU[w] = false;
        addCandidate(w, wParent,
                     wDist); // restore w as a candidate for later siblings
    }
}

void ConnectedInducedSubgraphEnumerator::extendFiltered(
    int s, ConditionalPredicate &predicate) {
    const int u = U.back();
    const int du = dist[u];

    // See extend()'s identical snapshot for why this is a reused per-depth
    // buffer rather than a fresh std::vector<int>(C.begin(), C.end()).
    std::vector<int> &siblings = siblingBuf[U.size()];
    siblings.clear();
    for (int v = candNext[0]; v != 0; v = candNext[v])
        siblings.push_back(v);

    for (int w : siblings) {
        const int dw = dist[w];
        const bool validChild = (w > s) && (dw > du || (dw == du && w > u));
        if (!validChild)
            continue;

        const int wDist = dw;
        const int wParent = parentOf[w];

        detachToU(w);
        U.push_back(w);
        inU[w] = true;

        std::vector<int> &introduced = introducedBuf[U.size()];
        introduced.clear();
        for (int x : adj[w])
            if (x > s && !inU[x] && !inC[x]) {
                addCandidate(x, w, wDist + 1);
                introduced.push_back(x);
            }

        if (predicate.tryAdd(w)) { // local, incremental check
            (*report)(U);
            extendFiltered(s, predicate); // only descend on a pass
            predicate.undo(w); // reverse tryAdd(w) -- see class contract
        }
        // else: tryAdd made no net change (transactional contract), so
        // there's nothing to undo -- just prune.

        for (int x : introduced)
            removeCandidate(x);
        U.pop_back();
        inU[w] = false;
        addCandidate(w, wParent, wDist);
    }
}
