//
//  enumerate_cis.cpp
//

#include "enumerate_cis.h"

void ConnectedInducedSubgraphEnumerator::enumerateFromRoot(
    int s, const std::function<void(const std::vector<int> &)> &visit) {
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
