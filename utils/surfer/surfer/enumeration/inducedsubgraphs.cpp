//
//  enumerate_cis.cpp
//
//  Created by John Teague on 07/21/2026.
//

#include "surfer/enumeration/inducedsubgraphs.h"

#include <algorithm>
#include <cassert>
#include <queue>
#include <set>
#include <stdexcept>
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

ConnectedInducedSubgraphEnumerator::Outcome
ConnectedInducedSubgraphEnumerator::enumerateFromRootFiltered(
    int s, const std::function<void(const std::vector<int> &)> &visit,
    ConditionalPredicate &predicate, Position *position) {
    // Resume from the recorded position, if any, and record afresh.
    Position resumeFrom;
    const Position *resume = nullptr;
    if (position && !position->empty()) {
        resumeFrom = std::move(*position);
        resume = &resumeFrom;
    }
    if (position) {
        position->levels.clear();
        position->deepest = 0;
    }
    // A rejected root: a budget refusal suspends the root before it starts
    // (it resumes from scratch); a refused re-add while resuming can only be
    // an external stop.
    auto rootRefused = [&] {
        if (resume)
            return Outcome::stopped;
        return position && predicate.budgetExhausted() ? Outcome::suspended
                                                       : Outcome::completed;
    };
    // A stopped pass never got back down to where it was to carry on, so it
    // leaves the position as it came in: the root can still resume from
    // there (a SearchFrontier records it).
    auto finish = [&](Outcome outcome) {
        if (outcome == Outcome::stopped && position && resume)
            *position = std::move(resumeFrom);
        return outcome;
    };

    if (isSeeded_) {
        if (maxSize_ && U.size() >= maxSize_)
            return Outcome::completed; // a cap of 0 added faces: not even the root
        report = &visit;
        const int w = s;
        const int prev = candPrev[w], next = candNext[w];

        detachToU(w);
        U.push_back(w);
        inU[w] = true;

        Outcome outcome;
        // Introduce w's neighbours only once w passes; see extendFiltered().
        if (predicate.tryAdd(w)) {
            std::vector<int> &introduced = introducedBuf[U.size()];
            introduced.clear();
            for (int x : adj[w])
                if (x > 1 && !inU[x] && !inC[x]) {
                    addCandidate(x, w, dist[w] + 1);
                    introduced.push_back(x);
                }

            if (!resume)
                (*report)(U); // a resumed root was reported by an earlier pass
            outcome = extendFiltered(1, predicate, resume, position);
            predicate.undo(w);

            for (int x : introduced)
                removeCandidate(x);
        } else {
            outcome = rootRefused();
        }

        U.pop_back();
        inU[w] = false;
        // Back where it was, so the list is canonical again for the next
        // root; see relink().
        relink(w, prev, next);

        return finish(outcome);
    }

    report = &visit;

    U.push_back(s);
    inU[s] = true;
    dist[s] = 0;

    Outcome outcome;
    if (predicate.tryAdd(s)) {
        if (!resume)
            (*report)(U); // output {s}

        for (int w : adj[s])
            if (w > s && !inU[w] && !inC[w])
                addCandidate(w, s, 1);

        outcome = extendFiltered(s, predicate, resume, position);

        while (!listEmpty())
            removeCandidate(listFront()); // reset before next root
        predicate.undo(s);
    } else {
        outcome = rootRefused();
    }

    U.pop_back();
    inU[s] = false;
    dist[s] = -1;
    return finish(outcome);
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

    // A neighbour of the seed that fails with the seed alone is pruned for
    // good: the predicate is anti-monotonic, and every set this enumerator
    // will visit contains the seed, so none containing that neighbour can
    // pass. Leave it out of the list, but keep it in C, so no later vertex
    // re-introduces it. On the profiled row, 880 of 1,019 seed neighbours:
    // 84% of every scan, and 91% of the failing tryAdd() calls.
    //
    // The list's order right now is canonical: every later root restores
    // membership but rotates the order, so this is the state
    // resetCandidateOrder() returns to.
    canonicalCandidates_ = roots_;
    resetCandidateOrder();
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

ConnectedInducedSubgraphEnumerator::Outcome
ConnectedInducedSubgraphEnumerator::extendFiltered(
    int s, ConditionalPredicate &predicate, const Position *resume,
    Position *record) {
    const size_t level = U.size();
    if (maxSize_ && level >= maxSize_)
        return Outcome::completed; // at the cap; see setMaxSize()
    const int u = U.back();
    const int du = dist[u];

    // See extend()'s identical snapshot for why this is a reused per-depth
    // buffer rather than a fresh std::vector<int>(C.begin(), C.end()).
    std::vector<int> &siblings = siblingBuf[level];
    siblings.clear();
    for (int v = candNext[0]; v != 0; v = candNext[v])
        siblings.push_back(v);
    std::vector<std::array<int, 3>> &pruned = prunedBuf[level];
    pruned.clear();

    // Resuming: this node's snapshot is the one the suspended pass saw
    // (the list is a function of the path; see relink()), so put its pruned
    // children aside again, in the same order, and pick up the loop at the
    // recorded child. Above the deepest level that child is on the recorded
    // path: re-added without being reported, and descended into. At the
    // deepest level it is the one whose attempt the budget refused: tried
    // afresh.
    const bool resumeHere = resume && level <= resume->deepest &&
                            level < resume->levels.size();
    size_t start = 0;
    if (resumeHere) {
        const Position::Level &at = resume->levels[level];
        for (int p : at.pruned) {
            const int prev = candPrev[p], next = candNext[p];
            listErase(p); // stays in C
            pruned.push_back({p, prev, next});
        }
        start = static_cast<size_t>(at.childIndex);
    }

    // Back in place, last pruned first, so the list leaves this node exactly
    // as it came in; see relink().
    auto restorePruned = [&] {
        for (auto it = pruned.rbegin(); it != pruned.rend(); ++it)
            relink((*it)[0], (*it)[1], (*it)[2]);
    };
    auto recordHere = [&](size_t index, int child) {
        if (record->levels.size() <= level)
            record->levels.resize(level + 1);
        Position::Level &at = record->levels[level];
        at.childIndex = static_cast<long>(index);
        at.child = child;
        at.pruned.clear();
        for (const auto &p : pruned)
            at.pruned.push_back(p[0]);
    };

    for (size_t i = start; i < siblings.size(); ++i) {
        const int w = siblings[i];
        const bool rebuilding = resumeHere && i == start &&
                                level < resume->deepest;
        const int dw = dist[w];
        const bool validChild = (w > s) && (dw > du || (dw == du && w > u));
        // A rebuilt node must offer exactly the child the suspended pass was
        // exploring; anything else means the rebuild went wrong, and carrying
        // on would silently skip or repeat part of the tree.
        if (rebuilding && (w != resume->levels[level].child || !validChild))
            throw std::logic_error(
                "ConnectedInducedSubgraphEnumerator: a resumed pass did not "
                "rebuild the snapshot it recorded");
        if (!validChild)
            continue;

        const int wDist = dw;
        // w's list neighbours, to put it back exactly there; see relink().
        const int prev = candPrev[w], next = candNext[w];

        detachToU(w);
        U.push_back(w);
        inU[w] = true;

        // w's neighbours join C only once w passes. Nearly every child fails
        // (97% on a profiled row), and introducing its neighbours only to
        // remove them again cost that row ~30% of its wall time. Doing it
        // after the check changes nothing observable: no predicate reads C,
        // and a vertex outside U and C always has dist -1 and parentOf 0, so
        // introducing and then removing candidates restores C exactly.
        if (predicate.tryAdd(w)) { // local, incremental check
            std::vector<int> &introduced = introducedBuf[U.size()];
            introduced.clear();
            for (int x : adj[w])
                if (x > s && !inU[x] && !inC[x]) {
                    addCandidate(x, w, wDist + 1);
                    introduced.push_back(x);
                }

            if (!rebuilding)
                (*report)(U); // a rebuilt node was reported by an earlier pass
            // Only descend on a pass.
            const Outcome below = extendFiltered(
                s, predicate, rebuilding ? resume : nullptr, record);
            predicate.undo(w); // reverse tryAdd(w) -- see class contract

            for (int x : introduced)
                removeCandidate(x);

            U.pop_back();
            inU[w] = false;
            relink(w, prev, next); // in place, not at the tail
            if (below != Outcome::completed) {
                if (below == Outcome::suspended)
                    recordHere(i, w); // this node is on the path
                restorePruned();
                return below;
            }
            continue;
        }

        U.pop_back();
        inU[w] = false;
        if (rebuilding) {
            // It passed when first added, so only an external stop refuses
            // it now.
            relink(w, prev, next);
            restorePruned();
            return Outcome::stopped;
        }
        if (record && predicate.budgetExhausted()) {
            // Suspend: w was not judged, so it is tried first next time.
            relink(w, prev, next);
            recordHere(i, w);
            record->deepest = level;
            restorePruned();
            return Outcome::suspended;
        }

        // tryAdd made no net change (transactional contract), so there is
        // nothing to undo. w is pruned at this node, and so throughout its
        // subtree: every set there containing w contains U + w, and the
        // predicate is anti-monotonic (see enumerateFiltered()). So keep it
        // out of the list for the rest of this loop -- no later sibling's
        // subtree scans or tries it again -- but in C, so none re-introduces
        // it. (Without a `record`, a budget refusal lands here too, which is
        // harmless: nothing after it in this pass is reported.)
        inC[w] = true;
        pruned.push_back({w, prev, next});
    }
    restorePruned();
    return Outcome::completed;
}
