#ifndef ENUMERATE_CIS_H

#define ENUMERATE_CIS_H

// enumerate_cis.cpp
//
// C++ implementation of a modified algorithm from:
//   M. Alokshiya, S. Salem, F. Abed,
//   "A linear delay algorithm for enumerating all connected induced subgraphs",
//   BMC Bioinformatics 20(Suppl 12):319, 2019.

#include <functional>
#include <list>
#include <vector>

// A stateful, backtrackable predicate -- for use when checking whether a
// candidate vertex can be added is naturally incremental (only needs local
// information about the vertex being added, not the whole current set).
//
// Contract:
//  - tryAdd(v) is called exactly when v has just become the newest member
//    of U (i.e. U.back() == v). It must be TRANSACTIONAL: either it commits
//    whatever internal state change corresponds to v's addition and
//    returns true, or it makes NO net change to internal state and returns
//    false. If checking v requires several elementary sub-steps and one of
//    them fails partway through, tryAdd must roll back the sub-steps it
//    already performed within this same call before returning false.
//  - undo(v) is called if and only if the matching tryAdd(v) returned true,
//    immediately after the recursive exploration of v's subtree completes.
//    It must exactly reverse the changes tryAdd(v) made.
//  - Calls nest LIFO, exactly like U's own push_back/pop_back.
class ConditionalPredicate {
  public:
    virtual bool tryAdd(int v) = 0;
    virtual void undo(int v) = 0;
    virtual ~ConditionalPredicate() = default;
};

class ConnectedInducedSubgraphEnumerator {
  public:
    ConnectedInducedSubgraphEnumerator(int n,
                                       const std::vector<std::vector<int>> &adj)
        : n(n), adj(adj), dist(n + 1, -1), parentOf(n + 1, 0),
          inU(n + 1, false), inC(n + 1, false), itOf(n + 1) {}

    // Calls visit(U) once for every non-empty connected induced subgraph of
    // the graph, where U is given as a vertex list (in discovery order).
    // Sequential entry point -- just walks every root.
    void enumerate(const std::function<void(const std::vector<int> &)> &visit) {
        for (int s = 1; s <= n; ++s)
            enumerateFromRoot(s, visit);
    }

    // Enumerates every connected induced subgraph anchored at vertex s (i.e.
    // s is the smallest-id vertex of every subgraph produced). Subtrees for
    // different s are fully independent -- this is the unit of parallelism:
    // give each thread its own enumerator instance (own dist/parentOf/inU/
    // inC/U/C state) over a shared, read-only adjacency list, and hand out
    // values of s to threads dynamically.
    void enumerateFromRoot(
        int s, const std::function<void(const std::vector<int> &)> &visit);

    // --- Filtered variant -----------------------------------------------
    //
    // In addition to connectivity, only descends into (and reports) vertex
    // sets satisfying a predicate.
    //
    // CORRECTNESS REQUIREMENT: the predicate must be *hereditary* -- for
    // every connected U* satisfying it, every connected subset of U* must
    // also satisfy it. Given that, these find every connected induced
    // subgraph satisfying the predicate, with no duplicates, while visiting
    // only the predicate-satisfying nodes plus a thin "boundary" of their
    // failing children -- not the whole search space. If the predicate is
    // NOT hereditary, this can silently miss results.
    void enumerateFiltered(
        const std::function<void(const std::vector<int> &)> &visit,
        const std::function<bool(const std::vector<int> &)> &predicate) {
        StatelessPredicate adapter(U, predicate);
        enumerateFiltered(visit, adapter);
    }

    void enumerateFromRootFiltered(
        int s, const std::function<void(const std::vector<int> &)> &visit,
        const std::function<bool(const std::vector<int> &)> &predicate) {
        StatelessPredicate adapter(U, predicate);
        enumerateFromRootFiltered(s, visit, adapter);
    }

    void enumerateFiltered(
        const std::function<void(const std::vector<int> &)> &visit,
        ConditionalPredicate &predicate) {
        for (int s = 1; s <= n; ++s)
            enumerateFromRootFiltered(s, visit, predicate);
    }

    void enumerateFromRootFiltered(
        int s, const std::function<void(const std::vector<int> &)> &visit,
        ConditionalPredicate &predicate);

  private:
    int n;
    const std::vector<std::vector<int>> &adj;

    std::vector<int> U; // current connected vertex set (insertion order)
    std::list<int> C;   // current candidate vertices (extension frontier)
    std::vector<int>
        dist; // dist[v]: distance from anchor(U); valid for v in U or C
    std::vector<int>
        parentOf; // parentOf[v]: vertex of U that first discovered v
    std::vector<bool> inU, inC;
    std::vector<std::list<int>::iterator> itOf;
    const std::function<void(const std::vector<int> &)> *report = nullptr;

    // Adapts a stateless whole-U predicate into the ConditionalPredicate
    // interface, so the stateless
    // enumerateFiltered()/enumerateFromRootFiltered() overloads above can just
    // forward into the ConditionalPredicate ones instead of duplicating the
    // traversal. tryAdd re-evaluates the predicate over the current U (which,
    // by the time tryAdd(v) is called, already includes v -- U.push_back
    // happens before the check); undo is a no-op since there's no persistent
    // state of its own to roll back.
    class StatelessPredicate : public ConditionalPredicate {
        const std::vector<int> &U;
        const std::function<bool(const std::vector<int> &)> &predicate;

      public:
        StatelessPredicate(
            const std::vector<int> &U_,
            const std::function<bool(const std::vector<int> &)> &predicate_)
            : U(U_), predicate(predicate_) {}
        bool tryAdd(int) override { return predicate(U); }
        void undo(int) override {}
    };

    void addCandidate(int v, int par, int d);

    // Fully discard v from candidate consideration: used for candidates that
    // never actually enter U (e.g. temporary "introduced" candidates that
    // get discarded again on backtrack).
    void removeCandidate(int v);

    // Detach v from the candidate list because it is being moved into U.
    // dist[v] is deliberately preserved (NOT reset): while v sits in U it
    // still needs to report its distance from anchor(U), since utmost(U) is
    // always U.back() and its distance is read straight out of dist[].
    void detachToU(int v);

    // Try to extend the current connected set U (with fixed anchor s) using
    // every vertex currently in the candidate set C.
    void extend(int s);

    // Traversal for the filtered variant: only reports/descends when
    // predicate.tryAdd(w) passes -- see enumerateFiltered() above for the
    // correctness requirement this relies on. Used directly by the
    // ConditionalPredicate overloads, and indirectly by the stateless
    // overloads via StatelessPredicate.
    void extendFiltered(int s, ConditionalPredicate &predicate);
};

#endif // ENUMERATE_CIS_H
