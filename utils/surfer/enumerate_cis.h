#ifndef ENUMERATE_CIS_H

#define ENUMERATE_CIS_H

// enumerate_cis.h
//
// C++ implementation of a modified algorithm from:
//   M. Alokshiya, S. Salem, F. Abed,
//   "A linear delay algorithm for enumerating all connected induced subgraphs",
//   BMC Bioinformatics 20(Suppl 12):319, 2019.

#include <atomic>
#include <cstdint>
#include <functional>
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

// Decorates another ConditionalPredicate with an external stop signal:
// once *stopRequested is set, every tryAdd() is rejected outright (inner is
// never even consulted), so a caller enumerating with this predicate prunes
// everywhere and unwinds -- an early, cooperative exit from an otherwise
// unbounded DFS. Knows nothing about why the flag might be set (a signal
// handler, another thread, a timeout, ...) or what `inner` actually checks;
// it only adds interruptibility on top.
class InterruptiblePredicate : public ConditionalPredicate {
    ConditionalPredicate &inner_;
    const std::atomic<bool> &stopRequested_;

  public:
    InterruptiblePredicate(ConditionalPredicate &inner,
                           const std::atomic<bool> &stopRequested)
        : inner_(inner), stopRequested_(stopRequested) {}

    bool tryAdd(int v) override {
        if (stopRequested_.load(std::memory_order_relaxed))
            return false;
        return inner_.tryAdd(v);
    }

    // Only called when the matching tryAdd(v) returned true (per
    // ConditionalPredicate's contract), which only happens by delegating to
    // inner_ above -- so this can unconditionally delegate too.
    void undo(int v) override { inner_.undo(v); }
};

class ConnectedInducedSubgraphEnumerator {
  public:
    // A graph obtained by contracting a connected seed set of some original
    // graph into a single vertex, always given the globally-smallest id
    // (1) -- see contractSeed() and the seeded constructors below, which
    // rely on this convention: enumerateFromRoot(s, ...) finds every
    // connected set anchored at s (i.e. whose smallest-id member is s), so
    // if s is the smallest id in the whole graph, "anchored at s" and
    // "contains s" coincide, and enumerateFromRoot(1, ...) on a SeededGraph
    // finds exactly the connected supersets of the original seed.
    struct SeededGraph {
        int n; // vertex count of the contracted graph
        std::vector<std::vector<int>> adj; // 1-indexed; vertex 1 = seed
        // originalOf[i], for i >= 2: original graph id of contracted vertex
        // i. originalOf[1] is unused (-1): vertex 1 is the whole seed, not
        // a single original vertex.
        std::vector<int> originalOf;
    };

    // Contracts `seed` (which must be a connected induced subgraph of
    // (origN, origAdj) -- checked via assert) into a single vertex,
    // producing a SeededGraph as described above.
    static SeededGraph
    contractSeed(int origN, const std::vector<std::vector<int>> &origAdj,
                const std::vector<int> &seed);

    ConnectedInducedSubgraphEnumerator(int n,
                                       const std::vector<std::vector<int>> &adj)
        : n(n), adj(adj), candNext(n + 1, 0), candPrev(n + 1, 0),
          dist(n + 1, -1), parentOf(n + 1, 0), inU(n + 1, false),
          inC(n + 1, false), siblingBuf(n + 1), introducedBuf(n + 1) {
        initUnseededRoots_();
    }

    // Seeded variant: treats vertex 1 as a pre-contracted seed (see
    // contractSeed() above -- the seed always ends up as the
    // globally-smallest id, vertex 1, precisely so that "anchored at 1" and
    // "contains the seed" coincide) and fast-forwards this instance's state
    // to exactly where
    // enumerateFromRootFiltered(1, ...) would be right before its first call
    // to extendFiltered(): U = {1}, with the seed already committed via
    // predicate.tryAdd(1) and the candidate set populated with every
    // neighbor of 1. roots_ is then populated with every such neighbor that
    // individually passes predicate.tryAdd/undo (tested and immediately
    // reverted -- these are independent, mutually exclusive branches, not a
    // commitment), ready to be handed out as independent units of work (see
    // enumerateFromRoot/enumerateFromRootFiltered below, which are seed-aware
    // when isSeeded is set).
    ConnectedInducedSubgraphEnumerator(int n,
                                       const std::vector<std::vector<int>> &adj,
                                       bool isSeeded,
                                       ConditionalPredicate &predicate)
        : n(n), adj(adj), candNext(n + 1, 0), candPrev(n + 1, 0),
          dist(n + 1, -1), parentOf(n + 1, 0), inU(n + 1, false),
          inC(n + 1, false), siblingBuf(n + 1), introducedBuf(n + 1),
          isSeeded_(isSeeded) {
        if (isSeeded_)
            seedFastForward_(&predicate);
        else
            initUnseededRoots_();
    }

    // Same, but with no predicate: roots_ is every neighbor of the seed,
    // unfiltered.
    ConnectedInducedSubgraphEnumerator(int n,
                                       const std::vector<std::vector<int>> &adj,
                                       bool isSeeded)
        : n(n), adj(adj), candNext(n + 1, 0), candPrev(n + 1, 0),
          dist(n + 1, -1), parentOf(n + 1, 0), inU(n + 1, false),
          inC(n + 1, false), siblingBuf(n + 1), introducedBuf(n + 1),
          isSeeded_(isSeeded) {
        if (isSeeded_)
            seedFastForward_(nullptr);
        else
            initUnseededRoots_();
    }

    // Elements of roots_ are exactly the values it is valid to call
    // enumerateFromRoot()/enumerateFromRootFiltered() on: every graph vertex
    // when unseeded, or every seed-neighbor surviving the predicate (if one
    // was given at construction) when seeded.
    const std::vector<int> &getRoots() const { return roots_; }

    // Calls visit(U) once for every non-empty connected induced subgraph of
    // the graph, where U is given as a vertex list (in discovery order).
    // Sequential entry point -- just walks every root.
    void enumerate(const std::function<void(const std::vector<int> &)> &visit) {
        for (int s : roots_)
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
        for (int s : roots_)
            enumerateFromRootFiltered(s, visit, predicate);
    }

    void enumerateFromRootFiltered(
        int s, const std::function<void(const std::vector<int> &)> &visit,
        ConditionalPredicate &predicate);

  private:
    int n;
    const std::vector<std::vector<int>> &adj;

    std::vector<int> U; // current connected vertex set (insertion order)

    // Candidate set C (extension frontier), stored as an intrusive doubly
    // linked list over preallocated arrays keyed by vertex id, rather than
    // std::list<int>: candidates are added/removed extremely often (once per
    // DFS node visited, of which there can be astronomically many), and a
    // node-based std::list heap-allocates/frees on every single insert/erase.
    // Index 0 is a sentinel: candNext[0]/candPrev[0] are the list's
    // head/tail, so an empty list is candNext[0] == 0.
    std::vector<int> candNext, candPrev;
    void listPushBack(int v) {
        int tail = candPrev[0];
        candNext[tail] = v;
        candPrev[v] = tail;
        candNext[v] = 0;
        candPrev[0] = v;
    }
    void listErase(int v) {
        int p = candPrev[v], nx = candNext[v];
        candNext[p] = nx;
        candPrev[nx] = p;
    }
    int listFront() const { return candNext[0]; }
    bool listEmpty() const { return candNext[0] == 0; }

    std::vector<int>
        dist; // dist[v]: distance from anchor(U); valid for v in U or C
    std::vector<int>
        parentOf; // parentOf[v]: vertex of U that first discovered v
    // Plain byte flags, not std::vector<bool>: these are read twice per
    // candidate neighbor in extend()/extendFiltered()'s innermost loop (the
    // busiest loop in the whole search), and std::vector<bool>'s bit-packed
    // storage turns each read into a shift/mask/test sequence instead of a
    // single byte load -- measurable under profiling as one of the largest
    // remaining self-time contributors. The extra memory (bytes instead of
    // bits, per graph node) is negligible.
    std::vector<uint8_t> inU, inC;
    const std::function<void(const std::vector<int> &)> *report = nullptr;

    // See the seeded constructors above. When isSeeded_ is set, vertex 1 is
    // a pre-contracted seed, U/dist/inU/candidates are already
    // fast-forwarded to reflect it, and enumerateFromRoot/
    // enumerateFromRootFiltered reinterpret their vertex argument as "which
    // already-populated candidate (sibling of the seed) to descend into",
    // always keeping vertex 1 -- not the sibling -- as the true anchor
    // passed to extend()/extendFiltered(). When false, everything behaves
    // exactly as it always has.
    bool isSeeded_ = false;

    // Elements are exactly the values valid to call enumerateFromRoot()/
    // enumerateFromRootFiltered() on -- see getRoots() above.
    std::vector<int> roots_;

    // Populates roots_ with every vertex 1..n (the unseeded case).
    void initUnseededRoots_() {
        roots_.resize(n);
        for (int s = 1; s <= n; ++s)
            roots_[s - 1] = s;
    }

    // Fast-forwards state to treat vertex 1 as a seed already anchored in U
    // (see the seeded constructors above for the full contract). predicate
    // may be nullptr, in which case the seed is not validated/committed
    // against anything and roots_ is every neighbor of 1, unfiltered.
    void seedFastForward_(ConditionalPredicate *predicate);

    // Reusable per-depth scratch buffers for extend()/extendFiltered(),
    // indexed by the current recursion depth (U.size()) -- replaces two
    // heap allocations that used to happen on every single recursive call
    // (the "siblings" snapshot of C, and the "introduced" list per
    // candidate). Depth never repeats concurrently within one enumerator
    // instance (recursion is strictly sequential and depth strictly
    // increases on descent/decreases on backtrack), so each depth's buffer
    // is safe to clear and reuse across calls instead of reallocating.
    std::vector<std::vector<int>> siblingBuf, introducedBuf;

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
