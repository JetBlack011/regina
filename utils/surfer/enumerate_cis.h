#ifndef ENUMERATE_CIS_H

#define ENUMERATE_CIS_H

//
//  enumerate_cis.h
//
//  C++ implementation of a modified algorithm from:
//    M. Alokshiya, S. Salem, F. Abed,
//    "A linear delay algorithm for enumerating all connected induced
//    subgraphs", BMC Bioinformatics 20(Suppl 12):319, 2019.
//

#include <atomic>
#include <cstdint>
#include <functional>
#include <vector>

/*! \file utils/surfer/enumerate_cis.h
 *  \brief Enumerates connected induced subgraphs of a graph, with optional
 *  filtering and seeding.
 */

/**
 * A stateful, backtrackable predicate over the current vertex set U, for
 * checks that are naturally incremental (need only local information
 * about the vertex being added, not all of U).
 *
 * tryAdd(v) is called exactly when v becomes the newest member of U (i.e.
 * `U.back() == v`), and must be transactional: either it commits and
 * returns \c true, or it makes no net change and returns \c false. undo(v)
 * is called if and only if the matching tryAdd(v) returned \c true,
 * immediately after v's subtree has been fully explored, and must exactly
 * reverse it. Calls nest LIFO, exactly like U's own push_back()/pop_back().
 */
class ConditionalPredicate {
  public:
    /** Attempts to commit v's addition; see the class documentation. */
    virtual bool tryAdd(int v) = 0;

    /** Reverses a prior successful tryAdd(v); see the class documentation. */
    virtual void undo(int v) = 0;

    virtual ~ConditionalPredicate() = default;
};

/**
 * Decorates another ConditionalPredicate with an external stop signal:
 * once `*stopRequested` is set, every tryAdd() is rejected outright
 * without consulting \a inner, so a caller enumerating with this predicate
 * prunes everywhere and unwinds -- an early, cooperative exit from an
 * otherwise unbounded DFS.
 */
class InterruptiblePredicate : public ConditionalPredicate {
    ConditionalPredicate &inner_; /**< The predicate being decorated. */
    const std::atomic<bool> &stopRequested_;
        /**< External stop signal, checked on every tryAdd(). */

  public:
    /** Wraps `inner`, checking `stopRequested` on every tryAdd(). */
    InterruptiblePredicate(ConditionalPredicate &inner,
                           const std::atomic<bool> &stopRequested)
        : inner_(inner), stopRequested_(stopRequested) {}

    bool tryAdd(int v) override {
        if (stopRequested_.load(std::memory_order_relaxed))
            return false;
        return inner_.tryAdd(v);
    }

    // A successful tryAdd(v) always delegates to inner_, so undo(v) --
    // called only after such a tryAdd() -- can delegate unconditionally.
    void undo(int v) override { inner_.undo(v); }
};

/**
 * Enumerates every connected induced subgraph of a fixed graph, optionally
 * restricted to those satisfying a hereditary predicate, and optionally
 * seeded from a pre-selected connected subset.
 */
class ConnectedInducedSubgraphEnumerator {
  public:
    /**
     * A graph obtained by contracting a connected seed set of some
     * original graph into a single vertex, always given the
     * globally-smallest id (1).
     *
     * This convention is why seeding works: enumerateFromRoot(s, ...)
     * finds every connected set anchored at s (i.e. whose smallest-id
     * member is s), so with the seed at id 1, "anchored at 1" and
     * "contains the seed" coincide -- meaning enumerateFromRoot(1, ...)
     * on a SeededGraph finds exactly the connected supersets of the
     * original seed.
     */
    struct SeededGraph {
        int n; /**< The vertex count of the contracted graph. */
        std::vector<std::vector<int>> adj;
            /**< 1-indexed adjacency list; vertex 1 is the seed. */
        std::vector<int> originalOf;
            /**< originalOf[i], for i >= 2, is the original graph's id for
                 contracted vertex i. originalOf[1] is unused (-1), since
                 vertex 1 represents the whole seed, not a single original
                 vertex. */
    };

    /**
     * Contracts `seed` into a single vertex, producing a SeededGraph as
     * described above.
     *
     * \pre `seed` is a connected induced subgraph of (origN, origAdj)
     * (checked via assert).
     */
    static SeededGraph
    contractSeed(int origN, const std::vector<std::vector<int>> &origAdj,
                const std::vector<int> &seed);

    /**
     * Creates an unseeded enumerator over the graph (n, adj).
     */
    ConnectedInducedSubgraphEnumerator(int n,
                                       const std::vector<std::vector<int>> &adj)
        : n(n), adj(adj), candNext(n + 1, 0), candPrev(n + 1, 0),
          dist(n + 1, -1), parentOf(n + 1, 0), inU(n + 1, false),
          inC(n + 1, false), siblingBuf(n + 1), introducedBuf(n + 1) {
        initUnseededRoots_();
    }

    /**
     * Seeded variant: treats vertex 1 as a pre-contracted seed (see
     * SeededGraph) and fast-forwards to the state
     * enumerateFromRootFiltered(1, ...) would reach just before its first
     * descent: U = {1}, committed via `predicate.tryAdd(1)`, with roots_
     * set to every neighbor of 1 that individually passes
     * tryAdd()/undo() (tested and reverted, not committed) -- ready to
     * hand out as independent units of work.
     */
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

    /**
     * As above, but with no predicate: roots_ becomes every neighbor of
     * the seed, unfiltered.
     */
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

    /**
     * Returns the values it is valid to call enumerateFromRoot()/
     * enumerateFromRootFiltered() on: every graph vertex when unseeded, or
     * every seed-neighbor surviving the predicate (if one was given at
     * construction) when seeded.
     */
    const std::vector<int> &getRoots() const { return roots_; }

    /**
     * Calls `visit(U)` once for every non-empty connected induced
     * subgraph of the graph, with U given as a vertex list in discovery
     * order. Sequential entry point: walks every root in turn.
     */
    void enumerate(const std::function<void(const std::vector<int> &)> &visit) {
        for (int s : roots_)
            enumerateFromRoot(s, visit);
    }

    /**
     * Enumerates every connected induced subgraph anchored at vertex `s`
     * (i.e. `s` is the smallest-id vertex of every subgraph produced).
     *
     * Subtrees for different `s` are fully independent -- this is the
     * unit of parallelism: give each thread its own enumerator instance
     * over a shared, read-only adjacency list, and hand out values of `s`
     * to threads dynamically.
     */
    void enumerateFromRoot(
        int s, const std::function<void(const std::vector<int> &)> &visit);

    /**
     * As enumerate(), but only descends into (and reports) vertex sets
     * satisfying `predicate`, in addition to connectivity.
     *
     * \pre `predicate` is hereditary: for every connected U* satisfying
     * it, every connected subset of U* also satisfies it. Given that,
     * this finds every connected induced subgraph satisfying the
     * predicate, with no duplicates, visiting only the
     * predicate-satisfying nodes plus a thin "boundary" of their failing
     * children -- not the whole search space.
     *
     * \warning If `predicate` is not hereditary, this can silently miss
     * results.
     */
    void enumerateFiltered(
        const std::function<void(const std::vector<int> &)> &visit,
        const std::function<bool(const std::vector<int> &)> &predicate) {
        StatelessPredicate adapter(U, predicate);
        enumerateFiltered(visit, adapter);
    }

    /** Per-root counterpart to enumerateFiltered() above. */
    void enumerateFromRootFiltered(
        int s, const std::function<void(const std::vector<int> &)> &visit,
        const std::function<bool(const std::vector<int> &)> &predicate) {
        StatelessPredicate adapter(U, predicate);
        enumerateFromRootFiltered(s, visit, adapter);
    }

    /**
     * As enumerateFiltered() above, but taking a stateful
     * ConditionalPredicate instead of a stateless one -- see
     * ConditionalPredicate for the incremental-checking contract this
     * enables.
     */
    void enumerateFiltered(
        const std::function<void(const std::vector<int> &)> &visit,
        ConditionalPredicate &predicate) {
        for (int s : roots_)
            enumerateFromRootFiltered(s, visit, predicate);
    }

    /** Per-root counterpart to the ConditionalPredicate overload above. */
    void enumerateFromRootFiltered(
        int s, const std::function<void(const std::vector<int> &)> &visit,
        ConditionalPredicate &predicate);

  private:
    int n; /**< The number of vertices in the graph (numbered 1..n). */
    const std::vector<std::vector<int>> &adj;
        /**< The graph's adjacency list, 1-indexed. */

    std::vector<int> U; /**< The current connected vertex set, in insertion order. */

    /**
     * The candidate set C (extension frontier), stored as an intrusive
     * doubly linked list over preallocated arrays keyed by vertex id,
     * rather than std::list<int> -- candidates are added/removed once per
     * DFS node visited (astronomically many), and a node-based std::list
     * heap-allocates on every insert/erase. Index 0 is a sentinel:
     * candNext[0]/candPrev[0] are the list's head/tail, so an empty list
     * is `candNext[0] == 0`.
     */
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

    std::vector<int> dist;
        /**< dist[v]: distance from anchor(U); valid for v in U or C. */
    std::vector<int> parentOf;
        /**< parentOf[v]: the vertex of U that first discovered v. */

    /**
     * Membership flags for U and C, as plain byte flags rather than
     * std::vector<bool>: these are read twice per candidate neighbor in
     * the innermost loop of extend()/extendFiltered() (the busiest loop
     * in the whole search), and vector<bool>'s bit-packed storage turns
     * each read into a shift/mask/test instead of a single byte load --
     * measurable under profiling as one of the largest remaining
     * self-time contributors. The extra memory is negligible.
     */
    std::vector<uint8_t> inU, inC;

    const std::function<void(const std::vector<int> &)> *report = nullptr;
        /**< The visit callback for the enumeration currently in progress. */

    /**
     * Whether vertex 1 is a pre-contracted seed (see the seeded
     * constructors above). When set, enumerateFromRoot()/
     * enumerateFromRootFiltered() reinterpret their vertex argument as
     * "which candidate sibling of the seed to descend into", always
     * keeping vertex 1 as the true anchor.
     */
    bool isSeeded_ = false;

    std::vector<int> roots_; /**< See getRoots(). */

    /** Populates roots_ with every vertex 1..n (the unseeded case). */
    void initUnseededRoots_() {
        roots_.resize(n);
        for (int s = 1; s <= n; ++s)
            roots_[s - 1] = s;
    }

    /**
     * Fast-forwards state to treat vertex 1 as a seed already anchored in
     * U; see the seeded constructors above for the full contract.
     *
     * \a predicate may be \c null, in which case the seed is not
     * validated/committed against anything and roots_ becomes every
     * neighbor of 1, unfiltered.
     */
    void seedFastForward_(ConditionalPredicate *predicate);

    /**
     * Reusable per-depth scratch buffers for extend()/extendFiltered(),
     * indexed by the current recursion depth (`U.size()`) -- replaces two
     * heap allocations that used to happen on every recursive call.
     * Recursion is strictly sequential, so a given depth's buffer is
     * always safe to clear and reuse rather than reallocate.
     */
    std::vector<std::vector<int>> siblingBuf, introducedBuf;

    /**
     * Adapts a stateless whole-U predicate into the ConditionalPredicate
     * interface, so the stateless enumerateFiltered()/
     * enumerateFromRootFiltered() overloads above can forward into the
     * ConditionalPredicate ones instead of duplicating the traversal.
     */
    class StatelessPredicate : public ConditionalPredicate {
        const std::vector<int> &U;
        const std::function<bool(const std::vector<int> &)> &predicate;

      public:
        StatelessPredicate(
            const std::vector<int> &U_,
            const std::function<bool(const std::vector<int> &)> &predicate_)
            : U(U_), predicate(predicate_) {}

        // U already includes v by the time tryAdd(v) is called
        // (U.push_back happens before the check); undo is a no-op since
        // there's no persistent state of its own to roll back.
        bool tryAdd(int) override { return predicate(U); }
        void undo(int) override {}
    };

    void addCandidate(int v, int par, int d);

    /**
     * Fully discards v from candidate consideration: used for candidates
     * that never actually enter U (e.g. temporary "introduced" candidates
     * discarded again on backtrack).
     */
    void removeCandidate(int v);

    /**
     * Detaches v from the candidate list because it is being moved into
     * U. dist[v] is deliberately preserved (not reset): while v sits in U
     * it still needs to report its distance from anchor(U), since
     * utmost(U) is always `U.back()`, read straight out of dist[].
     */
    void detachToU(int v);

    /** Tries to extend U (with fixed anchor `s`) using every candidate in C. */
    void extend(int s);

    /**
     * As extend(), but only reports/descends when `predicate.tryAdd(w)`
     * passes -- see enumerateFiltered() for the correctness requirement
     * this relies on.
     */
    void extendFiltered(int s, ConditionalPredicate &predicate);
};

#endif // ENUMERATE_CIS_H
