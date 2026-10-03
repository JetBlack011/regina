/**
 * @file searchfrontier.h
 * @brief How far one search got, exactly: its breadth, and a place to resume.
 *
 * A search's traversal is a function of its search graph, its roots and its
 * schedule alone (EmbeddingSearch::runSearch_: each root's walk depends on
 * that root only, and budget passes carry on from a recorded
 * ConnectedInducedSubgraphEnumerator::Position). So where a search stopped
 * is fully described by: the iterative-deepening round it was in, and for
 * each root of that round, whether it finished, how many budget passes it
 * had, and where its last pass stopped. Rounds before it ran to the end.
 *
 * A SearchFrontier records exactly that, plus a fingerprint of everything
 * that fixes the traversal and what it accepts (see
 * EmbeddingSearch::frontierFingerprint_). Handed back to a search with the
 * same fingerprint, it resumes there: the two runs together report what one
 * uninterrupted run reports, nothing visited twice and nothing skipped
 * (submanifoldsearch_test's frontier tests). A stop suspends each in-flight
 * root at its exact position (BudgetedPredicate::suspendOnStop), so a
 * frontier never overclaims.
 *
 * In the same ORDER only when passes are unbudgeted. Within a round, a root
 * that spends its budget without finishing goes to the back of the root
 * queue (runSearch_'s worker, `rootQueue.push_back(idx)`), while a resumed
 * round rebuilds the queue in index order (runRound). So with root budgets
 * a resumed run visits a permutation of what one run visits: a search
 * stopped at a surface target and resumed reaches its next target through
 * other surfaces than one run would have (8_8 at the production shape,
 * 2026-10-02: 25,255 of 120,000 differ over six resumed steps), while
 * running it to exhaustion gives exactly one run's surfaces. Preserving the
 * order would need the queue in the frontier: a format change.
 *
 * It is also the search's breadth, as data: which round, how many roots
 * finished, how far the rest got, and the cumulative counts over every run
 * that continued it.
 *
 * Soundness for a caller that skips a prefix: the prefix's surfaces were
 * reported by the run that recorded it, so only a frontier written after
 * that run's finds were durable may be resumed. A caller that keeps its
 * finds in a file records it in the frontier (`pending`), with its durable
 * length, and checks it before resuming (cobound: `sign` must have signed
 * that far).
 */

#ifndef SEARCHFRONTIER_H
#define SEARCHFRONTIER_H

#include "surfer/enumeration/inducedsubgraphs.h"

#include <iosfwd>
#include <optional>
#include <string>
#include <vector>

struct SearchFrontier {
    /**
     * The traversal's version, part of every fingerprint. Bump it with any
     * change to which candidates the enumeration visits, in what order, or
     * what a Position means -- a frontier from another version must never
     * resume. embeddingsearch_test's test_traversal_pinned pins a small
     * search's visit order to this number, so an unbumped change fails it.
     */
    static constexpr int kTraversalVersion = 1;

    std::string fingerprint; /**< See EmbeddingSearch::frontierFingerprint_. */
    unsigned round = 1;      /**< The IDDFS round in progress, 1-based. */
    unsigned rounds = 1;     /**< Rounds in all: the capped ones plus the final. */
    std::optional<long long> cap; /**< That round's face cap; none if uncapped. */
    long long suppressBelow = 0;  /**< The previous round's cap (its finds are not re-reported). */
    std::optional<long long> deepestExhausted; /**< As SearchStats::deepestExhaustedCap. */
    bool complete = false;   /**< Every round ran to the end: nothing is left. */

    /** One root of the round in progress, by its index in the sorted root list. */
    struct Root {
        bool done = false;  /**< Enumerated to the end within the round's cap. */
        unsigned level = 0; /**< Budget passes had; the next pass's ration is start * growth^level. */
        ConnectedInducedSubgraphEnumerator::Position position; /**< Where its walk stopped. */
    };
    std::vector<Root> roots; /**< Every root of the round (all done when complete). */

    /**
     * Where the caller recorded this search's finds, when it keeps them in a
     * file (cobound: the search's pending file), and that file's fsynced
     * byte length when the frontier was taken. Format 2 writes it (save():
     * relative to the frontier's directory); a format-1 frontier has none.
     */
    struct Pending {
        std::string path;
        long long bytes = 0;
    };
    std::optional<Pending> pending;

    // Cumulative over the run that recorded this and every run it resumed.
    unsigned runs = 0;
    long long satisfying = 0; /**< Candidates satisfying the boundary condition and accepted. */
    long long found = 0;      /**< Candidates reported. */
    long long attempts = 0;   /**< Enumeration attempts (tryAdd() calls). */

    size_t rootsDone() const;
    /** Roots neither done nor untouched: part-walked. */
    size_t rootsStarted() const;
    /** The largest per-root pass count. */
    unsigned maxLevel() const;

    /**
     * One line: "round r/R cap c, roots d done + s part-walked of N, max
     * level L; runs n, satisfying S, found F, attempts A[, complete]".
     */
    std::string summary() const;

    void write(std::ostream &out) const;
    /** \throws std::runtime_error on anything but write()'s format. */
    static SearchFrontier read(std::istream &in);

    /**
     * Writes `path` atomically (a temporary beside it, then a rename). The
     * pending file is recorded relative to `path`'s directory (its own name
     * when no relative path exists), so a work tree that is packed, synced
     * or copied keeps its resume.
     */
    void save(const std::string &path) const;
    /**
     * Reads `path`; nullopt if it does not exist. A relative pending record
     * is resolved against `path`'s directory (an absolute one is kept), and
     * when that file does not exist but one of its name lies beside the
     * frontier, that one is named: `pending->path` is absolute either way.
     * \throws as read().
     */
    static std::optional<SearchFrontier> load(const std::string &path);
};

#endif // SEARCHFRONTIER_H
