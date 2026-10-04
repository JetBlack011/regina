//
//  censusnaming.h
//
//  Naming a complement by the census.
//

#ifndef SURFER_LINKNAMING_CENSUSNAMING_H
#define SURFER_LINKNAMING_CENSUSNAMING_H

#include <atomic>
#include <utility>
#include <chrono>
#include <mutex>
#include <optional>
#include <string>

#include <triangulation/dim3.h>

#include "linknaming/complement/linkcomplement.h"

/*! \file utils/surfer/linknaming/census/censusnaming.h
 *  \brief Names an EdgeComplement/Link by its complement: the unknot and
 *  unlinks (unlinknaming.h), then the local census, regina::Census::lookup()
 *  and a bounded Pachner search. Answers are cached in complementcache.h.
 *
 *  See linkcomplement.h for the pure edge-set/complement representation
 *  this operates on.
 */

namespace census {

/**
 * A mutex guarding only regina::Census::lookup() itself: Regina's census
 * databases are not safe to query from more than one thread at a time --
 * concurrent calls have been observed to crash outright (Tokyo Cabinet
 * reports a "threading error"), not just contend. Building or simplifying a
 * complement is unaffected, and stays parallel across callers. The
 * memoized cache of complement answers is guarded separately, so a cache
 * hit never blocks behind an in-flight lookup on another thread.
 */
extern std::mutex censusLookupMutex;

/**
 * Names `e` by its complement, without printing: the census name if there
 * is one, "Unknot" if a genus-1 handlebody, or else the bare isoSig as a
 * fallback identifier.
 */
std::string nameComplement(const EdgeComplement &e);

/**
 * Test-only. When set, nameComplement() appends a fresh suffix to every
 * name it returns, so no two names it gives ever agree -- exactly the
 * non-determinism a census entry number already has, pushed to the limit.
 * Which surfaces a search accepts must not change (tests/
 * name_independence_test.sh): names are recorded, never gated on.
 * cobound sets it from the SURFER_TEST_PERTURB_NAMES environment variable.
 */
extern std::atomic<bool> perturbNamesForTesting;

/** `name`, perturbed as nameComplement() perturbs its own when
 *  perturbNamesForTesting is set; for other namers (outgoing::OutgoingNamer)
 *  to honour the same test. */
std::string perturbedForTesting(std::string name);

/**
 * Names `l` by its complement, without printing (all of `l`'s components
 * drilled together, i.e. Link::buildComplement()): `"<n>-component
 * unlink"` if `l`'s complement is proven split (see
 * unlinknaming.h's groupProvesUnlink() -- a fast, sound
 * generalization of nameComplement(const EdgeComplement&)'s genus-1/"Unknot"
 * check to n > 1 components, via a free-group check rather than a
 * handlebody genus check: unlike a single unknotted curve, a split
 * multi-component unlink's complement has multiple torus boundary
 * components, so it is never itself a handlebody), the census name if
 * there is one, "Unknot" if `l` is a single unknotted curve, or else the
 * bare isoSig as a fallback identifier.
 *
 * A genuine overload, not virtual dispatch: since Link publicly inherits
 * EdgeComplement, nameComplement(const EdgeComplement&) would already run and
 * build the right (possibly multi-cusp) complement if called with a Link,
 * and (for n == 1) already returns exactly the same answer this does. This
 * overload is picked automatically by ordinary overload resolution
 * whenever the argument is a Link -- it only changes behavior for n > 1,
 * by trying the split-unlink fast path first.
 */
std::string nameComplement(const Link &l);

/**
 * Prints and returns whether `e`'s complement is named: either as a
 * genus-1 handlebody (the complement of a single unknotted component), or
 * as a census hit.
 */
bool reportComplement(const EdgeComplement &e);

/** Prints whether each component of `l`'s complement is named; see reportComplement(const EdgeComplement&). */
void reportComplement(const Link &l);

} // namespace census

namespace census {

/**
 * Whether resolveAnswer() should attempt retriangulateAndLookup() on
 * a local-census/real-Census::lookup() miss. Off by default -- materially
 * more expensive than a plain lookup, so existing callers (surfer.cpp)
 * must opt in via --retriangulate-on-miss; cobound's run without a goal
 * defaults it on (retriangulate_on_miss), since resolving non-hyperbolic
 * outgoing ends is core to its purpose.
 *
 * Set once, before any search worker thread is spawned, same contract as
 * linkcomplement.h's simplifyComplements.
 */
extern std::atomic<bool> retriangulateOnMiss;

/**
 * Whether retriangulateOnMiss also applies to multi-component (link)
 * complements. Off by default: a link's name, read from its complement
 * alone, never bears a slice-genus bound during a search, and the atlas's
 * per-cobordism naming names every outgoing link afterwards, so the
 * search spending up to retriangulateTimeBudgetSeconds per link complement
 * on it bought nothing. Same set-once contract as retriangulateOnMiss.
 */
extern std::atomic<bool> retriangulateLinks;

/**
 * Parameters resolveAnswer() passes to retriangulateAndLookup() on
 * every retriangulateOnMiss attempt -- retriangulateAndLookup()'s own
 * defaults (height=2, candidateBudget=8000, timeBudget=20s) exist for
 * direct/test callers only; production naming always goes through
 * these, so they're the actual knobs for trading match rate against
 * per-curve worst-case cost. retriangulate()'s candidate count grows
 * roughly exponentially in height, so height=1 (rather than the default
 * 2) is usually the single biggest lever if boundary-processing throughput
 * matters more than catching every non-canonically-triangulated match.
 * Same set-once-before-any-worker-thread-spawns contract as
 * retriangulateOnMiss.
 */
extern std::atomic<int> retriangulateHeight;
extern std::atomic<size_t> retriangulateCandidateBudget;
extern std::atomic<long long> retriangulateTimeBudgetSeconds;

/**
 * Looks up `sig` (an isoSig, e.g. from EdgeComplement::buildComplement()'s
 * result) in the local SQLite census -- a copy of the 3 census databases
 * that can ever match a cusped boundary complement, plus any
 * SnapPy-named isoSigs Regina's own census misses, plus anything
 * inserted directly via insertCensusEntry() (cobound: a knot target's
 * complement, and Pachner hits).
 * Unlike the real regina::Census::lookup(), needs no mutex: each thread
 * lazily opens its own read-only connection, and SQLite supports many
 * concurrent readers natively.
 *
 * Returns nullopt on a miss (isoSig not in the census, or the file isn't
 * present at all -- e.g. a fresh checkout that hasn't run the generator
 * yet); callers should fall through to the real Census::lookup() in that
 * case, exactly as if the file didn't exist.
 */
std::optional<std::string> localCensusLookup(const std::string &sig);

/**
 * Sets the path localCensusLookup() opens. Every thread's cached
 * connection is invalidated and lazily reopened against the new path on
 * its next lookup. Returns whether a file currently exists at `path`, so
 * callers can report the census as loaded/skipped without a second
 * existence check.
 *
 * Set once, before any search worker thread is spawned (see surfer.cpp's
 * --census-db flag, defaulting to the SURFER_CENSUS_PATH compile
 * definition from CMakeLists.txt) -- same contract as
 * linkcomplement.h's simplifyComplements/complement::cacheLimit.
 */
bool setCensusPath(const std::string &path);

/**
 * Points localCensusLookup() at a path guaranteed not to exist, so it
 * always misses. Test-only: production code should only ever call
 * setCensusPath() once, at startup.
 */
void resetCensusForTesting();

/**
 * Whether the local census is ever written (the plan's census_updates): off,
 * insertCensusEntry() writes nothing, whoever calls it -- a search's
 * end-of-search insert of its own knot's complement and the Pachner-hit
 * insert inside naming alike. On by default; a goal-directed run switches it
 * off. Set once, before any search worker thread is spawned, like
 * retriangulateOnMiss.
 */
extern std::atomic<bool> censusUpdates;

/**
 * Inserts (isoSig -> name) into the local census at the current path
 * (creating the file/table if missing), tagged with `source` for
 * provenance. INSERT OR IGNORE: never overwrites an existing entry. Bumps
 * the same generation counter setCensusPath() does, so other threads'
 * cached read connections pick up the change on their next lookup.
 * Thread-safe (its own mutex-guarded read-write connection, opened
 * lazily, WAL + busy_timeout enabled for robustness against another
 * concurrently-running cobound instance sharing the same
 * file). No-op (returns false) if the path can't be opened for writing, and
 * with censusUpdates off (counted neither ok nor failed: nothing was tried).
 */
bool insertCensusEntry(const std::string &isoSig, const std::string &name,
                       const std::string &source = "verifyslicegenus");

/**
 * How many insertCensusEntry() calls have succeeded and failed in this
 * process, (ok, failed). An INSERT OR IGNORE of an existing key counts as
 * ok. Production callers ignore the return value, so this is the only
 * place a census that cannot be written to shows up.
 */
std::pair<long long, long long> insertCounts();

/**
 * When a direct census lookup on `complement` misses, searches up to
 * `height` Pachner moves' worth of alternate triangulations of the same
 * manifold (regina::Triangulation<3>::retriangulate(), single-threaded to
 * avoid oversubscribing an already-multithreaded search) for one whose
 * isoSig hits the census, stopping at the first match or once
 * `candidateBudget`/`timeBudget` is exhausted.
 *
 * This exists because neither the local census nor the real
 * Census::lookup() can match two different triangulations of the same
 * NON-hyperbolic manifold if they happened to simplify() to
 * non-isomorphic results (unlike hyperbolic manifolds, where SnapPea's
 * canonical cell decomposition guarantees the same triangulation
 * regardless of construction path) -- this performs a bounded
 * canonicalization-by-search instead.
 *
 * On a hit, additionally inserts (complement's ORIGINAL isoSig -> matched
 * name) into the census, self-reinforcing so the identical non-canonical
 * triangulation resolves directly next time without re-searching.
 *
 * \return the matched name, or nullopt if the budget is exhausted with no
 * match.
 */
std::optional<std::string> retriangulateAndLookup(
    const regina::Triangulation<3> &complement, int height = 2,
    size_t candidateBudget = 8000,
    std::chrono::seconds timeBudget = std::chrono::seconds(20));

} // namespace census

#endif // SURFER_LINKNAMING_CENSUSNAMING_H
