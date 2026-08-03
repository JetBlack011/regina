//
//  identifycomplement.h
//
//  Created by John Teague on 07/16/2026.
//

#ifndef IDENTIFYCOMPLEMENT_H

#define IDENTIFYCOMPLEMENT_H

#include <atomic>
#include <chrono>
#include <functional>
#include <mutex>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include <triangulation/dim3.h>

#include "linkcomplement.h"

/*! \file utils/surfer/identifycomplement.h
 *  \brief Recognizes and names an EdgeComplement/Link's complement: genus
 *  checks, the local census, and regina::Census::lookup().
 *
 *  See linkcomplement.h for the pure edge-set/complement representation
 *  this operates on.
 */

namespace identify {

/**
 * The outcome of recognizing a complement, memoized by isomorphism
 * signature (see identifycomplement.cpp's recognitionCache).
 */
struct RecognitionResult {
    /**
     * Triangulation<3>::recogniseHandlebody()'s result, if it has been
     * computed for this isoSig: 1 means the unknot (a genus-1 handlebody),
     * -1 means "not any handlebody" (the only case eligible for a census
     * lookup), and any other value is the pre-existing "almost definitely a
     * bug" case (see recognizeComplement(const Link&)). nullopt means not
     * yet computed.
     */
    std::optional<ssize_t> genus;

    /**
     * Whether Census::lookup() has been run for this isoSig. Only ever
     * attempted when genus == -1.
     */
    bool censusChecked = false;

    /** The census hit's name, if censusChecked and a hit was found. */
    std::optional<std::string> censusName;

    /**
     * Whether census::retriangulateAndLookup() has already been tried for
     * this isoSig and failed to find a match -- so a recurring bare isoSig
     * within one run doesn't repeatedly burn its search budget. Only
     * meaningful once censusChecked is true and censusName is still
     * nullopt.
     */
    bool retriangulateAttempted = false;
};

/**
 * Counters for how much recomputation the recognition cache is actually
 * avoiding. See recognitionCacheStats().
 */
struct RecognitionCacheStats {
    long long genusChecks = 0;
    long long genusCacheHits = 0;
    long long censusChecks = 0;
    long long censusCacheHits = 0;

    /**
     * On a genus cache miss, how the genus was actually resolved: via the
     * fast fundamental-group-is-Z check (proves genus == 1), the fast
     * SnapPea hyperbolicity check (proves genus == -1), or the full
     * recogniseHandlebody() fallback (neither fast check applied, e.g. a
     * torus/satellite knot). These three always sum to
     * genusChecks - genusCacheHits.
     */
    long long groupFastPathHits = 0;
    long long snapPeaFastPathHits = 0;
    long long recogniseHandlebodyFallbacks = 0;

    /**
     * How many times the local census (census::localCensusLookup()) was
     * checked on a census cache miss, and how many of those were hits --
     * a local hit skips the mutex-guarded, on-disk-database-reopening
     * real Census::lookup() entirely. Always <= censusChecks - censusCacheHits.
     */
    long long localCensusChecks = 0;
    long long localCensusHits = 0;

    /** How many times recognitionCache has been fully cleared after exceeding recognitionCacheLimit. */
    long long cacheResets = 0;
};

/** A snapshot of the recognition cache's current hit/miss counters. */
RecognitionCacheStats recognitionCacheStats();

/** The recognition cache's current entry count (distinct isoSigs seen). */
size_t recognitionCacheSize();

/**
 * Clears recognitionCache and its stats outright, ignoring
 * recognitionCacheLimit. Test-only: production code should only ever see
 * this cache clear itself automatically via recognitionCacheLimit.
 */
void resetRecognitionCacheForTesting();

/**
 * Entry-count threshold past which recognitionCache (identifycomplement.cpp)
 * clears itself entirely before admitting the next new isoSig -- keyed by
 * content (an isoSig string), not a recyclable integer id, so a lookup
 * racing a clear is simply a clean miss, never a wrong hit against an
 * unrelated key (unlike PetalCache, no epoch-tagging is needed here).
 *
 * Set once, before any search worker thread is spawned, same contract as
 * linkcomplement.h's simplifyComplements.
 */
extern std::atomic<size_t> recognitionCacheLimit;

/**
 * A mutex guarding only regina::Census::lookup() itself: Regina's census
 * databases are not safe to query from more than one thread at a time --
 * concurrent calls have been observed to crash outright (Tokyo Cabinet
 * reports a "threading error"), not just contend. Building or simplifying a
 * complement is unaffected, and stays parallel across callers. The
 * memoized cache of recognition results is guarded separately, so a cache
 * hit never blocks behind an in-flight lookup on another thread.
 */
extern std::mutex censusLookupMutex;

/**
 * Non-printing identification of `e`'s complement: the census name if
 * recognized, "Unknot" if a genus-1 handlebody, or else the bare isoSig as
 * a fallback identifier.
 */
std::string identify(const EdgeComplement &e);

/**
 * Non-printing identification of `l`'s complement (all of `l`'s components
 * drilled together, i.e. Link::buildComplement()): `"<n>-component
 * unlink"` if `l`'s complement is proven split (see
 * identifycomplement.cpp's groupProvesUnlink() -- a fast, sound
 * generalization of identify(const EdgeComplement&)'s genus-1/"Unknot"
 * check to n > 1 components, via free-group recognition rather than a
 * handlebody genus check: unlike a single unknotted curve, a split
 * multi-component unlink's complement has multiple torus boundary
 * components, so it is never itself a handlebody), the census name if
 * recognized, "Unknot" if `l` is a single unknotted curve, or else the
 * bare isoSig as a fallback identifier.
 *
 * A genuine overload, not virtual dispatch: since Link publicly inherits
 * EdgeComplement, identify(const EdgeComplement&) would already run and
 * build the right (possibly multi-cusp) complement if called with a Link,
 * and (for n == 1) already returns exactly the same answer this does. This
 * overload is picked automatically by ordinary overload resolution
 * whenever the argument is a Link -- it only changes behavior for n > 1,
 * by trying the split-unlink fast path first.
 */
std::string identify(const Link &l);

/**
 * Whether `name` (as returned by identify(const Link&)) is safe to use for
 * a genus deduction regardless of which of the several oriented variants
 * of a link the drilled curves actually were: true for `"Unknot"` and any
 * `"<n>-component unlink"` result (a split unlink's components don't
 * interact, so each bounds its own disk in B^4 regardless of orientation),
 * false for any other (genuinely linked) multi-component result -- two
 * differently-oriented variants of the same link can have very different
 * true slice genus while sharing one complement, so a bare complement name
 * alone never determines which oriented link was actually found. Single-
 * curve names (from identify(const EdgeComplement&)) need no such check:
 * a single component has no orientation ambiguity to begin with.
 *
 * Kept here, next to where these strings are actually produced, rather
 * than pattern-matched elsewhere, so the two stay in sync if the format
 * ever changes.
 */
bool isOrientationSafeName(const std::string &name);

/**
 * Prints and returns whether `e`'s complement is recognized: either as a
 * genus-1 handlebody (the complement of a single unknotted component), or
 * as a census hit.
 */
bool recognizeComplement(const EdgeComplement &e);

/**
 * Cheap test for whether `e`'s complement is a genus-1 handlebody (a solid
 * torus) -- unlike recognizeComplement()/identify(), this never falls back
 * to the slower Census::lookup().
 *
 * \return \c true if and only if the complement is a genus-1 handlebody.
 */
bool isUnknot(const EdgeComplement &e);

/** Prints whether each component of `l`'s complement is recognized; see recognizeComplement(const EdgeComplement&). */
void recognizeComplement(const Link &l);

/** Counters for how much recomputation BoundarySignatureCache is actually avoiding. */
struct BoundarySignatureCacheStats {
    long long checks = 0;
    long long hits = 0;
    long long cacheResets = 0; /**< How many times this cache has been fully cleared after exceeding its clear threshold. */
};

/**
 * Memoizes identify()'s results for marked edge sets in one fixed ambient
 * boundary-component triangulation, canonicalized against that
 * triangulation's own automorphism group -- so a boundary curve already
 * seen (exactly, or up to a symmetry of the boundary component)
 * short-circuits before buildComplement()/simplify()/isoSig() ever runs,
 * rather than only being caught by the isoSig-keyed recognitionCache
 * *after* paying for that triangulation (see identifycomplement.cpp).
 *
 * recognitionCache is still needed on top of this: two curves can be the
 * same knot/link type without being combinatorially related by any
 * automorphism of the boundary component (e.g. two non-isomorphic edge
 * paths that happen to drill out to the same manifold), in which case only
 * recognitionCache -- keyed on the post-simplify isoSig -- catches the
 * duplicate. This cache is a cheaper pre-filter in front of that one, not a
 * replacement for it.
 *
 * The ambient boundary-component triangulation this tracks must stay fixed
 * for this cache's lifetime, and must be the exact triangulation
 * KnottedSurface::boundaryLinks() draws its edges from (i.e. built the same
 * way, via BoundaryComponent::build()) -- see SurfaceSearch, which owns one
 * instance per boundary component, shared across every worker thread's
 * KnottedSurface, the same way it shares one PetalCache.
 *
 * Thread-safe: the (lazily built, then read-only) automorphism group is
 * guarded by a std::once_flag, and the memoization table by its own mutex
 * -- so a cache hit never blocks behind another thread still computing the
 * group, and computing the group happens exactly once regardless of how
 * many threads race to trigger it.
 */
class BoundarySignatureCache {
  public:
    /** Default entry-count threshold; see the constructor's `clearThreshold` parameter. */
    static constexpr size_t DEFAULT_CLEAR_THRESHOLD = 200'000;

    /**
     * Tracks `boundary`, which must outlive this cache. Once this cache
     * holds `clearThreshold` distinct canonical signatures, it clears
     * itself entirely before admitting the next new one -- unlike
     * PetalCache, no epoch-tagging is needed here, since a cache_ entry is
     * keyed by its canonical signature (content), not a recyclable integer
     * id: a lookup racing a clear simply sees a clean miss, never a wrong
     * hit against an unrelated key.
     */
    explicit BoundarySignatureCache(
        const regina::Triangulation<3> &boundary,
        size_t clearThreshold = DEFAULT_CLEAR_THRESHOLD);

    /**
     * Returns the identify() result for `edgeIndices` -- this boundary
     * component's marked edges, sorted, e.g. from
     * EdgeComplement::edgeIndices() -- computing it via `compute` on a
     * cache miss and memoizing the result. `compute` is only invoked on a
     * miss; typically `[&]{ return identify::identify(curve); }`.
     */
    std::string identifyCached(const std::vector<size_t> &edgeIndices,
                                const std::function<std::string()> &compute);

    /**
     * Returns a snapshot of this cache's current hit/miss counters; see
     * BoundarySignatureCacheStats. Returned by value (not by reference)
     * since another thread may be concurrently updating the live counters.
     */
    BoundarySignatureCacheStats stats() const;

    /** This cache's current entry count (distinct canonical signatures seen). */
    size_t size() const;

  private:
    const regina::Triangulation<3> *boundary_;

    std::once_flag groupOnce_;
    std::vector<std::vector<size_t>> automorphismEdgeImages_;
        /**< Lazily built by ensureAutomorphismGroup_():
             automorphismEdgeImages_[g][e] is the image of edge `e` of
             *boundary_ under the g-th automorphism found. Read-only once
             built, so safe to read from any thread without further
             synchronization (std::once_flag already establishes
             happens-before for every caller, not just the one that ran
             the initializer). */

    size_t clearThreshold_;

    mutable std::mutex cacheMutex_; /**< Guards cache_ and stats_. */
    std::unordered_map<std::string, std::string> cache_;
    mutable BoundarySignatureCacheStats stats_;

    /** Builds automorphismEdgeImages_, exactly once, on first use. */
    void ensureAutomorphismGroup_();

    /**
     * Canonicalizes `edgeIndices` by minimizing its sorted image over
     * every automorphism in automorphismEdgeImages_ (triggering
     * ensureAutomorphismGroup_() first), returning the result as a
     * comma-joined string suitable as a cache_ key.
     */
    std::string canonicalKey_(const std::vector<size_t> &edgeIndices);
};

} // namespace identify

namespace census {

/**
 * Whether resolveRecognition() should attempt retriangulateAndLookup() on
 * a local-census/real-Census::lookup() miss. Off by default -- materially
 * more expensive than a plain lookup, so existing callers (surfer.cpp)
 * must opt in via --retriangulate-on-miss; verifyslicegenus defaults it
 * on, since resolving non-hyperbolic far-ends is core to its purpose.
 *
 * Set once, before any search worker thread is spawned, same contract as
 * linkcomplement.h's simplifyComplements.
 */
extern std::atomic<bool> retriangulateOnMiss;

/**
 * Parameters resolveRecognition() passes to retriangulateAndLookup() on
 * every retriangulateOnMiss attempt -- retriangulateAndLookup()'s own
 * defaults (height=2, candidateBudget=8000, timeBudget=20s) exist for
 * direct/test callers only; production recognition always goes through
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
 * SnapPy-identified isoSigs Regina's own census misses, plus anything
 * inserted directly via insertCensusEntry() (see tools/gen_census.py).
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
 * linkcomplement.h's simplifyComplements/identify::recognitionCacheLimit.
 */
bool setCensusPath(const std::string &path);

/**
 * Points localCensusLookup() at a path guaranteed not to exist, so it
 * always misses. Test-only: production code should only ever call
 * setCensusPath() once, at startup.
 */
void resetCensusForTesting();

/**
 * Inserts (isoSig -> name) into the local census at the current path
 * (creating the file/table if missing), tagged with `source` for
 * provenance. INSERT OR IGNORE: never overwrites an existing entry. Bumps
 * the same generation counter setCensusPath() does, so other threads'
 * cached read connections pick up the change on their next lookup.
 * Thread-safe (its own mutex-guarded read-write connection, opened
 * lazily, WAL + busy_timeout enabled for robustness against another
 * concurrently-running surfer/verifyslicegenus instance sharing the same
 * file). No-op (returns false) if the path can't be opened for writing.
 */
bool insertCensusEntry(const std::string &isoSig, const std::string &name,
                       const std::string &source = "verifyslicegenus");

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

#endif // IDENTIFYCOMPLEMENT_H
