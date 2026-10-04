//
//  complementcache.h
//
//  The one cache of answers about drilled complements.
//

#ifndef SURFER_LINKNAMING_COMPLEMENTCACHE_H
#define SURFER_LINKNAMING_COMPLEMENTCACHE_H

#include <atomic>
#include <functional>
#include <optional>
#include <string>
#include <sys/types.h>

/*! \file utils/surfer/linknaming/complement/complementcache.h
 *  \brief The one cache of answers about drilled complements, keyed by the
 *  complement's isomorphism signature: its handlebody genus (unlinknaming.h)
 *  and its census name (censusnaming.h), merged per isoSig, with one mutex,
 *  one set of counters and one entry limit.
 *
 *  The same boundary complement tends to recur across many found surfaces,
 *  and neither answer is free -- recogniseHandlebody() is normal-surface
 *  theory, and Regina's Census::lookup() reopens its on-disk databases under
 *  a mutex -- so caching turns "one recognition per surface" into "one
 *  recognition per distinct complement". The cache's own mutex is separate
 *  from the census lookup's, so a cache hit -- the common case once a search
 *  has run a while, and the only case isUnknot()'s hot path ever takes --
 *  never blocks behind a slow lookup on another thread.
 */

namespace complement {

/**
 * The outcome of recognizing a complement, memoized by isomorphism
 * signature (see the cache below).
 */
struct ComplementAnswer {
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
struct ComplementCacheStats {
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

    /** One kind of complement's Pachner-search (retriangulateAndLookup())
     * cost and yield. */
    struct PachnerCounts {
        long long attempts = 0;
        long long successes = 0;
        long long milliseconds = 0; /**< Wall time, summed over threads. */
    };
    PachnerCounts pachnerKnots; /**< One-cusp complements. */
    PachnerCounts pachnerLinks; /**< Multi-cusp complements (see
                                     census::retriangulateLinks). */
};

/** A snapshot of the recognition cache's current hit/miss counters. */
ComplementCacheStats cacheStats();

/** The recognition cache's current entry count (distinct isoSigs seen). */
size_t cacheSize();

/**
 * Clears recognitionCache and its stats outright, ignoring
 * recognitionCacheLimit. Test-only: production code should only ever see
 * this cache clear itself automatically via recognitionCacheLimit.
 */
void resetCacheForTesting();

/**
 * Entry-count threshold past which recognitionCache (complementcache.cpp)
 * clears itself entirely before admitting the next new isoSig -- keyed by
 * content (an isoSig string), not a recyclable integer id, so a lookup
 * racing a clear is simply a clean miss, never a wrong hit against an
 * unrelated key (unlike PetalCache, no epoch-tagging is needed here).
 *
 * Set once, before any search worker thread is spawned, same contract as
 * linkcomplement.h's simplifyComplements.
 */
extern std::atomic<size_t> cacheLimit;

/**
 * Returns a snapshot of sig's current cache entry, or nullopt if unseen.
 * By value, not by reference: another thread may concurrently complete
 * this same entry (e.g. filling in the census fields after this snapshot
 * was taken), so a reference into the map would be a data race on the
 * struct's fields even though the map itself never erases (and hence never
 * invalidates references to existing elements).
 */
std::optional<ComplementAnswer> lookupAnswer(const std::string &sig);

/**
 * Merges `update` into sig's entry monotonically -- genus, once computed,
 * is a deterministic function of sig, so first-write-wins is safe there;
 * censusChecked/retriangulateAttempted only ever move false -> true. This
 * is what keeps a thread racing to complete an entry from clobbering
 * another thread's already-finished result. Returns the post-merge
 * snapshot. Clears the whole cache first when sig is new and the cache is
 * at recognitionCacheLimit.
 */
ComplementAnswer cacheAnswer(const std::string &sig,
                                   const ComplementAnswer &update);

/**
 * One counted lookup, under the cache's mutex: adds one to the counter
 * `checks`, and returns a snapshot of sig's entry -- adding one to `hits`
 * -- when there is one and `answers` says it holds what the caller needs;
 * nullopt otherwise.
 */
std::optional<ComplementAnswer>
checkAnswer(const std::string &sig, long long ComplementCacheStats::*checks,
                 long long ComplementCacheStats::*hits,
                 const std::function<bool(const ComplementAnswer &)> &answers);

/** Applies `update` to the cache's counters, under its mutex. */
void countInCache(const std::function<void(ComplementCacheStats &)> &update);

} // namespace complement

#endif // SURFER_LINKNAMING_COMPLEMENTCACHE_H
