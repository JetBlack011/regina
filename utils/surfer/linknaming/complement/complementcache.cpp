//
//  complementcache.cpp
//

#include "linknaming/complement/complementcache.h"

#include <mutex>
#include <unordered_map>

std::atomic<size_t> complement::cacheLimit{200'000};

namespace {

// Memoizes complement answers (both recogniseHandlebody()'s genus and, for
// non-handlebody complements, Census::lookup()'s name) by isomorphism
// signature. The same boundary complement (e.g. a hyperbolic knot
// guaranteed to be a census hit) tends to recur across many found surfaces
// -- every real Census::lookup() reopens six on-disk census databases from
// scratch under censusLookupMutex, and recogniseHandlebody() itself is not
// free either, so caching turns "one answer per surface" into "one
// answer per distinct complement". Guarded by its own
// cachedAnswersMutex, separate from censusLookupMutex, so that a cache
// hit -- the common case once a search has been running a while, and the
// only case isUnknot()'s hot path ever takes -- never blocks behind a slow
// in-flight Census::lookup() on another thread.
std::mutex cachedAnswersMutex;
std::unordered_map<std::string, complement::ComplementAnswer> cachedAnswers;
complement::ComplementCacheStats answerStats;

} // namespace

namespace complement {

std::optional<ComplementAnswer> lookupAnswer(const std::string &sig) {
    std::lock_guard<std::mutex> lock(cachedAnswersMutex);
    auto it = cachedAnswers.find(sig);
    if (it == cachedAnswers.end())
        return std::nullopt;
    return it->second;
}

ComplementAnswer cacheAnswer(const std::string &sig,
                                   const ComplementAnswer &update) {
    std::lock_guard<std::mutex> lock(cachedAnswersMutex);
    if (cachedAnswers.find(sig) == cachedAnswers.end() &&
            cachedAnswers.size() >=
                complement::cacheLimit.load(
                    std::memory_order_relaxed)) {
        cachedAnswers.clear();
        ++answerStats.cacheResets;
    }
    ComplementAnswer &entry = cachedAnswers[sig];
    if (update.genus && !entry.genus)
        entry.genus = update.genus;
    if (update.censusChecked && !entry.censusChecked) {
        entry.censusChecked = true;
        entry.censusName = update.censusName;
    }
    if (update.retriangulateAttempted && !entry.retriangulateAttempted) {
        entry.retriangulateAttempted = true;
        if (update.censusName && !entry.censusName)
            entry.censusName = update.censusName;
    }
    return entry;
}

std::optional<ComplementAnswer>
checkAnswer(const std::string &sig, long long ComplementCacheStats::*checks,
                 long long ComplementCacheStats::*hits,
                 const std::function<bool(const ComplementAnswer &)> &answers) {
    std::lock_guard<std::mutex> lock(cachedAnswersMutex);
    ++(answerStats.*checks);
    auto it = cachedAnswers.find(sig);
    if (it != cachedAnswers.end() && answers(it->second)) {
        ++(answerStats.*hits);
        return it->second;
    }
    return std::nullopt;
}

void countInCache(const std::function<void(ComplementCacheStats &)> &update) {
    std::lock_guard<std::mutex> lock(cachedAnswersMutex);
    update(answerStats);
}

ComplementCacheStats cacheStats() {
    std::lock_guard<std::mutex> lock(cachedAnswersMutex);
    return answerStats;
}

size_t cacheSize() {
    std::lock_guard<std::mutex> lock(cachedAnswersMutex);
    return cachedAnswers.size();
}

void resetCacheForTesting() {
    std::lock_guard<std::mutex> lock(cachedAnswersMutex);
    cachedAnswers.clear();
    answerStats = ComplementCacheStats{};
}

} // namespace complement
