//
//  complementcache.cpp
//

#include "linknaming/complement/complementcache.h"

#include <mutex>
#include <unordered_map>

std::atomic<size_t> identify::recognitionCacheLimit{200'000};

namespace {

// Memoizes recognition results (both recogniseHandlebody()'s genus and, for
// non-handlebody complements, Census::lookup()'s name) by isomorphism
// signature. The same boundary complement (e.g. a hyperbolic knot
// guaranteed to be a census hit) tends to recur across many found surfaces
// -- every real Census::lookup() reopens six on-disk census databases from
// scratch under censusLookupMutex, and recogniseHandlebody() itself is not
// free either, so caching turns "one recognition per surface" into "one
// recognition per distinct complement". Guarded by its own
// recognitionCacheMutex, separate from censusLookupMutex, so that a cache
// hit -- the common case once a search has been running a while, and the
// only case isUnknot()'s hot path ever takes -- never blocks behind a slow
// in-flight Census::lookup() on another thread.
std::mutex recognitionCacheMutex;
std::unordered_map<std::string, identify::RecognitionResult> recognitionCache;
identify::RecognitionCacheStats recognitionStats;

} // namespace

namespace identify {

std::optional<RecognitionResult> lookupRecognition(const std::string &sig) {
    std::lock_guard<std::mutex> lock(recognitionCacheMutex);
    auto it = recognitionCache.find(sig);
    if (it == recognitionCache.end())
        return std::nullopt;
    return it->second;
}

RecognitionResult storeRecognition(const std::string &sig,
                                   const RecognitionResult &update) {
    std::lock_guard<std::mutex> lock(recognitionCacheMutex);
    if (recognitionCache.find(sig) == recognitionCache.end() &&
            recognitionCache.size() >=
                identify::recognitionCacheLimit.load(
                    std::memory_order_relaxed)) {
        recognitionCache.clear();
        ++recognitionStats.cacheResets;
    }
    RecognitionResult &entry = recognitionCache[sig];
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

std::optional<RecognitionResult>
checkRecognition(const std::string &sig, long long RecognitionCacheStats::*checks,
                 long long RecognitionCacheStats::*hits,
                 const std::function<bool(const RecognitionResult &)> &answers) {
    std::lock_guard<std::mutex> lock(recognitionCacheMutex);
    ++(recognitionStats.*checks);
    auto it = recognitionCache.find(sig);
    if (it != recognitionCache.end() && answers(it->second)) {
        ++(recognitionStats.*hits);
        return it->second;
    }
    return std::nullopt;
}

void countRecognition(const std::function<void(RecognitionCacheStats &)> &update) {
    std::lock_guard<std::mutex> lock(recognitionCacheMutex);
    update(recognitionStats);
}

RecognitionCacheStats recognitionCacheStats() {
    std::lock_guard<std::mutex> lock(recognitionCacheMutex);
    return recognitionStats;
}

size_t recognitionCacheSize() {
    std::lock_guard<std::mutex> lock(recognitionCacheMutex);
    return recognitionCache.size();
}

void resetRecognitionCacheForTesting() {
    std::lock_guard<std::mutex> lock(recognitionCacheMutex);
    recognitionCache.clear();
    recognitionStats = RecognitionCacheStats{};
}

} // namespace identify
