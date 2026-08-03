//
//  identifycomplement.cpp
//
//  Created by John Teague on 07/16/2026.
//

#include "identifycomplement.h"

#include <cstdint>
#include <fstream>
#include <iostream>
#include <list>
#include <sstream>

#include <sqlite3.h>

#include <census/census.h>
#include <snappea/snappeatriangulation.h>

#include "pairsig.h"
#include "linknames.h"

std::mutex identify::censusLookupMutex;
std::atomic<size_t> identify::recognitionCacheLimit{200'000};
std::atomic<bool> census::retriangulateOnMiss{false};
std::atomic<int> census::retriangulateHeight{2};
std::atomic<size_t> census::retriangulateCandidateBudget{8000};
std::atomic<long long> census::retriangulateTimeBudgetSeconds{20};

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

// Returns a snapshot of sig's current cache entry, or nullopt if unseen.
// By value, not by reference: another thread may concurrently complete
// this same entry (e.g. filling in the census fields after this snapshot
// was taken), so a reference into the map would be a data race on the
// struct's fields even though the map itself never erases (and hence never
// invalidates references to existing elements).
std::optional<identify::RecognitionResult>
lookupRecognition(const std::string &sig) {
    std::lock_guard<std::mutex> lock(recognitionCacheMutex);
    auto it = recognitionCache.find(sig);
    if (it == recognitionCache.end())
        return std::nullopt;
    return it->second;
}

// Merges `update` into sig's entry monotonically -- genus, once computed,
// is a deterministic function of sig, so first-write-wins is safe there;
// censusChecked/retriangulateAttempted only ever move false -> true. This
// is what keeps a thread racing to complete an entry from clobbering
// another thread's already-finished result. Returns the post-merge
// snapshot.
identify::RecognitionResult
storeRecognition(const std::string &sig,
                 const identify::RecognitionResult &update) {
    std::lock_guard<std::mutex> lock(recognitionCacheMutex);
    if (recognitionCache.find(sig) == recognitionCache.end() &&
            recognitionCache.size() >=
                identify::recognitionCacheLimit.load(
                    std::memory_order_relaxed)) {
        recognitionCache.clear();
        ++recognitionStats.cacheResets;
    }
    identify::RecognitionResult &entry = recognitionCache[sig];
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

// The actual (uncached) Census::lookup() call, serialized under
// censusLookupMutex per its documented contract. No memoization here --
// that's entirely recognitionCache's job now.
std::optional<std::string> censusLookupName(
    const regina::Triangulation<3> &complement) {
    std::lock_guard<std::mutex> lock(identify::censusLookupMutex);
    std::list<regina::CensusHit> hits = regina::Census::lookup(complement);
    if (hits.empty())
        return std::nullopt;

    std::string raw = hits.front().name();
    if (auto classical = linknames::name(raw))
        return *classical + " (" + raw + ")";
    return raw;
}

// Guards censusPathOverride_ -- touched only by census::setCensusPath()
// (effectively write-once, before any search worker thread is spawned; see
// its header doc comment) and by CensusConnection_::ensureCurrent() below
// when a thread lazily (re)opens its connection, which happens far less
// often than census lookups themselves.
std::mutex censusPathConfigMutex_;
std::string censusPathOverride_;

// Bumped by census::setCensusPath()/census::insertCensusEntry(), so every
// thread's already-open CensusConnection_ notices its path (or the file's
// content) may be stale and reopens lazily on its next lookup.
std::atomic<uint64_t> censusGeneration_{0};

// One read-only SQLite connection (plus its one prepared statement) per
// thread, opened lazily on first use -- unlike censusLookupMutex's single
// serialized regina::Census::lookup(), SQLite's read-only mode supports
// many concurrent readers natively, so no cross-thread locking is needed
// here at all.
struct CensusConnection_ {
    uint64_t generation = static_cast<uint64_t>(-1);
    sqlite3 *conn = nullptr;
    sqlite3_stmt *stmt = nullptr;

    ~CensusConnection_() { close(); }

    void close() {
        if (stmt) {
            sqlite3_finalize(stmt);
            stmt = nullptr;
        }
        if (conn) {
            sqlite3_close(conn);
            conn = nullptr;
        }
    }

    std::string currentPath() const {
        std::lock_guard<std::mutex> lock(censusPathConfigMutex_);
        return censusPathOverride_.empty() ? SURFER_CENSUS_PATH
                                           : censusPathOverride_;
    }

    // Reopens against the current path override (or the compiled-in
    // default) if `gen` is newer than what this connection last opened
    // against. A failed open (e.g. the census hasn't been generated yet)
    // just leaves stmt null -- localCensusLookup() treats that as a
    // permanent miss until the next generation bump, not an error.
    void ensureCurrent(uint64_t gen) {
        if (generation == gen)
            return;
        close();
        generation = gen;

        std::string path = currentPath();

        if (sqlite3_open_v2(path.c_str(), &conn, SQLITE_OPEN_READONLY,
                            nullptr) != SQLITE_OK) {
            close();
            return;
        }
        if (sqlite3_prepare_v2(conn,
                "SELECT name, source FROM census WHERE isosig = ?1", -1,
                &stmt, nullptr) != SQLITE_OK) {
            close();
        }
    }
};

// Fast, sound, one-sided proof that `t`'s genus is 1 (the unknot): its
// fundamental group is Z if and only if it's the unknot (Dehn's lemma).
// group() already tries to simplify the presentation internally (same
// idiom as engine/triangulation/dim3/knot.cpp's Poincare-conjecture
// fast path), so a presentation of exactly one generator and no relations
// -- i.e. <a|>, which *is* Z by construction -- is conclusive. A "false"
// here is only ever inconclusive (simplify() didn't collapse it that far),
// never a wrong answer: this can only shorten the path to genus == 1, so
// it's always safe to fall back to recogniseHandlebody() when it fails.
bool groupProvesUnknot(const regina::Triangulation<3> &t) {
    const regina::GroupPresentation &g = t.group();
    return g.countGenerators() == 1 && g.countRelations() == 0;
}

// Fast, sound, one-sided proof that `t` (a LINK complement, possibly
// multiple components) is split -- i.e. `t` is identify(const Link&)'s
// n-component-unlink case -- generalizing groupProvesUnknot() above from
// n == 1 to any n. A presentation with zero relations is free by
// construction (same "sound regardless of how simplify() got there"
// argument as groupProvesUnknot()), and a free fundamental group forces a
// link to be split: by Milnor's prime decomposition theorem, an orientable
// 3-manifold's prime decomposition realizes the Grushko free-product
// decomposition of its fundamental group, so a free pi_1 of rank k means
// `t` splits along essential spheres into k pieces, each with a single
// torus boundary component and pi_1 == Z -- which, by the exact same
// Dehn's-lemma argument groupProvesUnknot() itself relies on, forces each
// piece to be a solid torus. So `t` is the complement of k unknotted,
// pairwise split components, i.e. the k-component unlink. (No need to
// separately check k against the link's actual component count: a link
// complement's H_1 always has rank == component count via Alexander
// duality/meridians, regardless of link type, so if the group is ALSO
// free, its rank -- which must match H_1's rank -- is automatically the
// component count.) As with groupProvesUnknot(), a "false" here is only
// ever inconclusive (simplify() didn't collapse the presentation that
// far), never wrong -- always safe to fall through to the normal
// resolveRecognition() path when it fails.
bool groupProvesUnlink(const regina::Triangulation<3> &t) {
    return t.group().countRelations() == 0;
}

// Fast, sound, one-sided proof that `t`'s genus is -1 (not any
// handlebody, of any genus): a genuine hyperbolic structure rules out
// every handlebody by geometrization (no handlebody -- solid tori
// included -- admits a complete hyperbolic structure, since its boundary
// is compressible). Symmetric to groupProvesUnknot() above: a "false"
// here is inconclusive, never wrong, so it's always safe to fall back to
// recogniseHandlebody(). SnapPeaTriangulation's construction and
// solutionType() are not documented thread-safe, but each caller here
// builds one from its own local, independently-owned Triangulation<3> --
// no state is shared across threads -- and this was stress-tested under
// concurrent load (8 threads, thousands of calls) with no crashes or
// inconsistent results.
bool hyperbolicityProvesNotHandlebody(const regina::Triangulation<3> &t) {
    regina::SnapPeaTriangulation snappea(t);
    return snappea.solutionType() ==
           regina::SnapPeaTriangulation::Solution::Geometric;
}

// Cached recogniseHandlebody(): never touches censusLookupMutex, so this is
// safe to call from identify::isUnknot()'s hot, highly-parallel path
// without risking contention with an in-flight Census::lookup(). Tries the
// two fast, sound one-sided checks above before falling back to the
// expensive normal-surface-theory path -- both are cheap regardless of
// outcome, and either resolving conclusively skips recogniseHandlebody()
// entirely. Neither helps for knots that are neither the unknot nor
// hyperbolic (torus/satellite knots), which still fall through to
// recogniseHandlebody() same as before.
ssize_t cachedGenus(const regina::Triangulation<3> &complement,
                    const std::string &sig) {
    {
        std::lock_guard<std::mutex> lock(recognitionCacheMutex);
        ++recognitionStats.genusChecks;
        auto it = recognitionCache.find(sig);
        if (it != recognitionCache.end() && it->second.genus) {
            ++recognitionStats.genusCacheHits;
            return *it->second.genus;
        }
    }

    ssize_t genus;
    enum class Path { Group, SnapPea, Fallback } path;
    if (groupProvesUnknot(complement)) {
        genus = 1;
        path = Path::Group;
    } else if (hyperbolicityProvesNotHandlebody(complement)) {
        genus = -1;
        path = Path::SnapPea;
    } else {
        genus = complement.recogniseHandlebody();
        path = Path::Fallback;
    }

    {
        std::lock_guard<std::mutex> lock(recognitionCacheMutex);
        if (path == Path::Group)
            ++recognitionStats.groupFastPathHits;
        else if (path == Path::SnapPea)
            ++recognitionStats.snapPeaFastPathHits;
        else
            ++recognitionStats.recogniseHandlebodyFallbacks;
    }
    return *storeRecognition(sig, identify::RecognitionResult{.genus = genus})
                .genus;
}

// Full resolution: genus, and (only if genus == -1) a census check
// (local census, then the real Census::lookup(), then -- if
// retriangulateOnMiss -- census::retriangulateAndLookup()).
identify::RecognitionResult
resolveRecognition(const regina::Triangulation<3> &complement,
                   const std::string &sig) {
    ssize_t genus = cachedGenus(complement, sig);
    if (genus != -1)
        return *lookupRecognition(sig); // fully resolved; census never applies

    {
        std::lock_guard<std::mutex> lock(recognitionCacheMutex);
        ++recognitionStats.censusChecks;
        auto it = recognitionCache.find(sig);
        if (it != recognitionCache.end() && it->second.censusChecked &&
                (it->second.censusName || it->second.retriangulateAttempted)) {
            ++recognitionStats.censusCacheHits;
            return it->second;
        }
    }

    {
        std::lock_guard<std::mutex> lock(recognitionCacheMutex);
        ++recognitionStats.localCensusChecks;
    }
    auto name = census::localCensusLookup(sig);
    if (name) {
        std::lock_guard<std::mutex> lock(recognitionCacheMutex);
        ++recognitionStats.localCensusHits;
    } else {
        name = censusLookupName(complement);
    }

    bool retriangulateAttempted = false;
    if (!name && census::retriangulateOnMiss.load(std::memory_order_relaxed)) {
        name = census::retriangulateAndLookup(
            complement,
            census::retriangulateHeight.load(std::memory_order_relaxed),
            census::retriangulateCandidateBudget.load(
                std::memory_order_relaxed),
            std::chrono::seconds(census::retriangulateTimeBudgetSeconds.load(
                std::memory_order_relaxed)));
        retriangulateAttempted = true;
    }

    return storeRecognition(
        sig, identify::RecognitionResult{
                 .genus = -1,
                 .censusChecked = true,
                 .censusName = name,
                 .retriangulateAttempted = retriangulateAttempted});
}

// Shared tail of identify(const EdgeComplement&)/identify(const Link&),
// once a complement is already built and resolveRecognition() has already
// run for it: turns the result into identify()'s final answer. Factored
// out (rather than having identify(const Link&) just call identify(const
// EdgeComplement&)) specifically so neither caller ever builds the same
// complement twice -- buildComplement() is not cheap.
std::string nameFromRecognition(const identify::RecognitionResult &result,
                                const std::string &sig) {
    if (result.genus == 1)
        return "Unknot";
    if (result.censusName)
        return *result.censusName;
    return sig;
}

} // namespace

namespace identify {

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

std::string identify(const EdgeComplement &e) {
    auto complement = e.buildComplement();
    std::string sig = complement.isoSig();
    RecognitionResult result = resolveRecognition(complement, sig);
    return nameFromRecognition(result, sig);
}

std::string identify(const Link &l) {
    auto complement = l.buildComplement();

    if (l.countComponents() > 1 && groupProvesUnlink(complement))
        return std::to_string(l.countComponents()) + "-component unlink";

    std::string sig = complement.isoSig();
    RecognitionResult result = resolveRecognition(complement, sig);
    return nameFromRecognition(result, sig);
}

bool isOrientationSafeName(const std::string &name) {
    return name == "Unknot" || name.ends_with("-component unlink");
}

bool recognizeComplement(const EdgeComplement &e) {
    auto complement = e.buildComplement();
    std::string sig = complement.isoSig();
    RecognitionResult result = resolveRecognition(complement, sig);

    if (result.genus == 1) {
        std::cout << "      unknot, " << sig << "\n";
        return true;
    }
    if (result.censusName) {
        std::cout << "      recognized as " << *result.censusName << ", "
                  << sig << "\n";
        return true;
    }
    return false;
}

bool isUnknot(const EdgeComplement &e) {
    auto complement = e.buildComplement();
    return cachedGenus(complement, complement.isoSig()) == 1;
}

void recognizeComplement(const Link &l) {
    if (recognizeComplement(static_cast<const EdgeComplement &>(l))) {
        return;
    }

    // buildComplement()/simplify() runs again here (the whole-link
    // recognizeComplement() above already built the same complement) --
    // a pre-existing inefficiency this cache doesn't address, since
    // simplify() has no isoSig to key off of until after it's run. The
    // cache does at least make this second recogniseHandlebody() free.
    auto complement = l.buildComplement();
    std::string sig = complement.isoSig();
    ssize_t genus = cachedGenus(complement, sig);
    int numComponents = l.countComponents();

    if (genus != -1 && numComponents == 1) {
        std::cout << "[!] WARNING! Recognized as a genus " << genus
                  << " handlebody, " << sig << "\n";
        std::cout << "[!] This is almost definitely a bug, please "
                     "report it!\n";
    } else if (numComponents == 1) {
        std::cout << "      NOT unknot, " << sig << "\n";
    } else if (numComponents > 1) {
        for (int i = 0; i < numComponents; ++i) {
            regina::Triangulation<3> compI = l.buildComplement(i);
            std::string sigI = compI.isoSig();
            RecognitionResult result = resolveRecognition(compI, sigI);

            std::cout << "    Component " << i + 1 << ": ";
            if (result.genus == 1) {
                std::cout << "unknot, " << sigI << "\n";
            } else if (result.genus && *result.genus != -1) {
                std::cout << "\n[!] WARNING! Recognized as a genus "
                          << *result.genus << " handlebody, ";
                std::cout << "\n[!] This is almost definitely a bug, "
                             "please "
                             "report it!\n";
            } else if (result.censusName) {
                std::cout << "      recognized as " << *result.censusName
                          << ", " << sigI << "\n";
            } else {
                std::cout << "NOT unknot, " << sigI << "\n";
            }
        }
    }
}

BoundarySignatureCache::BoundarySignatureCache(
    const regina::Triangulation<3> &boundary, size_t clearThreshold)
    : boundary_(&boundary), clearThreshold_(clearThreshold) {}

void BoundarySignatureCache::ensureAutomorphismGroup_() {
    std::call_once(groupOnce_, [this] {
        size_t numEdges = boundary_->countEdges();
        boundary_->findAllIsomorphisms(*boundary_,
            [&](const regina::Isomorphism<3> &alpha) {
                std::vector<size_t> image(numEdges);
                for (size_t e = 0; e < numEdges; ++e)
                    image[e] = resolveFaceIndex<3, 1>(
                        *boundary_,
                        applyIsomorphism(alpha,
                            faceDescriptor<3, 1>(*boundary_,
                                static_cast<int>(e))));
                automorphismEdgeImages_.push_back(std::move(image));
                return false; // keep enumerating every automorphism
            });
    });
}

std::string BoundarySignatureCache::canonicalKey_(
        const std::vector<size_t> &edgeIndices) {
    ensureAutomorphismGroup_();

    // Starting from edgeIndices itself (already sorted, per this method's
    // precondition) means correctness never depends on the identity
    // automorphism actually being among automorphismEdgeImages_ -- only on
    // finding it, or something at least as good, being harmless.
    std::vector<size_t> best = edgeIndices;
    std::vector<size_t> candidate;
    for (const auto &perm : automorphismEdgeImages_) {
        candidate.clear();
        candidate.reserve(edgeIndices.size());
        for (size_t e : edgeIndices)
            candidate.push_back(perm[e]);
        std::ranges::sort(candidate);
        if (candidate < best)
            best = candidate;
    }

    std::ostringstream out;
    for (size_t i = 0; i < best.size(); ++i) {
        if (i)
            out << ',';
        out << best[i];
    }
    return out.str();
}

std::string BoundarySignatureCache::identifyCached(
        const std::vector<size_t> &edgeIndices,
        const std::function<std::string()> &compute) {
    std::string key = canonicalKey_(edgeIndices);

    {
        std::lock_guard<std::mutex> lock(cacheMutex_);
        ++stats_.checks;
        auto it = cache_.find(key);
        if (it != cache_.end()) {
            ++stats_.hits;
            return it->second;
        }
    }

    std::string result = compute();

    {
        std::lock_guard<std::mutex> lock(cacheMutex_);
        if (cache_.find(key) == cache_.end() &&
                cache_.size() >= clearThreshold_) {
            cache_.clear();
            ++stats_.cacheResets;
        }
        cache_.emplace(std::move(key), result);
    }
    return result;
}

BoundarySignatureCacheStats BoundarySignatureCache::stats() const {
    std::lock_guard<std::mutex> lock(cacheMutex_);
    return stats_;
}

size_t BoundarySignatureCache::size() const {
    std::lock_guard<std::mutex> lock(cacheMutex_);
    return cache_.size();
}

} // namespace identify

namespace census {

bool setCensusPath(const std::string &path) {
    {
        std::lock_guard<std::mutex> lock(censusPathConfigMutex_);
        censusPathOverride_ = path;
    }
    censusGeneration_.fetch_add(1, std::memory_order_relaxed);
    return std::ifstream(path).good();
}

void resetCensusForTesting() {
    setCensusPath("/nonexistent-census-for-testing.sqlite");
}

std::optional<std::string> localCensusLookup(const std::string &sig) {
    thread_local CensusConnection_ mc;
    mc.ensureCurrent(censusGeneration_.load(std::memory_order_relaxed));
    if (!mc.stmt)
        return std::nullopt;

    sqlite3_reset(mc.stmt);
    sqlite3_clear_bindings(mc.stmt);
    sqlite3_bind_text(mc.stmt, 1, sig.c_str(), static_cast<int>(sig.size()),
                      SQLITE_TRANSIENT);

    if (sqlite3_step(mc.stmt) != SQLITE_ROW)
        return std::nullopt;

    std::string rawName(
        reinterpret_cast<const char *>(sqlite3_column_text(mc.stmt, 0)));
    std::string source(
        reinterpret_cast<const char *>(sqlite3_column_text(mc.stmt, 1)));

    // Regina-sourced rows store the raw census hit name, same as
    // censusLookupName() gets from CensusHit::name() -- format it
    // identically for byte-identical output. SnapPy/verifyslicegenus/
    // retriangulate-sourced rows already store a best-effort pretty name,
    // not a raw census name, so linknames::name() would just miss on those
    // -- return as-is.
    if (source == "regina") {
        if (auto classical = linknames::name(rawName))
            return *classical + " (" + rawName + ")";
        return rawName;
    }
    return rawName;
}

namespace {

// One lazily-opened, mutex-guarded read-write connection for
// insertCensusEntry() -- distinct from CensusConnection_'s per-thread
// read-only connections above. A single shared connection (rather than
// one per thread) is fine here: writes are low-volume (at most one per
// knot processed by verifyslicegenus, not a hot search-path operation),
// so serializing them costs nothing measurable.
std::mutex writeConnMutex_;
sqlite3 *writeConn_ = nullptr;
std::string writeConnPath_;

// Opens (or reopens, if the configured path has changed since the last
// call) the write connection, creating the file/table if missing and
// enabling WAL + a busy timeout so concurrent readers (this process's own
// CensusConnection_s, or another process's) are never blocked by a brief
// write. Caller must hold writeConnMutex_. Returns false (leaving
// writeConn_ null) if the path can't be opened for writing.
bool ensureWriteConn_() {
    std::string path;
    {
        std::lock_guard<std::mutex> lock(censusPathConfigMutex_);
        path = censusPathOverride_.empty() ? SURFER_CENSUS_PATH
                                           : censusPathOverride_;
    }
    if (writeConn_ && writeConnPath_ == path)
        return true;

    if (writeConn_) {
        sqlite3_close(writeConn_);
        writeConn_ = nullptr;
    }
    writeConnPath_ = path;

    if (sqlite3_open_v2(path.c_str(), &writeConn_,
                        SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE,
                        nullptr) != SQLITE_OK) {
        if (writeConn_)
            sqlite3_close(writeConn_);
        writeConn_ = nullptr;
        return false;
    }

    sqlite3_exec(writeConn_, "PRAGMA journal_mode=WAL;", nullptr, nullptr,
                nullptr);
    sqlite3_exec(writeConn_, "PRAGMA busy_timeout=5000;", nullptr, nullptr,
                nullptr);
    sqlite3_exec(writeConn_,
        "CREATE TABLE IF NOT EXISTS census ("
        "  isosig TEXT PRIMARY KEY, name TEXT NOT NULL, source TEXT NOT NULL"
        ");",
        nullptr, nullptr, nullptr);
    return true;
}

} // namespace

bool insertCensusEntry(const std::string &isoSig, const std::string &name,
                       const std::string &source) {
    std::lock_guard<std::mutex> lock(writeConnMutex_);
    if (!ensureWriteConn_())
        return false;

    sqlite3_stmt *stmt = nullptr;
    if (sqlite3_prepare_v2(writeConn_,
            "INSERT OR IGNORE INTO census(isosig, name, source) "
            "VALUES (?1, ?2, ?3)",
            -1, &stmt, nullptr) != SQLITE_OK)
        return false;

    sqlite3_bind_text(stmt, 1, isoSig.c_str(), static_cast<int>(isoSig.size()),
                      SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 2, name.c_str(), static_cast<int>(name.size()),
                      SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 3, source.c_str(), static_cast<int>(source.size()),
                      SQLITE_TRANSIENT);
    bool ok = sqlite3_step(stmt) == SQLITE_DONE;
    sqlite3_finalize(stmt);

    if (ok)
        censusGeneration_.fetch_add(1, std::memory_order_relaxed);
    return ok;
}

std::optional<std::string> retriangulateAndLookup(
    const regina::Triangulation<3> &complement, int height,
    size_t candidateBudget, std::chrono::seconds timeBudget) {
    auto start = std::chrono::steady_clock::now();
    size_t candidatesTried = 0;
    std::optional<std::string> found;

    // Single-threaded (the `threads` argument below): this already runs
    // from within an already-multithreaded search, and retriangulate()'s
    // own action callback is mutex-serialized internally regardless, so
    // nested parallelism here would only add contention, not throughput.
    complement.retriangulate(height, /*threads=*/1, /*tracker=*/nullptr,
        [&](regina::Triangulation<3> &&t) -> bool {
            if (++candidatesTried > candidateBudget)
                return true; // stop -- budget exhausted
            if (std::chrono::steady_clock::now() - start >= timeBudget)
                return true; // stop -- time budget exhausted

            std::string sig = t.isoSig();
            if (auto name = localCensusLookup(sig)) {
                found = name;
                return true; // stop -- found a match
            }
            return false; // keep searching
        });

    if (found)
        insertCensusEntry(complement.isoSig(), *found, "retriangulate");
    return found;
}

} // namespace census
