//
//  censusnaming.cpp
//

#include "linknaming/census/censusnaming.h"

#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <list>
#include <sstream>

#include <sqlite3.h>

#include <census/census.h>

#include "linknaming/census/censusnames.h"
#include "linknaming/complement/complementcache.h"
#include "linknaming/complement/unlinknaming.h"

std::mutex census::censusLookupMutex;
std::atomic<bool> census::censusUpdates{true};
std::atomic<bool> census::retriangulateOnMiss{false};
std::atomic<bool> census::retriangulateLinks{false};
std::atomic<int> census::retriangulateHeight{2};
std::atomic<size_t> census::retriangulateCandidateBudget{8000};
std::atomic<long long> census::retriangulateTimeBudgetSeconds{20};

namespace {

using complement::cachedGenus;
using complement::checkAnswer;
using complement::countInCache;
using complement::lookupAnswer;
using complement::ComplementCacheStats;
using complement::cacheAnswer;

// The actual (uncached) Census::lookup() call, serialized under
// censusLookupMutex per its documented contract. No memoization here --
// that's entirely complement cache's job now.
std::optional<std::string> censusLookupName(
    const regina::Triangulation<3> &complement) {
    std::lock_guard<std::mutex> lock(census::censusLookupMutex);
    std::list<regina::CensusHit> hits = regina::Census::lookup(complement);
    if (hits.empty())
        return std::nullopt;

    std::string raw = hits.front().name();
    if (auto classical = census::name(raw))
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

// Full resolution: genus, and (only if genus == -1) a census check
// (local census, then the real Census::lookup(), then -- if
// retriangulateOnMiss, and for a link complement only if
// retriangulateLinks too -- census::retriangulateAndLookup()).
complement::ComplementAnswer
resolveAnswer(const regina::Triangulation<3> &complement,
                   const std::string &sig) {
    ssize_t genus = cachedGenus(complement, sig);
    if (genus != -1) {
        // Fully resolved; census never applies. Another thread may have
        // cleared the cache since cachedGenus() stored the genus.
        if (auto cached = lookupAnswer(sig))
            return *cached;
        return complement::ComplementAnswer{.genus = genus};
    }

    // The Pachner search is the expensive rung. A knot's name can bear a
    // slice-genus bound, so it is worth it there; a link's name, taken from
    // its complement alone, never can, so for links it is skipped unless
    // asked for.
    const bool multiComponent = complement.countBoundaryComponents() > 1;
    const bool mayRetriangulate =
        census::retriangulateOnMiss.load(std::memory_order_relaxed) &&
        (!multiComponent ||
         census::retriangulateLinks.load(std::memory_order_relaxed));

    // An entry is final once the census was checked and either named it,
    // or the Pachner search was tried, or is not wanted for it at all --
    // otherwise every repeat would redo the census lookups.
    if (auto hit = checkAnswer(
            sig, &ComplementCacheStats::censusChecks,
            &ComplementCacheStats::censusCacheHits,
            [mayRetriangulate](const complement::ComplementAnswer &r) {
                return r.censusChecked &&
                       (r.censusName || r.retriangulateAttempted ||
                        !mayRetriangulate);
            }))
        return *hit;

    countInCache([](ComplementCacheStats &s) { ++s.localCensusChecks; });
    auto name = census::localCensusLookup(sig);
    if (name) {
        countInCache([](ComplementCacheStats &s) { ++s.localCensusHits; });
    } else {
        name = censusLookupName(complement);
    }

    bool retriangulateAttempted = false;
    if (!name && mayRetriangulate) {
        const auto pachnerStart = std::chrono::steady_clock::now();
        name = census::retriangulateAndLookup(
            complement,
            census::retriangulateHeight.load(std::memory_order_relaxed),
            census::retriangulateCandidateBudget.load(
                std::memory_order_relaxed),
            std::chrono::seconds(census::retriangulateTimeBudgetSeconds.load(
                std::memory_order_relaxed)));
        retriangulateAttempted = true;
        const long long ms =
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - pachnerStart)
                .count();
        countInCache([&](ComplementCacheStats &s) {
            auto &counts = multiComponent ? s.pachnerLinks : s.pachnerKnots;
            ++counts.attempts;
            counts.successes += name ? 1 : 0;
            counts.milliseconds += ms;
        });
    }

    return cacheAnswer(
        sig, complement::ComplementAnswer{
                 .genus = -1,
                 .censusChecked = true,
                 .censusName = name,
                 .retriangulateAttempted = retriangulateAttempted});
}

// Shared tail of nameComplement(const EdgeComplement&)/nameComplement(const Link&),
// once a complement is already built and resolveAnswer() has already
// run for it: turns the result into nameComplement()'s final answer.
// Factored out (rather than having nameComplement(const Link&) just call
// nameComplement(const EdgeComplement&)) specifically so neither caller ever builds the same
// complement twice -- buildComplement() is not cheap.
std::string nameFromAnswer(const complement::ComplementAnswer &result,
                                const std::string &sig) {
    if (result.genus == 1)
        return complement::unlinkName(1);
    if (result.censusName)
        return *result.censusName;
    return sig;
}

} // namespace

namespace census {

namespace {
// See perturbNamesForTesting: every name nameComplement() returns gets a fresh
// suffix, so any decision that depends on a name -- rather than on geometry
// -- changes, and the name-independence test sees it.
std::string perturbed(std::string name) {
    if (!perturbNamesForTesting.load(std::memory_order_relaxed))
        return name;
    static std::atomic<unsigned long long> counter{0};
    return name + " [~" + std::to_string(counter.fetch_add(1) + 1) + "]";
}
} // namespace

std::atomic<bool> perturbNamesForTesting{false};

std::string perturbedForTesting(std::string name) { return perturbed(std::move(name)); }

namespace {
// nameComplement()'s answer for a built complement: the genus check, the census,
// else the isoSig.
std::string nameBuiltComplement(const regina::Triangulation<3> &complement) {
    std::string sig = complement.isoSig();
    complement::ComplementAnswer result = resolveAnswer(complement, sig);
    return perturbed(nameFromAnswer(result, sig));
}
} // namespace

std::string nameComplement(const EdgeComplement &e) {
    return nameBuiltComplement(e.buildComplement());
}

std::string nameComplement(const Link &l) {
    auto complement = l.buildComplement();

    if (l.countComponents() > 1 && complement::groupProvesUnlink(complement))
        return perturbed(complement::unlinkName(l.countComponents()));

    return nameBuiltComplement(complement);
}

bool reportComplement(const EdgeComplement &e) {
    auto complement = e.buildComplement();
    std::string sig = complement.isoSig();
    complement::ComplementAnswer result = resolveAnswer(complement, sig);

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

void reportComplement(const Link &l) {
    if (reportComplement(static_cast<const EdgeComplement &>(l))) {
        return;
    }

    // buildComplement()/simplify() runs again here (the whole-link
    // reportComplement() above already built the same complement) --
    // a pre-existing inefficiency this cache doesn't address, since
    // simplify() has no isoSig to key off of until after it's run. The
    // cache does at least make this second recogniseHandlebody() free.
    auto complement = l.buildComplement();
    std::string sig = complement.isoSig();
    ssize_t genus = complement::cachedGenus(complement, sig);
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
            complement::ComplementAnswer result = resolveAnswer(compI, sigI);

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

} // namespace census

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

    // Reset on every way out. A statement left stepped on a hit keeps this
    // thread's read transaction -- and its SHARED lock -- open until the
    // thread's next lookup. With a dozen threads, some always hold one, so
    // insertCensusEntry() could never commit: each insert waited out the
    // 5 s busy_timeout and failed, silently. Every Pachner success paid
    // those 5 s, and nothing it named ever reached the census for a later
    // search (measured 2026-09-26: c3 added no census row in 694 searches).
    struct ResetOnExit {
        sqlite3_stmt *stmt;
        ~ResetOnExit() { sqlite3_reset(stmt); }
    } resetOnExit{mc.stmt};

    if (sqlite3_step(mc.stmt) != SQLITE_ROW)
        return std::nullopt;

    std::string rawName(
        reinterpret_cast<const char *>(sqlite3_column_text(mc.stmt, 0)));
    std::string source(
        reinterpret_cast<const char *>(sqlite3_column_text(mc.stmt, 1)));

    // Regina-sourced rows store the raw census hit name, same as
    // censusLookupName() gets from CensusHit::name() -- format it
    // identically for byte-identical output. SnapPy-, search-
    // retriangulate-sourced rows already store a best-effort pretty name,
    // not a raw census name, so census::name() would just miss on those
    // -- return as-is.
    if (source == "regina") {
        if (auto classical = census::name(rawName))
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
// knot a search processed, not a hot search-path operation),
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

std::atomic<long long> insertsOk_{0};
std::atomic<long long> insertsFailed_{0};

std::pair<long long, long long> insertCounts() {
    return {insertsOk_.load(std::memory_order_relaxed),
            insertsFailed_.load(std::memory_order_relaxed)};
}

bool insertCensusEntry(const std::string &isoSig, const std::string &name,
                       const std::string &source) {
    if (!censusUpdates.load(std::memory_order_relaxed))
        return false;
    std::lock_guard<std::mutex> lock(writeConnMutex_);
    if (!ensureWriteConn_()) {
        insertsFailed_.fetch_add(1, std::memory_order_relaxed);
        return false;
    }

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

    (ok ? insertsOk_ : insertsFailed_).fetch_add(1, std::memory_order_relaxed);
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
