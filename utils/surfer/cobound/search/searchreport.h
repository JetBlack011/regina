//
//  searchreport.h
//
//  What one search reports as it runs and when it ends.
//

#ifndef SURFER_COBOUND_SEARCHREPORT_H
#define SURFER_COBOUND_SEARCHREPORT_H

#include <atomic>
#include <chrono>
#include <filesystem>
#include <map>
#include <mutex>
#include <optional>
#include <string>

#include "cobound/search/search.h"
#include "surfer/enumeration/surfacesearch.h"
#include "surfer/report/progress.h"

/*! \file utils/surfer/cobound/search/searchreport.h
 *  \brief One search's report: its live progress block on stderr, the
 *  per-search files it appends to (--surface-stats,
 *  --self-intersection-census), and its log lines. Every format here is
 *  frozen (the atlas's tools parse them: plan, Hard constraint 2).
 */

namespace rowsearch {

/**
 * The status block redrawn in place (surfer/report/progress): the process's
 * one. One is enough: the live-DFS progress (callbacks.onProgress) and the
 * post-search boundary-processing progress (callbacks.onBoundaryProcessing*)
 * never run at the same time (the latter only starts once the former's DFS
 * phase has fully joined), so both share the same redraw region.
 */
extern report::RollingReport progressBlock;

/** Fired from callbacks.onProgress once per second while a search runs. */
void printProgress(const SearchStats &stats, SurfaceSearch &e);

/**
 * Fired from callbacks.onBoundaryProcessingProgress once per second during
 * the post-search boundary-identification phase
 * (processRemainingSurfaceBoundaries). `resolvedGenus`, if set, is the
 * literature target once the search has found a constructive witness for it
 * (resolution is all-or-nothing, so there is no partial genus to show, only
 * "not yet" vs. the exact target once hit), so it's visible at a glance how
 * the current genus stands against the literature target even while this
 * phase is still churning through whatever else was queued.
 */
void printBoundaryProgress(size_t processed, size_t total,
                           std::chrono::steady_clock::duration elapsed,
                           std::optional<int> resolvedGenus, int target);

/**
 * Aggregate distribution of the surfaces a search found: how they spread
 * across the face budget, and across homeomorphism type.
 *
 * The face histogram answers whether --max-faces is actually binding. IDDFS
 * reports each surface once, at the round where it first fits (see
 * suppressBelow/prevCap in embeddingsearch.cpp), so a surface's `triangles`
 * is its true face count and the population sitting at exactly --max-faces
 * is precisely what the cap is truncating.
 *
 * Deliberately an AGGREGATE rather than a per-surface log. --surface-log
 * already writes a row per surface, but a single search can find millions of
 * them, it is overwritten each knot, and it forces a pairSig() (an isoSig
 * computation) per surface. The distinct keys here are bounded by
 * faces x genus x punctures -- a few hundred -- so the whole distribution
 * costs one map lookup per surface and a few hundred CSV lines per row.
 */
struct SurfaceStatsKey {
  long long triangles;
  bool orientable;
  int genus;
  int punctures;
  int tubedGenus;
  int closedComponents;
  bool connected;

  auto operator<=>(const SurfaceStatsKey &) const = default;
};

/**
 * Thread-safe counts keyed by SurfaceStatsKey.
 *
 * record() is called from onSurfaceBoundaryProcessed, which runs on the
 * boundary-identification worker threads. That phase is bounded by
 * identification throughput (a few hundred surfaces a second), so a single
 * mutex is far below the noise floor and buys simplicity over the
 * thread-local-then-merge dance SurfaceTypeTally needs on the hot DFS path.
 */
class SurfaceStatsTally {
public:
  void record(const SurfaceStatsKey &key) {
    std::lock_guard<std::mutex> lock(mutex_);
    ++counts_[key];
  }

  /** Returns the accumulated counts and resets, ready for the next row. */
  std::map<SurfaceStatsKey, long long> take() {
    std::lock_guard<std::mutex> lock(mutex_);
    std::map<SurfaceStatsKey, long long> out;
    out.swap(counts_);
    return out;
  }

private:
  mutable std::mutex mutex_;
  std::map<SurfaceStatsKey, long long> counts_;
};

/**
 * Appends one row's distribution to `path`, creating it (with a header) if
 * it does not yet exist.
 *
 * Appended per row rather than written once at the end so that an
 * interrupted sweep keeps the statistics of every row that did finish --
 * the same reasoning that makes witness checkpoints a per-row operation.
 */
void appendSurfaceStats(const std::filesystem::path &path,
                        const std::string &rowName, long long maxFaces,
                        const std::map<SurfaceStatsKey, long long> &counts);

/**
 * Appends one row of --self-intersection-census output (see
 * SelfIntersectionCensus) for `rowName`, with the search's own resolved
 * count alongside. Per row, like appendSurfaceStats(), so an interrupted
 * run keeps what it measured.
 */
void appendSelfIntersectionCensus(const std::filesystem::path &path,
                                  const std::string &rowName,
                                  long long maxFaces, bool resolveUnlinked,
                                  const SearchStats &stats,
                                  SelfIntersectionCensus &census);

/** How many rejected surfaces per reason per row --rejection-sample-log keeps. */
constexpr int REJECTION_SAMPLES_PER_REASON = 20;

/*
 * verifyslicegenus's per-row lines, from what the search returned
 * (cascade::HopRun) and nothing else. Each writes to `out` exactly as the
 * row loop did, stream state included (the breadth line leaves `out` in
 * std::fixed's precision, as it always has).
 */

/** `breadth:` (the atlas's search_breadth.py parses it): the recorded
 *  frontier, whether the frontier the search was offered was carried on
 *  from, and the frontier's cost. */
void printSweepBreadth(std::ostream &out, const std::string &name,
                       const cascade::HopRun &run);

/** `N new witnesses, outcome X` and `accounting:` (dispatch.py's RE_OUTCOME
 *  and RE_ACCOUNTING). */
void printOutcome(std::ostream &out, const std::string &name, const cascade::HopRun &run);

/** `identification:` (the census and recognition counters as the search ended, against
 *  the search's start), `diagram naming:` and their warnings. */
void printIdentification(std::ostream &out, const std::string &name,
                         const cascade::HopRun &run);

/** `search profile:` (bench_search.sh parses it), with the linking audit
 *  when it is on. */
void printSearchProfile(std::ostream &out, const std::string &name,
                        const cascade::HopRun &run);

} // namespace rowsearch

#endif // SURFER_COBOUND_SEARCHREPORT_H
