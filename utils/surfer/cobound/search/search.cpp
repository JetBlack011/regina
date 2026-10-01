//
//  search.cpp
//

#include "cobound/search/search.h"

#include <algorithm>
#include <map>
#include <stdexcept>
#include <unordered_set>

#include <sys/resource.h>

#include "cobound/driver/timers.h"
#include "cobound/search/incoming.h"
#include "cobound/search/preconditions.h"
#include "surfer/enumeration/surfacesearch.h"

namespace rowsearch {

BoundaryCondition conditionFor(BoundaryConditionMode mode, int componentCount) {
    switch (mode) {
    case BoundaryConditionMode::proper:
        return BoundaryCondition::proper;
    case BoundaryConditionMode::connected:
    case BoundaryConditionMode::automatic:
    default:
        return componentCount == 1 ? BoundaryCondition::connected
                                   : BoundaryCondition::proper;
    }
}

RowWatchdog::RowWatchdog(WatchdogLimits limits,
                         std::function<void(const char *)> endRow)
    : limits_(std::move(limits)), endRow_(std::move(endRow)) {
    if (!limits_.any())
        return;
    thread_ = std::thread([this] {
        const auto rowDeadline =
            std::chrono::steady_clock::now() +
            std::chrono::duration<double>(limits_.rowSeconds.value_or(0));
        while (!done_.load(std::memory_order_relaxed)) {
            {
                std::unique_lock<std::mutex> lock(wakeMutex_);
                wake_.wait_for(lock, std::chrono::milliseconds(200), [this] {
                    return done_.load(std::memory_order_relaxed);
                });
            }
            if (done_.load(std::memory_order_relaxed))
                break;
            // Checked before the clocks: see the class comment.
            if (limits_.surfaceTarget &&
                satisfying_.load(std::memory_order_relaxed) >=
                    *limits_.surfaceTarget) {
                endRow_("surface-target");
                break;
            }
            const auto now = std::chrono::steady_clock::now();
            if (limits_.rowSeconds && now >= rowDeadline) {
                endRow_("timeout");
                break;
            }
            if (limits_.sweepSeconds &&
                now - limits_.sweepStart >
                    std::chrono::duration<double>(*limits_.sweepSeconds)) {
                endRow_("timeout");
                break;
            }
            // Quiescence: this row has stopped teaching us anything new.
            // Meaningful during the drain too, which is where witnesses are
            // actually identified.
            if (limits_.quiescenceSeconds &&
                limits_.idleMillis() >
                    static_cast<long long>(*limits_.quiescenceSeconds * 1000)) {
                endRow_("quiescent");
                break;
            }
        }
    });
}

void RowWatchdog::stop() {
    {
        std::lock_guard<std::mutex> lock(wakeMutex_);
        done_.store(true, std::memory_order_relaxed);
    }
    wake_.notify_all();
    if (thread_.joinable())
        thread_.join();
}

RowWatchdog::~RowWatchdog() { stop(); }

} // namespace rowsearch

namespace cascade {

HopSearcher::HopSearcher(const farside::SignatureTable &signatures,
                         const exactnaming::ExactTables *exact, HopShape shape,
                         unsigned threads,
                         std::shared_ptr<exactnaming::TableCaches> exactCaches)
    : signatures_(signatures), exact_(exact), shape_(shape), threads_(threads),
      exactCaches_(std::move(exactCaches)) {
  if (exact_ && !exactCaches_)
    exactCaches_ = std::make_shared<exactnaming::TableCaches>(*exact_);
}

HopRun HopSearcher::run(const farside::WitnessRedrawer &row,
                        const std::string &rowName, long long surfaceTarget,
                        double seconds,
                        const std::function<bool(const KeptSurface &)> &stop,
                        const SearchFrontier *resume) const {
  const auto wall0 = std::chrono::steady_clock::now();
  const double cpu0 = timers::processCpuSeconds();
  const rowsearch::RowBuild &rb = row.rowBuild();
  if (rb.seedFaces.empty())
    throw std::runtime_error("hop: the row has no collar seed");

  // Declared before the search, which holds a pointer to it.
  farside::DiagramNamer namer(rb.link.tri, rb.pdcode.size(), *rb.cob, signatures_);
  if (exact_) namer.enableExactNames(*exact_, exactCaches_);

  SurfaceSearch e(rb.tri, rb.seedFaces, rb.searchSideBC);
  SurfaceSearchLimits limits;
  limits.pendingSurfaceCap = shape_.pendingSurfaceCap;
  limits.petalCacheLimit = shape_.petalCacheLimit;
  limits.boundarySignatureCacheLimit = shape_.boundarySignatureCacheLimit;
  // Only a multi-curve component's curve COUNT is ever used, as in
  // verifyslicegenus; and no pair signatures (faces are kept instead).
  limits.nameLinkCurves = false;
  limits.capturePairSig = false;
  e.configureLimits(limits);
  e.setBoundaryNamer(namer);
  // Always recorded (it costs one fingerprint): a later hop from this node
  // carries on from it instead of searching this prefix again.
  e.setRecordFrontier(true);
  e.setResumeFrontier(resume);
  if (const size_t touching = e.countSearchableFacesTouching(rb.searchSideBC))
    throw std::runtime_error("hop: " + std::to_string(touching) +
                             " searchable non-seed triangles touch the search "
                             "side, so found surfaces could change it");
  e.primeBoundaryName(rb.searchSideBC, rb.searchEdges, rowName);
  e.configureSelfIntersections({.resolveUnlinked = shape_.resolveUnlinked,
                                .census = nullptr});

  std::string outcome = "exhausted";
  std::mutex outcomeMutex;
  auto noteStop = [&](const char *why) {
    std::lock_guard<std::mutex> lock(outcomeMutex);
    if (outcome == "exhausted") outcome = why;
  };

  rowsearch::RowAccounting acct;
  std::mutex keptMutex;
  std::unordered_set<std::string> keys;
  HopRun out;
  std::atomic<bool> stopped{false};

  std::optional<rowsearch::RowWatchdog> watchdog;
  SurfaceSearchCallbacks callbacks;
  callbacks.onInterrupted = [&] { noteStop("interrupted"); };
  callbacks.onProgress = [&](const SearchStats &stats) {
    if (watchdog) watchdog->publishSatisfying(stats.satisfyingCount);
  };
  callbacks.surfaceTarget = surfaceTarget;
  callbacks.onSurfaceTarget = [&] { noteStop("surface-target"); };
  callbacks.onBoundaryProcessingStarted = [&](size_t total, unsigned) { out.drainTail = total; };
  callbacks.onBoundaryProcessingComplete = [&](size_t, std::chrono::steady_clock::duration d) {
    out.drainTailSeconds = std::chrono::duration<double>(d).count();
  };

  callbacks.onSurfaceBoundaryProcessed = [&](const SurfaceBoundaryInfo &info) {
    acct.described.fetch_add(1, std::memory_order_relaxed);
    const rowsearch::GatedSurface g = rowsearch::gateSurface(info, rb);
    if (!g.accepted()) {
      acct.reject(g.gate);
      return;
    }
    auto link = farside::orientedOutgoingLink(g.orientedLinks, g.surfaceOf,
                                              row.outgoing(), row.row(),
                                              rb.searchSideBC);
    if (!link) {
      // incomingFlips() fails exactly where classifyRowOrientation() does,
      // which the gate has just passed: impossible.
      acct.reject(rowsearch::Gate::orientationBroken);
      return;
    }
    cobordismgraph::Witness w;
    w.subject = rowName;
    w.subjectComponents = rb.componentCount;
    w.genus = info.tubedGenus;
    w.tubed = !info.connected;
    w.resolvedVertices = info.resolvedVertices;
    std::string farName;
    if (g.split.otherSides.empty()) {
      w.kind = cobordismgraph::WitnessKind::direct;
    } else {
      w.kind = cobordismgraph::WitnessKind::cobordism;
      farName = rowsearch::farSideName(g, rb, &namer);
      w.other = farName;
      w.otherComponents = g.split.otherSides.front().components;
    }
    std::string key = rowsearch::keptKey(w, *link, row);
    {
      std::lock_guard<std::mutex> lock(keptMutex);
      if (!keys.insert(key).second) {
        acct.duplicate.fetch_add(1, std::memory_order_relaxed);
        return;
      }
    }
    acct.recorded.fetch_add(1, std::memory_order_relaxed);
    KeptSurface k{.link = std::move(*link),
                  .genus = info.tubedGenus,
                  .resolvedVertices = info.resolvedVertices,
                  .farName = std::move(farName),
                  .faces = info.captureFaces(),
                  .key = std::move(key),
                  .witness = std::move(w)};
    std::lock_guard<std::mutex> lock(keptMutex);
    if (stop && !stopped.load() && stop(k)) {
      stopped.store(true);
      noteStop("stopped");
      e.requestStop();
      e.skipRemainingBoundaryProcessing();
    }
    out.kept.push_back(std::move(k));
  };

  watchdog.emplace(rowsearch::WatchdogLimits{.surfaceTarget = surfaceTarget,
                                             .rowSeconds = seconds},
                   [&](const char *why) {
                     noteStop(why);
                     e.requestStop();
                   });
  const auto searchStart = std::chrono::steady_clock::now();
  out.setup = std::chrono::duration<double>(searchStart - wall0).count();
  const SearchStats stats = e.search(
      threads_, BoundaryCondition::proper, callbacks, shape_.iddfsIterations,
      shape_.iddfsStep, shape_.iddfsStart, std::nullopt,
      /*orientableOnly=*/true, shape_.maxFaces, shape_.rootBudgetStart,
      shape_.rootBudgetGrowth);
  out.search = std::chrono::duration<double>(std::chrono::steady_clock::now() - searchStart).count();
  for (auto r : stats.profile.rounds) out.rounds.push_back(std::chrono::duration<double>(r).count());
  watchdog->stop();

  const bool drainSkipped = e.boundaryProcessingSkipped();
  {
    const farside::NamingStats &ns = namer.stats();
    out.naming = ns.summary();
    out.namingDiagramSeconds = ns.microsDiagram / 1e6;
    out.namingFallbackSeconds = ns.microsFallback / 1e6;
    out.namingExactSeconds = ns.microsExact / 1e6;
    out.namingSlowestSeconds = ns.slowestMicros() / 1e6;
  }
  out.accepted = stats.satisfyingCount;
  out.accounting = acct.summary(out.accepted, drainSkipped);
  out.accountingFailure = acct.failure(out.accepted, e.rebuildFailures(), drainSkipped);
  out.outcome = outcome;
  out.resumed = e.resumedFrontier();
  out.resumeRefusal = e.resumeRefusal();
  // Only a prefix whose every surface was examined may be skipped later.
  if (out.accountingFailure.empty() && !drainSkipped)
    out.frontier = e.frontier();
  out.wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();
  out.cpu = timers::processCpuSeconds() - cpu0;
  return out;
}

} // namespace cascade
