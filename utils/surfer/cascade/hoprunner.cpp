// hoprunner.cpp

#include "hoprunner.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <map>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <thread>
#include <unordered_set>

#include <sys/resource.h>

#include "cobordismgraph.h"
#include "pairsig.h"
#include "rowsearch.h"
#include "surfacesearch.h"

namespace cascade {

namespace {

double cpuSeconds() {
  rusage u{};
  getrusage(RUSAGE_SELF, &u);
  return static_cast<double>(u.ru_utime.tv_sec + u.ru_stime.tv_sec) +
         static_cast<double>(u.ru_utime.tv_usec + u.ru_stime.tv_usec) / 1e6;
}

// Which row components and how many far-side curves each surface component
// carries, as a canonical string: surface components are unlabelled, so the
// per-component entries are sorted.
std::string groupingOf(const farside::OutgoingLink &link,
                       const farside::WitnessRedrawer &row) {
  std::map<size_t, std::pair<std::vector<size_t>, int>> bySurface;
  for (size_t i = 0; i < link.incomingFirstEdge.size(); ++i)
    bySurface[link.incomingSurfaceComponent[i]].first.push_back(
        row.rowComponentOf(link.incomingFirstEdge[i]));
  for (size_t sc : link.surfaceComponent) ++bySurface[sc].second;
  std::vector<std::string> parts;
  for (auto &[sc, entry] : bySurface) {
    std::sort(entry.first.begin(), entry.first.end());
    std::string s;
    for (size_t rc : entry.first) s += std::to_string(rc) + '.';
    parts.push_back(s + ':' + std::to_string(entry.second));
  }
  std::sort(parts.begin(), parts.end());
  std::string out;
  for (const std::string &p : parts) out += p + '|';
  return out;
}

} // namespace

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
  const double cpu0 = cpuSeconds();
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
  e.setBoundaryNamer(&namer);
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
    std::string key = cobordismgraph::witnessIdentity(w) + '\x1f' + groupingOf(*link, row);
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
  out.cpu = cpuSeconds() - cpu0;
  return out;
}

std::string pairSigOf(const regina::Triangulation<4> &thickening,
                      const std::vector<int> &faces) {
  return pairSig<4, 2>(thickening, faces);
}

std::vector<std::string> pairSigsOf(const std::vector<SignRequest> &requests,
                                    unsigned threads, const std::string &cacheDir) {
  std::map<std::pair<std::string, int>, std::vector<size_t>> byRow;
  for (size_t i = 0; i < requests.size(); ++i)
    byRow[{requests[i].rowPD, requests[i].layers}].push_back(i);
  std::vector<const std::pair<const std::pair<std::string, int>, std::vector<size_t>> *> rows;
  for (const auto &entry : byRow) rows.push_back(&entry);

  std::vector<std::string> out(requests.size());
  std::atomic<size_t> next{0};
  std::mutex errorMutex;
  std::string error;
  auto work = [&] {
    for (size_t r; (r = next.fetch_add(1)) < rows.size();) {
      try {
        const auto &[row, indices] = *rows[r];
        rowsearch::RowBuild rb;
        rowsearch::buildRow(row.first, row.second, row.second, /*useCone=*/false, rb);
        const std::unique_ptr<PairSigContext<4, 2>> context =
            cacheDir.empty() ? std::make_unique<PairSigContext<4, 2>>(rb.tri)
                             : PairSigContext<4, 2>::cached(rb.tri, cacheDir);
        for (size_t i : indices) out[i] = context->sig(requests[i].faces);
      } catch (const std::exception &e) {
        std::lock_guard<std::mutex> lock(errorMutex);
        error = e.what();
      }
    }
  };
  std::vector<std::thread> pool;
  const size_t n = std::min<size_t>(std::max(threads, 1u), rows.size());
  for (size_t t = 0; t < n; ++t) pool.emplace_back(work);
  for (auto &t : pool) t.join();
  if (!error.empty()) throw std::runtime_error("pairSigsOf: " + error);
  return out;
}

} // namespace cascade
