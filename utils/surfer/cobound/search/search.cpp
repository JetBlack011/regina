//
//  search.cpp
//

#include "cobound/search/search.h"

#include <algorithm>
#include <condition_variable>
#include <filesystem>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

#include <sys/resource.h>

#include "cobound/cobordisms/pending.h"
#include "cobound/driver/signals.h"
#include "cobound/driver/timers.h"
#include "cobound/search/incoming.h"
#include "cobound/search/preconditions.h"
#include "cobound/search/searchreport.h"
#include "cobound/frozen.h"
#include "linknaming/census/censusnaming.h"
#include "linknaming/complement/linkcomplement.h"
#include "surfer/submanifold/linkingnumber.h"
#include "surfer/enumeration/surfacesearch.h"
#include "surfer/report/csvwriter.h"

namespace search {

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
            if (limits_.rowSeconds && std::chrono::steady_clock::now() >= rowDeadline) {
                endRow_("timeout");
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

} // namespace search

namespace search {

SearchShape searchShape(const RunShape &shape) {
  if (!shape.layers || !shape.maxFaces || !shape.iddfsIterations || !shape.iddfsStart ||
      !shape.iddfsStep || !shape.rootBudgetStart || !shape.rootBudgetGrowth)
    throw std::logic_error("searchShape(): the hop shape leaves a field that decides what a "
                           "search visits unset (the config sets every one)");
  SearchShape s;
  s.condition = BoundaryCondition::proper;
  s.iddfsIterations = *shape.iddfsIterations;
  s.iddfsStep = *shape.iddfsStep;
  s.iddfsStart = *shape.iddfsStart;
  s.maxFaces = *shape.maxFaces;
  s.rootBudgetStart = *shape.rootBudgetStart;
  s.rootBudgetGrowth = *shape.rootBudgetGrowth;
  s.resolveUnlinked = shape.resolveUnlinked;
  if (shape.pendingSurfaceCap) s.limits.pendingSurfaceCap = *shape.pendingSurfaceCap;
  if (shape.petalCacheLimit) s.limits.petalCacheLimit = *shape.petalCacheLimit;
  if (shape.boundarySignatureCacheLimit)
    s.limits.boundarySignatureCacheLimit = *shape.boundarySignatureCacheLimit;
  // Only a multi-curve component's curve COUNT is ever used, as in
  // verifyslicegenus; and no pair signatures (faces are kept instead).
  s.limits.nameLinkCurves = false;
  s.limits.capturePairSig = false;
  return s;
}

SeedInvariantFailure::SeedInvariantFailure(size_t touching)
    : std::runtime_error("hop: " + std::to_string(touching) +
                         " searchable non-seed triangles touch the search "
                         "side, so found surfaces could change it"),
      touching(touching) {}

Searcher::Searcher(const linknaming::SignatureTable &signatures,
                         const linknaming::ExactTables *exact, RunShape shape,
                         unsigned threads,
                         std::shared_ptr<linknaming::TableCaches> exactCaches)
    : signatures_(&signatures), exact_(exact), shape_(shape), threads_(threads),
      exactCaches_(std::move(exactCaches)) {
  if (exact_ && !exactCaches_)
    exactCaches_ = std::make_shared<linknaming::TableCaches>(*exact_);
}

Searcher::Searcher(const linknaming::SignatureTable *signatures,
                         const linknaming::ExactTables *exact, unsigned threads)
    : signatures_(signatures), exact_(exact), threads_(threads) {
  if (exact_)
    exactCaches_ = std::make_shared<linknaming::TableCaches>(*exact_);
}

SearchResult Searcher::run(const outgoing::OutgoingReader &row,
                        const std::string &rowName, long long surfaceTarget,
                        double seconds,
                        const std::function<bool(const KeptSurface &)> &stop,
                        const SearchFrontier *resume,
                        std::optional<std::string> censusName) const {
  SearchRequest request = requestFor(row, rowName, surfaceTarget, seconds);
  request.resume = resume;
  request.stop = stop;
  request.censusName = std::move(censusName);
  return run(row.rowBuild(), request);
}

SearchRequest Searcher::requestFor(const outgoing::OutgoingReader &row,
                                      const std::string &rowName, long long surfaceTarget,
                                      double seconds) const {
  SearchRequest request;
  request.name = rowName;
  request.row = &row;
  request.shape = searchShape(shape_);
  request.surfaceTarget = surfaceTarget;
  request.seconds = seconds;
  // Always recorded (it costs one fingerprint): a later hop from this node
  // carries on from it instead of searching this prefix again.
  request.recordFrontier = true;
  request.layers = *shape_.layers;
  return request;
}

SearchResult Searcher::run(const search::RowBuild &rb,
                        const SearchRequest &request) const {
  const auto wall0 = std::chrono::steady_clock::now();
  const double cpu0 = timers::processCpuSeconds();
  const SearchShape &shape = request.shape;
  if (!shape.resolveUnlinked)
    throw std::logic_error("HopSearcher::run(): the search shape does not say whether "
                           "resolvable surfaces count (resolve_unlinked has no default)");
  const bool resolveUnlinked = *shape.resolveUnlinked;
  const LiteratureInterval &literature = request.literature;
  const SearchOutputs &outputs = request.outputs;
  if (rb.seedFaces.empty())
    throw SearchRefused("hop: the row has no collar seed");
  if (!request.row)
    throw std::logic_error("HopSearcher::run(): the request has no row to read its finds on");

  // Declared before the search, which holds pointers to them. Every
  // boundary is named by its complement unless the row draws its far sides.
  const outgoing::ComplementNamer complementNamer{};
  std::optional<outgoing::DiagramNamer> namer;
  SurfaceSearch e(rb.tri, rb.seedFaces, rb.incomingBC);
  e.configureLimits(shape.limits);
  // The process owns SIGINT and SIGTERM (driver/signals.h), not the library.
  e.setSigintHandling(false);
  // The search's frontier: carried on from, and recorded. A frontier whose
  // finds may not be signed yet is refused first (divergence 1): skipping
  // its prefix would leave them nowhere but a pending file nobody signs.
  std::string pendingRefusal;
  if (request.resume && request.resume->pending) {
    const SearchFrontier::Pending &p = *request.resume->pending;
    const auto under = [](const std::string &file, const std::string &dir) {
      const std::string f = std::filesystem::weakly_canonical(file).string();
      const std::string d = std::filesystem::weakly_canonical(dir).string();
      return f.size() > d.size() && f.compare(0, d.size(), d) == 0 && f[d.size()] == '/';
    };
    const bool ours = request.runDirectory && under(p.path, *request.runDirectory);
    const long long signedTo = cobordisms::signedThrough(p.path);
    if (!ours && signedTo < p.bytes)
      pendingRefusal = "its pending file " + p.path + " is signed through " +
                       std::to_string(signedTo) + " of the " + std::to_string(p.bytes) +
                       " bytes its prefix needs (sign it first)";
  }
  if (request.resume && pendingRefusal.empty())
    e.setResumeFrontier(request.resume);
  e.setRecordFrontier(request.recordFrontier);
  if (request.pairSigCacheDir)
    e.setPairSigCacheDir(*request.pairSigCacheDir);
  e.setBoundaryNamer(complementNamer);
  if (signatures_) {
    // A row whose far sides cannot be drawn is refused (divergence 2): its
    // T does not read back, and naming it some other way would hide that.
    try {
      namer.emplace(rb.link.tri, rb.pdcode.size(), *rb.cob, *signatures_);
      if (exact_) namer->enableExactNames(*exact_, exactCaches_);
    } catch (const std::exception &ex) {
      throw SearchRefused(ex.what());
    }
    e.setBoundaryNamer(*namer);
  }

  {
    // The invariant that makes the search side fixed: no searchable
    // triangle other than the seed has an edge on it. Checked once here
    // rather than re-derived for every surface found.
    if (const size_t touching = e.countSearchableFacesTouching(rb.incomingBC))
      throw SeedInvariantFailure(touching);
    // Its name is known by construction; never identify it.
    e.primeBoundaryName(rb.incomingBC, rb.incomingEdges, request.name);
  }

  // Fresh per search, so each census line describes one search.
  std::optional<SelfIntersectionCensus> selfIntersectionCensus;
  if (outputs.selfIntersectionCensus) {
    selfIntersectionCensus.emplace();
    selfIntersectionCensus->incomingBoundary = static_cast<long>(rb.incomingBC);
  }
  e.configureSelfIntersections(
      {.resolveUnlinked = resolveUnlinked,
       .census = selfIntersectionCensus ? &*selfIntersectionCensus : nullptr});

  std::optional<search::SurfaceStatsTally> surfaceStats;
  if (outputs.surfaceStats)
    surfaceStats.emplace();

  std::optional<CsvWriter> surfaceLog;
  if (outputs.surfaceLog)
    surfaceLog.emplace(*outputs.surfaceLog,
                       "orientable,genus,tubed_genus,punctures,triangles,"
                       "pairsig",
                       threads_);

  // Why the search stopped. Set by whoever requests the stop; "exhausted"
  // means nobody did and the search ran out of candidates on its own, which
  // is only possible with a face cap (see EmbeddingSearch::search()'s
  // hardFaceCap).
  std::string outcome = "exhausted";
  std::mutex outcomeMutex;
  auto noteStop = [&](const char *why) {
    std::lock_guard<std::mutex> lock(outcomeMutex);
    if (outcome == "exhausted") outcome = why;
  };

  SearchResult out;
  out.recognitionBefore = complement::recognitionCacheStats();
  out.censusWritesBefore = census::insertCounts();

  // Every surface the drain describes lands in exactly one of its buckets
  // (see the accounting after the search).
  search::RowAccounting acct;
  // The fixtures for divergence 2's test: with SURFER_TEST_UNACCOUNTED
  // naming this search, its first described surface is dropped from every
  // bucket, as a surface lost between the drain and the record would be;
  // with SURFER_TEST_IMPOSSIBLE, it is counted in an impossible bucket.
  std::atomic<bool> dropOne{false}, impossibleOne{false};
  if (const char *v = std::getenv("SURFER_TEST_UNACCOUNTED"); v && request.name == v) {
    dropOne.store(true);
    std::cerr << "[!] SURFER_TEST_UNACCOUNTED: one surface of " << request.name
              << " is dropped from the accounting (test mode)\n";
  }
  if (const char *v = std::getenv("SURFER_TEST_IMPOSSIBLE"); v && request.name == v) {
    impossibleOne.store(true);
    std::cerr << "[!] SURFER_TEST_IMPOSSIBLE: one surface of " << request.name
              << " is counted in a state that cannot occur (test mode)\n";
  }

  // --rejection-sample-log: the pair signatures of the first few surfaces
  // each rejection reason turns away in this search, so a misbehaving gate
  // can be audited offline without re-searching.
  std::mutex rejectionSampleMutex;
  std::map<std::string, int> rejectionSamplesTaken;
  auto sampleRejection = [&](const char *reason, const SurfaceBoundaryInfo &info) {
    if (!outputs.rejectionSamples)
      return;
    {
      std::lock_guard<std::mutex> lock(rejectionSampleMutex);
      if (rejectionSamplesTaken[reason]++ >= search::REJECTION_SAMPLES_PER_REASON)
        return;
    }
    std::string sig = info.capturePairSig ? info.capturePairSig() : std::string{};
    std::lock_guard<std::mutex> lock(rejectionSampleMutex);
    *outputs.rejectionSamples << csvField(request.name) << ',' << reason << ','
                              << info.tubedGenus << ','
                              << (info.connected ? "true" : "false") << ','
                              << csvField(info.boundaryDescription) << ','
                              << csvField(sig) << '\n';
    outputs.rejectionSamples->flush();
  };

  // One kept surface per keptKey() (divergence 7), none whose identity the
  // database a run loaded already holds.
  std::mutex keptMutex;
  std::unordered_set<std::string> keys;
  // Distinct identities among them (under keptMutex).
  std::unordered_set<std::string> newIdentities;
  // The search's pending file (divergence 7): every kept surface, appended
  // and fsynced as the search runs.
  std::optional<cobordisms::PendingWriter> pending;
  if (request.pending)
    pending.emplace(*request.pending);
  std::atomic<bool> stopped{false};
  // request.judge: whether the graph has judged a find constructive.
  std::atomic<bool> constructive{false};
  std::mutex fatalMutex;
  // request.judge: the cobordism graph's judgement of each find, on its own
  // thread, finds queued in the order they were kept (the drain never waits
  // for the graph). Started just before the search; drained after it.
  struct JudgeQueue {
    std::mutex mutex;
    std::condition_variable wake;
    std::deque<KeptSurface> finds;
    bool closing = false;
    std::thread thread;
  } judging;
  auto judgeFind = [&](const KeptSurface &k) {
    if (!request.judge) return;
    {
      std::lock_guard<std::mutex> lock(judging.mutex);
      judging.finds.push_back(k);
    }
    judging.wake.notify_one();
  };
  // Started just before the search (below). It polls at 200ms but has no
  // access to SearchStats; onProgress is the only place the live count is
  // handed to us, so it publishes it there.
  std::optional<search::RowWatchdog> watchdog;
  SurfaceSearchCallbacks callbacks;
  // A stop nobody else noted (a signal: divergence 8) is not running out of
  // candidates.
  callbacks.onInterrupted = [&] { noteStop("interrupted"); };
  callbacks.onProgress = [&](const SearchStats &stats) {
    // The process's signal policy (driver/signals.h): a first SIGINT or
    // SIGTERM stops this search within a tick; its drain still finishes.
    if (runsignals::interrupted()) e.requestStop();
    if (watchdog) watchdog->publishSatisfying(stats.satisfyingCount);
    if (outputs.progress) search::printProgress(stats, e);
    // The pending file's 60 s writes run through the enumeration too, not
    // only the post-search drain (below): the drain describes surfaces while
    // the search runs, so a kill during enumeration loses at most a minute.
    if (pending) pending->checkpoint(/*force=*/false);
  };
  // The equalising rule, checked where each surface is counted, so a search
  // stops at the target rather than a progress tick later (see
  // SearchCallbacks::surfaceTarget). The watchdog below still checks it
  // too, as a backstop. As there, the drain is let finish.
  if (request.surfaceTarget) {
    callbacks.surfaceTarget = *request.surfaceTarget;
    callbacks.onSurfaceTarget = [&] { noteStop("surface-target"); };
  }
  // For the `search profile:` line: petal-cache counters as root filtering
  // ends.
  callbacks.onRootsReady = [&] { out.petalsAtRoots = e.petalCacheStats(); };
  callbacks.onBoundaryProcessingStarted = [&](size_t total, unsigned threads) {
    out.drainTail = total;
    if (outputs.progress) {
      search::progressBlock.forget();
      std::cerr << "[+] boundary processing: " << total << " queued surfaces, "
                << threads << " threads\n";
    }
  };
  if (outputs.progress || pending)
    callbacks.onBoundaryProcessingProgress =
        [&](size_t processed, size_t total,
            std::chrono::steady_clock::duration elapsed) {
          if (outputs.progress)
            search::printBoundaryProgress(
                processed, total, elapsed,
                constructive.load() ? std::optional<int>(literature.literatureLo)
                                    : std::nullopt,
                literature.literatureHi);
          // The drain is the long pole of a search and the phase most likely
          // to be interrupted, so checkpoint from here: a kill loses at most
          // a minute of found surfaces.
          if (pending)
            pending->checkpoint(/*force=*/false);
        };
  callbacks.onBoundaryProcessingComplete =
      [&](size_t total, std::chrono::steady_clock::duration elapsed) {
        out.drainTailTime = elapsed;
        out.drainTailSeconds = std::chrono::duration<double>(elapsed).count();
        if (outputs.progress) {
          search::progressBlock.forget();
          std::cerr << "[+] boundary processing: done (" << total
                    << " processed in " << formatElapsed(elapsed) << ")\n";
        }
      };

  callbacks.onSurfaceBoundaryProcessed = [&](const SurfaceBoundaryInfo &info) {
    // Recorded before any of the filtering below: the question this answers
    // is what the SEARCH found, not what survived the checks that decide
    // whether a surface bounds this particular row.
    if (surfaceStats)
      surfaceStats->record(search::SurfaceStatsKey{
          .triangles = info.triangleCount,
          .orientable = info.orientable,
          .genus = info.genus,
          .punctures = info.punctures,
          .tubedGenus = info.tubedGenus,
          .closedComponents = info.closedComponents,
          .connected = info.connected});
    if (surfaceLog) {
      std::string pairSig = info.capturePairSig ? info.capturePairSig() : std::string{};
      std::ostringstream line;
      line << (info.orientable ? "true" : "false") << ',' << info.genus << ','
           << info.tubedGenus << ',' << info.punctures << ',' << info.triangleCount
           << ',' << csvField(pairSig);
      surfaceLog->writeRow(line.str());
    }

    acct.described.fetch_add(1, std::memory_order_relaxed);
    if (dropOne.exchange(false))
      return;
    if (impossibleOne.exchange(false)) {
      acct.reject(search::Gate::orientationBroken);
      return;
    }

    // Orientable, search side intact, the row's own oriented variant, one
    // far side (search::gateSurface()). The orientation check needs the
    // search side's oriented curves, and an exact oriented far-side name the
    // rest, so the gate captures them once.
    const search::GatedSurface g = search::gateSurface(info, rb);
    if (!g.accepted()) {
      acct.reject(g.gate);
      sampleRejection(search::gateReason(g.gate), info);
      return;
    }
    // The oriented outgoing link, kept and judged by: oriented by the gate's
    // own judgement of the incoming curves (g.flips: the row is rb's,
    // request.row->row() is *rb.orientation).
    std::optional<outgoing::OutgoingLink> link;
    {
      link = outgoing::orientedOutgoingLink(g.orientedLinks, g.surfaceOf,
                                           request.row->outgoing(), g.flips,
                                           rb.incomingBC);
      if (!link) {
        // Only a surface component off the row, which the gate's flips (one
        // per component meeting the row) rule out: impossible.
        acct.reject(search::Gate::orientationBroken);
        return;
      }
    }

    // A disconnected find is NOT discarded. Its components tube into a
    // single connected surface with the same boundary and genus exactly
    // info.tubedGenus (see KnottedSurface::tubedSurfaceType), so it
    // witnesses precisely what a connected find of that genus would. This
    // is what makes multi-component links tractable at all: their seeded
    // collar starts as one disjoint annulus per component, and nothing
    // forces the DFS to ever bridge them.
    cobordisms::Cobordism w;
    w.subject = request.name;
    w.subjectComponents = rb.componentCount;
    w.genus = info.tubedGenus;
    w.tubed = !info.connected;
    w.resolvedVertices = info.resolvedVertices;
    // The cobordism as the database will record it: its provenance now.
    w.sourceRow = request.name;
    w.thickenLayers = request.layers;
    w.maxFaces = shape.maxFaces.value_or(0);
    std::string outgoingName;
    if (g.split.otherSides.empty()) {
      w.kind = cobordisms::CobordismKind::direct;
    } else {
      // Exactly one: the gate turns away more (multi-far-side). A
      // genuinely-linked far side is recorded but, unless it is a knot or a
      // proven unlink, it will not carry a bound: a complement does not
      // determine a link. It is kept because the observation is real and is
      // exactly what a later per-witness naming needs as input; the solver's
      // farSideBearsBound() is what declines it.
      const search::BoundarySide &outgoingSide = g.split.otherSides.front();
      // Normalized, so a census hit and the table name are one graph node,
      // and oriented where exact names are on (search::farSideName()).
      outgoingName = search::nameOutgoing(g, namer ? &*namer : nullptr);
      w.kind = cobordisms::CobordismKind::cobordism;
      w.other = outgoingName;
      w.otherComponents = outgoingSide.components;
      // other_candidates is derived from the name when it is signed, as every
      // reader of the field re-derives it.
    }

    // The dedup matters a lot, since every search harvests: a single search
    // reports thousands of near-identical surfaces. One kept surface per
    // keptKey(), and none the loaded database holds already.
    std::string identity = cobordisms::cobordismIdentity(w);
    if (request.knownIdentities && request.knownIdentities->count(identity)) {
      acct.duplicate.fetch_add(1, std::memory_order_relaxed);
      return;
    }
    std::string key = search::keptKey(w, *link, *request.row);
    {
      std::lock_guard<std::mutex> lock(keptMutex);
      if (!keys.insert(key).second) {
        acct.duplicate.fetch_add(1, std::memory_order_relaxed);
        return;
      }
      newIdentities.insert(std::move(identity));
    }
    acct.recorded.fetch_add(1, std::memory_order_relaxed);
    KeptSurface k{.link = std::move(*link),
                  .genus = info.tubedGenus,
                  .resolvedVertices = info.resolvedVertices,
                  .outgoingName = std::move(outgoingName),
                  .faces = info.captureFaces(),
                  .key = std::move(key),
                  .cobordism = std::move(w)};
    std::lock_guard<std::mutex> lock(keptMutex);
    if (pending)
      pending->add(cobordisms::PendingCobordism{k.cobordism, request.rowPD, request.layers, k.faces});
    if (request.stop && !stopped.load() && request.stop(k)) {
      stopped.store(true);
      noteStop("stopped");
      e.requestStop();
      e.skipRemainingBoundaryProcessing();
    }
    judgeFind(k);
    out.kept.push_back(std::move(k));
  };

  // Stops the DFS, and by default LETS THE BOUNDARY DRAIN FINISH.
  //
  // The time limit bounds the *search*, not the row. Identification is the
  // product: an unidentified surface says nothing whatever about a slice
  // genus, so a surface found and then discarded unexamined is pure waste.
  // Cutting the drain short measured 1% identification on a run that
  // produced 1,216,027 qualifying surfaces -- 1.2 million boundaries thrown
  // away unlooked-at, which is why that run's "found nothing" meant nothing.
  //
  // The cost is that a search overruns its nominal limit by however long the
  // queue takes to drain, which can be minutes. That is the right trade:
  // waiting is cheaper than searching a region and then refusing to look at
  // what it found.
  watchdog.emplace(
      search::WatchdogLimits{.surfaceTarget = request.surfaceTarget,
                                .rowSeconds = request.seconds},
      [&](const char *why) {
        noteStop(why);
        e.requestStop();
      });
  if (request.judge)
    judging.thread = std::thread([&] {
      for (;;) {
        KeptSurface k;
        {
          std::unique_lock<std::mutex> lock(judging.mutex);
          judging.wake.wait(lock, [&] { return judging.closing || !judging.finds.empty(); });
          if (judging.finds.empty()) return;
          k = std::move(judging.finds.front());
          judging.finds.pop_front();
        }
        FindJudgement j;
        try {
          j = request.judge(k);
        } catch (const std::exception &ex) {
          j.contradiction = std::string("the cobordism graph failed: ") + ex.what();
        }
        // Only a CONSTRUCTIVE result settles a row: an assisted one is a
        // correct deduction but not an independent verification. Announced
        // on stdout, once, the first time: otherwise the only sign is "--
        // ACHIEVED" in the redrawn stderr progress block, invisible to any
        // log filter, and a row could sit verified-but-unwritten for hours
        // with nothing in the log to say so. Checkpointed at once too: this
        // is the single most valuable moment in a row, and every search
        // harvests, so the row may keep running for hours afterwards.
        if (j.constructive && !constructive.exchange(true, std::memory_order_relaxed)) {
          std::cout << "[+] " << request.name << kFrozenConstructiveWitnessFound
                    << *j.constructive
                    << " (literature [" << literature.literatureLo << ", " << literature.literatureHi
                    << "]). Checkpointing now.\n"
                    << std::flush;
          if (pending)
            pending->checkpoint(/*force=*/true);
        }
        if (j.contradiction.empty()) continue;
        // The graph's gates: something the search or the naming computed is
        // wrong. Nothing more is judged; the search stops at once.
        {
          std::lock_guard<std::mutex> lock(fatalMutex);
          if (out.fatal.empty())
            out.fatal = request.name + ": the cobordism graph's gates: " + j.contradiction;
        }
        noteStop("fatal-bug");
        e.requestStop();
        e.skipRemainingBoundaryProcessing();
        std::lock_guard<std::mutex> lock(judging.mutex);
        judging.finds.clear();
        judging.closing = true;
        return;
      }
    });
  // The judge's last finds are judged before anything reads this search's
  // outcome; on every way out, its thread is joined.
  auto finishJudging = [&] {
    if (!judging.thread.joinable()) return;
    {
      std::lock_guard<std::mutex> lock(judging.mutex);
      judging.closing = true;
    }
    judging.wake.notify_all();
    judging.thread.join();
  };
  struct JudgeJoin {
    std::function<void()> f;
    ~JudgeJoin() { f(); }
  } judgeJoin{finishJudging};

  const auto searchStart = std::chrono::steady_clock::now();
  out.setup = std::chrono::duration<double>(searchStart - wall0).count();
  const SearchStats stats = e.search(
      threads_, shape.condition, callbacks, shape.iddfsIterations,
      shape.iddfsStep, shape.iddfsStart, std::nullopt,
      /*orientableOnly=*/true, shape.maxFaces, shape.rootBudgetStart,
      shape.rootBudgetGrowth);
  out.search = std::chrono::duration<double>(std::chrono::steady_clock::now() - searchStart).count();
  for (auto r : stats.profile.rounds) out.rounds.push_back(std::chrono::duration<double>(r).count());
  finishJudging();
  // Every kept surface is described and judged: the rest of the pending
  // file, fsynced, before anything reads this search's record. A failure
  // here throws: what the search kept would otherwise never be signed.
  if (pending) {
    pending->flush();
    out.pendingPath = pending->path().string();
    out.pendingBytes = pending->syncedBytes();
  }
  watchdog->stop();
  if (surfaceLog)
    surfaceLog->finalize();

  // Accounting: every surface the search accepted must have been described
  // by the drain (unless the drain was deliberately cut short), and every
  // described surface must sit in exactly one bucket. Anything else means
  // surfaces vanished unexamined -- the failure mode that once emptied whole
  // rows without a trace.
  const bool drainSkipped = e.boundaryProcessingSkipped();
  if (namer) {
    const linknaming::NamingStats &ns = namer->stats();
    out.naming = ns.summary();
    out.namingDiagramSeconds = ns.microsDiagram / 1e6;
    out.namingFallbackSeconds = ns.microsFallback / 1e6;
    out.namingExactSeconds = ns.microsExact / 1e6;
    out.namingSlowestSeconds = ns.slowestMicros() / 1e6;
    out.diagramNamed = true;
    out.nonPlanar = ns.nonPlanar;
  }
  out.accepted = stats.satisfyingCount;
  out.accounting = acct.summary(out.accepted, drainSkipped);
  out.accountingFailure = acct.failure(out.accepted, e.rebuildFailures(), drainSkipped);
  out.described = acct.described.load();
  out.recorded = acct.recorded.load();
  out.newCobordisms = static_cast<long long>(newIdentities.size());
  out.otherOrientation = acct.orientation.load();
  out.drainSkipped = drainSkipped;
  // Surfaces were accepted, yet not one reached the record. That can be
  // genuine (every one witnesses another oriented variant), but it is also
  // exactly what a broken gate looks like, so it never licenses a negative.
  out.nothingExamined = acct.nothingExamined();
  out.impossible = acct.impossible();

  if (surfaceStats)
    search::appendSurfaceStats(*outputs.surfaceStats, request.name,
                                  shape.maxFaces.value_or(0), surfaceStats->take());
  if (selfIntersectionCensus)
    search::appendSelfIntersectionCensus(*outputs.selfIntersectionCensus,
                                            request.name, shape.maxFaces.value_or(0),
                                            resolveUnlinked, stats,
                                            *selfIntersectionCensus);
  if (outputs.progress)
    search::progressBlock.forget();

  // Divergence 2: a search that cannot account for its surfaces vouches for
  // no negative: its outcome says so, and it records no frontier (below) and
  // claims no exhaustion (the drivers read accountingFailure).
  if (!out.accountingFailure.empty() || out.nothingExamined)
    noteStop("unaccounted");
  out.outcome = outcome;
  out.resumed = e.resumedFrontier();
  out.resumeRefusal = pendingRefusal.empty() ? e.resumeRefusal() : pendingRefusal;
  if (request.resume)
    out.resumeOfferedRuns = request.resume->runs;
  out.recordedFrontier = e.frontier();
  // It records the pending file and its fsynced length (above: flushed).
  if (out.recordedFrontier && pending)
    out.recordedFrontier->pending = SearchFrontier::Pending{
        std::filesystem::absolute(pending->path()).string(), out.pendingBytes};
  out.frontierSeconds = e.frontierSeconds();
  // Only a prefix whose every surface was examined may be skipped later
  // (divergence 1): the accounting balanced, something was examined, the
  // drain completed, and the pending file is fsynced (it was, above, or the
  // search threw).
  if (out.accountingFailure.empty() && !drainSkipped && !out.nothingExamined)
    out.frontier = out.recordedFrontier;
  out.stats = stats;
  out.petals = e.petalCacheStats();
  out.boundaryCache = e.boundarySignatureCacheStats();
  out.recognitionAfter = complement::recognitionCacheStats();
  out.censusWritesAfter = census::insertCounts();
  out.linkingAudit = linkingnumber::auditLinkingNumbers.load();

  // The incoming knot's complement into the census under its table name,
  // after the counters above were read (it never counted in the row's own
  // identification line). Knots only. A link's complement is shared by
  // infinitely many links (Rolfsen twisting), so writing `isoSig -> L6a3`
  // into a shared cache would assert, permanently and for every future far
  // side landing on that isoSig, an identification the complement cannot
  // support -- exactly the claim linknames.h forbids adding "from a
  // complement match alone". For a knot the same entry is sound by
  // Gordon-Luecke.
  if (request.censusName && rb.componentCount == 1 &&
      census::censusUpdates.load(std::memory_order_relaxed)) {
    Link incoming(rb.link.tri, rb.link.edges);
    census::insertCensusEntry(incoming.buildComplement().isoSig(), *request.censusName);
  }
  out.wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();
  out.cpu = timers::processCpuSeconds() - cpu0;
  return out;
}

} // namespace search
