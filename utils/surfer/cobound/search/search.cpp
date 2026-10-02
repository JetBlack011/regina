//
//  search.cpp
//

#include "cobound/search/search.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

#include <sys/resource.h>

#include "cobound/cobordisms/pairsigner.h"
#include "cobound/driver/timers.h"
#include "cobound/search/incoming.h"
#include "cobound/search/preconditions.h"
#include "cobound/search/searchreport.h"
#include "linknaming/census/censusnaming.h"
#include "surfer/submanifold/linkingnumber.h"
#include "surfer/enumeration/surfacesearch.h"
#include "surfer/report/csvwriter.h"

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

SearchShape searchShape(const HopShape &shape) {
  SearchShape s;
  s.condition = BoundaryCondition::proper;
  s.iddfsIterations = shape.iddfsIterations;
  s.iddfsStep = shape.iddfsStep;
  s.iddfsStart = shape.iddfsStart;
  s.iddfsFinalThreads = std::nullopt;
  s.maxFaces = shape.maxFaces;
  s.rootBudgetStart = shape.rootBudgetStart;
  s.rootBudgetGrowth = shape.rootBudgetGrowth;
  s.resolveUnlinked = shape.resolveUnlinked;
  s.limits.pendingSurfaceCap = shape.pendingSurfaceCap;
  s.limits.petalCacheLimit = shape.petalCacheLimit;
  s.limits.boundarySignatureCacheLimit = shape.boundarySignatureCacheLimit;
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

HopSearcher::HopSearcher(const farside::SignatureTable &signatures,
                         const exactnaming::ExactTables *exact, HopShape shape,
                         unsigned threads,
                         std::shared_ptr<exactnaming::TableCaches> exactCaches)
    : signatures_(&signatures), exact_(exact), shape_(shape), threads_(threads),
      exactCaches_(std::move(exactCaches)) {
  if (exact_ && !exactCaches_)
    exactCaches_ = std::make_shared<exactnaming::TableCaches>(*exact_);
}

HopSearcher::HopSearcher(const farside::SignatureTable *signatures,
                         const exactnaming::ExactTables *exact, SearchPolicy policy,
                         unsigned threads)
    : signatures_(signatures), exact_(exact), policy_(policy), threads_(threads) {
  if (exact_)
    exactCaches_ = std::make_shared<exactnaming::TableCaches>(*exact_);
}

HopRun HopSearcher::run(const farside::WitnessRedrawer &row,
                        const std::string &rowName, long long surfaceTarget,
                        double seconds,
                        const std::function<bool(const KeptSurface &)> &stop,
                        const SearchFrontier *resume) const {
  SearchRequest request;
  request.name = rowName;
  request.row = &row;
  request.shape = searchShape(shape_);
  request.surfaceTarget = surfaceTarget;
  request.seconds = seconds;
  // Always recorded (it costs one fingerprint): a later hop from this node
  // carries on from it instead of searching this prefix again.
  request.recordFrontier = true;
  request.resume = resume;
  request.stop = stop;
  return run(row.rowBuild(), request);
}

HopRun HopSearcher::run(const rowsearch::RowBuild &rb,
                        const SearchRequest &request) const {
  const auto wall0 = std::chrono::steady_clock::now();
  const double cpu0 = timers::processCpuSeconds();
  const SearchShape &shape = request.shape;
  const SweepInputs &sweep = request.sweep;
  const SearchOutputs &outputs = request.outputs;
  // Divergence 7: a witness per identity across the database, signed during
  // the search (verifyslicegenus); else a kept surface per keptKey().
  const bool signDuringSearch =
      policy_.signing == SearchPolicy::Signing::duringSearch;
  if (rb.seedFaces.empty() && !request.unseeded)
    throw std::runtime_error("hop: the row has no collar seed");
  if (signDuringSearch ? !sweep.record || !sweep.names : !request.row)
    throw std::logic_error("HopSearcher::run(): the request lacks what its "
                           "signing needs");
  if (policy_.judgeInSearch && (!sweep.bounds || !sweep.names))
    throw std::logic_error("HopSearcher::run(): judging needs bounds and names");

  // Declared before the search, which holds pointers to them. Every
  // boundary is named by its complement unless the row draws its far sides.
  const farside::ComplementNamer complementNamer{};
  std::optional<farside::DiagramNamer> namer;
  std::optional<SurfaceSearch> eOpt;
  if (rb.seedFaces.empty())
    eOpt.emplace(rb.tri);
  else
    eOpt.emplace(rb.tri, rb.seedFaces, rb.searchSideBC);
  SurfaceSearch &e = *eOpt;
  e.configureLimits(shape.limits);
  // The search's frontier: carried on from, and recorded.
  if (request.resume)
    e.setResumeFrontier(request.resume);
  e.setRecordFrontier(request.recordFrontier);
  if (request.pairSigCacheDir)
    e.setPairSigCacheDir(*request.pairSigCacheDir);
  e.setBoundaryNamer(complementNamer);
  if (signatures_ && request.diagramNaming) {
    if (policy_.failures == SearchPolicy::Failures::halt) {
      // Divergence 2: a row whose namer cannot be built stays on the
      // complement route.
      try {
        namer.emplace(rb.link.tri, rb.pdcode.size(), *rb.cob, *signatures_);
        if (exact_) namer->enableExactNames(*exact_, exactCaches_);
        e.setBoundaryNamer(*namer);
      } catch (const std::exception &ex) {
        std::cerr << "[!] " << request.name
                  << ": diagram naming off for this row (" << ex.what()
                  << ")\n";
      }
    } else {
      namer.emplace(rb.link.tri, rb.pdcode.size(), *rb.cob, *signatures_);
      if (exact_) namer->enableExactNames(*exact_, exactCaches_);
      e.setBoundaryNamer(*namer);
    }
  }

  if (!rb.seedFaces.empty()) {
    // The invariant that makes the search side fixed: no searchable
    // triangle other than the seed has an edge on it. Checked once here
    // rather than re-derived for every surface found.
    if (const size_t touching = e.countSearchableFacesTouching(rb.searchSideBC))
      throw SeedInvariantFailure(touching);
    // Its name is known by construction; never identify it.
    e.primeBoundaryName(rb.searchSideBC, rb.searchEdges, request.name);
  }

  // Fresh per search, so each census line describes one search.
  std::optional<SelfIntersectionCensus> selfIntersectionCensus;
  if (outputs.selfIntersectionCensus) {
    selfIntersectionCensus.emplace();
    selfIntersectionCensus->searchSideBoundary = static_cast<long>(rb.searchSideBC);
  }
  e.configureSelfIntersections(
      {.resolveUnlinked = shape.resolveUnlinked,
       .census = selfIntersectionCensus ? &*selfIntersectionCensus : nullptr});

  std::optional<rowsearch::SurfaceStatsTally> surfaceStats;
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

  HopRun out;
  out.recognitionBefore = identify::recognitionCacheStats();
  out.censusWritesBefore = census::insertCounts();

  // Every surface the drain describes lands in exactly one of its buckets
  // (see the accounting after the search).
  rowsearch::RowAccounting acct;

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
      if (rejectionSamplesTaken[reason]++ >= rowsearch::REJECTION_SAMPLES_PER_REASON)
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

  // Signing::deferred: one kept surface per keptKey().
  std::mutex keptMutex;
  std::unordered_set<std::string> keys;
  std::atomic<bool> stopped{false};
  // judgeInSearch: whether this search has found a constructive witness.
  std::atomic<bool> constructive{false};
  std::mutex fatalMutex;
  // Signing::duringSearch: the search's new witnesses are signed off the
  // drain threads, the pair-signature context built by the signer at the
  // first of them (WitnessSigner). Started just before the search.
  std::optional<WitnessSigner> signer;
  if (signDuringSearch)
    sweep.record->markActivity();

  // Started just before the search (below). It polls at 200ms but has no
  // access to SearchStats; onProgress is the only place the live count is
  // handed to us, so it publishes it there.
  std::optional<rowsearch::RowWatchdog> watchdog;
  SurfaceSearchCallbacks callbacks;
  // A stop nobody else noted (SIGINT) is not running out of candidates.
  callbacks.onInterrupted = [&] { noteStop("interrupted"); };
  callbacks.onProgress = [&](const SearchStats &stats) {
    if (watchdog) watchdog->publishSatisfying(stats.satisfyingCount);
    if (outputs.progress) rowsearch::printProgress(stats, e);
  };
  // The equalising rule, checked where each surface is counted, so a search
  // stops at the target rather than a progress tick later (see
  // SearchCallbacks::surfaceTarget). The watchdog below still checks it
  // too, as a backstop. As there, the drain is let finish unless
  // skipDrainOnTimeout.
  if (request.surfaceTarget) {
    callbacks.surfaceTarget = *request.surfaceTarget;
    callbacks.onSurfaceTarget = [&] {
      noteStop("surface-target");
      if (request.skipDrainOnTimeout)
        e.skipRemainingBoundaryProcessing();
    };
  }
  // For the `search profile:` line: petal-cache counters as root filtering
  // ends.
  callbacks.onRootsReady = [&] { out.petalsAtRoots = e.petalCacheStats(); };
  callbacks.onBoundaryProcessingStarted = [&](size_t total, unsigned threads) {
    out.drainTail = total;
    if (outputs.progress) {
      rowsearch::progressBlock.forget();
      std::cerr << "[+] boundary processing: " << total << " queued surfaces, "
                << threads << " threads\n";
    }
  };
  if (outputs.progress || signDuringSearch)
    callbacks.onBoundaryProcessingProgress =
        [&](size_t processed, size_t total,
            std::chrono::steady_clock::duration elapsed) {
          if (outputs.progress)
            rowsearch::printBoundaryProgress(
                processed, total, elapsed,
                constructive.load() ? std::optional<int>(sweep.literatureLo)
                                    : std::nullopt,
                sweep.literatureHi);
          // The drain is the long pole of a search and the phase most likely
          // to be interrupted, so checkpoint from here.
          if (signDuringSearch)
            sweep.record->checkpoint(/*force=*/false);
        };
  callbacks.onBoundaryProcessingComplete =
      [&](size_t total, std::chrono::steady_clock::duration elapsed) {
        out.drainTailTime = elapsed;
        out.drainTailSeconds = std::chrono::duration<double>(elapsed).count();
        if (outputs.progress) {
          rowsearch::progressBlock.forget();
          std::cerr << "[+] boundary processing: done (" << total
                    << " processed in " << formatElapsed(elapsed) << ")\n";
        }
      };

  callbacks.onSurfaceBoundaryProcessed = [&](const SurfaceBoundaryInfo &info) {
    // Recorded before any of the filtering below: the question this answers
    // is what the SEARCH found, not what survived the checks that decide
    // whether a surface bounds this particular row.
    if (surfaceStats)
      surfaceStats->record(rowsearch::SurfaceStatsKey{
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

    // Orientable, search side intact, the row's own oriented variant, one
    // far side (rowsearch::gateSurface()). The orientation check needs the
    // search side's oriented curves, and an exact oriented far-side name the
    // rest, so the gate captures them once.
    const rowsearch::GatedSurface g = rowsearch::gateSurface(info, rb);
    if (!g.accepted()) {
      acct.reject(g.gate);
      sampleRejection(rowsearch::gateReason(g.gate), info);
      return;
    }
    // Signing::deferred keeps the oriented outgoing link: oriented by the
    // gate's own judgement of the incoming curves (g.flips: the row is rb's,
    // request.row->row() is *rb.orientation).
    std::optional<farside::OutgoingLink> link;
    if (!signDuringSearch) {
      link = farside::orientedOutgoingLink(g.orientedLinks, g.surfaceOf,
                                           request.row->outgoing(), g.flips,
                                           rb.searchSideBC);
      if (!link) {
        // Only a surface component off the row, which the gate's flips (one
        // per component meeting the row) rule out: impossible.
        acct.reject(rowsearch::Gate::orientationBroken);
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
    cobordismgraph::Witness w;
    w.subject = request.name;
    w.subjectComponents = rb.componentCount;
    w.genus = info.tubedGenus;
    w.tubed = !info.connected;
    w.resolvedVertices = info.resolvedVertices;
    if (signDuringSearch) {
      // The witness as the database records it: its provenance now.
      w.sourceRow = request.name;
      w.thickenLayers = sweep.thickenLayers;
      w.maxFaces = shape.maxFaces.value_or(0);
    }
    std::string farName;
    if (g.split.otherSides.empty()) {
      w.kind = cobordismgraph::WitnessKind::direct;
    } else {
      // Exactly one: the gate turns away more (multi-far-side). A
      // genuinely-linked far side is recorded but, unless it is a knot or a
      // proven unlink, it will not carry a bound: a complement does not
      // determine a link. It is kept because the observation is real and is
      // exactly what a later per-witness naming needs as input; the solver's
      // farSideBearsBound() is what declines it.
      const cobordismgraph::BoundarySide &far = g.split.otherSides.front();
      // Normalized, so a census hit and the table name are one graph node,
      // and oriented where exact names are on (rowsearch::farSideName()).
      farName = rowsearch::farSideName(g, namer ? &*namer : nullptr);
      w.kind = cobordismgraph::WitnessKind::cobordism;
      w.other = farName;
      w.otherComponents = far.components;
      // Re-derived from the name, as every reader of other_candidates does.
      if (signDuringSearch)
        w.otherCandidates = sweep.names->candidates(farName, far.components);
    }

    if (!signDuringSearch) {
      std::string key = rowsearch::keptKey(w, *link, *request.row);
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
      if (request.stop && !stopped.load() && request.stop(k)) {
        stopped.store(true);
        noteStop("stopped");
        e.requestStop();
        e.skipRemainingBoundaryProcessing();
      }
      out.kept.push_back(std::move(k));
      return;
    }

    // Signing::duringSearch. The dedup matters a lot, since every search
    // harvests: a single search reports thousands of near-identical
    // surfaces, and signing each would dominate the run. The signer takes the witness's
    // faces and signs it off this (drain) thread; its identity is claimed
    // now, so the dedup stays exact.
    if (!sweep.record->claim(w)) {
      acct.duplicate.fetch_add(1, std::memory_order_relaxed);
      return;
    }
    signer->add(w, info.captureFaces());
    acct.recorded.fetch_add(1, std::memory_order_relaxed);
    if (!policy_.judgeInSearch)
      return;

    // Does this witness alone already settle the row? Checked cheaply
    // against the bounds as of the last full solve rather than by
    // re-running the solver here: the far side's own bound can only have
    // improved since, so this under-reports at worst, and the next solve
    // picks up anything it missed. Whether the bound leans on a literature
    // value decides whether the row may be considered settled: stopping the
    // search on an assisted one would forfeit the chance of finding the
    // surface that upgrades it.
    const cobordismgraph::UpperBoundVia via =
        cobordismgraph::upperBoundVia(w, *sweep.bounds, *sweep.names);
    const int implied = via.genus;
    if (implied != cobordismgraph::NO_UPPER_BOUND && implied < sweep.literatureLo) {
      std::ostringstream msg;
      msg << request.name << ": witness implies an upper bound of " << implied
          << ", BELOW the literature lower bound " << sweep.literatureLo << ".";
      {
        std::lock_guard<std::mutex> lock(fatalMutex);
        if (out.fatal.empty())
          out.fatal = msg.str();
      }
      noteStop("fatal-bug");
      e.requestStop();
      e.skipRemainingBoundaryProcessing();
      return;
    }
    if (implied != cobordismgraph::NO_UPPER_BOUND && implied <= sweep.literatureLo &&
        !via.assisted) {
      // Only a CONSTRUCTIVE result settles a row. An assisted one is a
      // correct deduction but not an independent verification.
      // Announce on stdout, once, the first time this row resolves.
      // Previously the only sign was "-- ACHIEVED" appearing in the
      // redrawn stderr progress block, which is invisible to any log filter
      // and vanishes as soon as the next block overwrites it -- so a row
      // could sit verified-but-unwritten for hours with nothing in the log
      // to say so. Checkpoint immediately too: this is the single most
      // valuable moment in a row, and every search harvests, so the row may
      // keep running for hours afterwards.
      if (!constructive.exchange(true, std::memory_order_relaxed)) {
        std::cout << "[+] " << request.name
                  << ": CONSTRUCTIVE witness found -- reaches genus " << implied
                  << " (literature [" << sweep.literatureLo << ", "
                  << sweep.literatureHi << "]). Checkpointing now.\n"
                  << std::flush;
        sweep.record->checkpoint(/*force=*/true);
      }
    }
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
  // what it found. skipDrainOnTimeout restores the old behaviour for when
  // throughput genuinely matters more.
  {
    rowsearch::WatchdogLimits watchdogLimits{
        .surfaceTarget = request.surfaceTarget,
        .rowSeconds = request.seconds,
        .sweepSeconds = request.sweepSeconds,
        .sweepStart = request.sweepStart,
        // Quiescence: this search has stopped teaching us anything new, so
        // spending the rest of its budget enumerating more of the same is
        // worse than moving on to one we know nothing about.
        .quiescenceSeconds = request.quiescenceSeconds};
    if (sweep.record)
      watchdogLimits.idleMillis = [&] { return sweep.record->millisSinceNew(); };
    watchdog.emplace(std::move(watchdogLimits), [&](const char *why) {
      noteStop(why);
      e.requestStop();
      if (request.skipDrainOnTimeout)
        e.skipRemainingBoundaryProcessing();
    });
  }
  if (signDuringSearch)
    signer.emplace(
        [&e]() -> const PairSigContext<4, 2> & { return e.pairSigContext(); },
        [&](cobordismgraph::Witness &&w) { sweep.record->publish(std::move(w)); });

  const auto searchStart = std::chrono::steady_clock::now();
  out.setup = std::chrono::duration<double>(searchStart - wall0).count();
  const SearchStats stats = e.search(
      threads_, shape.condition, callbacks, shape.iddfsIterations,
      shape.iddfsStep, shape.iddfsStart, shape.iddfsFinalThreads,
      /*orientableOnly=*/true, shape.maxFaces, shape.rootBudgetStart,
      shape.rootBudgetGrowth);
  out.search = std::chrono::duration<double>(std::chrono::steady_clock::now() - searchStart).count();
  for (auto r : stats.profile.rounds) out.rounds.push_back(std::chrono::duration<double>(r).count());
  if (signer) {
    // Every surface is described; sign what is still queued, with every
    // thread, before anything reads or writes this search's witnesses.
    signer->finish(threads_);
    out.pairSigsSigned = signer->signedCount();
    out.pairSigMillis = signer->signMillis();
    out.pairSigContextSeconds = signer->contextSeconds();
    out.pairSigContextLoaded = e.pairSigContextLoaded();
    out.pairSigFinishSeconds = signer->finishSeconds();
    rowsearch::printPairSignatures(std::cout, request.name, out);
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
    const farside::NamingStats &ns = namer->stats();
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
  out.otherOrientation = acct.orientation.load();
  out.drainSkipped = drainSkipped;
  // Surfaces were accepted, yet not one reached the record. That can be
  // genuine (every one witnesses another oriented variant), but it is also
  // exactly what a broken gate looks like, so it never licenses a negative.
  out.nothingExamined = acct.nothingExamined();

  if (surfaceStats)
    rowsearch::appendSurfaceStats(*outputs.surfaceStats, request.name,
                                  shape.maxFaces.value_or(0), surfaceStats->take());
  if (selfIntersectionCensus)
    rowsearch::appendSelfIntersectionCensus(*outputs.selfIntersectionCensus,
                                            request.name, shape.maxFaces.value_or(0),
                                            shape.resolveUnlinked, stats,
                                            *selfIntersectionCensus);
  if (outputs.progress)
    rowsearch::progressBlock.forget();

  // Divergence 2: a search that cannot account for its surfaces vouches for
  // no negative.
  if (policy_.failures == SearchPolicy::Failures::halt &&
      (!out.accountingFailure.empty() || out.nothingExamined))
    noteStop("unaccounted");
  out.outcome = outcome;
  out.resumed = e.resumedFrontier();
  out.resumeRefusal = e.resumeRefusal();
  if (request.resume)
    out.resumeOfferedRuns = request.resume->runs;
  out.recordedFrontier = e.frontier();
  out.frontierSeconds = e.frontierSeconds();
  // Only a prefix whose every surface was examined may be skipped later
  // (divergence 1: and, by verifyslicegenus's rule, only when something was).
  if (out.accountingFailure.empty() && !drainSkipped &&
      !(policy_.frontierNeedsExamined && out.nothingExamined))
    out.frontier = e.frontier();
  out.stats = stats;
  out.petals = e.petalCacheStats();
  out.boundaryCache = e.boundarySignatureCacheStats();
  out.recognitionAfter = identify::recognitionCacheStats();
  out.censusWritesAfter = census::insertCounts();
  out.linkingAudit = linkingnumber::auditLinkingNumbers.load();
  out.wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();
  out.cpu = timers::processCpuSeconds() - cpu0;
  return out;
}

} // namespace cascade
