//
//  searchreport.cpp
//

#include "cobound/search/searchreport.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "linknaming/census/censusnaming.h"
#include "linknaming/complement/complementcache.h"
#include "surfer/report/csvwriter.h"
#include "surfer/submanifold/linkingnumber.h"

namespace rowsearch {

report::RollingReport progressBlock;

void appendSurfaceStats(const std::filesystem::path &path,
                        const std::string &rowName, long long maxFaces,
                        const std::map<SurfaceStatsKey, long long> &counts) {
  if (counts.empty())
    return;
  const bool needHeader = !std::filesystem::exists(path);
  std::ofstream out(path, std::ios::app);
  if (!out) {
    std::cerr << "[!] could not open " << path << " for --surface-stats\n";
    return;
  }
  if (needHeader)
    out << "row,max_faces,triangles,orientable,genus,punctures,tubed_genus,"
           "closed_components,connected,count\n";
  for (const auto &[key, n] : counts)
    out << csvField(rowName) << ',' << maxFaces << ',' << key.triangles << ','
        << (key.orientable ? "true" : "false") << ',' << key.genus << ','
        << key.punctures << ',' << key.tubedGenus << ','
        << key.closedComponents << ',' << (key.connected ? "true" : "false")
        << ',' << n << '\n';
}

void appendSelfIntersectionCensus(const std::filesystem::path &path,
                                  const std::string &rowName,
                                  long long maxFaces, bool resolveUnlinked,
                                  const SearchStats &stats,
                                  SelfIntersectionCensus &census) {
  const bool needHeader = !std::filesystem::exists(path);
  std::ofstream out(path, std::ios::app);
  if (!out) {
    std::cerr << "[!] could not open " << path
              << " for --self-intersection-census\n";
    return;
  }
  if (needHeader)
    out << "row,max_faces,resolve_unlinked,satisfying,embedded,resolved,"
           "singular,interior_unlinked,interior_uncertified,"
           "boundary_unlinked,boundary_uncertified,multi_open,"
           "multi_open_search_side,multi_open_far,multi_open_far_clean,"
           "multi_open_far_simple,far_configs,far_clean_configs,"
           "configs_saturated,audited,audit_knotted,knotted_pairsigs\n";
  size_t farConfigs, farCleanConfigs;
  bool saturated;
  {
    std::lock_guard<std::mutex> lock(census.configsMutex);
    farConfigs = census.farConfigs.size();
    farCleanConfigs = census.farCleanConfigs.size();
    saturated = census.configsSaturated;
  }
  std::string hits;
  {
    std::lock_guard<std::mutex> lock(census.hitsMutex);
    for (size_t i = 0; i < census.knottedPairSigs.size(); ++i)
      hits += (i ? ";" : "") + census.knottedPairSigs[i];
  }
  auto get = [](const std::atomic<long long> &a) {
    return a.load(std::memory_order_relaxed);
  };
  out << csvField(rowName) << ',' << maxFaces << ','
      << (resolveUnlinked ? "true" : "false") << ',' << stats.satisfyingCount
      << ',' << stats.embeddedCount << ',' << stats.resolvedCount << ','
      << get(census.singular) << ',' << get(census.interiorUnlinked) << ','
      << get(census.interiorUncertified) << ','
      << get(census.boundaryUnlinked) << ','
      << get(census.boundaryUncertified) << ',' << get(census.multiOpen)
      << ',' << get(census.multiOpenSearchSide) << ','
      << get(census.multiOpenFar) << ',' << get(census.multiOpenFarClean)
      << ',' << get(census.multiOpenFarSimple) << ',' << farConfigs << ','
      << farCleanConfigs << ',' << (saturated ? "true" : "false") << ','
      << get(census.audited) << ',' << get(census.auditKnotted) << ','
      << csvField(hits) << '\n';
}

void printProgress(const SearchStats &stats, SurfaceSearch &e) {
  std::ostringstream report;
  report << "[+] elapsed: " << formatElapsed(stats.elapsed)
         << " | candidates examined: " << stats.foundCount
         << " | embedded surfaces found: " << stats.embeddedCount
         << " | satisfying boundary condition: " << stats.satisfyingCount
         << "\n";
  // Which iterative-deepening round we are in, and how deep finds have
  // actually gone. Without this a run looks like it is exploring up to
  // --max-faces when it may never have finished round 1 -- the calibration
  // row on L6a1{1} spent its whole 600s budget inside round 1 (cap 5) and
  // never started rounds 2-4, which is invisible from candidate counts alone.
  report << "[+] iddfs round " << stats.iddfsRound << "/"
         << stats.iddfsTotalRounds;
  if (stats.iddfsCapped)
    report << " (cap " << stats.iddfsCap << " faces)";
  else
    report << " (final, uncapped)";
  report << " | deepest satisfying find: " << stats.largestSatisfying
         << " faces";
  // Both numbers matter: the total says how big the surface is, the
  // seed-relative count says how far the search actually reached, and the
  // collar makes those differ by a couple of orders of magnitude.
  if (stats.seedFaces > 0 && stats.largestSatisfying >= stats.seedFaces)
    report << " (+" << (stats.largestSatisfying - stats.seedFaces)
           << " beyond the " << stats.seedFaces << "-face seed)";
  if (stats.rootBudget > 0)
    report << " | root budget " << stats.rootBudget;
  report << "\n";
  // rootsExhausted, not rootsCompleted: with per-root budgets a root is
  // re-walked once per pass, so rootsCompleted counts visits and can exceed
  // the root count. rootsExhausted counts each root at most once, so this is
  // the round's real progress.
  report << "[+] roots exhausted this round: " << stats.rootsExhausted << "/"
         << stats.rootsPerPass << " (visits " << stats.rootsCompleted << ")\n";
  report << "[+] surface homeomorphism types found so far: "
         << e.surfaceTypeTally().summary() << "\n";
  progressBlock.draw(report.str());
}

void printBoundaryProgress(size_t processed, size_t total,
                           std::chrono::steady_clock::duration elapsed,
                           std::optional<int> resolvedGenus, int target) {
  double elapsedSec = std::chrono::duration<double>(elapsed).count();
  double rate = (elapsedSec > 0.0) ? static_cast<double>(processed) / elapsedSec
                                   : 0.0;

  std::ostringstream report;
  report << "[+] boundary processing: elapsed " << formatElapsed(elapsed)
         << " | processed " << processed << "/" << total;
  if (rate > 0.0 && processed < total) {
    auto eta = std::chrono::duration_cast<std::chrono::steady_clock::duration>(
        std::chrono::duration<double>(
            static_cast<double>(total - processed) / rate));
    report << " | ETA " << formatElapsed(eta);
  }
  report << "\n";
  report << "[+] slice genus: "
         << (resolvedGenus ? std::to_string(*resolvedGenus) : "?") << "/"
         << target << " (literature target)"
         << (resolvedGenus ? " -- ACHIEVED" : "") << "\n";
  progressBlock.draw(report.str());
}

void printPairSignatures(std::ostream &out, const std::string &name,
                         const cascade::WitnessSigner &signer, bool contextLoaded) {
  out << "[+] " << name << ": pair signatures: " << signer.signedCount()
      << " signed in " << std::fixed << std::setprecision(1)
      << signer.signMillis() / 1000.0 << "s of thread time; context ";
  if (signer.signedCount() > 0)
    out << (contextLoaded ? "loaded from the cache, " : "built, ") << "ready "
        << signer.contextSeconds() << "s after the search began";
  else
    out << "not needed";
  out << "; " << signer.finishSeconds() << "s after the drain finishing\n"
      << std::defaultfloat;
}

void printSweepBreadth(std::ostream &out, const std::string &name,
                       const cascade::HopRun &run, const SearchFrontier *resumeFrom) {
  out << "[+] " << name << ": breadth: ";
  if (const auto &f = run.recordedFrontier)
    out << f->summary() << "; fingerprint " << f->fingerprint.substr(0, 12);
  else
    out << "not recorded";
  out << "; resumed "
      << (!resumeFrom ? std::string("none")
          : run.resumed
              ? std::string("yes (") + std::to_string(resumeFrom->runs) + " runs before)"
              : "no: " + run.resumeRefusal)
      << "; frontier " << std::fixed << std::setprecision(2) << run.frontierSeconds
      << "s, replayed " << run.stats.profile.replayed << " re-adds\n"
      << std::defaultfloat;
}

void printOutcome(std::ostream &out, const std::string &name, const cascade::HopRun &run) {
  out << "[+] " << name << ": " << run.recorded << " new witnesses, outcome "
      << run.outcome;
  if (run.otherOrientation > 0)
    out << ", " << run.otherOrientation << " surfaces rejected on orientation mismatch";
  out << "\n";
  // Its own line, parsed by tools/orchestrate/dispatch.py (RE_ACCOUNTING);
  // the summary line above stays exactly as RE_OUTCOME expects.
  out << "[+] " << name << ": accounting: " << run.accounting << "\n";
}

void printIdentification(std::ostream &out, const std::string &name,
                         const cascade::HopRun &run) {
  const identify::RecognitionCacheStats r = identify::recognitionCacheStats();
  const identify::RecognitionCacheStats &before = run.recognitionBefore;
  const identify::BoundarySignatureCacheStats &b = run.boundaryCache;
  auto secs = [](long long ms) {
    std::ostringstream o;
    o << std::fixed << std::setprecision(1) << ms / 1000.0;
    return o.str();
  };
  const auto censusWritesNow = census::insertCounts();
  const long long censusOk = censusWritesNow.first - run.censusWritesBefore.first;
  const long long censusFailed = censusWritesNow.second - run.censusWritesBefore.second;
  out << "[+] " << name << ": identification: boundary cache " << b.hits << "/"
      << b.checks << " hits, census checks " << (r.censusChecks - before.censusChecks)
      << " (local hits " << (r.localCensusHits - before.localCensusHits)
      << "), Pachner knots " << (r.pachnerKnots.attempts - before.pachnerKnots.attempts)
      << " tried/" << (r.pachnerKnots.successes - before.pachnerKnots.successes)
      << " named/"
      << secs(r.pachnerKnots.milliseconds - before.pachnerKnots.milliseconds)
      << "s, links " << (r.pachnerLinks.attempts - before.pachnerLinks.attempts) << "/"
      << (r.pachnerLinks.successes - before.pachnerLinks.successes) << "/"
      << secs(r.pachnerLinks.milliseconds - before.pachnerLinks.milliseconds)
      << "s, pairsigs " << run.pairSigsSigned << "/" << secs(run.pairSigMillis)
      << "s, census writes " << censusOk << " ok/" << censusFailed << " failed\n";
  if (run.diagramNamed) {
    out << "[+] " << name << ": diagram naming: " << run.naming << "\n";
    // Harmless to the names (each went to the complement route), but each
    // is a drawer defect that must be found.
    if (run.nonPlanar > 0)
      out << "[!] " << name << ": WARNING: " << run.nonPlanar
          << " far-side drawings were not planar diagrams (drawer "
             "defect; named by the complement route instead)\n";
  }
  // A census that cannot be written to costs nothing in correctness, but
  // every name it fails to keep is recomputed by every later row.
  if (censusFailed > 0)
    out << "[!] " << name << ": WARNING: " << censusFailed
        << " census writes failed (names found here will not "
           "reach later rows)\n";
}

void printSearchProfile(std::ostream &out, const std::string &name,
                        const cascade::HopRun &run) {
  // Where the search's time went. Measurement only; parsed by
  // cobound/tests/bench_search.sh.
  const SearchStats::Profile &p = run.stats.profile;
  const PetalCache::Stats &petals = run.petals;
  const PetalCache::Stats atRoots = run.petalsAtRoots.value_or(petals);
  auto dsecs = [](std::chrono::steady_clock::duration d) {
    std::ostringstream o;
    o << std::fixed << std::setprecision(1) << std::chrono::duration<double>(d).count();
    return o.str();
  };
  auto nsecs = [](long long nanos) {
    std::ostringstream o;
    o << std::fixed << std::setprecision(1) << nanos / 1e9;
    return o.str();
  };
  auto unknotMisses = [](const PetalCache::Stats &s) {
    return s.unknotChecks - s.unknotCacheHits;
  };
  auto linkingMisses = [](const PetalCache::Stats &s) {
    return s.linkingChecks - s.linkingCacheHits;
  };
  out << "[+] " << name << ": search profile: prototype " << dsecs(p.prototype)
      << "s (unknot misses " << unknotMisses(atRoots) << " in "
      << nsecs(atRoots.unknotMissNanos) << "s, linking misses " << linkingMisses(atRoots)
      << " in " << nsecs(atRoots.linkingMissNanos) << "s); rounds";
  for (auto round : p.rounds)
    out << " " << dsecs(round) << "s";
  out << "; drain tail " << run.drainTail << " surfaces in " << dsecs(run.drainTailTime)
      << "s; nodes " << p.nodes << ", attempts " << p.attempts << ", evaluated "
      << p.evaluated << ", charged " << p.charged << ", replayed " << p.replayed
      << "; petal misses: unknot " << unknotMisses(petals) << " in "
      << nsecs(petals.unknotMissNanos) << "s, linking " << linkingMisses(petals) << " in "
      << nsecs(petals.linkingMissNanos) << "s (cochains " << petals.linkingFast
      << ", fallbacks " << petals.linkingFallbacks << ")";
  if (linkingnumber::auditLinkingNumbers.load())
    out << "; linking audit: " << petals.linkingAudited << " checked ("
        << petals.linkingAuditNonzero << " linked), " << petals.linkingDisagreements
        << " disagree, drilling route " << nsecs(petals.linkingAuditOldNanos) << "s";
  out << "\n";
}

} // namespace rowsearch
