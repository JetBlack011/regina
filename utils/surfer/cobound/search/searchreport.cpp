//
//  searchreport.cpp
//

#include "cobound/search/searchreport.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include "surfer/report/csvwriter.h"

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

} // namespace rowsearch
