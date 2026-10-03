//
//  scheduler.h
//
//  Which search next: a goal-directed run's scheduler.
//

#ifndef SURFER_COBOUND_SCHEDULER_H
#define SURFER_COBOUND_SCHEDULER_H

#include <string>
#include <vector>

#include "cobound/search/search.h"

/*! \file utils/surfer/cobound/driver/scheduler.h
 *  \brief Which search next, with a goal: from the target outwards, the
 *  links that could still move the target's bound, each loaded with the
 *  database's cobordisms first (bounds/databasecobordisms.h), then
 *  searched, until the goal has a proof (written as a certificate,
 *  bounds/certificate.h) or the limits are spent. Its run directory's
 *  records are driver/runrecords.h's.
 *
 *  Exit codes: 0 met, 1 not met (a first SIGINT or SIGTERM included), 2
 *  halted (an impossible state), 3 a contradiction -- even when the goal is
 *  met (plan divergences 2, 6, 8).
 */

namespace cascade {

/// A goal run's options.
struct GoalOptions {
  std::string targetPD, targetName, work, verify, knotTable, linkTable, knotSymmetry, censusDb;
  int goalGenus = 0;
  bool goalDisjoint = false; // goal on the singleton partition (disjoint discs)
  long hopSurfaces = 100000, maxHopSurfaces = 1000000;
  int threads = 10;
  int maxExpansions = 20;
  double cpuBudget = 3600.0 * 4;
  size_t maxCrossings = 24;
  std::string strategy = "best";
  bool literature = true;
  /// The atlas's master cobordisms.csv (read only): a node that IS a table
  /// entry the atlas searched gets that row's witnesses as free edges.
  std::string masterWitnesses;
  /// "process": each hop searched in this process (search.h); "child": a
  /// verifyslicegenus row per hop, its witnesses read back from their pair
  /// signatures. Both search the same thickening with the same shape.
  std::string hopMode = "process";
  /// Each hop's search shape: the campaign's (hosts.conf) unless the
  /// --hop-* options change it. Either hop mode uses it.
  HopShape hopShape;
  bool verbose = false;
  /// The atlas-format witness store every kept surface is recorded in
  /// (pending.h); empty for none. Read-only stores to deduplicate against
  /// (the master, say), and the run's name for cascade: subjects.
  std::string witnessStore, runName;
  std::string pairSigCache; ///< --pair-sig-cache: stored pair-signature contexts
  /// --read-back-cache: master witnesses' read-backs kept across runs
  /// (RowReadBacks); empty for none.
  std::string readBackCache;
  std::vector<std::string> dedupeAgainst;
  /// Write lower_report.jsonl at the end; the sources file (atlas
  /// data/lower_bound_sources.csv: name,...,special) marks which literature
  /// bounds are not Lipschitz.
  bool lowerReport = false;
  std::string lowerSources;
  /// --goal-lower G: also stop once lower(target, goal partition) >= G
  /// (README.md, "Lower-bound mode"); -1 for no lower goal. Needs
  /// --lower-sources.
  int goalLower = -1;
  /// A node kept only for the lower goal is never expanded above this many
  /// crossings: the chain must come back to a table entry.
  size_t lowerMaxCrossings = 16;
  /// --master-loads lazy|eager: a node's stored rows loaded when it is
  /// about to be expanded (the default since 2026-09-30), or as soon as it
  /// is met and useful (every table node the graph reaches).
  bool lazyMasterLoads = true;
  /// Hub breadth (John, 2026-09-29): a node with at least hubDegree witness
  /// edges is expanded once at hubSurfaces or more. 0: off.
  size_t hubDegree = 0;
  long hubSurfaces = 0;
};

/**
 * Runs `options`' goal-directed run in its work directory and returns its
 * exit code (above). Writes, under the work directory: hop_<k>_n<node>/
 * (each search's log, frontier and pending file), cascade.jsonl,
 * profiles.jsonl, node_bounds.jsonl, nodes.csv (with a witness store),
 * lower_report.jsonl (lowerReport), and certificate.json or
 * lower_certificate.json when the goal is met.
 * \throws std::exception for a run that cannot start (a split target, a
 *         target PD that is not its named entry, unreadable tables).
 */
int runToGoal(const GoalOptions &options);

} // namespace cascade

#endif // SURFER_COBOUND_SCHEDULER_H
