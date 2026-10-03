//
//  scheduler.h
//
//  Which search next: a goal-directed run's scheduler.
//

#ifndef SURFER_COBOUND_SCHEDULER_H
#define SURFER_COBOUND_SCHEDULER_H

#include <string>
#include <vector>

#include "cobound/driver/config.h"
#include "cobound/search/search.h"

/*! \file utils/surfer/cobound/driver/scheduler.h
 *  \brief Which search next. A run's targets are searched by one scheduler:
 *  without a goal each target once, in order (driver/sweep.h: the old sweep,
 *  at depth 0); with a goal, from the target outwards (here): the links that
 *  could still move the target's bound, each loaded with the database's
 *  cobordisms first (bounds/databasecobordisms.h), then searched, until the
 *  goal has a proof (written as a certificate, bounds/certificate.h) or the
 *  limits are spent. Its run directory's records are driver/runrecords.h's.
 *
 *  Exit codes, with a goal: 0 met, 1 not met (a first SIGINT or SIGTERM
 *  included), 2 halted (an impossible state), 3 a contradiction -- even
 *  when the goal is met (plan divergences 2, 6, 8).
 */

namespace cascade {

/// A goal run's options: every value from the config (goalOptions()).
struct GoalOptions {
  std::string targetPD, targetName, work, knotTable, linkTable, knotSymmetry, censusDb;
  /// Whether the census opened (setup::applyRunSettings()).
  bool censusLoaded = false;
  int goalGenus = 0;
  bool goalDisjoint = false; // goal on the singleton partition (disjoint discs)
  long hopSurfaces = 0, maxHopSurfaces = 0;
  int threads = 1;
  int maxExpansions = 0;
  double cpuBudget = 0;
  /// Each search's wall-clock backstop.
  double searchSeconds = 0;
  size_t maxCrossings = 0;
  std::string strategy;
  bool literature = true;
  /// A read-only database (the atlas's master): a node that IS a table
  /// entry the database holds rows of gets those cobordisms as free edges.
  std::string masterWitnesses;
  /// Each search's shape.
  HopShape hopShape;
  /// The complement cache's limit, as the profile line prints it (it is
  /// process-wide: setup::applyRunSettings() sets it).
  size_t complementCacheLimit = 0;
  /// The database every kept surface is signed into at the run's end
  /// (pending.h); empty for none. Read-only stores to deduplicate against
  /// (the master, say), and the run's name for cascade: subjects.
  std::string witnessStore, runName;
  std::string pairSigCache; ///< stored pair-signature contexts
  /// The database cobordisms' read-backs kept across runs (RowReadBacks);
  /// empty for none.
  std::string readBackCache;
  std::vector<std::string> dedupeAgainst;
  /// Write lower_report.jsonl at the end; the sources file (atlas
  /// data/lower_bound_sources.csv: name,...,special) marks which literature
  /// bounds are not Lipschitz.
  bool lowerReport = false;
  std::string lowerSources;
  /// Also stop once lower(target, goal partition) >= goalLower (README.md,
  /// "Lower-bound mode"); -1 for no lower goal. Needs lowerSources.
  int goalLower = -1;
  /// A node kept only for the lower goal is never expanded above this many
  /// crossings: the chain must come back to a table entry.
  size_t lowerMaxCrossings = 0;
  /// Hub breadth (John, 2026-09-29): a node with at least hubDegree witness
  /// edges is expanded once at hubSurfaces or more. 0: off.
  size_t hubDegree = 0;
  long hubSurfaces = 0;
};

/// A goal run's options from its config (context goal).
/// \exception config::Error a configuration a goal run cannot run with.
GoalOptions goalOptions(const config::Config &cfg);

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
