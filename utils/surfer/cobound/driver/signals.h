//
//  signals.h
//
//  The process's one policy for SIGINT and SIGTERM.
//

#ifndef SURFER_COBOUND_SIGNALS_H
#define SURFER_COBOUND_SIGNALS_H

/*! \file utils/surfer/cobound/driver/signals.h
 *  \brief How a run that searches answers SIGINT and SIGTERM (plan divergence
 *  8): one process-level policy, the same for a run with a goal and without.
 *
 *  - The first signal ends the search that is running, cleanly: it stops
 *    within a second (its search polls interrupted()), its drain finishes,
 *    its pending cobordisms are fsynced, its frontier is ruled on as any
 *    stopped search's, and it records the outcome `interrupted`. The driver
 *    then starts no further search, finishes its run as it would at any
 *    other end, and exits with what the run was for: 0 without a goal (the
 *    search is thinner, not failed), 1 with one unless it is met.
 *  - A second signal ends the process at once: _Exit(128 + the signal).
 *
 *  SIGTERM is what stopping a scope, `timeout` or a manual kill sends. A run
 *  over remote_run.sh's MemoryMax is killed by the cgroup's OOM killer with
 *  SIGKILL, which nothing catches: the pending file's fsyncs are what bound
 *  that loss. The library's own per-search SIGINT handling (SigintScope) is
 *  turned off for every search a cobound driver runs
 *  (SubmanifoldSearch::setSigintHandling()), so this is the only handler.
 */

namespace runsignals {

/** Installs the policy for SIGINT and SIGTERM. Call once, from main, before
 *  any search. */
void install();

/** Whether a first signal has arrived. Safe from any thread. */
bool interrupted();

/** Which signal it was (SIGINT or SIGTERM); 0 before one arrives. */
int which();

/** "SIGINT", "SIGTERM" or "" (for log lines). */
const char *name();

} // namespace runsignals

#endif // SURFER_COBOUND_SIGNALS_H
