//
//  sweep.h
//
//  A run without a goal: each target searched once, in order.
//

#ifndef SURFER_COBOUND_SWEEP_H
#define SURFER_COBOUND_SWEEP_H

#include "cobound/driver/config.h"

/*! \file utils/surfer/cobound/driver/sweep.h
 *  \brief The scheduler without a goal (depth 0, the old sweep): each
 *  target -- the rows of `targets` within max_crossings, in crossing order,
 *  or the one diagram `target_pd` -- searched once, each by the one search
 *  (search/search.h), judged by its own cobordism graph (SearchJudge), its
 *  finds written to its pending file (<work>/hop_<k>_n0/kept.csv) and all
 *  of them signed into `cobordisms` at the run's end. Each searched row's
 *  search record goes to the verdicts file; its status and bounds are
 *  `solve`'s (plan divergence 10). The log lines are verifyslicegenus's,
 *  every one (frozen).
 */

namespace sweep {

/**
 * Runs a run without a goal from its config (context run). Returns its exit
 * code: 0 (a first SIGINT or SIGTERM included: the run searches no further
 * row, records `interrupted` and ends as usual), 1 when its tables or
 * symmetry types cannot be read, or 2 when any search failed its accounting
 * (after the run's other rows) or the run halted (fatal.h).
 * \throws config::Error for a configuration a search cannot run with.
 */
int run(const config::Config &cfg);

} // namespace sweep

#endif // SURFER_COBOUND_SWEEP_H
