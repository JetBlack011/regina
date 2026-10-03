//
//  setup.h
//
//  What a cobound run sets once per process, before its first naming.
//

#ifndef SURFER_COBOUND_SETUP_H
#define SURFER_COBOUND_SETUP_H

#include "cobound/driver/config.h"

/*! \file utils/surfer/cobound/driver/setup.h
 *  \brief Once per process (plan, "Startup per process"): the settings that
 *  are process-wide -- the census and whether it is written, Pachner
 *  searches on a census miss and their time budget, the complement cache's
 *  size, the linking audit -- set from the config before the first naming,
 *  so a run never inherits another driver's (the retired cascadesearch set
 *  its own in its run(); verifyslicegenus in its main).
 *
 *  The tables are loaded once: a run's exact tables (the depth-0 graph's,
 *  or a goal run's namer's), with diagram naming's signature table derived
 *  from them (linknaming::SignatureTable::fromTables(const ExactTables &)).
 *  What is NOT loaded unless asked for: a database's cobordisms as a goal
 *  run's free edges (master_witnesses, a goal key: never without a goal),
 *  and link classes (computed lazily per base by a goal run's namer).
 */

namespace setup {

/// Applies a run's (with a goal or without) process-wide settings from its
/// config. Returns whether the census opened.
bool applyRunSettings(const config::Config &cfg);

} // namespace setup

#endif // SURFER_COBOUND_SETUP_H
