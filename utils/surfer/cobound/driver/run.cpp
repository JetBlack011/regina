//
//  run.cpp
//
//  cobound run: a search from each target, scheduled by the goal.
//

#include <filesystem>
#include <iostream>
#include <optional>

#include "cobound/driver/commands.h"
#include "cobound/driver/config.h"
#include "cobound/driver/scheduler.h"
#include "cobound/driver/setup.h"
#include "cobound/driver/signals.h"
#include "cobound/driver/sweep.h"
#include "surfer/report/atomicwrite.h"

// `run` (was verifyslicegenus's search runs and cascadesearch): the config
// (driver/config.h) names the targets, the work directory and, optionally, a
// goal. Without one each target is searched once (sweep.h); with one the
// scheduler searches outwards from the target until the goal has a proof or
// the limits are spent (scheduler.h). The configuration it ran with is
// written to <work>/cobound.conf first: atomically, but not fsynced. It is a
// record of the run, rewritten by every run, and nothing reads it back to
// resume; an fsync before the search can stall for a second or more while
// the system writes back other data (phase 5's C2 measurement).
int commands::run(const std::vector<std::string> &args) {
  try {
    const config::Config cfg = config::forCommand("run", std::nullopt, args);
    const std::string work = cfg.text("work");
    std::filesystem::create_directories(work);
    report::atomicWrite(
        work + "/cobound.conf", [&](std::ostream &out) { cfg.writeEffective(out, "run"); },
        report::Durability::cache);
    if (cfg.context() == config::Context::run)
      return sweep::run(cfg);
    cascade::GoalOptions options = cascade::goalOptions(cfg);
    // The process-wide settings, from the config, before the first naming.
    options.censusLoaded = setup::applyRunSettings(cfg);
    // One policy for SIGINT and SIGTERM (divergence 8): the first ends the
    // running search cleanly and the run searches no further; a second ends
    // the process.
    runsignals::install();
    return cascade::runToGoal(options);
  } catch (const std::exception &e) {
    std::cerr << "cobound run: " << e.what() << "\n";
    return 2;
  }
}
