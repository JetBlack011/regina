//
//  setup.cpp
//

#include "cobound/driver/setup.h"

#include <cstdlib>
#include <iostream>
#include <string>

#include "linknaming/census/censusnaming.h"
#include "linknaming/complement/complementcache.h"
#include "surfer/submanifold/linkingnumber.h"

namespace setup {

bool applyRunSettings(const config::Config &cfg) {
  const bool goal = cfg.context() == config::Context::goal;
  complement::recognitionCacheLimit.store(
      static_cast<size_t>(cfg.integer("complement_cache_limit")), std::memory_order_relaxed);
  census::censusUpdates.store(cfg.flag("census_updates"), std::memory_order_relaxed);
  census::retriangulateOnMiss.store(cfg.flag("retriangulate_on_miss"),
                                    std::memory_order_relaxed);
  census::retriangulateTimeBudgetSeconds.store(cfg.integer("retriangulate_time_budget"),
                                               std::memory_order_relaxed);
  if (!goal) {
    linkingnumber::auditLinkingNumbers.store(cfg.flag("audit_linking"));
    // The test hook of name_independence_test (a run without a goal, as
    // verifyslicegenus honoured it).
    if (const char *perturb = std::getenv("SURFER_TEST_PERTURB_NAMES");
        perturb && *perturb && std::string(perturb) != "0") {
      census::perturbNamesForTesting.store(true);
      std::cerr << "[!] SURFER_TEST_PERTURB_NAMES: every identified name is "
                   "perturbed (test mode)\n";
    }
  }
  return census::setCensusPath(cfg.text("census"));
}

} // namespace setup
