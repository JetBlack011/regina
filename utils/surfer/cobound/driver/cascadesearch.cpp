// cascadesearch.cpp
//
// Goal-directed chained searches for one target knot or link (the
// scheduler: scheduler.h). See README.md.
//
// Usage:
//   cascadesearch --target-pd '<PD>' --target-name NAME --work DIR
//                 --knot-table CSV --link-table CSV --census-db PATH
//                 [--knot-symmetry CSV] [--master-witnesses CSV]
//                 [--hop-mode process|child] [--verifyslicegenus PATH]
//                 [--goal-genus G] [--goal disjoint|connected]
//                 [--hop-surfaces N] [--max-hop-surfaces N] [--threads T]
//                 [--max-expansions K] [--cpu-budget SECONDS]
//                 [--max-crossings C] [--strategy best|dfs|bfs]
//                 [--literature | --constructive]
//                 (--resolve-unlinked | --no-resolve-unlinked)
//                 [--witness-store CSV --run-name NAME [--dedupe-against CSV]...]
//   cascadesearch --sign-only --work DIR --witness-store CSV
//                 --knot-table CSV --link-table CSV [--dedupe-against CSV]...
//
// With --witness-store, every surface a hop keeps is also recorded for the
// atlas (pending.h): its hop appends it to <hop dir>/kept.csv at once, and
// the run's end signs the ones whose witness is new (to the store and to
// every --dedupe-against file) and appends them to the store, as the sweep
// would have recorded them. A hop's subject is the target's own name, a
// node's proved table name, or cascade:<run name>/<target>/n<node>, so no two
// runs can ever record different links under one name. --sign-only does the
// end-of-run step alone, for a run that was killed.

#include <iostream>
#include <stdexcept>
#include <string>

#include "cobound/cobordisms/pending.h"
#include "cobound/driver/scheduler.h"
#include "cobound/driver/signals.h"
#include "cobound/solver/literature.h"

using namespace cascade;

int main(int argc, char **argv) {
  GoalOptions c;
  bool signOnly = false;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) throw std::invalid_argument("missing value for " + a);
      return argv[++i];
    };
    if (a == "--target-pd") c.targetPD = next();
    else if (a == "--target-name") c.targetName = next();
    else if (a == "--work") c.work = next();
    else if (a == "--verifyslicegenus") c.verify = next();
    else if (a == "--knot-table") c.knotTable = next();
    else if (a == "--link-table") c.linkTable = next();
    else if (a == "--knot-symmetry") c.knotSymmetry = next();
    else if (a == "--census-db") c.censusDb = next();
    else if (a == "--goal-genus") c.goalGenus = std::stoi(next());
    else if (a == "--goal") c.goalDisjoint = next() == "disjoint";
    else if (a == "--hop-surfaces") c.hopSurfaces = std::stol(next());
    else if (a == "--max-hop-surfaces") c.maxHopSurfaces = std::stol(next());
    else if (a == "--threads") c.threads = std::stoi(next());
    else if (a == "--max-expansions") c.maxExpansions = std::stoi(next());
    else if (a == "--cpu-budget") c.cpuBudget = std::stod(next());
    else if (a == "--max-crossings") c.maxCrossings = std::stoul(next());
    else if (a == "--strategy") c.strategy = next();
    else if (a == "--literature") c.literature = true;
    else if (a == "--constructive") c.literature = false;
    else if (a == "--master-witnesses") c.masterWitnesses = next();
    else if (a == "--hop-mode") c.hopMode = next();
    else if (a == "--hop-max-faces") c.hopShape.maxFaces = std::stoll(next());
    else if (a == "--hop-iddfs-start") c.hopShape.iddfsStart = std::stoll(next());
    else if (a == "--hop-iddfs-iterations") c.hopShape.iddfsIterations = std::stoul(next());
    else if (a == "--hop-root-budget") c.hopShape.rootBudgetStart = std::stoll(next());
    else if (a == "--hop-iddfs-step") c.hopShape.iddfsStep = std::stoll(next());
    else if (a == "--hop-root-growth") c.hopShape.rootBudgetGrowth = std::stoll(next());
    else if (a == "--hop-pending-cap") c.hopShape.pendingSurfaceCap = std::stoull(next());
    else if (a == "--hop-petal-cache") c.hopShape.petalCacheLimit = std::stoull(next());
    else if (a == "--hop-boundary-cache")
      c.hopShape.boundarySignatureCacheLimit = std::stoull(next());
    else if (a == "--hop-recognition-cache")
      c.hopShape.recognitionCacheLimit = std::stoull(next());
    else if (a == "--resolve-unlinked") c.hopShape.resolveUnlinked = true;
    else if (a == "--no-resolve-unlinked") c.hopShape.resolveUnlinked = false;
    else if (a == "--verbose") c.verbose = true;
    else if (a == "--witness-store") c.witnessStore = next();
    else if (a == "--pair-sig-cache") c.pairSigCache = next();
    else if (a == "--read-back-cache") c.readBackCache = next();
    else if (a == "--run-name") c.runName = next();
    else if (a == "--dedupe-against") c.dedupeAgainst.push_back(next());
    else if (a == "--sign-only") signOnly = true;
    else if (a == "--hub-degree") c.hubDegree = std::stoul(next());
    else if (a == "--hub-surfaces") c.hubSurfaces = std::stol(next());
    else if (a == "--lower-report") c.lowerReport = true;
    else if (a == "--lower-sources") c.lowerSources = next();
    else if (a == "--goal-lower") c.goalLower = std::stoi(next());
    else if (a == "--master-loads") {
      const std::string v = next();
      if (v != "lazy" && v != "eager") throw std::invalid_argument("--master-loads lazy|eager");
      c.lazyMasterLoads = v == "lazy";
    }
    else if (a == "--lower-max-crossings") c.lowerMaxCrossings = std::stoul(next());
    else {
      std::cerr << "unknown argument " << a << "\n";
      return 2;
    }
  }
  if (signOnly) {
    // The end-of-run store step alone, for a run that was killed after its
    // hops wrote kept.csv (the subjects are already in those lines).
    if (c.work.empty() || c.witnessStore.empty() || c.knotTable.empty() ||
        c.linkTable.empty()) {
      std::cerr << "cascadesearch --sign-only: --work, --witness-store, --knot-table and "
                   "--link-table are required\n";
      return 2;
    }
    try {
      const cobordismgraph::NameTable names =
          witnessstore::loadTableNames(c.knotTable, c.linkTable, "");
      const StoreResult s = signPending(c.work, c.witnessStore, c.dedupeAgainst,
                                        names, static_cast<unsigned>(c.threads),
                                        c.pairSigCache);
      std::cout << "[+] witness store: " << s.kept << " kept, " << s.fresh << " new, "
                << s.appended << " appended to " << c.witnessStore << "\n";
      return 0;
    } catch (const std::exception &e) {
      std::cerr << "cascadesearch: " << e.what() << "\n";
      return 2;
    }
  }
  if (!c.hopShape.resolveUnlinked) {
    // Plan divergence 3: no default. It changes which surfaces a hop's
    // surface target counts, and every frontier's fingerprint.
    std::cerr << "cascadesearch: --resolve-unlinked or --no-resolve-unlinked is required "
                 "(resolve_unlinked has no default)\n";
    return 2;
  }
  if (!c.witnessStore.empty() && (c.runName.empty() || c.hopMode != "process")) {
    // A child hop writes its own witnesses (hop dir cob.csv), under its own
    // row name; only in-process hops are recorded through pending.h.
    std::cerr << "cascadesearch: --witness-store needs --run-name, and in-process hops\n";
    return 2;
  }
  if (c.goalLower >= 0 && c.lowerSources.empty()) {
    std::cerr << "cascadesearch: --goal-lower needs --lower-sources (the special sources' "
                 "largest bound is what makes the lower gate prune)\n";
    return 2;
  }
  if (c.targetPD.empty() || c.work.empty() || c.knotTable.empty() ||
      c.linkTable.empty() || c.censusDb.empty() ||
      (c.hopMode == "child" && c.verify.empty())) {
    std::cerr << "cascadesearch: --target-pd, --work, --knot-table, --link-table and "
                 "--census-db are required, and --verifyslicegenus with --hop-mode child\n";
    return 2;
  }
  if (c.targetName.empty()) c.targetName = "target";
  try {
    // One policy for SIGINT and SIGTERM (divergence 8): the first ends the
    // running search cleanly and the run searches no further; a second ends
    // the process.
    runsignals::install();
    return runToGoal(c);
  } catch (const std::exception &e) {
    std::cerr << "cascadesearch: " << e.what() << "\n";
    return 2;
  }
}
