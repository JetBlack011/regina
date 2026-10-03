// config_test.cpp: cobound's one configuration schema (../driver/config.h).
//
// Every key's default reproduces the reference behaviour (plan, phase 5):
// without a goal verifyslicegenus's (9eb3cec4c: its option defaults and the
// library defaults it left in place), with a goal cascadesearch's (its
// Config and HopShape), and each command's its tool's. One check per key and
// context, against the reference value; the library's own constants where
// the reference took them from the library. Then the parser: files, --set,
// `none`, required keys, types, unknown and inapplicable keys, aliases, the
// run's context, and the effective configuration read back.

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <unistd.h>

#include "cobound/driver/config.h"
#include "linknaming/census/censusnaming.h"
#include "linknaming/complement/complementcache.h"
#include "linknaming/linknamer.h"
#include "linknaming/tests/check.h"
#include "surfer/enumeration/surfacesearch.h"

using config::Assignment;
using config::Config;
using config::Context;

namespace {

// The keys a context requires, set to something, so its defaults can be read.
std::vector<Assignment> minimal(Context c) {
  std::vector<Assignment> a;
  auto set = [&](const char *k, const char *v) { a.push_back({k, v, "test"}); };
  switch (c) {
  case Context::run:
    set("resolve_unlinked", "0");
    set("work", "/w");
    set("verdicts", "/v.csv");
    set("targets", "/t.csv");
    set("knot_table", "/k.csv");
    set("link_table", "/l.csv");
    break;
  case Context::goal:
    set("goal_genus", "0");
    set("resolve_unlinked", "1");
    set("work", "/w");
    set("target_pd", "[[1;5;2;4];[3;1;4;6];[5;3;6;2]]");
    set("knot_table", "/k.csv");
    set("link_table", "/l.csv");
    set("census", "/c.sqlite");
    break;
  case Context::solve:
    set("verdicts", "/v.csv");
    set("targets", "/t.csv");
    set("knot_table", "/k.csv");
    set("link_table", "/l.csv");
    break;
  case Context::sign:
    set("work", "/w");
    set("cobordisms", "/db.csv");
    set("knot_table", "/k.csv");
    set("link_table", "/l.csv");
    break;
  case Context::name:
    set("knot_table", "/k.csv");
    set("link_table", "/l.csv");
    break;
  case Context::draw:
    break;
  }
  return a;
}

std::string value(const Config &c, const std::string &key) {
  const config::Key *k = config::findKey(key);
  if (!c.has(key)) return "none";
  switch (k->type) {
  case config::Type::flag: return c.flag(key) ? "1" : "0";
  case config::Type::integer: return std::to_string(c.integer(key));
  case config::Type::real: {
    std::ostringstream o;
    o << c.real(key);
    return o.str();
  }
  case config::Type::threads: return std::to_string(c.threads(key));
  case config::Type::paths: {
    std::string s;
    for (const std::string &p : c.paths(key)) s += (s.empty() ? "" : ",") + p;
    return s;
  }
  default: return c.text(key);
  }
}

const std::string hardware = std::to_string(
    std::thread::hardware_concurrency() == 0 ? 1u : std::thread::hardware_concurrency());

// The reference's value of `key` in context `ctx`, and where it came from.
struct Reference {
  const char *key;
  Context ctx;
  std::string value;
  const char *source;
};

std::vector<Reference> references() {
  const SurfaceSearchLimits library;   // verifyslicegenus left these in place
  const linknaming::NamerLimits namer; // farsidename's defaults
  using C = Context;
  return {
      // the tables and files
      {"knot_symmetry", C::run, "none", "verifyslicegenus knotSymmetryPath"},
      {"knot_symmetry", C::goal, "none", "cascadesearch Config::knotSymmetry"},
      {"knot_symmetry", C::solve, "none", "verifyslicegenus knotSymmetryPath"},
      {"knot_symmetry", C::name, "none", "farsidename symmetry"},
      {"target_pd", C::run, "none", "(a run without a goal searches the targets' rows)"},
      {"target_name", C::run, "target", "cascadesearch: an unnamed target is `target`"},
      {"target_name", C::goal, "target", "cascadesearch: an unnamed target is `target`"},
      {"cobordisms", C::run, "cobordisms.csv", "verifyslicegenus cobordismsPath"},
      {"cobordisms", C::solve, "cobordisms.csv", "verifyslicegenus cobordismsPath"},
      {"cobordisms", C::goal, "none", "cascadesearch Config::witnessStore"},
      {"census", C::run, SURFER_CENSUS_PATH, "verifyslicegenus censusPath"},
      {"census", C::solve, SURFER_CENSUS_PATH, "verifyslicegenus censusPath"},
      {"dedupe_against", C::goal, "none", "cascadesearch Config::dedupeAgainst"},
      {"dedupe_against", C::sign, "none", "cascadesearch Config::dedupeAgainst"},
      {"pair_sig_cache", C::run, "none", "verifyslicegenus pairSigCacheDir"},
      {"pair_sig_cache", C::goal, "none", "cascadesearch Config::pairSigCache"},
      {"pair_sig_cache", C::sign, "none", "cascadesearch Config::pairSigCache"},
      {"pair_sig_cache", C::draw, "none", "farsidediagram sigCache"},
      // each search
      {"threads", C::run, hardware, "verifyslicegenus: hardware_concurrency()"},
      {"threads", C::goal, "10", "cascadesearch Config::threads"},
      {"threads", C::sign, "10", "cascadesearch --sign-only Config::threads"},
      {"threads", C::name, hardware, "farsidename: the namer's hardware_concurrency() pool"},
      {"layers", C::run, "2", "verifyslicegenus thickenLayers = collarLayers"},
      {"layers", C::goal, "2", "HopShape::layers"},
      {"layers", C::draw, "2", "farsidediagram layers"},
      {"boundary_condition", C::run, "proper", "verifyslicegenus boundaryConditionMode"},
      {"max_faces", C::run, "none", "verifyslicegenus maxFaces"},
      {"max_faces", C::goal, "5", "HopShape::maxFaces"},
      {"iddfs_iterations", C::run, "0", "verifyslicegenus iddfsIterations"},
      {"iddfs_iterations", C::goal, "2", "HopShape::iddfsIterations"},
      {"iddfs_start", C::run, "none", "verifyslicegenus iddfsStart"},
      {"iddfs_start", C::goal, "4", "HopShape::iddfsStart"},
      {"iddfs_step", C::run, "0", "verifyslicegenus iddfsStep"},
      {"iddfs_step", C::goal, "1", "HopShape::iddfsStep"},
      {"root_budget_start", C::run, "0", "verifyslicegenus rootBudgetStart"},
      {"root_budget_start", C::goal, "840", "HopShape::rootBudgetStart"},
      {"root_budget_growth", C::run, "2", "verifyslicegenus rootBudgetGrowth"},
      {"root_budget_growth", C::goal, "2", "HopShape::rootBudgetGrowth"},
      {"surface_target", C::run, "none", "verifyslicegenus surfaceTarget"},
      {"surface_target", C::goal, "100000", "cascadesearch Config::hopSurfaces"},
      {"max_surface_target", C::goal, "1000000", "cascadesearch Config::maxHopSurfaces"},
      {"search_seconds", C::run, "none", "verifyslicegenus perKnotTimeLimit"},
      {"search_seconds", C::goal, "7200", "cascadesearch: each hop's fixed 7200 s"},
      {"pending_surface_cap", C::run, std::to_string(library.pendingSurfaceCap),
       "SurfaceSearchLimits::pendingSurfaceCap"},
      {"pending_surface_cap", C::goal, "20000000", "HopShape::pendingSurfaceCap"},
      {"petal_cache_limit", C::run, std::to_string(library.petalCacheLimit),
       "SurfaceSearchLimits::petalCacheLimit"},
      {"petal_cache_limit", C::goal, "12000000", "HopShape::petalCacheLimit"},
      {"boundary_signature_cache_limit", C::run,
       std::to_string(library.boundarySignatureCacheLimit),
       "SurfaceSearchLimits::boundarySignatureCacheLimit"},
      {"boundary_signature_cache_limit", C::goal, "1000000",
       "HopShape::boundarySignatureCacheLimit"},
      {"complement_cache_limit", C::run, std::to_string(complement::recognitionCacheLimit.load()),
       "identify::recognitionCacheLimit"},
      {"complement_cache_limit", C::goal, "1500000", "HopShape::recognitionCacheLimit"},
      {"exact_far_side_names", C::run, "0", "verifyslicegenus exactFarSideNames"},
      {"exact_far_side_names", C::goal, "1", "cascadesearch: every hop names exactly"},
      {"census_updates", C::run, census::censusUpdates.load() ? "1" : "0",
       "census::censusUpdates (verifyslicegenus left it on)"},
      {"census_updates", C::goal, "0", "cascadesearch Cascade::run(): censusUpdates off"},
      {"retriangulate_on_miss", C::run, "1", "verifyslicegenus retriangulateOnMissArg"},
      {"retriangulate_on_miss", C::goal, "0", "cascadesearch Cascade::run(): off"},
      {"retriangulate_time_budget", C::run,
       std::to_string(census::retriangulateTimeBudgetSeconds.load()),
       "census::retriangulateTimeBudgetSeconds"},
      {"retriangulate_time_budget", C::goal,
       std::to_string(census::retriangulateTimeBudgetSeconds.load()),
       "census::retriangulateTimeBudgetSeconds (the cascade left it)"},
      {"max_crossings", C::run, "13", "verifyslicegenus maxCrossings"},
      {"max_crossings", C::solve, "13", "verifyslicegenus maxCrossings"},
      {"max_crossings", C::goal, "24", "cascadesearch Config::maxCrossings"},
      // what a search without a goal writes
      {"frontier_dir", C::run, "none", "verifyslicegenus frontierDir"},
      {"resume_frontier_dir", C::run, "none", "verifyslicegenus resumeFrontierDir"},
      {"surface_log", C::run, "none", "verifyslicegenus surfaceLogPath"},
      {"surface_stats", C::run, "none", "verifyslicegenus surfaceStatsPath"},
      {"self_intersection_census", C::run, "none", "verifyslicegenus"},
      {"rejection_sample_log", C::run, "none", "verifyslicegenus"},
      {"audit_linking", C::run, "0", "linkingnumber::auditLinkingNumbers"},
      // a goal
      {"goal_genus", C::goal, "0", "cascadesearch Config::goalGenus"},
      {"goal_lower", C::goal, "none", "cascadesearch Config::goalLower = -1"},
      {"goal_partition", C::goal, "connected", "cascadesearch Config::goalDisjoint = false"},
      {"literature", C::goal, "1", "cascadesearch Config::literature"},
      {"max_searches", C::goal, "20", "cascadesearch Config::maxExpansions"},
      {"cpu_budget", C::goal, "14400", "cascadesearch Config::cpuBudget"},
      {"strategy", C::goal, "best", "cascadesearch Config::strategy"},
      {"master_witnesses", C::goal, "none", "cascadesearch Config::masterWitnesses"},
      {"read_back_cache", C::goal, "none", "cascadesearch Config::readBackCache"},
      {"run_name", C::goal, "none", "cascadesearch Config::runName"},
      {"hub_degree", C::goal, "0", "cascadesearch Config::hubDegree"},
      {"hub_surfaces", C::goal, "0", "cascadesearch Config::hubSurfaces"},
      {"lower_report", C::goal, "0", "cascadesearch Config::lowerReport"},
      {"lower_sources", C::goal, "none", "cascadesearch Config::lowerSources"},
      {"lower_max_crossings", C::goal, "16", "cascadesearch Config::lowerMaxCrossings"},
      // solve
      {"name_aliases", C::solve, "none", "verifyslicegenus nameAliasPath"},
      {"far_side_resolutions", C::solve, "none", "verifyslicegenus"},
      {"far_side_exact", C::solve, "none", "verifyslicegenus"},
      {"link_classes", C::solve, "none", "verifyslicegenus"},
      {"cascade_proofs", C::solve, "none", "verifyslicegenus"},
      {"sum_rules", C::solve, "0", "verifyslicegenus sumRules"},
      // name
      {"namer_search_height", C::name, std::to_string(namer.searchHeight), "NamerLimits"},
      {"namer_search_visits", C::name, std::to_string(namer.searchVisits), "NamerLimits"},
      {"namer_simplify_tries", C::name, std::to_string(namer.simplifyTries), "NamerLimits"},
      {"namer_exhaustive_height", C::name, std::to_string(namer.exhaustiveHeight), "NamerLimits"},
      {"namer_max_search_crossings", C::name, std::to_string(namer.maxSearchCrossings),
       "NamerLimits"},
      {"namer_deep_height", C::name, std::to_string(namer.deepHeight), "NamerLimits"},
      {"namer_deep_visits", C::name, std::to_string(namer.deepVisits), "NamerLimits"},
      {"namer_max_deep_crossings", C::name, std::to_string(namer.maxDeepCrossings),
       "NamerLimits"},
      {"name_profile", C::name, "0", "farsidename profile"},
      {"name_reference", C::name, "0", "farsidename reference"},
      // draw
      {"draw_gauss", C::draw, "0", "farsidediagram gauss"},
      {"draw_faces", C::draw, "0", "farsidediagram facesInput"},
      {"draw_pairsig", C::draw, "0", "farsidediagram pairsigOut"},
  };
}

// Every key of every context it applies to is checked against a reference,
// once (a key added to the schema without one fails here).
void testDefaults() {
  const std::vector<Reference> refs = references();
  for (const Reference &r : refs) {
    const Config c(r.ctx, minimal(r.ctx));
    CHECK_EQ(value(c, r.key), r.value,
             std::string(r.key) + " in " + config::contextName(r.ctx) + " = " + r.source);
  }
  // Required keys: the reference had no default for them (or, `work`, the
  // plan takes its default away).
  struct Required {
    const char *key;
    Context ctx;
  };
  for (const Required &q : std::vector<Required>{
           {"resolve_unlinked", Context::run},  {"resolve_unlinked", Context::goal},
           {"work", Context::run},              {"work", Context::goal},
           {"work", Context::sign},             {"verdicts", Context::run},
           {"verdicts", Context::solve},        {"target_pd", Context::goal},
           {"knot_table", Context::run},        {"link_table", Context::run},
           {"knot_table", Context::solve},      {"link_table", Context::solve},
           {"targets", Context::run},           {"targets", Context::solve},
           {"knot_table", Context::goal},       {"link_table", Context::goal},
           {"knot_table", Context::sign},       {"link_table", Context::sign},
           {"knot_table", Context::name},       {"link_table", Context::name},
           {"census", Context::goal},           {"cobordisms", Context::sign}}) {
    std::vector<Assignment> a;
    for (const Assignment &m : minimal(q.ctx))
      if (m.key != q.key) a.push_back(m);
    bool threw = false;
    try {
      Config(q.ctx, a);
    } catch (const config::Error &e) {
      threw = std::string(e.what()).find(q.key) != std::string::npos;
    }
    CHECK(threw, std::string(q.key) + " is required in " + config::contextName(q.ctx));
  }
  // targets is required without a goal unless target_pd names the one
  // diagram the run searches instead.
  {
    std::vector<Assignment> a;
    for (const Assignment &m : minimal(Context::run))
      if (m.key != "targets") a.push_back(m);
    a.push_back({"target_pd", "[[1;5;2;4];[3;1;4;6];[5;3;6;2]]", "test"});
    bool built = false, hasTargets = true;
    try {
      const Config c(Context::run, a);
      built = true;
      hasTargets = c.has("targets");
    } catch (const config::Error &) {
    }
    CHECK(built && !hasTargets, "targets: not required by a run given target_pd");
  }
  // goal_genus's own default (0, cascadesearch's), with a run made
  // goal-directed by goal_lower alone.
  {
    std::vector<Assignment> a;
    for (const Assignment &m : minimal(Context::goal))
      if (m.key != "goal_genus") a.push_back(m);
    a.push_back({"goal_lower", "1", "test"});
    a.push_back({"lower_sources", "/s.csv", "test"});
    CHECK_EQ(config::runContext(a) == Context::goal, true, "goal_lower alone: a goal");
    CHECK_EQ(Config(Context::goal, a).integer("goal_genus"), 0LL,
             "goal_genus in run with a goal = cascadesearch Config::goalGenus");
  }
  // Each (key, context) the schema has is in the table above or required.
  size_t pairs = 0;
  for (const config::Key &k : config::schema())
    for (const auto &[ctx, rule] : k.rules) {
      ++pairs;
      bool found = rule.kind == config::Rule::Kind::required;
      for (const Reference &r : refs)
        found = found || (r.ctx == ctx && k.name == r.key);
      CHECK(found, k.name + " in " + config::contextName(ctx) + " has a reference check");
    }
  CHECK(pairs > 100, "the schema covers every command");
}

std::string tempFile(const std::string &text) {
  char tmpl[] = "/tmp/config_test_XXXXXX";
  const int fd = mkstemp(tmpl);
  close(fd);
  std::ofstream(tmpl) << text;
  return tmpl;
}

void testParsing() {
  const std::string file = tempFile("# a comment\n\nwork = /from/file\n"
                                    "resolve_unlinked = yes\n"
                                    "target_pd = PD[X[4; 1; 3; 2]; X[2; 3; 1; 4]]\n"
                                    "max_faces = 3\nverdicts = /v.csv\n"
                                    "knot_table = /k.csv\nlink_table = /l.csv\n");
  config::CommandLine cl = config::parseCommandLine(
      {"--config", file, "--set", "max_faces=4", "--set", "work = /from/set", "positional"});
  CHECK_EQ(cl.positional.size(), size_t(1), "one positional argument");
  CHECK_EQ(config::runContext(cl.assignments) == Context::run, true,
           "no goal key: a run without a goal");
  const Config c(Context::run, cl.assignments);
  CHECK_EQ(c.integer("max_faces"), 4LL, "--set wins over the file");
  CHECK_EQ(c.text("work"), std::string("/from/set"), "values are trimmed");
  CHECK_EQ(c.flag("resolve_unlinked"), true, "yes is a flag");
  CHECK_EQ(c.text("target_pd"), std::string("PD[X[4; 1; 3; 2]; X[2; 3; 1; 4]]"),
           "a value keeps its inner spaces");

  // none unsets a key without a default; a key with one cannot be unset.
  std::vector<Assignment> a = minimal(Context::run);
  a.push_back({"max_faces", "none", "test"});
  CHECK_EQ(Config(Context::run, a).has("max_faces"), false, "none unsets max_faces");
  a.push_back({"layers", "none", "test"});
  bool threw = false;
  try {
    Config(Context::run, a);
  } catch (const config::Error &) {
    threw = true;
  }
  CHECK(threw, "layers (default 2) cannot be none");

  auto refused = [](const std::string &key, const std::string &v) {
    std::vector<Assignment> b = minimal(Context::run);
    b.push_back({key, v, "test"});
    try {
      Config(Context::run, b);
    } catch (const config::Error &) {
      return true;
    }
    return false;
  };
  CHECK(refused("max_faces", "five"), "a word is no whole number");
  CHECK(refused("max_faces", "5x"), "nor is 5x");
  CHECK(refused("resolve_unlinked", "maybe"), "maybe is no flag");
  CHECK(refused("boundary_condition", "loose"), "a choice outside its list");
  CHECK(refused("threads", "0"), "threads must be > 0");
  CHECK(refused("no_such_key", "1"), "an unknown key");
  CHECK(!refused("threads", "auto"), "threads may be auto");

  // A key of another command is ignored, and said to be.
  std::vector<Assignment> g = minimal(Context::run);
  g.push_back({"max_searches", "4", "test"});
  const Config ig(Context::run, g);
  CHECK_EQ(ig.ignored().size(), size_t(1), "max_searches is no key of a run without a goal");

  // A goal key makes the run goal-directed (`none` does not).
  CHECK_EQ(config::runContext({{"goal_lower", "1", "t"}}) == Context::goal, true,
           "goal_lower: a goal");
  CHECK_EQ(config::runContext({{"goal_genus", "none", "t"}}) == Context::run, true,
           "goal_genus = none: no goal");

  // An unknown option is refused; an alias sets its key.
  threw = false;
  try {
    config::parseCommandLine({"--max-faces", "3"});
  } catch (const config::Error &) {
    threw = true;
  }
  CHECK(threw, "the retired options are refused");
  const config::CommandLine al = config::parseCommandLine(
      {"--knots", "/k", "--profile", "--links", "/l"},
      {{"--knots", "knot_table"}, {"--links", "link_table"}, {"--profile", "name_profile", false, "1"}});
  const Config nm(Context::name, al.assignments);
  CHECK_EQ(nm.text("knot_table"), std::string("/k"), "an alias with a value");
  CHECK_EQ(nm.flag("name_profile"), true, "an alias without one");
  std::filesystem::remove(file);
}

// The effective configuration, read back, is the same configuration.
void testEffective() {
  std::vector<Assignment> a = minimal(Context::goal);
  a.push_back({"target_pd", "PD[X[4; 1; 3; 2]; X[2; 3; 1; 4]]", "--set"});
  a.push_back({"max_faces", "3", "--set"});
  a.push_back({"threads", "auto", "--set"});
  a.push_back({"dedupe_against", "/a.csv, /b.csv", "--set"});
  const Config c(Context::goal, a);
  std::ostringstream out;
  c.writeEffective(out, "run");
  const std::string file = tempFile(out.str());
  const config::CommandLine back = config::parseCommandLine({"--config", file});
  CHECK_EQ(config::runContext(back.assignments) == Context::goal, true,
           "read back: still goal-directed");
  const Config d(Context::goal, back.assignments);
  bool same = true;
  for (const config::Key &k : config::schema())
    if (k.rule(Context::goal)) same = same && value(c, k.name) == value(d, k.name);
  CHECK(same, "every key reads back to the value it had");
  CHECK_EQ(d.paths("dedupe_against").size(), size_t(2), "a list reads back");
  CHECK(out.str().find("threads = " + hardware) != std::string::npos,
        "auto threads are written as the number used");
  std::filesystem::remove(file);
}

} // namespace

int main() {
  testDefaults();
  testParsing();
  testEffective();
  return checks::finish("config_test");
}
