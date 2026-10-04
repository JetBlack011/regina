//
//  config.cpp
//

#include "cobound/driver/config.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <thread>

namespace config {

namespace {

using C = Context;

Rule def(std::string v) { return {Rule::Kind::value, std::move(v)}; }
const Rule unset{Rule::Kind::unset, ""};
const Rule required{Rule::Kind::required, ""};

Key key(std::string name, Type type, std::vector<std::pair<Context, Rule>> rules,
        std::string replaces, std::string doc, std::vector<std::string> choices = {}) {
    Key k;
    k.name = std::move(name);
    k.type = type;
    k.rules = std::move(rules);
    k.replaces = std::move(replaces);
    k.doc = std::move(doc);
    k.choices = std::move(choices);
    return k;
}

// A key renamed in phase 6, with the spelling it had before.
Key formerly(Key k, std::string oldName) {
    k.oldNames.push_back(std::move(oldName));
    return k;
}

std::string trim(const std::string &s) {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return s.substr(a, b - a);
}

std::optional<long long> parseInteger(const std::string &text) {
    if (text.empty()) return std::nullopt;
    size_t used = 0;
    try {
        const long long v = std::stoll(text, &used);
        if (used != text.size()) return std::nullopt;
        return v;
    } catch (const std::exception &) {
        return std::nullopt;
    }
}

std::optional<double> parseReal(const std::string &text) {
    if (text.empty()) return std::nullopt;
    size_t used = 0;
    try {
        const double v = std::stod(text, &used);
        if (used != text.size()) return std::nullopt;
        return v;
    } catch (const std::exception &) {
        return std::nullopt;
    }
}

std::vector<std::string> splitPaths(const std::string &text) {
    std::vector<std::string> out;
    size_t start = 0;
    while (start <= text.size()) {
        const size_t comma = text.find(',', start);
        const std::string item =
            trim(text.substr(start, comma == std::string::npos ? std::string::npos : comma - start));
        if (!item.empty()) out.push_back(item);
        if (comma == std::string::npos) break;
        start = comma + 1;
    }
    return out;
}

// `value` checked against `k`'s type, in its stored spelling; throws Error.
std::string canonical(const Key &k, const std::string &value, const std::string &source) {
    auto bad = [&](const std::string &what) -> Error {
        return Error(source + ": " + k.name + " = '" + value + "': " + what);
    };
    switch (k.type) {
    case Type::flag: {
        auto f = parseFlag(value);
        if (!f) throw bad("not a flag (0 or 1)");
        return *f ? "1" : "0";
    }
    case Type::integer: {
        auto v = parseInteger(value);
        if (!v) throw bad("not a whole number");
        return std::to_string(*v);
    }
    case Type::real:
        if (!parseReal(value)) throw bad("not a number");
        return value;
    case Type::threads: {
        if (value == "auto") return value;
        auto v = parseInteger(value);
        if (!v || *v <= 0) throw bad("not `auto` or a whole number > 0");
        return std::to_string(*v);
    }
    case Type::path:
        if (value.empty()) throw bad("an empty path");
        return value;
    case Type::paths:
        if (splitPaths(value).empty()) throw bad("no paths");
        return value;
    case Type::choice:
        if (std::find(k.choices.begin(), k.choices.end(), value) == k.choices.end()) {
            std::string list;
            for (const std::string &c : k.choices) list += (list.empty() ? "" : "|") + c;
            throw bad("not one of " + list);
        }
        return value;
    case Type::text:
        return value;
    }
    return value;
}

} // namespace

const char *contextName(Context c) {
    switch (c) {
    case C::run: return "run";
    case C::goal: return "run with a goal";
    case C::solve: return "solve";
    case C::sign: return "sign";
    case C::name: return "name";
    case C::draw: return "draw";
    }
    return "?";
}

const Rule *Key::rule(Context c) const {
    for (const auto &[ctx, r] : rules)
        if (ctx == c) return &r;
    return nullptr;
}

std::optional<bool> parseFlag(const std::string &text) {
    if (text == "1" || text == "true" || text == "yes" || text == "on") return true;
    if (text == "0" || text == "false" || text == "no" || text == "off") return false;
    return std::nullopt;
}

const std::vector<Key> &schema() {
    static const std::vector<Key> keys = {
        // ---- the tables and the files a command reads or writes ----
        key("knot_table", Type::path,
            {{C::run, required}, {C::goal, required}, {C::solve, required},
             {C::sign, required}, {C::name, required}},
            "--knot-table (verifyslicegenus, cascadesearch), --knots (farsidename)",
            "The knot table (Name,PD Notation,Genus-4D): literature bounds, table PDs, names. "
            "Required everywhere: verifyslicegenus's default, the file's name in the working "
            "directory, would read whatever CSV of that name happens to be there."),
        key("link_table", Type::path,
            {{C::run, required}, {C::goal, required}, {C::solve, required},
             {C::sign, required}, {C::name, required}},
            "--link-table (verifyslicegenus, cascadesearch), --links (farsidename)",
            "The link table, as knot_table."),
        key("knot_symmetry", Type::path,
            {{C::run, unset}, {C::goal, unset}, {C::solve, unset}, {C::name, unset}},
            "--knot-symmetry (verifyslicegenus, cascadesearch), --symmetry (farsidename)",
            "Knot symmetry types (data/knot_symmetry.csv): slice composites beyond 3_1#m3_1 "
            "and 4_1#4_1 are anchors only with them."),
        key("targets", Type::path, {{C::run, required}, {C::solve, required}},
            "--input (verifyslicegenus)",
            "Without a goal: the table rows to search, each once, those within max_crossings, "
            "in crossing order. For solve: the rows the verdicts list first. Required (no "
            "working-directory default), except by a run whose target_pd names its one "
            "diagram."),
        key("target_pd", Type::text, {{C::run, unset}, {C::goal, required}},
            "--target-pd (cascadesearch)",
            "One diagram to search from, as a PD code. Without a goal it is searched once, "
            "instead of the rows of targets."),
        key("target_name", Type::text, {{C::run, def("target")}, {C::goal, def("target")}},
            "--target-name (cascadesearch)",
            "target_pd's name: its table name, or any name (its cobordisms are recorded under "
            "it)."),
        key("verdicts", Type::path, {{C::run, required}, {C::solve, required}},
            "--output (verifyslicegenus)",
            "The verdicts file (verify_genus_v2.csv's columns). A run writes each searched "
            "row's search record; solve re-derives every status and bound."),
        key("cobordisms", Type::path,
            {{C::run, def("cobordisms.csv")}, {C::goal, unset}, {C::solve, def("cobordisms.csv")},
             {C::sign, required}},
            "--cobordisms (verifyslicegenus), --witness-store (cascadesearch)",
            "The cobordism database. Without a goal, loaded (a search keeps no cobordism it "
            "holds) and signed into at the run's end; with a goal, signed into at the run's "
            "end when set (with run_name); sign signs into it; solve reads it."),
        key("work", Type::path, {{C::run, required}, {C::goal, required}, {C::sign, required}},
            "--work (verifyslicegenus, whose default <cobordisms>.pending retires; "
            "cascadesearch)",
            "The run directory: each search's hop_<k>_n<node>/ (its pending cobordisms, "
            "kept.csv), and with a goal the run's records. sign signs every pending file "
            "in it."),
        key("census", Type::path,
            {{C::run, def(SURFER_CENSUS_PATH)}, {C::goal, required}, {C::solve, def(SURFER_CENSUS_PATH)}},
            "--census-db (verifyslicegenus, cascadesearch)",
            "The local census database naming complements."),
        key("dedupe_against", Type::paths, {{C::goal, unset}, {C::sign, unset}},
            "--dedupe-against (cascadesearch)",
            "Read-only databases whose cobordisms are not recorded again when signing."),
        key("pair_sig_cache", Type::path,
            {{C::run, unset}, {C::goal, unset}, {C::sign, unset}, {C::draw, unset}},
            "--pair-sig-cache (verifyslicegenus, cascadesearch), --sig-cache (farsidediagram)",
            "Where pair-signature contexts (.pairsigctx) are kept between runs."),

        // ---- how each search runs ----
        key("threads", Type::threads,
            {{C::run, def("auto")}, {C::goal, def("10")}, {C::sign, def("10")}, {C::name, def("auto")}},
            "--threads (verifyslicegenus, cascadesearch)",
            "Threads for searching, naming and signing (`auto`: the machine's). At one thread "
            "the end-of-run signing is not overlapped with the search, so a one-thread run "
            "pays its pair-signature contexts in full after it (R1s x1.7, accepted on "
            "2026-10-03)."),
        key("layers", Type::integer, {{C::run, def("2")}, {C::goal, def("2")}, {C::draw, def("2")}},
            "--thicken-layers and --collar-layers (verifyslicegenus; collared through every "
            "layer), farsidediagram --layers",
            "Layers of the thickening a search runs in (collared through all of them)."),
        key("boundary_condition", Type::choice, {{C::run, def("proper")}},
            "--boundary-condition (verifyslicegenus)",
            "Which surfaces a search without a goal accepts: `proper` (every row), `connected` "
            "or `auto` (connected for a knot). A goal run always searches `proper`.",
            {"auto", "connected", "proper"}),
        key("resolve_unlinked", Type::flag, {{C::run, required}, {C::goal, required}},
            "--resolve-unlinked / --no-resolve-unlinked",
            "Also accept surfaces whose only self-intersections are unlinked. No default: it "
            "changes what counts toward a surface target and every frontier's fingerprint."),
        key("max_faces", Type::integer, {{C::run, unset}, {C::goal, def("5")}},
            "--max-faces (verifyslicegenus), --hop-max-faces (cascadesearch)",
            "Faces a search may add to the collar seed (unset: no cap)."),
        key("iddfs_iterations", Type::integer, {{C::run, def("0")}, {C::goal, def("2")}},
            "--iddfs-iterations, --hop-iddfs-iterations", "Iterative-deepening rounds."),
        key("iddfs_start", Type::integer, {{C::run, unset}, {C::goal, def("4")}},
            "--iddfs-start, --hop-iddfs-start", "The first round's face cap."),
        key("iddfs_step", Type::integer, {{C::run, def("0")}, {C::goal, def("1")}},
            "--iddfs-step, --hop-iddfs-step", "How much each round raises the cap."),
        key("root_budget_start", Type::integer, {{C::run, def("0")}, {C::goal, def("840")}},
            "--root-budget-start (verifyslicegenus), --hop-root-budget (cascadesearch)",
            "Each root's first pass's ration of attempts (0: unbudgeted)."),
        key("root_budget_growth", Type::integer, {{C::run, def("2")}, {C::goal, def("2")}},
            "--root-budget-growth, --hop-root-growth", "How the ration grows between passes."),
        key("surface_target", Type::integer, {{C::run, unset}, {C::goal, def("100000")}},
            "--surface-target (verifyslicegenus), --hop-surfaces (cascadesearch)",
            "Stop a search once this many surfaces satisfy its boundary condition. With a "
            "goal, each search's first target, doubled up to max_surface_target once nothing "
            "useful is left at it."),
        key("max_surface_target", Type::integer, {{C::goal, def("1000000")}},
            "--max-hop-surfaces (cascadesearch)", "The largest surface target a goal run uses."),
        key("search_seconds", Type::real, {{C::run, unset}, {C::goal, def("7200")}},
            "--per-knot-time-limit (verifyslicegenus); the cascade's fixed 7200 s per hop",
            "A wall-clock backstop per search; its drain still finishes."),
        key("pending_surface_cap", Type::integer,
            {{C::run, def("500000")}, {C::goal, def("20000000")}},
            "--pending-surface-cap, --hop-pending-cap",
            "Surfaces queued for the drain before the search helps drain them."),
        key("petal_cache_limit", Type::integer,
            {{C::run, def("2000000")}, {C::goal, def("12000000")}},
            "--petal-cache-limit, --hop-petal-cache", "Petal cache entries before it is cleared."),
        key("boundary_signature_cache_limit", Type::integer,
            {{C::run, def("200000")}, {C::goal, def("1000000")}},
            "--boundary-signature-cache-limit, --hop-boundary-cache",
            "Outgoing names cached per search before the cache is cleared."),
        key("complement_cache_limit", Type::integer,
            {{C::run, def("200000")}, {C::goal, def("1500000")}},
            "--recognition-cache-limit (verifyslicegenus), --hop-recognition-cache "
            "(cascadesearch)",
            "Complement answers cached per process before the cache is cleared."),
        formerly(key("outgoing_names", Type::flag, {{C::run, def("1")}, {C::goal, def("1")}},
                     "--exact-far-side-names (verifyslicegenus; a goal run's searches always did)",
                     "Always 1: outgoing links are always named, an oriented link per "
                     "surface. Kept so configurations that pass it still parse; 0 is "
                     "refused."),
                 "exact_far_side_names"),
        key("census_updates", Type::flag, {{C::run, def("1")}, {C::goal, def("0")}},
            "--no-census-updates (verifyslicegenus; the cascade never wrote)",
            "Write the census: a knot row's complement after its search, and Pachner hits."),
        key("retriangulate_on_miss", Type::flag, {{C::run, def("1")}, {C::goal, def("0")}},
            "--no-retriangulate-on-miss (verifyslicegenus; the cascade's searches never did)",
            "Run a Pachner search for a knot complement the census misses."),
        key("retriangulate_time_budget", Type::integer, {{C::run, def("20")}, {C::goal, def("20")}},
            "--retriangulate-time-budget (verifyslicegenus)",
            "Seconds per Pachner search."),
        key("max_crossings", Type::integer,
            {{C::run, def("13")}, {C::goal, def("24")}, {C::solve, def("13")}},
            "--max-crossings (verifyslicegenus: rows above it are skipped; cascadesearch: no "
            "link above it is searched)",
            "Without a goal (and for solve), rows above it are not searched (skipped); with "
            "one, no link above it is searched."),

        // ---- what a search without a goal also writes ----
        key("frontier_dir", Type::path, {{C::run, unset}}, "--frontier-dir (verifyslicegenus)",
            "Record each row's frontier as <dir>/<row>.frontier."),
        key("resume_frontier_dir", Type::path, {{C::run, unset}},
            "--resume-frontier-dir (verifyslicegenus)",
            "Carry each row on from <dir>/<row>.frontier."),
        key("surface_log", Type::path, {{C::run, unset}}, "--surface-log (verifyslicegenus)",
            "Every surface each search finds (rewritten per row): G-surfaces' instrument."),
        key("surface_stats", Type::path, {{C::run, unset}}, "--surface-stats (verifyslicegenus)",
            "Appends each searched row's surface statistics (surface_stats.csv)."),
        key("self_intersection_census", Type::path, {{C::run, unset}},
            "--self-intersection-census (verifyslicegenus)",
            "Measurement only: appends each row's self-intersection census."),
        key("rejection_sample_log", Type::path, {{C::run, unset}},
            "--rejection-sample-log (verifyslicegenus)",
            "The first surfaces each gate turns away, with their pair signatures."),
        key("audit_linking", Type::flag, {{C::run, def("0")}}, "--audit-linking (verifyslicegenus)",
            "Validation only: every petal linking number computed twice; halt on a "
            "disagreement."),

        // ---- a goal ----
        key("goal_genus", Type::integer, {{C::goal, def("0")}}, "--goal-genus (cascadesearch)",
            "The goal: the target's genus at most this (setting it makes the run "
            "goal-directed)."),
        key("goal_lower", Type::integer, {{C::goal, unset}}, "--goal-lower (cascadesearch)",
            "Also stop once the target's lower bound reaches this (needs lower_sources)."),
        key("goal_partition", Type::choice, {{C::goal, def("connected")}},
            "--goal (cascadesearch)", "The goal's surface: one connected piece, or disjoint "
            "pieces.", {"connected", "disjoint"}),
        key("literature", Type::flag, {{C::goal, def("1")}},
            "--literature / --constructive (cascadesearch)",
            "Whether literature upper bounds may enter a proof."),
        key("max_searches", Type::integer, {{C::goal, def("20")}},
            "--max-expansions (cascadesearch)", "Searches a goal run may make."),
        key("cpu_budget", Type::real, {{C::goal, def("14400")}}, "--cpu-budget (cascadesearch)",
            "CPU seconds the searches may spend."),
        key("strategy", Type::choice, {{C::goal, def("best")}}, "--strategy (cascadesearch)",
            "Which link to search next.", {"best", "dfs", "bfs"}),
        formerly(key("master_cobordisms", Type::path, {{C::goal, unset}},
                     "--master-witnesses (cascadesearch)",
                     "A read-only database whose cobordisms of each link about to be searched "
                     "enter the graph first (never loaded without a goal)."),
                 "master_witnesses"),
        key("read_back_cache", Type::path, {{C::goal, unset}},
            "--read-back-cache (cascadesearch)",
            "Where those cobordisms' outgoing links are kept between runs."),
        key("run_name", Type::text, {{C::goal, unset}}, "--run-name (cascadesearch)",
            "The run's name in cascade:<run>/<target>/n<node> subjects (needed with "
            "cobordisms)."),
        key("hub_degree", Type::integer, {{C::goal, def("0")}}, "--hub-degree (cascadesearch)",
            "A link with this many cobordisms is searched once at hub_surfaces (0: off)."),
        key("hub_surfaces", Type::integer, {{C::goal, def("0")}}, "--hub-surfaces (cascadesearch)",
            "A hub's surface target."),
        key("lower_report", Type::flag, {{C::goal, def("0")}}, "--lower-report (cascadesearch)",
            "Write lower_report.jsonl at the run's end."),
        key("lower_sources", Type::path, {{C::goal, unset}}, "--lower-sources (cascadesearch)",
            "Which literature lower bounds are special (data/lower_bound_sources.csv)."),
        key("lower_max_crossings", Type::integer, {{C::goal, def("16")}},
            "--lower-max-crossings (cascadesearch)",
            "A link kept only for the lower goal is not searched above this."),

        // ---- solve ----
        key("name_aliases", Type::path, {{C::solve, unset}}, "--name-aliases (verifyslicegenus)",
            "Observed -> classical outgoing names, applied when solving."),
        formerly(key("outgoing_resolutions", Type::path, {{C::solve, unset}},
                     "--far-side-resolutions (verifyslicegenus)",
                     "Per-cobordism proved outgoing names."),
                 "far_side_resolutions"),
        formerly(key("outgoing_names_file", Type::path, {{C::solve, unset}},
                     "--far-side-exact (verifyslicegenus)",
                     "Per-cobordism names of outgoing links (name)."),
                 "far_side_exact"),
        key("link_classes", Type::path, {{C::solve, unset}}, "--link-classes (verifyslicegenus)",
            "The table's link classes (tableclasses)."),
        formerly(key("certified_bounds", Type::path, {{C::solve, unset}},
                     "--cascade-proofs (verifyslicegenus)",
                     "Certified bounds (data/cascade_proofs.csv)."),
                 "cascade_proofs"),
        key("sum_rules", Type::flag, {{C::solve, def("0")}}, "--sum-rules (verifyslicegenus)",
            "Bound sums along components and splits with link factors from their pieces."),

        // ---- name ----
        key("namer_search_height", Type::integer, {{C::name, def("2")}},
            "--search-height (farsidename)", "linknaming::NamerLimits::searchHeight."),
        key("namer_search_visits", Type::integer, {{C::name, def("20000")}},
            "--search-visits (farsidename)", "NamerLimits::searchVisits."),
        key("namer_simplify_tries", Type::integer, {{C::name, def("24")}},
            "--simplify-tries (farsidename)", "NamerLimits::simplifyTries."),
        key("namer_exhaustive_height", Type::integer, {{C::name, def("1")}},
            "--exhaustive-height (farsidename)", "NamerLimits::exhaustiveHeight."),
        key("namer_max_search_crossings", Type::integer, {{C::name, def("16")}},
            "--max-search-crossings (farsidename)", "NamerLimits::maxSearchCrossings."),
        key("namer_deep_height", Type::integer, {{C::name, def("3")}},
            "--deep-height (farsidename)", "NamerLimits::deepHeight."),
        key("namer_deep_visits", Type::integer, {{C::name, def("200000")}},
            "--deep-visits (farsidename)", "NamerLimits::deepVisits."),
        key("namer_max_deep_crossings", Type::integer, {{C::name, def("12")}},
            "--max-deep-crossings (farsidename)", "NamerLimits::maxDeepCrossings."),
        key("name_profile", Type::flag, {{C::name, def("0")}}, "--profile (farsidename)",
            "Per-step times on stderr at the end."),
        key("name_reference", Type::flag, {{C::name, def("0")}}, "--reference (farsidename)",
            "Read each outgoing link by the reference (uncached) carry."),

        // ---- draw ----
        key("draw_gauss", Type::flag, {{C::draw, def("0")}}, "--gauss (farsidediagram)",
            "Append signed Gauss data (and the ROW line's build=)."),
        key("draw_faces", Type::flag, {{C::draw, def("0")}}, "--faces (farsidediagram)",
            "Read `<id> <f1,f2,...>` (faces of the row's thickening) instead of pair "
            "signatures."),
        key("draw_pairsig", Type::flag, {{C::draw, def("0")}}, "--pairsig (farsidediagram)",
            "With draw_faces: append each surface's pair signature."),
    };
    return keys;
}

const Key *findKey(const std::string &name) {
    for (const Key &k : schema())
        if (k.name == name ||
            std::find(k.oldNames.begin(), k.oldNames.end(), name) != k.oldNames.end())
            return &k;
    return nullptr;
}

CommandLine parseCommandLine(const std::vector<std::string> &args,
                             const std::vector<Alias> &aliases) {
    CommandLine out;
    std::vector<Assignment> sets;
    for (size_t i = 0; i < args.size(); ++i) {
        const std::string &a = args[i];
        auto next = [&]() -> const std::string & {
            if (i + 1 >= args.size()) throw Error(a + " needs a value");
            return args[++i];
        };
        if (a == "--config") {
            const std::string file = next();
            std::ifstream in(file);
            if (!in) throw Error("cannot read the config file " + file);
            std::string line;
            for (size_t n = 1; std::getline(in, line); ++n) {
                const std::string t = trim(line);
                if (t.empty() || t[0] == '#') continue;
                const size_t eq = t.find('=');
                if (eq == std::string::npos)
                    throw Error(file + ":" + std::to_string(n) + ": not `key = value`");
                out.assignments.push_back(
                    {trim(t.substr(0, eq)), trim(t.substr(eq + 1)), file + ":" + std::to_string(n)});
            }
        } else if (a == "--set") {
            const std::string kv = next();
            const size_t eq = kv.find('=');
            if (eq == std::string::npos) throw Error("--set " + kv + ": not key=value");
            sets.push_back({trim(kv.substr(0, eq)), trim(kv.substr(eq + 1)), "--set"});
        } else if (auto al = std::find_if(aliases.begin(), aliases.end(),
                                          [&](const Alias &x) { return x.flag == a; });
                   al != aliases.end()) {
            sets.push_back({al->key, al->takesValue ? next() : al->value, a});
        } else if (a.size() > 2 && a.compare(0, 2, "--") == 0) {
            throw Error("unknown option " + a + " (keys are set with --config FILE or --set "
                        "key=value; `cobound help keys` lists them)");
        } else {
            out.positional.push_back(a);
        }
    }
    out.assignments.insert(out.assignments.end(), sets.begin(), sets.end());
    return out;
}

Config forCommand(const std::string &command, std::optional<Context> context,
                  const std::vector<std::string> &args, const std::vector<Alias> &aliases,
                  std::vector<std::string> *positional) {
    CommandLine cl = parseCommandLine(args, aliases);
    if (positional)
        *positional = cl.positional;
    else if (!cl.positional.empty())
        throw Error("cobound " + command + ": unexpected argument '" + cl.positional.front() +
                    "' (keys are set with --config FILE or --set key=value)");
    Config cfg(context.value_or(runContext(cl.assignments)), cl.assignments);
    for (const Assignment &a : cfg.ignored())
        std::cerr << "[!] config: " << a.key << " is not a key of " << contextName(cfg.context())
                  << " (" << a.source << "); ignored\n";
    return cfg;
}

Context runContext(const std::vector<Assignment> &assignments) {
    std::map<std::string, std::string> last;
    for (const Assignment &a : assignments) last[a.key] = a.value;
    for (const char *k : {"goal_genus", "goal_lower"})
        if (auto it = last.find(k); it != last.end() && it->second != "none") return C::goal;
    return C::run;
}

Config::Config(Context context, const std::vector<Assignment> &assignments)
    : context_(context) {
    std::map<std::string, const Assignment *> last;
    for (const Assignment &a : assignments) {
        const Key *k = findKey(a.key);
        if (!k) throw Error(a.source + ": unknown key '" + a.key + "'");
        if (!k->rule(context)) {
            ignored_.push_back(a);
            continue;
        }
        last[k->name] = &a; // an old spelling assigns the key it names
    }
    for (const Key &k : schema()) {
        const Rule *r = k.rule(context);
        if (!r) continue;
        Value v;
        if (auto it = last.find(k.name); it != last.end()) {
            const Assignment &a = *it->second;
            v.source = a.source;
            if (a.key != k.name) v.spelling = a.key;
            if (a.value == "none") {
                if (r->kind != Rule::Kind::unset)
                    throw Error(a.source + ": " + k.name + " cannot be none in " +
                                contextName(context));
            } else {
                v.text = canonical(k, a.value, a.source);
            }
        } else {
            v.source = "default";
            // targets names the rows a run searches; a run given one diagram
            // (target_pd) searches that instead and needs none.
            const auto pd = last.find("target_pd");
            const bool oneDiagram = k.name == "targets" && context == Context::run &&
                                    pd != last.end() && pd->second->value != "none";
            if (r->kind == Rule::Kind::required && !oneDiagram)
                throw Error(k.name + " is required for " + std::string(contextName(context)) +
                            " (it has no default)");
            if (r->kind == Rule::Kind::value) v.text = r->value;
        }
        values_.emplace(k.name, std::move(v));
    }
}

const Config::Value &Config::value_(const std::string &key, Type expected) const {
    auto it = values_.find(key);
    if (it == values_.end())
        throw std::logic_error("config: " + key + " is not a key of " + contextName(context_));
    const Key *k = findKey(key);
    if (k->type != expected && !(expected == Type::text && k->type == Type::path))
        throw std::logic_error("config: " + key + " read as the wrong type");
    return it->second;
}

bool Config::has(const std::string &key) const {
    auto it = values_.find(key);
    if (it == values_.end())
        throw std::logic_error("config: " + key + " is not a key of " + contextName(context_));
    return it->second.text.has_value();
}

bool Config::flag(const std::string &key) const {
    const Value &v = value_(key, Type::flag);
    return v.text && *v.text == "1";
}

long long Config::integer(const std::string &key) const {
    auto v = optionalInteger(key);
    if (!v) throw std::logic_error("config: " + key + " is unset");
    return *v;
}

std::optional<long long> Config::optionalInteger(const std::string &key) const {
    const Value &v = value_(key, Type::integer);
    if (!v.text) return std::nullopt;
    return std::stoll(*v.text);
}

double Config::real(const std::string &key) const {
    auto v = optionalReal(key);
    if (!v) throw std::logic_error("config: " + key + " is unset");
    return *v;
}

std::optional<double> Config::optionalReal(const std::string &key) const {
    const Value &v = value_(key, Type::real);
    if (!v.text) return std::nullopt;
    return std::stod(*v.text);
}

std::string Config::text(const std::string &key) const {
    return optionalText(key).value_or("");
}

std::optional<std::string> Config::optionalText(const std::string &key) const {
    const Key *k = findKey(key);
    const Value &v = value_(key, k && k->type == Type::choice ? Type::choice : Type::text);
    return v.text;
}

std::vector<std::string> Config::paths(const std::string &key) const {
    const Value &v = value_(key, Type::paths);
    return v.text ? splitPaths(*v.text) : std::vector<std::string>{};
}

unsigned Config::threads(const std::string &key) const {
    const Value &v = value_(key, Type::threads);
    if (!v.text || *v.text == "auto") {
        const unsigned n = std::thread::hardware_concurrency();
        return n == 0 ? 1 : n;
    }
    return static_cast<unsigned>(std::stoul(*v.text));
}

void Config::writeEffective(std::ostream &out, const std::string &command) const {
    out << "# cobound " << command << " (" << contextName(context_)
        << "): the configuration it ran with, every key it reads.\n"
           "# A value not the default is preceded by where it came from. Read back with "
           "--config, it runs the same.\n";
    for (const Key &k : schema()) {
        auto it = values_.find(k.name);
        if (it == values_.end()) continue;
        std::string value = it->second.text.value_or("none");
        if (k.type == Type::threads && value == "auto") {
            value = std::to_string(threads(k.name));
            out << "# " << it->second.source << ": auto\n";
        } else if (it->second.source != "default") {
            out << "# " << it->second.source
                << (it->second.spelling.empty() ? "" : " (as " + it->second.spelling + ")")
                << "\n";
        }
        out << k.name << " = " << value << "\n";
    }
    for (const Assignment &a : ignored_)
        out << "# ignored (not a key of " << contextName(context_) << ", " << a.source
            << "): " << a.key << " = " << a.value << "\n";
}

} // namespace config
