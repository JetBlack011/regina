# cobound: bounding g₄ by cobordisms

One program, `cobound`, that searches a link's thickening for surfaces
(`surfer`), names their outgoing links (`linknaming`), keeps every cobordism
it finds in an append-only database, and bounds the smooth 4-genus from
cobordisms and the literature. A run has targets, an optional goal and limits:
without a goal it searches each target once; with one it searches outwards
from the target, link by link, until the goal has a proof. Paper: §6
(`sec:cobordisms`) and §7 (`sec:compresults`): `def:search`,
`def:cobordism-record`, `sec:bounds` and `sec:partition-genera`.

`cobound` is the only part where search meets naming. It keeps two bound
engines apart on purpose: the **cobordism graph** (`bounds/`), which a goal run
and each search's judge build, and the **atlas solver** (`solver/`), which
`cobound solve` runs and the atlas's `tools/frontier.py` reimplements
independently. They share no code, so that `frontier.py --check` (which must
say AGREE after any change to either solver) remains a check.

## The code

| folder | what |
|---|---|
| `main.cpp`, `driver/` | the commands (`run.cpp`, `solve.cpp`, `sign.cpp`, `draw.cpp`, `name.cpp`, `meridians.cpp`); `config`, the one configuration schema; `setup`, the process-wide settings; `targets`; `scheduler`, a goal run; `runrecords`, a goal run's files; `signals`; `fatal`, halting on an impossible state; `timers.h` |
| `search/` | `search`: THE search (`Searcher::run()`) and its watchdog; `incoming`: the incoming link's thickening, orientation and certification; `preconditions`: what a found surface must satisfy, the accounting, the dedupe key; `searchreport`: the progress block, per-search files and log lines |
| `cobordisms/` | `cobordism`: one record and its identity; `database`: the file; `pending`: a search's pending file and signing; `pairsigner`; `cobordismkey`; `appendonly`: locked, fsynced appends |
| `outgoing/` | `outgoinglink`: a surface's outgoing curves on T, oriented against the incoming link; `outgoingnamer`: the search's `BoundaryNamer`; `fromdatabase`: a stored cobordism's outgoing link, read back |
| `bounds/` | the cobordism graph: `partition`, `partitiongenera`, `cobordismgraph`, `links` (the registry), `axioms`, `searchcobordisms` (a search's finds into the graph), `databasecobordisms`, `searchjudge` (a search's own graph), `certificate` |
| `solver/` | the atlas solver: `literature` (`NameTable`), `solver`, `solverinputs`, `verdicts` |
| `frozen.h` | every output token that something reads back or parses and that still spells a retired term (below, "Frozen formats") |
| `json.h`, `parallelfor.h` | the one JSON writer; the one thread loop |

## Commands

```
cobound run       [--config FILE]... [--set key=value]...
cobound solve     [--config FILE]... [--set key=value]...
cobound sign      [--config FILE]... [--set key=value]...
cobound draw      [--layers N] [--gauss] [--faces [--pairsig [--sig-cache DIR]]] '<incoming PD>' < input
cobound name      --knots CSV --links CSV [--symmetry CSV] [namer limits] [--profile] [--reference] < input
cobound meridians sig | dump | dump-link | dump-subset | slope < input
cobound help [keys]
```

| command | does |
|---|---|
| `run` | searches from each target: without a goal each once; with a goal (`goal_genus` or `goal_lower` set) outwards from the target until the goal has a proof or the limits are spent. Writes the configuration it ran with to `<work>/cobound.conf`, which, read back with `--config`, runs the same search |
| `solve` | re-derives every status and bound into the verdicts file from the database, the tables and the solver's other inputs. Never writes the database, never searches |
| `sign` | signs a work directory's pending cobordisms into the database: what a run's own end does, for a run that was killed |
| `draw` | stored cobordisms' outgoing links as oriented diagrams, from their pair signatures (or, with `--faces`, from triangle lists) on one incoming diagram's thickening. Its first line is `ROW components=<n> lk=<matrix>`, then one `W <id> ok ...` or `W <id> FAILED <reason>` per cobordism; the format is in `driver/draw.cpp`'s header |
| `name` | stored cobordisms' outgoing links named by the full namer (`linknaming::LinkNamer::name()`), from `<id> <incoming> <layers> <pair signature>` lines; one tab-separated line each (`driver/name.cpp`'s header) |
| `meridians` | drilling with signed meridians, for the atlas's SnapPy pipeline: `dump` (from pair signatures), `dump-link` (from PD codes), `dump-subset`, `slope`, `sig`; the record format is in `driver/meridians.cpp`'s header. Each input line is answered on its own |

`draw` and `name` also take their flags as config keys (`--config`, `--set`);
`meridians` takes only its subcommand. Every command is a stdin/stdout tool
except `run`, `solve` and `sign`.

**Exit codes.**

| command | 0 | 1 | 2 | 3 |
|---|---|---|---|---|
| `run` without a goal | done (a first SIGINT or SIGTERM included) | the tables or symmetry types cannot be read | a configuration error; a search whose accounting failed or whose output write failed (after the run's other targets); a halt: a contradiction, an impossible state, a broken seed invariant, a linking-audit disagreement | · |
| `run` with a goal | the goal met | not met (limits spent, nothing useful left, or a first signal) | a configuration error, a table target whose given diagram does not certify, a run that cannot start; a halt: an impossible state, a broken seed invariant, a failed output write | a contradiction, even when the goal is met |
| `solve` | done | an input cannot be read | a configuration error; a contradiction (after writing the verdicts) | · |
| `sign` | done | · | any failure | · |
| `draw`, `name` | done (a cobordism that cannot be read is a `FAILED` line) | · | a usage or configuration error | · |
| `meridians` | done (`FAILED` records) | a usage error | · | · |

During a run, a second SIGINT or SIGTERM ends the process at once, `_Exit(128 + signal)` (below, "Signals").

## Configuration

`run`, `solve` and `sign` (and `draw` and `name`, besides their flags) read
`key = value` lines from `--config FILE` (any number, in order), then
`--set key=value` (any number; the last wins). One assignment per line; the
value runs to the end of the line, trimmed (a PD code may hold spaces); blank
lines and lines starting with `#` are skipped. `none` unsets a key that may be
unset; `paths` are comma-separated; a flag is `0`/`1` (or `true`/`false`,
`yes`/`no`, `on`/`off`).

One schema serves every command (`driver/config.h`). A key applies in some
contexts (`run`, `run` with a goal, `solve`, `sign`, `draw`, `name`), and in
each has a default, no default (unset), or none at all (required). A run is
goal-directed when `goal_genus` or `goal_lower` is assigned. An unknown key, a
value of the wrong type or a missing required key is an error; a known key
that does not apply is reported on stderr and ignored. `cobound help keys`
prints the schema, from which this table is made (`·`: the key does not
apply):

| key | type | run | with a goal | solve | sign | draw | name | what |
|---|---|---|---|---|---|---|---|---|
| `knot_table` | path | **required** | **required** | **required** | **required** | · | **required** | The knot table (`Name,PD Notation,Genus-4D`): literature bounds, table PDs, names. |
| `link_table` | path | **required** | **required** | **required** | **required** | · | **required** | The link table, likewise. |
| `knot_symmetry` | path | unset | unset | unset | · | · | unset | Knot symmetry types (the atlas's `knot_symmetry.csv`); without them only `3_1#m3_1` and `4_1#4_1` are slice-composite anchors. |
| `targets` | path | **required** | · | **required** | · | · | · | Without a goal: the table rows to search, each once, in crossing order, up to `max_crossings`. For `solve`: the rows the verdicts list first. Not needed when `target_pd` names the one diagram. |
| `target_pd` | text | unset | **required** | · | · | · | · | One diagram to search from, as a PD code (without a goal, instead of `targets`). |
| `target_name` | text | `target` | `target` | · | · | · | · | `target_pd`'s name: its table name, or any name; its cobordisms are recorded under it. |
| `verdicts` | path | **required** | · | **required** | · | · | · | The verdicts file. A run writes each searched table row's search record; `solve` re-derives every status and bound. |
| `cobordisms` | path | `cobordisms.csv` | unset | `cobordisms.csv` | **required** | · | · | The database. Without a goal: loaded (a search keeps no cobordism it holds) and signed into at the run's end. With a goal: signed into at the end when set (with `run_name`). `sign` signs into it; `solve` reads it. |
| `work` | path | **required** | **required** | · | **required** | · | · | The run directory: each search's pending file, and with a goal the run's records. `sign` signs every pending file in it. |
| `census` | path | `SURFER_CENSUS_PATH` | **required** | `SURFER_CENSUS_PATH` | · | · | · | The local census. |
| `dedupe_against` | paths | · | unset | · | unset | · | · | Read-only databases whose cobordisms are not recorded again when signing (comma-separated). |
| `pair_sig_cache` | path | unset | unset | · | unset | unset | · | Where pair-signature contexts (`.pairsigctx`) are kept between runs. |
| `threads` | threads | `auto` | `10` | · | `10` | · | `auto` | Threads for searching, naming and signing (`auto`: the machine's). At one thread, signing is not overlapped with the search. |
| `layers` | integer | `2` | `2` | · | · | `2` | · | Layers of the thickening a search runs in, collared through all of them. |
| `boundary_condition` | auto \| connected \| proper | `proper` | · | · | · | · | · | Which surfaces a search without a goal accepts: `proper`, or `connected` (`auto`: connected for a knot). A goal run always searches `proper`. |
| `resolve_unlinked` | flag | **required** | **required** | · | · | · | · | Also accept surfaces whose only self-intersections are unlinked. No default: it changes what counts toward a surface target and every frontier's fingerprint. |
| `max_faces` | integer | unset | `5` | · | · | · | · | Triangles a search may add to the collar (unset: no cap). |
| `iddfs_iterations` | integer | `0` | `2` | · | · | · | · | Iterative-deepening rounds. |
| `iddfs_start` | integer | unset | `4` | · | · | · | · | The first round's cap. |
| `iddfs_step` | integer | `0` | `1` | · | · | · | · | How much each round raises the cap. |
| `root_budget_start` | integer | `0` | `840` | · | · | · | · | Each root's first pass's ration of attempts (0: unbudgeted); see `surfer/README.md`, "Performance". |
| `root_budget_growth` | integer | `2` | `2` | · | · | · | · | How the ration grows between passes (at least 2). |
| `surface_target` | integer | unset | `100000` | · | · | · | · | Stop a search once this many surfaces satisfy its boundary condition. With a goal, the first target, doubled up to `max_surface_target` once nothing useful is left. |
| `max_surface_target` | integer | · | `1000000` | · | · | · | · | The largest surface target a goal run uses. |
| `search_seconds` | real | unset | `7200` | · | · | · | · | A wall-clock backstop per search; its drain still finishes. |
| `pending_surface_cap` | integer | `500000` | `20000000` | · | · | · | · | Surfaces queued for the drain before the search pauses to drain them. |
| `petal_cache_limit` | integer | `2000000` | `12000000` | · | · | · | · | Petal cache entries before it is cleared. |
| `boundary_signature_cache_limit` | integer | `200000` | `1000000` | · | · | · | · | Names cached per boundary component per search before the cache is cleared. |
| `complement_cache_limit` | integer | `200000` | `1500000` | · | · | · | · | Complement answers cached per process before the cache is cleared. |
| `outgoing_names` (formerly `exact_far_side_names`) | flag | `1` | `1` | · | · | · | · | Always 1: outgoing links are always named. Kept so configurations that pass it parse; 0 is refused. |
| `census_updates` | flag | `1` | `0` | · | · | · | · | Write the census: a knot target's complement after its search, and Pachner hits. |
| `retriangulate_on_miss` | flag | `1` | `0` | · | · | · | · | A Pachner search for a knot complement the census misses. |
| `retriangulate_time_budget` | integer | `20` | `20` | · | · | · | · | Seconds per Pachner search. |
| `max_crossings` | integer | `13` | `24` | `13` | · | · | · | Without a goal (and for `solve`), table rows above it are not searched; with one, no link above it is searched. |
| `frontier_dir` | path | unset | · | · | · | · | · | Record each target's frontier as `<dir>/<target>.frontier`. |
| `resume_frontier_dir` | path | unset | · | · | · | · | · | Carry each target on from `<dir>/<target>.frontier`. |
| `surface_log` | path | unset | · | · | · | · | · | Every surface each search describes, rewritten per search: the equivalence check's instrument. |
| `surface_stats` | path | unset | · | · | · | · | · | Appends each search's surface statistics (`surface_stats.csv`). |
| `self_intersection_census` | path | unset | · | · | · | · | · | Measurement only: appends each search's self-intersection census. |
| `rejection_sample_log` | path | unset | · | · | · | · | · | The first surfaces each gate turns away, with their pair signatures. |
| `audit_linking` | flag | `0` | · | · | · | · | · | Validation only: every petal linking number computed twice; halt on a disagreement. |
| `goal_genus` | integer | · | `0` | · | · | · | · | The goal: the target's genus at most this. Setting it, or `goal_lower`, makes the run goal-directed. |
| `goal_lower` | integer | · | unset | · | · | · | · | Also stop once the target's lower bound reaches this (needs `lower_sources`). |
| `goal_partition` | connected \| disjoint | · | `connected` | · | · | · | · | The goal's surface: one connected piece, or disjoint pieces. |
| `literature` | flag | · | `1` | · | · | · | · | Whether literature upper bounds may enter a proof. |
| `max_searches` | integer | · | `20` | · | · | · | · | Searches a goal run may make. |
| `cpu_budget` | real | · | `14400` | · | · | · | · | CPU seconds its searches may spend (checked between searches). |
| `strategy` | best \| dfs \| bfs | · | `best` | · | · | · | · | Which link to search next. |
| `master_cobordisms` (formerly `master_witnesses`) | path | · | unset | · | · | · | · | A read-only database whose cobordisms of each table link about to be searched enter the graph first. |
| `read_back_cache` | path | · | unset | · | · | · | · | Where those cobordisms' outgoing links are kept between runs. |
| `run_name` | text | · | unset | · | · | · | · | The run's name in untabulated links' subject names (needed with `cobordisms`). |
| `hub_degree` | integer | · | `0` | · | · | · | · | A link with this many cobordisms is searched once at `hub_surfaces` (0: off). |
| `hub_surfaces` | integer | · | `0` | · | · | · | · | A hub's surface target. |
| `lower_report` | flag | · | `0` | · | · | · | · | Write `lower_report.jsonl` at the run's end. |
| `lower_sources` | path | · | unset | · | · | · | · | Which literature lower bounds are special (the atlas's `lower_bound_sources.csv`). |
| `lower_max_crossings` | integer | · | `16` | · | · | · | · | A link kept only for the lower goal is not searched above this. |
| `name_aliases` | path | · | · | unset | · | · | · | Observed → table outgoing names, applied when solving. |
| `outgoing_resolutions` (formerly `far_side_resolutions`) | path | · | · | unset | · | · | · | Per-cobordism proved outgoing names (keyed by cobordism key). |
| `outgoing_names_file` (formerly `far_side_exact`) | path | · | · | unset | · | · | · | Per-cobordism outgoing names (the atlas's `far_side_exact.csv`). |
| `link_classes` | path | · | · | unset | · | · | · | The tables' link classes (`tableclasses`). |
| `certified_bounds` (formerly `cascade_proofs`) | path | · | · | unset | · | · | · | Bounds certified from goal runs (the atlas's `cascade_proofs.csv`). |
| `sum_rules` | flag | · | · | `0` | · | · | · | Bound sums along components, and splits with link factors, from their pieces. |
| `namer_search_height` | integer | · | · | · | · | · | `2` | `NamerLimits::searchHeight`. |
| `namer_search_visits` | integer | · | · | · | · | · | `20000` | `NamerLimits::searchVisits`. |
| `namer_simplify_tries` | integer | · | · | · | · | · | `24` | `NamerLimits::simplifyTries`. |
| `namer_exhaustive_height` | integer | · | · | · | · | · | `1` | `NamerLimits::exhaustiveHeight`. |
| `namer_max_search_crossings` | integer | · | · | · | · | · | `16` | `NamerLimits::maxSearchCrossings`. |
| `namer_deep_height` | integer | · | · | · | · | · | `3` | `NamerLimits::deepHeight`. |
| `namer_deep_visits` | integer | · | · | · | · | · | `200000` | `NamerLimits::deepVisits`. |
| `namer_max_deep_crossings` | integer | · | · | · | · | · | `12` | `NamerLimits::maxDeepCrossings`. |
| `name_profile` | flag | · | · | · | · | · | `0` | Per-step times on stderr at the end. |
| `name_reference` | flag | · | · | · | · | · | `0` | Read each outgoing link back by the uncached carry. |
| `draw_gauss` | flag | · | · | · | · | `0` | · | Append signed Gauss data (and the incoming line's `build=`). |
| `draw_faces` | flag | · | · | · | · | `0` | · | Read `<id> <f1,f2,...>` (triangles of the incoming diagram's thickening) instead of pair signatures. |
| `draw_pairsig` | flag | · | · | · | · | `0` | · | With `draw_faces`: append each surface's pair signature. |

`resolve_unlinked` and `work` have no default anywhere: `resolve_unlinked`
changes what a surface target counts and is part of every frontier's
fingerprint, and `work` is where a run's pending files go. Nor do the tables
and `targets`. A run without a goal also needs a stopping rule: `max_faces`,
`surface_target` or `search_seconds`.

**Old spellings.** Five keys were renamed; the old spelling is still accepted
everywhere a key is read, and `<work>/cobound.conf` writes the new name with
`# <source> (as <old>)`:

| key | formerly |
|---|---|
| `outgoing_names` | `exact_far_side_names` |
| `master_cobordisms` | `master_witnesses` |
| `outgoing_resolutions` | `far_side_resolutions` |
| `outgoing_names_file` | `far_side_exact` |
| `certified_bounds` | `cascade_proofs` |

A goal run, for example (`D` the atlas's `data/`):

```
target_pd = <PD>
target_name = 10_27
work = <dir>
knot_table = $D/4d_smooth_slice_genus_13_crossings_pd_codes.csv
link_table = $D/links_4d_smooth_slice_genus_11_crossings_pd_codes.csv
knot_symmetry = $D/knot_symmetry.csv
census = <census copy>
goal_genus = 1
literature = 0
resolve_unlinked = 1
surface_target = 50000
max_surface_target = 200000
max_searches = 8
```

**Once per process** (`driver/setup`), before the first naming: the census
path and whether it is written, Pachner searches on a census miss and their
time budget, the complement cache's limit, and (without a goal) the linking
audit. The tables are loaded once and shared by the search's namer and the
graph's. Nothing else is loaded unless asked for: a database's cobordisms as
free cobordisms of the graph only with a goal (`master_cobordisms`), and link classes
lazily, per base, only by a goal run's namer.

## A search's life

One search (`search::Searcher::run()`) is the same whoever runs it: a run
without a goal searches each target through it, and a goal run each link it
expands.

1. **The thickening.** T is built from the PD code exactly as given
   (`diagramtriangulation::parsePDCode()`: the same integers in the same
   order give the same T), thickened `layers` times and collared through
   every layer (`search::buildIncoming()`). The incoming link's edges, and
   its orientation on them, are pinned by the collar's own edges, never by an
   arbitrary isomorphism (T has automorphisms moving L).
2. **Certification.** T's incoming link, drawn back
   (`diagramtriangulation::DiagramDrawer`), must be isomorphic as a diagram to
   the PD's own, orientation kept and no mirror (`search::certifyIncoming()`).
   Every search does this. A table row that does not certify is refused, never
   searched on another diagram, since that would silently change T: without a
   goal it is recorded as `build-failed` and the run goes on; with a goal,
   for the target, the run exits 2. Only an untabulated goal target may fall
   back to its simplified diagram, and the run says so. Every table row
   certifies (17,153 of 17,153).
3. **Preconditions.** A search is refused (`SearchRefused`) without a collar
   seed or when its outgoing namer cannot be built. If any searchable triangle
   other than the seed's touches the incoming boundary, found surfaces could
   change the incoming link: an impossible state (`SeedInvariantFailure`), on
   which every run halts once what it found is written.
4. **The search** (`surfer`'s `SurfaceSearch`, its SIGINT handling off): under
   `proper` (or, without a goal and for a knot, `connected` if asked), with the
   shape's rounds, face cap and budgets, carried on from a frontier if one is
   offered and valid (step 8). It stops at the surface target (the search's
   breadth), when it runs out of candidates at its face cap, at
   `search_seconds` (the watchdog, checked after the surface target so that a
   tie records the surface target), on a signal, or on an output write that
   fails. Its drain finishes in every case.
5. **Gates and accounting** (`search/preconditions`). Each surface the drain
   describes is judged, in order: orientable; its boundary split into the
   incoming side and the others **by geometry** (the incoming boundary
   component; names are never compared, since census names are not
   canonical); the incoming side has the incoming link's component count; the
   incoming curves run with the incoming link's orientation, **judged per
   surface component** (each component may be reversed independently and the
   components still tube into one oriented surface); at most one outgoing
   side. Every described surface lands in exactly one bucket, and the search
   prints them:
   ```
   [+] <name>: accounting: accepted A, described D, recorded R, duplicate U,
       other-orientation O, search-side-elsewhere 0, impossible I, drain complete|skipped, ok|WARNING
   ```
   (one line). `other-orientation` counts cobordisms from another oriented
   variant of the incoming link. `impossible` counts states a correct build
   never reaches (non-orientable under orientability pruning, a broken
   incoming side, an incoherent curve, two outgoing sides, an unnamed side).
   `search-side-elsewhere` is a frozen bucket of a retired search mode and is
   always 0. The accounting **fails** when accepted ≠ described (unless the
   drain was cut short), described ≠ the buckets' sum, a surface could not be
   rebuilt, or `impossible` > 0. `WARNING` means surfaces were described but
   none recorded or deduplicated, which is what a broken gate looks like.
6. **Keeping.** An accepted surface's outgoing link is named
   (`outgoing::OutgoingNamer`: a knot by its edge set's name, a link by its
   oriented name per surface; see `linknaming/README.md`), and one surface is
   kept per key: the cobordism identity (kind, subject and its component
   count, outgoing name and its component count, genus, tubed, whether
   resolved) together with which incoming components and how many outgoing
   curves each surface component carries (`search::keptKey()`). A surface
   whose identity the loaded database already holds is a duplicate.
7. **Pending.** Each kept surface is appended to the search's pending file,
   `<work>/hop_<k>_n<link>/kept.csv` (the directory name is frozen), with its
   faces, incoming PD and layers. (A goal run without a database keeps its
   finds in its graph only.) The file is fsynced at most once a minute
   while the search and its drain run, at once at the first constructive find,
   and at the end, so a kill loses at most a minute of finds. A search signs
   nothing: pair signatures are made when the run signs its pending files
   (below, "The database and signing").
8. **The frontier** (`surfer/README.md`, "Frontiers"). A search's frontier is
   kept only if it can vouch for every surface in its prefix: the accounting
   balanced, something was examined (no `WARNING`), the drain completed, every
   output was written, and its pending file is fsynced. It records that file
   (relative to its own directory, so a moved work tree keeps its resume) and
   the file's fsynced length. A later search skips the prefix only if `sign`
   has signed that file at least that far (`<pending>.signed`), or the file is
   in its own run's work directory, whose sign step signs it; otherwise it
   refuses the frontier, says why, and searches from the start. A frontier
   that cannot be read (frontiers are written atomically but not fsynced) is
   treated as absent, with a warning. Without a goal, `frontier_dir` records
   `<dir>/<target>.frontier` and `resume_frontier_dir` carries each target on
   from one; a goal run keeps one frontier per link in its run.
9. **Judging** (without a goal). Each search has its own cobordism graph
   (`bounds::SearchJudge`): the searched link with its literature lower bound,
   and each find added as it is kept, on a thread of its own. A
   contradiction (a derived bound below a proved lower bound) ends the search
   and halts the run with 2 once its cobordisms are written. When a proof
   resting on no literature value reaches the searched link's literature lower
   bound, the search prints the frozen
   `[+] <name>: CONSTRUCTIVE witness found -- reaches genus G (literature [lo, hi]). Checkpointing now.`,
   shows `-- ACHIEVED` in its progress block, and fsyncs its pending file. The
   judge reads no database, so a link constructive only through another
   search's cobordisms is not announced; that is `solve`'s to find.
10. **Ending.** The outcome is `exhausted` (it ran out of candidates at its
    face cap), `surface-target`, `timeout`, `interrupted`, `stopped`,
    `unaccounted` (its accounting failed, or nothing was examined), `io-error`
    (an output write failed: no frontier, no exhaustion claim, what it kept is
    still signed), or `fatal-bug`; and, for a target that could not be built
    or was refused, `build-failed`. Only a search whose accounting balanced,
    with something examined, its drain complete and every write made, may
    claim exhaustion. Without a goal the run then prints
    `EXHAUSTIVE to N added faces` if it may, writes the target's search
    record into the verdicts file (`searched_faces`, `search_outcome`, and
    `exhausted_depth`, which is never lowered; status and bounds are left as
    they were), prints its `breadth:` line (when frontiers are kept), its
    outcome, accounting, `identification:`, `diagram naming:` and
    `search profile:` lines, and goes on to the next target. A goal run
    continues as below.

An accounting failure ends only its own search: without a goal the run
exits 2 after its other targets, since completeness is what such a run is for;
with a goal the search is marked `[!!] hop <k>: surface accounting failed --
...` and the run's code is unchanged, its certificate standing. An
impossible state halts any run (exit 2), after writing what was found, even
when the goal is met.

**Signals** (`driver/signals`). The first SIGINT or SIGTERM ends the running
search cleanly within a second: its drain finishes, its pending file is
fsynced, its frontier is ruled on like any stopped search's, and it records
`interrupted`. The run then starts no further search and ends as at any other
end: exit 0 without a goal (the search is thinner, not failed), 1 with one
unless the goal is met. An interrupted search is never taken as exhausted. A
second signal ends the process at once. SIGKILL (a cgroup's OOM killer) is
caught by nothing: the pending file's fsyncs bound that loss.

## The database and signing

`cobordisms.csv` is append-only and never rewritten. Its 13 columns, frozen:

```
kind,subject,subject_components,other,other_candidates,other_components,genus,tubed,pairsig,source_row,thicken_layers,max_faces,resolved_vertices
```

`kind` is `cobordism` or `direct` (no outgoing link: the surface bounds its
subject alone); `subject` is what the search searched; `other` the outgoing
link's name or description (`linknaming/README.md`) and `other_candidates`
the oriented variants a name could be; `genus` the tubed genus of the surface
(its components tubed together, closed ones discarded) and `tubed` whether it
was disconnected; `pairsig` the surface's pair signature over the searched
thickening; `source_row` the search's subject; `resolved_vertices` the
unlinked self-intersections, written empty when 0, so older lines read back
byte for byte. Loading keeps each line's byte offset instead of its pair
signature (~11 KB), reading one back only when needed; a torn last line is
ignored on load and cut before the next append; a file of another header is
refused for appending.

**Identity.** `cobordisms::cobordismIdentity()` is the one dedupe identity:
kind, subject and its component count, outgoing name and its component count,
genus, tubed, and whether resolved. The database holds one cobordism per
identity. An outgoing link that bears no bound is named, in the end, from
something several links can share, so two surfaces reaching different such
links can collapse to one identity; keying them by their own edges instead
made nearly every surface its own cobordism (19,690 from 23,672 surfaces of
one exhaustive search at cap 3, against 109), which no database can hold.

**Signing** (`cobordisms/pending`, `sign`). At a run's end (a halt included),
or by `cobound sign` after a kill, every pending file under `work` is read
from where it was last signed (a torn last line skipped). Each kept surface
whose identity is new to the database, to every `dedupe_against` database and
to the files read before it gets a pair signature (`cobordisms/pairsigner`:
each thickening's ambient part computed once, on the run's threads, kept in
`pair_sig_cache`), its `other_candidates` from the name tables, and is
appended and fsynced under an exclusive lock on `<database>.lock`, which
re-reads the database's identities, so runs sharing a database never record
one identity twice. Each file's `<file>.signed` records how far it was signed.
The run prints `[+] witness store: K kept, F new, A appended to <database>`.

**Which diagram.** A pair signature's ambient is the thickening of the diagram
searched, which for a goal run's untabulated or simplified link is not its
subject's table PD. So signing writes a line per appended cobordism in the
sidecar `<database>.rows.csv` (`witness,layers,row_pd`: its cobordism key,
layers and incoming PD; frozen), which goal runs read back and the atlas feeds
to `cobound name`. A run without a goal writes one only for a cobordism not
searched on its subject's table PD as written, so none for a table row.

**Subjects.** A goal run records a search's cobordisms under the target's own
name, a link's proved table name, or `cascade:<run>/<target>/n<link>` for an
untabulated link (frozen; link ids restart in every run), listed in
`<work>/nodes.csv`.

**Keys.** `cobordisms::cobordismKey()` is `sha1(pairsig)[:12]`: what the
atlas's per-cobordism files key on, computed there with Python's `hashlib`,
which `cobordismkey_test` checks it against.

## Goal runs

A run with a goal (`driver/scheduler`) searches from the target outwards. The
target is searched on its PD as given (certified once, step 2); every other
link on its registered diagram (simplified, reduced, every component in its
slot). Each search's surfaces enter the run's cobordism graph
(`bounds/searchcobordisms`), read straight off the search: no pair signature is
written and read back, and a kept surface keeps its triangles instead. The
loop, until a check ends it:

1. a halt (exit 2) or a contradiction (exit 3) ends the run, after signing
   what was found and writing the partition genera; the contradiction gate is
   checked **before** the goal, so a goal met through inconsistent facts is
   never certified;
2. the goal met (upper: the target's best surface on the goal partition has
   genus at most `goal_genus`; lower: `lower(target, goal partition)` at least
   `goal_lower`) ends the loop;
3. a first signal ends it (`interrupted`);
4. the target's own database cobordisms are loaded first, when it is a table
   entry `master_cobordisms` holds;
5. **choose** the next link (below); stop at `max_searches` or once the
   searches' CPU reaches `cpu_budget` (both checked between searches);
6. when nothing is useful at this surface target, double it and search every
   useful link again, deeper, carrying each on from its frontier; stop when the
   doubled target would pass `max_surface_target` (`nothing-useful`);
7. load the chosen link's database cobordisms first, if any (free cobordisms, which
   may close the proof without a search), then search it. A link with at least
   `hub_degree` cobordisms is searched once at `hub_surfaces`.

**Choosing.** A link is eligible when it has a diagram of at most
`max_crossings` crossings, was not refused, has not been searched at this
budget and is not searched out (its frontier complete). It is **useful** for
the upper goal if, given the best partition genera it could conceivably have
(every partition its linking numbers allow, at its proved lower bound), the
goal would follow over the cobordisms found so far. With a lower goal, a link
the upper gate rejects is kept if, given the best lower bounds it could ever
have, it would carry the lower goal to the target ("Lower goals"). `best`
orders by fewest crossings; then a link the upper gate keeps before one only
the lower gate keeps, and among those the most slack; then lower volume
(compared to 10⁻⁶, so runs are reproducible), shallower, older. `dfs` and
`bfs` order by depth first.

**Database cobordisms** (`master_cobordisms`, `bounds/databasecobordisms`). For
a table link about to be searched, the database's cobordisms of that subject,
and those of other subjects whose recorded outgoing link has the link's base
name (the reverse direction a forward search cannot find), are read back onto
the thickening of the diagram each was searched on (the `.rows.csv` sidecar,
else the subject's table PD; `outgoing/fromdatabase`), named and turned into
graph cobordisms. Names in the file are only a hint: every one is redrawn and
named exactly. A cobordism whose genus exceeds what any proof could spend is
skipped unread. Read-backs are kept across runs in `read_back_cache` (one file
per incoming diagram and layers, validated by the thickening's build digest;
a stale one is replaced).

**Other diagrams.** A search changes the surface only near the outgoing end of
one diagram's thickening, so a target that resists every search on one
diagram may fall at once on another diagram of the same knot: give it as
`target_pd` with its table name, and prove the identity separately (an
exterior isometry with the table PD's).

**The run's lines.** `[+] profile: ...` first (every setting that decides what
the run covers, as `key=value`), `[+] hop shape: ...`; per search
`[+] hop <k> <subject>: accounting: ...`, `: diagram naming: ...`,
`: breadth: ...` and `[+] hop <k>: node <id> (...)`; at the end `[+] done: ...
Target best: ...`, then `[+] GOAL MET: ...` or `[+] LOWER GOAL MET: ...` and
the outcome line a campaign parses,
`[+] <target>: N new witnesses, outcome met|expansion-limit|cpu-budget|nothing-useful|interrupted|contradiction|halted`.
All of these are frozen.

## The run directory

`work` holds, for every run: each search's `hop_<k>_n<link>/` (its pending
file `kept.csv`, when the run has a database to sign into; with a goal also `frontier.txt` and `log.txt`, whose first
line is the subject and the PD searched, in the `[[a;b;c;d];...]` spelling),
and `cobound.conf`. A goal run also writes (`driver/runrecords`,
`bounds/certificate`; names and formats frozen):

| file | what |
|---|---|
| `cascade.jsonl` | one JSON line per search (`{"hop":k,...}`: its link, crossings, surface target, wall and CPU, and where its time went: building and certifying the thickening, setup, search, rounds, drain tail, naming by route, adding its finds, naming new links, propagating; and the driver's time since the last record: choosing (with the gates' what-ifs), loading database cobordisms, the pending file and frontier), one per database load (`{"master":...}`, with its read-back, assembly, naming and propagation times and `skipped_genus`), one per refused link, and the run's own (`{"run":...}`: wall, CPU, cores, startup, loop, the searches' wall and CPU, the driver's totals, signing and the reports) |
| `profiles.jsonl` | every link's partition genera (written at every end, a halt included): name, table name, label, depth, crossings, whether searched, linking matrix, literature lower bound, Pareto entries (partition, genus, derivation), and for links of at most 5 components every partition with a positive lower bound or `forbidden` |
| `node_bounds.jsonl` | every link's identity, diagram, linking matrix, proved partition genera with whether each proof is constructive, and every lower bound with its reason's kind. Composites and links beyond the tables get bounds here that no table records |
| `lower_report.jsonl` | with `lower_report`: for every tabulated link Y that carries the target anything, what Y's literature lower bound carries to the target now (`carries`), and the most Y could ever carry (`could_carry`), each a what-if with every other lower bound forgotten; and whether Y is a special source |
| `nodes.csv` | the untabulated subjects recorded in the database, with their diagrams |
| `certificate.json` | when an upper goal is met: the proof, for the atlas's independent checker |
| `lower_certificate.json` | when a lower goal is met: the lower bound's proof as a tree of facts |

**Certificates.** An upper certificate holds the proof's derivations, children
before parents; for each cobordism, the searched diagram's PD, its layers,
the surface's triangles and the thickening's build digest (or, for a database
cobordism, its pair signature), the cobordism's shape and its maps; for each
outgoing piece, its match (method, component map, mirror and reversal); and
every link's diagram as signed Gauss data. An in-process cobordism's key is
`hop<k>#<i>`. A lower certificate holds the facts behind the lower bound, each
with its reason. The atlas's `tools/cascade_check.py` replays a certificate
without any of this code: it rebuilds each thickening from its PD and refuses
another digest, rebuilds each surface face by face with the search's own
checks (`cobound draw --faces --gauss`), finds every map with its own
diagram isomorphisms, re-proves every literature leaf's identity, and redoes
the arithmetic. Anything it cannot replay counts as a failure.

## Bounds from cobordisms

A goal run's knowledge, and each depth-0 search's judge: links, the
cobordisms, splits and sums that relate them, and immutable derivations of
every (partition, genus) bound (`bounds/cobordismgraph`). The sections below
are its theory; the paper states them in `sec:partition-genera`.

### Partition genera

A surface F in B⁴ bounded by a link (smooth or locally flat, oriented,
properly embedded, no closed components) has two numbers that matter here:
- its **partition**: which of the link's components bound the same piece;
- its **total genus**: the sum over the pieces.

A link's **partition genera** are the (partition, genus) pairs that surfaces
we can exhibit achieve (paper `def:partition-genera`). Tubing two pieces keeps
the total genus and merges their blocks (paper `lem:tubing`), so (P, g) makes
every (Q, h) with P refining Q and g ≤ h redundant, and the pairs form a
Pareto set (`PartitionGenera`). Special cases:
- the connected slice genus is the least genus over all entries, since every
  partition refines the one-block partition;
- "bounds disjoint discs" is the entry (singletons, 0);
- an unlink is all zeros.

**Why partitions, not g₄ alone.** A band move K → L is a pair of pants. If L
bounds two disjoint discs, K is slice; if L bounds only an annulus (g₄(L) = 0),
K bounds genus 1. Most genus-0 cobordisms from a knot in the atlas's database
go to links of 2 or 3 components (24,069 of 29,370), so a proof that a knot is
slice often passes through links bounding disjoint discs, which a
connected-genus rule cannot express.

### Gluing (`glue()`)

A cobordism C from link A to link B is recorded by its shape
(`CobordismShape`): its components, each boundary curve's component, and its
total genus (the database's tubed `genus`). Capping side B with a surface G
(partition Q of B's components, total genus h) gives a surface S = C ∪ G
bounded by A. Its pieces are the connected components K of the **gluing
graph**, whose vertices are C's components and G's pieces and whose edges are
B's curves (paper `lem:gluing-graph`). Each K has genus

    g(K) = Σ_{pieces p in K} g_p + b₁(K),    b₁(K) = E(K) − V(K) + 1.

Proof: gluing along circles adds nothing to χ, so

    χ(K) = Σ_p (2 − 2g_p − circles_p) = 2V − 2Σg_p − (a_K + 2E),

where a_K counts K's curves on A and each glued curve is a circle of two
pieces; also χ(K) = 2 − 2g(K) − a_K, and equating the two gives the formula.
Pieces with no curve on A close up and are discarded. `glue()` still counts
their genus, so its result is an upper bound, exact when nothing closes, and
never below its input's genus.

**Orientations.** Each component of a cobordism is oriented against the
incoming link, and each piece of G can be reversed on its own, so gluing needs
the outgoing link and the next graph link to agree only up to *global*
reversal (and mirror, which no bound here sees). Reversing a single component
is another link: `L7n1{0}` has g₄ 2 and `L7n1{1}` has 0.

### The linking condition (`linkingAllows()`)

Distinct pieces F_a, F_b are disjoint, so lk(∂F_a, ∂F_b) = F_a · F_b = 0: the
**total** linking number between any two blocks must vanish (paper
`lem:linking-condition`). This is not a per-pair condition: blocks {0} and
{1,2} with lk(0,1) = 1 and lk(0,2) = −1 are fine. It seeds lower bounds (a
forbidden partition, and each of its refinements, has no surface) and is a
contradiction gate: a derived partition that violates it means a component map
is wrong somewhere.

### Component maps

Partition genera are indexed by a link's components, so wherever one
diagram's components meet another's there must be an exact map. There are
three:
1. **Incoming to graph link.** A search's incoming curves come in T's
   component order (`DiagramDrawer::cyclesOf()`), not the PD's.
   `CobordismAssembler` draws the incoming link back from T and requires an
   orientation-preserving diagram isomorphism onto the searched diagram; its
   component map is the map (`search::certifyIncoming()`). A link without one
   is refused.
2. **Outgoing link to pieces.** `splitPieces()` keeps each component's
   origin, and `simplifyKeepingComponents()` refuses any simplification that
   changes a pairwise linking number. `simplify()` uses Reidemeister moves in
   place, and each piece keeps its own simplification, since `simplify()` is
   randomised.
3. **Piece to graph link.** `LinkRegistry::intern()` claims a match only from
   a diagram isomorphism (up to mirror and global reversal) or, for hyperbolic
   pieces, an isometry carrying meridians to meridians with one orientation
   sign (`linknaming::KernelLink`); each returns the map. A miss makes a new
   graph link, so a duplicate costs search time, never soundness.

**Split outgoing links** always get a fresh "whole" graph link joined to their
pieces by a split. Mirroring or reversing one piece changes a split link
(`K ⊔ −K` bounds an annulus, `K ⊔ K` need not), so wholes are never merged and
never searched.

### Splits

A split link is built from surfaces for its pieces placed in disjoint balls:
the whole's partition is the pieces' side by side, at the sum of their genera
(paper `lem:split-partitions`). A surface for the whole whose blocks each lie
in one piece restricts to a surface for that piece, of at most the whole's
genus. A surface mixing pieces restricts to nothing: `K ⊔ −K` bounds an
annulus, which says nothing about K.

### Sums along components

Most outgoing links that no table names are not beyond the tables but
composite: a knot summed into a link, or small links (chains of Hopf links)
summed along components. So every untabulated link is cut at its visible sum
spheres (`LinkNamer::decompose()`) into prime summands, each a graph link of
its own, joined to it by a **sum** (`CobordismGraph::addSum()`):
`pieceMap[k][c]` is the component of the whole that component c of summand k
becomes part of, a surjection whose sites must form a tree (a sphere
decomposition always does; a cycle is refused). Summands are named like any
link, so table pieces bring their literature facts.
- **Combine** (upper; paper `lem:sum-partitions`): surfaces for the summands,
  one each, give a surface for the whole whose genus is their sum and whose
  blocks are their blocks' images, merged wherever two summand components
  were summed together. A chain of Hopf links therefore bounds a connected
  genus-0 surface.
- **Lower** (paper `cor:sum-pieces`, connected bounds only):
  lower(whole) ≥ lower(piece_i) − Σ_{j≠i} (h_j + n_j − 1) over proved
  connected surfaces h_j of the other summands.

Certificates carry `sum-combine` derivations with the sum's maps, and the
checker cuts the whole's diagram at its visible spheres itself.

### Axioms

What a proof may rest on besides cobordisms (`bounds/axioms`; derivations of
kind leaf):
- the unknot's disc;
- a table entry's literature upper bound, when `literature` is on, and never
  for the target's own class (`mayUseLiteratureUpperBound()`): the target's
  own value proving the target would be circular, even through a duplicate
  graph link of it;
- a **direct** cobordism: a surface bounding the searched link alone;
- a **slice composite**: a knot the tables hold only summand by summand is
  named whole (each summand's chirality pinned), and when its summands cancel
  in concordance (`linknaming::isElementarySlice()`) it gets a genus-0 leaf
  `anchor <name>`, constructive like the unknot's (the ribbon disc is
  explicit). Never the target itself.

Every table link also gets its literature lower bound, the target's included:
lower bounds only feed the contradiction gates and the lower goal. Table names
are compared by link class (`LinkNamer::canonicalName()`), never by
`Tables::canonical()` alone, which joins by diagram only.

### The cobordism graph

The graph has cycles: every cobordism is used in **both directions** (read
backwards, a cobordism is a cobordism), an outgoing link can be a link met
before, and a search from an outgoing link can find a cobordism back to an
ancestor, improving it without searching it again. So bounds are a **fixed
point**, relaxed from a worklist whenever a cobordism, split, sum or leaf
arrives (`propagate()`).

**Derivations** are immutable and created only on a strict improvement, from
derivations that already exist, so every derivation's children have smaller
ids and every proof is a DAG. `glue()` never returns a genus below its
input's, so a cycle cannot justify a link by itself. No mutable "via" pointers:
once a child improves through its parent, they would loop.

**Contradiction gates.** Two things are contradictions: a genus below a proved
lower bound, and a partition the linking numbers forbid. They run in every
run, with literature lower bounds always loaded. A goal run judges each
search's finds in its graph once the search returns, and halts with 3, even
when the same find meets the goal. A search without a goal judges each find
in its own graph as it is kept (step 9) and halts with 2. Each means a bug or
a wrong naming somewhere. `recheck()` re-derives any derivation from its
children, a bookkeeping check (the certificate checker is the independent
one).

### Lower bounds

Reading a cobordism backwards transports lower bounds. The paper's
g₄(L₀) ≥ g₄(L₁) − g − n₀ + 1 charges n₀ − 1 for the worst gluing; capping with
a connected surface, the exact charge is n₀ − c, c the number of the
cobordism's pieces meeting L₀ (paper `lem:exact-charge`), so it vanishes when
each piece carries one component of L₀.

**The quantity.** `lower(n, Q)` bounds below the genus of every surface for
link n whose partition refines Q. A bound stored for Q applies to every
refinement of Q, and `lower()` reads it that way, so it is monotone by
construction. Tracked for links of at most 5 components (52 partitions).

**Seeds:** the literature lower bound, for every Q; the linking condition (a
forbidden Q, and each of its refinements, gets no surface).

**Across a cobordism, both directions** (paper `cor:lower-transport`): cap a
surface for X of partition P onto the cobordism. That gives a surface for the
other end Y of computable partition P′(P) and genus at most genus(F) plus what
the cap adds (glue() with genus 0), so lower(X, P) ≥ lower(Y, P′(P)) −
addition(P). By `lem:transport-monotone` this is monotone under refinement of
P, so the bound for all surfaces refining Q is the one at Q itself
(`transportedLower()`), and one cobordism suffices for every partition a
proof needs. Proof: refining the cap's partition adds vertices to the gluing
graph and no edges, so its components can only split and the addition
E − V + (components) can only fall; the induced partition only refines;
`lower(Y, ·)` is monotone and forbidden partitions are closed under
refinement.

**Splits:** a partition that never mixes pieces bounds the whole by the sum of
the pieces' bounds, a mixing one gets nothing (`K ⊔ −K`); a piece gets
lower(whole, Q ∪ P_B) − g_B for every proved surface of the other pieces.
**Sums:** as above.

**Gates.** The relaxation is capped at 200 passes: sound facts cannot pump a
bound past the truth, so non-convergence is reported as inconsistent input.
Afterwards every proved surface must satisfy every lower bound, else a
contradiction. A goal run relaxes after every search.

**Reasons.** Every raised lower bound remembers the fact that raised it
(`lowerWhy()`): a literature or linking leaf, a cobordism with the other end's
partition and the cap's addition, or a split or sum rule with the derivations
it subtracted. Raises are strict and no rule increases what it reads, so a
lower bound's proof is a tree.

**Why transported bounds rarely beat the tables.** An invariant f ≤ g₄ that
changes by at most a cobordism's charge across it satisfies
f(target) ≥ f(source) − charge, so a lower bound resting on such an f never
beats f computed at the target itself. That covers |σ_ω|/2, |τ|, |s|/2 and ν⁺
for knots, and the Murasugi–Tristram bound for links (its μ − 1 term is exactly
the slack splitting bands add), all of which the tables already record. So a
transported bound can close an entry only from a source whose bound is NOT of
this kind, and then only across a charge-0, genus-0 path. The atlas's
`data/lower_bound_sources.csv` marks such sources (`special`), which
`lower_sources` reads: of the 17,153 table entries, 2,360 knots whose lower
bound exceeds every such invariant, and 653 links above the Murasugi floor
(e.g. `L5a1`). The largest of their lower bounds is 4.

### Lower goals

`goal_lower` G keeps searching while lower(target, goal partition) < G, and
keeps a link the upper gate rejects if its best case could carry the lower goal
(`usefulLower()`): `lowerIf()` seeds the link in a copy of the graph that keeps
every bound it has, relaxes, and reads lower(target, goal). The seeds are
admissible: any transported bound is a literature seed minus non-negative
charges, so at most the largest special source's bound (4 in the tables); a
proved surface refining a partition caps it; the literature upper bound, a
connected surface, caps the coarsest. Nothing is seeded above what is known,
and a seed above a proved surface is refused. This is what makes the search
end: a chain from an open `[0;1]` target can spend at most that largest bound
less one in charge, and charge-0 searches are bounded by
`lower_max_crossings` (the chain must come back to a table entry), the search
and CPU limits and the surface-target ladder. The what-ifs for every candidate
the upper gate rejects run as one parallel batch per choice, cached per graph
version.

A met lower goal prints `[+] LOWER GOAL MET` and the chain of reasons, and
writes `lower_certificate.json`.

### Assumptions and their tests

| # | assumption | test |
|---|---|---|
| A1–A2 | partition normal form; refinement is a partial order; Bell numbers | `partitiongenera_test` `testPartition` |
| A3 | `glue()` matches an independent Euler-characteristic model on 20,000 random shapes, both sides, with closed pieces and cycles | `testGlueAgainstModel` |
| A4–A10 | the paper's rules as special cases (`lem:cobordism-inequality`, the unlink rule, a band and its reverse, the product, disconnected cobordisms, closed pieces) | `testPaperSpecialCases` |
| A11 | the linking condition is on block totals | `testLinking` |
| A12 | the Pareto set keeps exactly the non-implied entries | `testPartitionGenera` |
| A13 | `glue()` is monotone and never lowers genus, so cycles cannot self-improve | `testGlueMonotone` |
| B1 | a band to discs gives a slice knot; to an annulus, genus 1 | `cobordismgraph_test` `testBandToDisjointDiscs` |
| B2 | a cobordism found later improves an ancestor, with a well-founded proof | `testCycleImprovesAncestor` |
| B3–B4 | cycles without leaves derive nothing; cycles cannot self-improve | `testCycleWithoutLeafGivesNothing`, `testCycleCannotSelfImprove` |
| B5 | component maps matter: a permuted map gives a different, correct answer | `testComponentMapsMatter` |
| B6 | split combine and restrict; no false additivity | `testSplits` |
| B7 | the contradiction gates fire | `testContradictionGates` |
| B8 | non-bijective maps and boundaryless components are refused | `testSaturationAndMaps` |
| B9 | on 1,500 random graphs (over 300 with directed cycles) the fixed point is independent of arrival order, equals a naive closure, is saturated, and every derivation rechecks | `testRandomFixedPoints` |
| S | sums: a chain of Hopf links bounds a planar surface; a knot summed into a Hopf link (`lem:sum-along-components`(i)); the lower rule; sites forming a cycle, and pieces missing a component, refused | `testSums` |
| L1–L4 | lower bounds: transport along a concordance, the paper's reverse inequality (band and merging band), no penalty for disjoint annuli, the split rules | `testLowerConcordance`, `testLowerPaperCases`, `testLowerNoPenaltyForAnnuli`, `testLowerSplit` |
| L5 | on 600 random worlds (the upper closure taken as the truth, true minima as literature) no bound exceeds an achieved surface, and nothing contradicts | `testLowerSoundOnRandomWorlds` |
| L6 | the lower report's cleared what-if | `testLowerWhatIf` |
| L7 | transport is monotone under refinement, on 400 random worlds, every cobordism, both directions, every partition | `testLowerTransportMonotone` |
| L8 | `lowerIf()`, the what-if that keeps the graph's bounds | `testLowerIf` |
| L9 | every raised bound remembers its reason, down to a literature leaf | `testLowerWhy` |
| P1 | `profiles.jsonl`'s fields | `testPartitionGeneraFields` |
| C1–C5 | diagram isomorphisms: relabellings found with the right map; one reversed component refused; mirror and global reversal only when allowed; origins kept by `splitPieces()`; a required map realised exactly when it is a symmetry | `linknaming/tests/diagramiso_test` |
| D1, D6 | `simplify()` keeps origins and every linking number; nugatory crossings removed and nothing else | `linknaming/tests/simplification_test` |
| D2–D5 | registry: a relabelled diagram is a diagram hit with a correct map; scrambled diagrams of one link are one graph link, some found by isometry; orientation variants are different graph links; one unknot | `links_test` |
| D7 | a diagram with a nugatory crossing cannot be certified, and its reduced diagram can; a component over everything is lifted off and the pieces certify | `links_test` `testReducedDiagramsCertify`, `testLiftedDiagramsCertify` |
| E1–E2 | real searches (`10_3`, `L11n33{1}`) certify, every stored cobordism assembles, the incoming map is right, and `10_3` is proved slice constructively; the fast read-back equals the reference one | `searchcobordisms_test` |
| E3–E6 | an in-process search accounts for its surfaces as the canaries pin; every kept surface reads in process as it reads back from its pair signature, up to the incoming link's and the outgoing link's symmetries; keys distinct; a stop ends it; faces that are not a surface are refused | `search_test` |
| E7 | the build digest is the same for two builds of one diagram, and differs for another diagram or layer count | `search_test` `testBuildChecksum` |

Breaking `glue()` (dropping b₁, or off by one), checking linking pairwise, or
skipping dominance makes `partitiongenera_test` fail; dropping the genus
subtraction, or miscounting a split, makes `cobordismgraph_test` fail.

## The atlas solver (`cobound solve`)

`solve` reads the database (without pair signatures), the tables and its
other inputs, and re-derives every verdict (`solver/`). For a connected
genus-g cobordism from L₀ (n₀ components) to L₁ (n₁):

    g₄(L₀) ≤ g₄(L₁) + g + n₁ − 1        g₄(L₀) ≥ g₄(L₁) − g − n₀ + 1

relaxed over every cobordism in both directions to a fixed point
(`propagate()`; the derivation is in `solver/solver.h`), seeded with the
unknot, every n-component unlink mentioned (n discs tube into a connected
planar surface: no component penalty) and every slice composite mentioned, and
grounded by direct cobordisms and certified bounds. A disconnected find
counts: its `genus` is its tubed genus.

- **Which outgoing links bear a bound** (`outgoingBearsBound()`): a knot
  (Gordon–Luecke: the complement determines it up to mirror, which g₄ does not
  see), a proved unlink, or an outgoing link proved per cobordism
  (`outgoing_resolutions`). Any other link of 2 or more components bounds
  nothing: one complement belongs to infinitely many links, with different
  slice genera, and a max/min over a base name's orientation variants is not
  a bound over the real possibilities. Such cobordisms are still recorded; the
  solver declines them. A `complement:` description bears nothing.
- **The reverse direction** (bounding the outgoing link from the subject)
  needs the outgoing link's identity: a knot, or a name (`outgoing_names_file`).
- **Splits, composites and sums.** A split outgoing link's upper bound comes
  from its factors' surfaces placed side by side; its lower bound is not
  additive (`K ⊔ −K` bounds an annulus): g₄(A # B) − 1 ≤ g₄(A ⊔ B) ≤ g₄(A # B),
  split unknots drop out, and with f knot factors the lower bound is the
  composite rule's minus (f − 1). A composite knot's g₄ is at most the sum of
  its summands' and at least each summand's minus the others'. With
  `sum_rules`, sums along components and splits with link factors are bounded
  from their pieces.
- **Statuses** (`judge()`): `verified` (constructive, matching the literature),
  `verified-assisted` (resting on another name's literature value),
  `improved`, `pinned`, `bounded`, `unresolved`, and `skipped` for a table row
  above `max_crossings` with no bound. Derived bounds inconsistent with the
  literature, or with each other, are a contradiction: `solve` halts with 2
  after writing the verdicts.
- **Inputs.** `name_aliases`, `outgoing_resolutions`, `outgoing_names_file`,
  `link_classes`, `certified_bounds`, `knot_symmetry` and `sum_rules` are
  opt-in, each applied to a separate copy of the cobordisms (the database keeps
  what the search observed). `link_classes` is needed with
  `outgoing_names_file`: names are written as their class's canonical name, so
  bounds land on the wrong member's row without it.
- **The verdicts file** (`solver/verdicts`), one line per table row, columns
  frozen:
  `knot,resolved_genus,status,witness_kind,witness_pairsig,via_knot,via_edge_genus,depends_on,literature_lo,literature_hi,derived_lo,derived_hi,witness_basis,tubed,searched_faces,search_outcome,exhausted_depth`.
  `solve` rewrites every line it judges; a run writes only each searched row's
  search record (`searched_faces`, `search_outcome`, `exhausted_depth`, the
  last never lowered). `via_knot` and `depends_on` name a cobordism's other
  end, which is often its subject.

## Frozen formats

Everything the atlas reads or keys on keeps its bytes until the atlas changes
its readers with it. In `cobound/frozen.h` are the tokens that still spell a
retired term (`hop_`, `node`, `witness`, `far side`, `cascade:`, `row`, ...),
named for what they are; the rest are literals at the one place each is
written. Frozen:

- the database's 13 columns and its `.rows.csv` sidecar;
- the verdicts file's columns (the atlas's `verify_genus_v2.csv`, each
  campaign's verdicts shard);
- the run directory: `hop_<k>_n<link>/{kept.csv,frontier.txt,log.txt}`
  (`kept.csv`'s 16 columns, read by position), `cascade.jsonl`,
  `profiles.jsonl`, `node_bounds.jsonl`, `lower_report.jsonl`, `nodes.csv`;
- `certificate.json` and `lower_certificate.json`, every key;
- the log lines the atlas parses or keys on: the accounting line, `N new
  witnesses, outcome X`, `[+] exact far-side names: N table entries`,
  `breadth:`, `[+] hop <k> ...`, `[!!] hop <k>: surface accounting failed --
  ...`, `[+] profile:`, `hop shape:`, `EXHAUSTIVE to`, `boundary processing:`,
  `Target best:`, `GOAL MET`, `[+] Searching X`, the `diagram naming:` line's
  `far sides drawn` and counters, and the progress block;
- `surface_stats.csv`'s columns and the surface log's;
- `cobound draw`'s `ROW` and `W` lines, `cobound name`'s lines, and
  `cobound meridians`' records;
- every literal of the frontier fingerprint (`surfer/README.md`, "Frontiers").

## Tests

`ctest` in `build/utils/surfer/cobound_part`. The scripts take the `cobound`
binary and `tests/data/` (the canaries' rows, small tables through 6
crossings, and real cobordisms of two searches).

| test | pins |
|---|---|
| `canaries_test` | nine exhaustive searches at cap 3 (`tests/data/canaries.csv`: `3_1`, `6_1`, `8_8`, `8_20` and five links), census-free: each accounting line (recorded and duplicate summed) is `canaries.expected`'s, byte for byte; two goal runs (`3_1` goal 1, `L2a1{0}` goal 0): accounting, best genus and outcome are `goal_canaries.expected`'s |
| `search_defaults_test` | one set of defaults (2 layers, collared through both, `proper`): a config that states no shape runs the spelled-out one; `resolve_unlinked`, `work`, the tables and `targets` have no default; retired options are refused; `outgoing_names = 0` is refused under either spelling; `cobound.conf` read back runs the same search |
| `config_test` | every key's default in every context; files, `--set`, `none`, required keys, types, unknown and inapplicable keys, old spellings, the run's context, the effective configuration read back |
| `given_diagram_test` | a goal run searches a table target on its PD as written; a table row that does not certify is refused (`build-failed` without a goal, exit 2 with one); an untabulated target falls back to its simplified diagram, logged; what a target's searches record (`log.txt`, `kept.csv`, `.rows.csv`, `certificate.json`) is the given PD in the `[[a;b;c;d];...]` spelling |
| `frontier_pending_test` | a frontier records its pending file and length; a resume is refused, by name, until `sign` has signed that far, unless the file is the run's own; a moved work tree keeps its resume; a cut or garbage frontier is refused and the search starts afresh |
| `interrupted_outcome_test` | SIGINT and SIGTERM end the running search cleanly, at depth 0 and with a goal: `interrupted`, never `exhausted`; no further search; exit 0 without a goal, 1 with one |
| `unaccounted_search_test` | an imbalanced search (`SURFER_TEST_UNACCOUNTED`) ends only itself: `unaccounted`, no frontier, the run goes on, exit 2 at depth 0, the goal's own code with a goal; an impossible state (`SURFER_TEST_IMPOSSIBLE`) halts with 2 after writing what was found, even when the goal is met |
| `contradiction_halt_test` | the contradiction gates run in every run: with `L2a1{0}`'s 4-genus rewritten to 1, a depth-0 run halts with 2 after writing its cobordisms, and a goal run with 3 although the same find meets the goal |
| `io_failure_test` | a failed output write (a frontier path that is a directory, `/dev/full`, a file-size limit under the surface log, an unwritable search directory) ends its search as `io-error`: no frontier, no exhaustion claim, its finds signed, exit 2; never a signal |
| `name_independence_test` | the same exhaustive searches with every name and description perturbed (`SURFER_TEST_PERTURB_NAMES`) accept and account for exactly the same surfaces; only the recorded/duplicate split may move; the full naming route ran |
| `goal_layout_test` | a goal run explores the same links whatever the heap layout (five work-path lengths); needs the atlas's tables |
| `meridians_order_test` | `cobound meridians` answers each record on its own: the records reversed give the same output |
| `surface_log_searches_test` | the surface log over two searches on one thread logs both |
| `search_test` | see E3–E7 above; a request whose PD the thickening does not carry is refused |
| `searchcobordisms_test` | see E1–E2 above |
| `incomingmap_test` | the incoming map lands exactly on L × {0}, one closed directed curve per component, despite T's automorphisms; no searchable face but the seed's touches the incoming side; the bare collar's boundary matches per component; `certifyIncoming()` accepts the PD's own diagram and refuses a mirror and a nugatory kink (with a table as argument, every row) |
| `preconditions_test` | `conditionFor()`; the rejection names; every accounting bucket, failure message and the exact accounting body; the watchdog's order and reasons; `gateSurface()` on exhaustive cap-3 searches (`3_1` accepts 1,752; `L2a1{0}` 795 of 945, 150 the other orientation); a pair signature from captured faces equals the captured one |
| `outgoingnamer_test` | on real thickenings: the bare collar's outgoing link is the incoming link (a table knot by name, a table link per edge set as a variant and per surface as its own variant), straight from the diagram; a curve round one triangle is an unknot; the complement namers agree on an unlink with or without the census |
| `solver_test` | the solver: component-count terms, which outgoing links bear a bound, the reverse direction, splits, composites, anchors, cycles, support sets, statuses, `splitBoundary()`, per-component orientation, the cobordism identity, named outgoing links, `sum_rules` |
| `partitiongenera_test`, `cobordismgraph_test`, `links_test`, `axioms_test` | the assumptions above; `axioms_test` the literature leaf policy |
| `database_test` | a cobordism round-trips (commas quoted, `resolved_vertices` empty when 0); appends keep every byte and drop pair signatures; a torn last line is ignored, then cut; a 12-column file is refused for appending |
| `pending_test` | `kept.csv` round-trips; one database line per identity, each with the pair signature taken in the searched thickening; signing again, or against a database holding them, appends nothing; a torn last line is skipped |
| `fromdatabase_test` | read-backs serialise and parse exactly; a cache is read back whole by the next run; another build digest is replaced; a torn line is cut |
| `cobordismkey_test` | `cobordismKey()` against Python's `hashlib` |
| `thickening_rows_test` | the thickening's component count against each test-table row's name and Regina's count |
| `appendonly_test`, `json_test`, `parallelfor_test`, `timers_test` | whole appends under a lock with torn lines cut; the JSON writer's frozen output; `parallelFor()` calls each index once; the clocks |

Not in ctest: `bench_search.sh` and `compare_surface_sets.sh`
(`surfer/README.md`, "Performance").
