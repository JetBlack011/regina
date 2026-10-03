# cascade: goal-directed chained searches

`verifyslicegenus` searches every row the same way, and each search is
shallow: the collar `L × [0,2]` plus at most 4–5 triangles near `∂₊`. Depth
comes only from the solver chaining witnesses across other *table rows*.
The cascade works on one target instead. It searches the target, re-searches
only the far sides that could still improve the target's bound, and
re-diagrams and re-triangulates each far side from its simplified diagram
between hops. The plan is `~/.claude/plans/i-ve-put-you-into-idempotent-hamming.md`.

| file | what |
|---|---|
| `partition.{h,cpp}` | partitions of a link's components (restricted growth strings) |
| `profile.{h,cpp}` | cobordism shapes, `glue()`, the linking condition, Pareto profiles |
| `proofgraph.{h,cpp}` | nodes, witness and split edges, immutable proof records, the fixed point |
| `diagramiso.{h,cpp}` | diagram isomorphisms of signed Gauss data, with their component maps |
| `nodes.{h,cpp}` | the node registry: exact identity and component maps; `simplifyKeepingComponents()` |
| `hopedges.{h,cpp}` | a hop row's certificate, and each witness (`add()`, from its pair signature) or surface read in process (`addRead()`) turned into edges and nodes |
| `hoprunner.{h,cpp}` | a hop searched in process on `../rowsearch.h` (see "In-process hops") |
| `leaves.{h,cpp}` | literature leaves: table values, and which a proof may use |
| `cascadesearch.cpp` | the driver: hops in process (default) or as `verifyslicegenus` children (`--hop-mode child`) |
| `tools/cascade_check.py` | the independent checker; `tools/cascade_check_test.py` tests its own diagram code |
| `tests/` | one test per module; `data/` holds real witnesses from the Phase 0 hops and small tables |

The first four files have no Regina dependency.

Still to come:
- the SnapPy oracle (DG `PlausibleKnots`, `RibbonLinks`, HFK) for leaf facts
  beyond the tables;
- switching hops early, depth first, when a far side could beat the demand
  (`HopSearcher::run()`'s `stop` is the hook).

## Running it

One program, `cobound` (`main.cpp`, `driver/commands.h`), replaces
verifyslicegenus, cascadesearch, farsidename, farsidediagram and
peripheral_slopes (phase 5):

| command | was | does |
|---|---|---|
| `cobound run` | verifyslicegenus (a search run), cascadesearch | a search from each target: without a goal each once (`driver/sweep.h`), with one (`goal_genus` or `goal_lower` set) outwards from the target until the goal has a proof or the limits are spent (`driver/scheduler.h`) |
| `cobound solve` | verifyslicegenus --solve-only | every verdict re-derived from the database; never writes it |
| `cobound sign` | cascadesearch --sign-only | a work directory's pending cobordisms signed into the database |
| `cobound draw` | farsidediagram | stored cobordisms' outgoing links drawn; takes farsidediagram's flags |
| `cobound name` | farsidename | stored cobordisms' outgoing links named; takes farsidename's flags |
| `cobound meridians` | peripheral_slopes | `sig`, `dump`, `dump-link`, `dump-subset`, `slope`, byte-compatible |

`run`, `solve` and `sign` read a configuration (`driver/config.h`): `key =
value` lines from `--config FILE` (any number), then `--set key=value` (any
number; the last wins). One schema serves every command; each key's default
reproduces the retired tool's (without a goal verifyslicegenus's, with one
cascadesearch's; `config_test` checks every one), and `resolve_unlinked`
(with or without a goal) and `work` have none. `cobound help keys` lists
every key with its type, its default in each command and the option it
replaces. A run writes the configuration it ran with to `<work>/cobound.conf`;
read back with `--config`, it runs the same search.

A goal run, as Pool A ran it (`D` = the atlas's `data/`):

```
target_pd = <PD>
target_name = 10_27
work = <dir>
knot_table = $D/4d_smooth_slice_genus_13_crossings_pd_codes.csv
link_table = $D/links_4d_smooth_slice_genus_11_crossings_pd_codes.csv
knot_symmetry = $D/knot_symmetry.csv
census = <private census copy>
goal_genus = 1
literature = 0
threads = 14
resolve_unlinked = 1
surface_target = 50000
max_surface_target = 200000
max_searches = 8
```

| key (was) | what |
|---|---|
| `resolve_unlinked` (`--resolve-unlinked` / `--no-resolve-unlinked`) | required, with no default (plan divergence 3): whether surfaces whose only self-intersections are unlinked count (the campaign's: on). It changes what a surface target counts and every frontier's fingerprint |
| `goal_genus` (`--goal-genus`), `goal_partition` (`--goal connected\|disjoint`) | the bound to prove: connected g₄ ≤ g (the coarsest partition), or disjoint surfaces (singletons); setting `goal_genus` (or `goal_lower`) is what makes a run goal-directed |
| `goal_lower`, `lower_sources`, `lower_max_crossings` (`--goal-lower`, `--lower-sources`, `--lower-max-crossings`) | also stop once lower(target, goal partition) ≥ G ("Lower-bound mode" below); the sources file (the atlas's `data/lower_bound_sources.csv`) is required, and a node kept only for the lower goal is never expanded above n crossings (default 16) |
| `literature` (`--constructive` / `--literature`) | whether table values may be leaves (never the target's own) |
| `master_witnesses` (`--master-witnesses`) | a read-only database whose cobordisms are free edges for table nodes, loaded for a node just before it is searched (`bounds/databasecobordisms.h`; never without a goal) |
| `surface_target`, `max_surface_target` (`--hop-surfaces`, `--max-hop-surfaces`) | a search's surface target, and the most a revisit may raise it to |
| `max_searches`, `cpu_budget` (`--max-expansions`, `--cpu-budget`) | stop after this many searches, or this much search CPU (checked between searches) |
| `search_seconds` (the fixed 7200 s per hop) | each search's wall-clock backstop |
| `strategy` (`--strategy best\|dfs\|bfs`) | which node to search next |
| `max_crossings` (`--max-crossings`) | never search a node whose diagram has more than n crossings (default 24 with a goal) |
| `max_faces`, `iddfs_start`, `iddfs_iterations`, `iddfs_step`, `root_budget_start`, `root_budget_growth`, `layers` (`--hop-max-faces`, ...) | each search's shape; with a goal the defaults are the campaign's (cap 5, IDDFS 2 rounds from 4 step 1, root budget 840 doubling, 2 layers) |
| `pending_surface_cap`, `petal_cache_limit`, `boundary_signature_cache_limit`, `complement_cache_limit` (`--hop-pending-cap`, `--hop-petal-cache`, `--hop-boundary-cache`, `--hop-recognition-cache`) | hosts.conf's per-host limits (with a goal 20M, 12M, 1M, 1.5M) |
| `cobordisms`, `run_name` (`--witness-store`, `--run-name`) | record every kept surface for the atlas ("Witnesses for the atlas" below) |
| `dedupe_against` (`--dedupe-against`, comma-separated) | read-only databases whose identities the store must not repeat (the master) |
| `cobound sign` (`--sign-only`) | the store step alone, over `work`'s `hop_*/kept.csv` (a killed run) |
| `hub_degree`, `hub_surfaces` (`--hub-degree`, `--hub-surfaces`) | hub breadth (John, 2026-09-29): a node chosen for search with at least D witness edges is searched once at N surfaces (when above the current budget), for many more first-level far sides where many routes meet, as wide rows do; logged `[+] hub: ...` |
| `lower_report` (`--lower-report`) | write `lower_report.jsonl`: what each tabulated node's lower bound carries to the target ("Lower bounds" below) |

Retired with the old executables (plan, "Retired"): the child hop mode
(`--hop-mode child`, `--verifyslicegenus`), `--master-loads eager` (loads are
lazy), `--verbose`, and verifyslicegenus's `--cone`, `--sweep-time-limit`,
`--harvest-quiescence`, `--skip-drain-on-timeout`, `--no-diagram-naming`,
`--retriangulate-links`, `--retriangulate-height`,
`--retriangulate-candidate-budget`, `--no-simplify`, `--boundary-tally-cap`,
`--iddfs-final-threads`, `--rewrite-witnesses`, `--harvest` and
`--research-settled` (every run harvests and searches every target).

**Table classes.** A table name stands for its class: the variants of its
base that are one oriented link up to mirror and global reversal, joined by
coinciding diagrams or by an isometry carrying meridians to meridians with one
sign (`ExactNamer::canonicalName()`, the classes of
`data/table_link_classes.csv`). Nodes are named by that class, so every
comparison of a table name with a node's name goes through `classOf()`: the
target check (a target PD the namer proves to be another class is refused),
the literature guard (no leaf in the target's class, under any name) and
`--master-witnesses` (every class member's rows). `ExactTables::canonical()`
joins by diagram only and must not be compared with a node's name: until
`53e65b14b` it was, and 27 of `2026-09-open2`'s targets were refused as
their own class-mates (`L11a397{1;0}` "is" `L11a397{0;0}`). Each class-mate
was searched in its own right, so nothing was lost. `cascade_check.py`
applies the same rule independently: it refuses a literature leaf in the
target's class, by the class table or by its own isometry test, and a scan
of every certificate on both hosts (1,618) found none that used one.

It writes `cascade.jsonl` (one line per hop, with its phase timers),
`driver.log` and, when the goal is met, `certificate.json` for
`tools/cascade_check.py`. The driver log starts with a `[+] profile:` line
(every setting that decides what the run covers, as `key=value`), gives each
hop's accounting as `[+] hop <k> <subject>: accounting: ...` in
`verifyslicegenus`'s shape, then its far-side naming as `[+] hop <k>
<subject>: diagram naming: ...` (`NamingStats::summary()`, also in the hop's
`log.txt`: counts by outcome, time by route, and the slowest single name,
which is what can leave one thread draining alone; `cascade.jsonl` has the
times as `naming_diagram_s`, `naming_fallback_s`, `naming_exact_s` and
`naming_slowest_s`), and ends with the line a campaign parses,
`[+] <target>: <n> new witnesses, outcome <met|expansion-limit|cpu-budget|nothing-useful|contradiction>`
(n: the witnesses the store gained).

**`profiles.jsonl`** (since 2026-09-29, for the atlas page). The run's end
writes it, and so does a contradiction halt. It has one line per node of the
proof graph:
- `node`, `name`:
  - the target's name, a proved table name, `Unknot` or
    `<k>-component unlink` for a crossingless node;
  - otherwise the store's `cascade:<run>/<target>/n<id>`.
- `table`, `label`, `depth`, `crossings` (absent for a split far side's
  whole, which has no diagram of its own), and whether it was `searched`;
- `ProofGraph::profileFields()`:
  - its `linking` matrix and literature `genus_lower`;
  - its Pareto `entries` (partition, genus, record);
  - `lower`, for nodes of at most five components: every partition with a
    positive lower bound, or `forbidden` where the linking numbers rule it
    out.

A node's profile is its whole state in "Profiles" below. The certificates
keep only the partitions a proof used. `proofgraph_test`
`testProfileFields` pins the format.

## Witnesses for the atlas (`keptstore.h`, since 2026-09-28)

Every surface a hop keeps is a real cobordism, and one a later search's
propagation may use even when this run's goal is not met. With
`--witness-store`, each is recorded as `verifyslicegenus` would have
recorded it (13 columns, `witnessstore.h`, the same code the sweep appends
through):

- **When.** Each hop's search appends its kept surfaces to
  `hop_<k>_n<node>/kept.csv` as it runs (`cobordisms/pending.h`,
  PendingWriter: fsynced once a minute while the search and its drain
  run, at once at a first constructive find, and the rest when the search
  ends; plan
  divergence 7): every witness column but the pair signature, then the
  faces, row PD and layers. The run's end (every exit, a contradiction
  included) signs the ones whose witness identity is new and appends them to
  the store. `verifyslicegenus` does the same since divergence 7: each row's
  pending file is `<work>/hop_<k>_n0/kept.csv` (`--work`, default
  `<cobordisms>.pending`), signed into `--cobordisms` at the run's end.
  Signing is almost all the row's ambient (one `PairSigContext` per hop row,
  ~50 s for a 10-crossing row), so it is done once per row, rows in
  parallel, after the search. A killed run keeps its `kept.csv`, and
  `--sign-only` stores it later.
- **Which.** One per `cobordismgraph::witnessIdentity()`, the sweep's rule:
  identities already in the store or in a `--dedupe-against` file (the
  master) are skipped, before signing. The append holds `flock` on
  `<store>.lock` and re-reads the store's identities under it.
- **Under what name.** A hop's subject is the target's own name, a node's
  pinned table name (exact naming, a class up to mirror and global reversal,
  which g₄ does not see), or `cascade:<run name>/<target>/n<node>`, the
  store's scheme for untabulated nodes (`tools/cascade_record.py`). Node ids
  restart in every run, so the old `cascade_n<node>` would have put
  different links under one name. The `cascade:` subjects are listed in
  `<work>/nodes.csv` (the store's `nodes.csv` columns).
- **Other columns.** `other_candidates` comes from the name tables
  (`NameTable::candidates()`), as the sweep fills it; `source_row` is the
  subject; `thicken_layers` 2; `max_faces` the hop's cap.
- **Which row.** A pair signature's ambient is the hop row's thickening,
  built from the node's own diagram, which is not the table's PD even for a
  table-named node. So each appended witness also gets a line in
  `<store>.rows.csv`: its key (sha1(pairsig)[:12]), layers and row PD.
  `farsidename` accepts a PD in place of a row name. The atlas's
  `run_farsidename.py --row-pds` feeds it from the sidecar, and without the
  sidecar such witnesses cannot be redrawn.

Checked (2026-09-28): `keptstore_test` (a real exhaustive 3_1 hop: kept.csv
round-trips; one line per identity, each pair signature the one taken in the
searched thickening; storing again, or against a store holding them, appends
nothing). On a `10_27` run, all 24 stored witnesses were redrawn from their
pair signatures by `farsidename`: 19 give exactly the in-search name, and 5
a sharper one (raw isomorphism signatures that are `3_1#m5_2`, `3_1#4_1` and
`3_1#m3_1`, the census name `L108014` = `8_14`, and a `diagram:` name that is
the Hopf link), which is why the atlas names every merge exactly.

## Component maps

A profile is indexed by a link's components, so every place where one
diagram's components meet another's needs an exact map. There are three:
1. **Row to node.** A witness's incoming curves come in knotbuilder's
   component order (`DiagramDrawer::cyclesOf()`, ordered by first edge), not
   the PD's. `HopAssembler` draws the row's own link back from its
   triangulation. It requires an orientation-preserving diagram isomorphism
   onto the node's diagram, and uses that isomorphism's component map. This
   is also the per-node certificate that the searched triangulation carries
   the node's oriented link. A row without it is refused.
2. **Far side to pieces.** `exactnaming::splitPieces()` keeps each
   component's `origin`, and `simplifyKeepingComponents()` refuses any
   simplification that changes a pairwise linking number. Regina's
   `simplify()` uses r1, r2 and r3 in place; `simplifyExhaustive()` is never
   used. Each simplified piece stays with its own match, because
   `simplify()` is randomised.
3. **Piece to node.** `NodeRegistry::intern()` claims a match only from a
   diagram isomorphism (up to mirror and global reversal) or, for hyperbolic
   pieces, a SnapPea-kernel isometry carrying meridians to meridians with one
   orientation sign (`exactnaming::KernelLink`). Each returns the component
   map. A miss makes a new node, so a duplicate costs search time, never
   soundness.

**Split far sides** always get a fresh "whole" node joined to their pieces.
Mirroring or reversing one piece changes a split link: `K ⊔ −K` bounds an
annulus, `K ⊔ K` need not. So wholes are never merged, and never searched.

## Profiles

A surface F in B⁴ bounded by a link (smooth or locally flat, oriented,
proper, no closed components) has two numbers that matter here:
- its **partition**: which of the link's components bound the same piece;
- its **total genus**: the sum over pieces.

The **profile** of a link is the set of (partition, genus) pairs that
surfaces we can exhibit achieve. Tubing two pieces together keeps the total
genus and merges their blocks (paper lem:tubing), so a profile is a Pareto
set. (P, g) makes every (Q, h) with P refining Q and g ≤ h redundant.
Special cases:
- the connected slice genus is the least genus over all entries, since every
  partition refines the one-block partition;
- "bounds disjoint discs" is the entry (singletons, 0);
- an unlink is all zeros.

**Why partitions and not LinkInfo's g₄.** A band move K → L is a pair of
pants. If L bounds two disjoint discs, K is slice. If L bounds only an
annulus (g₄(L) = 0 in LinkInfo's sense), K bounds genus 1. On the master,
every genus-0 witness from a knot row goes to the knot itself or to a
2- or 3-component link, never to a different knot. So a sliceness chain must
pass through links bounding disjoint discs, which the solvers'
connected-genus rules cannot express.

## Gluing (`glue()`)

A cobordism C from link A to link B is recorded by its shape:
- its components;
- each boundary curve's component;
- its total genus (the witness file's tubed `genus`).

Capping side B with a surface G (partition Q of B's curves, total genus h)
gives a surface S = C ∪ G bounded by A. Its pieces are the connected
components K of the *gluing graph*:
- vertices are C's components and G's pieces;
- edges are B's curves.

Each K has genus

    g(K) = Σ_{pieces p in K} g_p + b₁(K),    b₁(K) = E(K) − V(K) + 1.

Proof. Gluing along circles adds nothing to χ. So

    χ(K) = Σ_p (2 − 2g_p − circles_p) = 2V − 2Σg_p − (a_K + 2E),

where a_K counts K's curves on A, and each glued curve is a circle of two
pieces. Also χ(K) = 2 − 2g(K) − a_K. Equating the two gives the formula.

Pieces with no curve on A close up and are discarded (lem:tubing). `glue()`
still counts their genus, so its result is an upper bound: exact when
nothing closes, and never below its input genus.

**Special cases.** `profile_test` checks each of these:
- the paper's `g₄(L₀) ≤ g₄(L₁) + g + n₁ − 1` (both surfaces connected);
- the unlink rule (no penalty for disjoint discs);
- product cobordisms as the identity;
- a band with discs versus with an annulus.

**Orientations.** Each component of a witness is oriented against the row
(`farside::incomingFlips`). Each piece of G can be reversed independently.
So gluing needs the far side and the next node to agree only up to *global*
reversal (and mirror, which every bound here is invariant under). Reversing a
single component is a different link: `L7n1{0}` has g₄ 2 and `L7n1{1}` has 0.

## The linking condition (`linkingAllows()`)

Distinct pieces F_a, F_b are disjoint, so lk(∂F_a, ∂F_b) = F_a·F_b = 0. The
*total* linking number between any two blocks must vanish. This is not a
per-pair condition: blocks {0} and {1,2} with lk(0,1) = 1 and lk(0,2) = −1
are fine. The proof graph uses it only as a contradiction gate: a derived
partition that violates it means a component map is wrong somewhere.

## Leaf facts

What a proof may rest on besides witnesses (`RecordKind::leaf`):
- the unknot's disc (the registry's unknot node);
- a table entry's literature upper bound (`leaves.h`), with `--literature`
  and never the target's own class, so its own value cannot prove it;
- a direct witness: a surface bounding the row alone;
- **a slice composite** (since 2026-09-30, John's observation on 10_99's
  node 68, which was `3_1#m3_1` and got a 335 CPU-s hop instead of a disc):
  a knot the prime-piece namer leaves untabulated is handed to the
  whole-diagram namer (`ExactNamer::name`: cut at the visible sum spheres,
  each summand named with its chirality pinned); an exact, pinned composite
  name is recorded, and when its summands cancel in concordance —
  `cobordismgraph::isElementarySlice`, the atlas solver's rule: the
  allowlist `3_1#m3_1`, `4_1#4_1`, or `K # m(Kʳ)` pairs by the summands'
  symmetry types from `--knot-symmetry` — the node gets a genus-0 leaf
  `anchor <name>`, constructive like the unknot's (the ribbon disc is
  explicit). Never the target itself. The checker replays it: cuts the
  diagram at its visible spheres (`connected_summands`), proves each
  summand to be its named knot with that chirality (fsid), and redoes the
  cancellation from the symmetry table (`elementary_slice`);
  `cascade_check_test.py` pins the cut on the square and granny knots and
  the rule on `8_17#m8_17` (not slice: `8_17` is not invertible) against
  `8_17#mr8_17`. Measured (yoga, 8 threads, 2026-09-30): `10_99`
  constructive at 2M-surface hops met its goal in **one hop** (224 s,
  1,926 CPU-s), node 18 = `3_1#m3_1` anchored as it was interned, against
  the 2026-09-28 retry's four hops (2M, 2M, 4M, 4M surfaces; ~3,880
  CPU-s), whose last hop proved the square knot slice by search; the
  certificate CERTIFIED, and renaming the anchor to `3_1#3_1`, `4_1#4_1`
  or `m3_1#m3_1` is refused (a summand's chirality is pinned by its Jones
  polynomial against the table diagram's, since an exterior isometry sees
  a knot only up to mirror).

Lower bounds have their own leaves: literature lower bounds (every node's,
the target's included, as contradiction gates) and the linking condition;
see "Lower bounds".

## Sums along components (since 2026-09-30)

The hub runs of 2026-09-30 showed why so few far sides got names: of
1,411 prime far-side pieces, 1,143 were untabulated, and 1,128 of those
had at most 11 crossings — they were not beyond the tables but *composite*:
507 with a visible knot summand, 610 chains of small links (Hopf links
summed along components), which the tables never list and which the
prime-piece namer (`ExactNamer::identify`) cannot describe. John: the
cascade must name and reason about these as the atlas solver does.

Now every untabulated node is cut at its visible sum spheres
(`ExactNamer::decompose`, exposed for this) into prime summands, each
interned as a node, and joined to it by a **sum edge** (`SumEdge`,
`ProofGraph::addSum`): `pieceMap[k][c]` is the component of the whole that
component c of summand k becomes part of, a surjection whose sites must
form a tree (a sphere decomposition always does; a cycle is refused). The
summand nodes are named like any other, so table pieces bring their
literature leaves. Rules:
- **Combine** (upper; paper `lem:sum-partitions`): surfaces for the
  summands, one each, give a surface for the whole whose genus is their
  sum and whose blocks are their blocks' images merged wherever two piece
  components were summed together (a boundary connected sum at each site).
  A chain of Hopf links therefore realizes a connected genus-0 surface, a
  knot summed into a link realizes the link's partition at the genera's sum.
- **Lower** (paper `cor:sum-pieces`, connected bounds only): `lower(whole)
  ≥ lower(piece_i) − Σ_{j≠i} (h_j + n_j − 1)` over proved connected surfaces
  `h_j` of the other summands.
`proofgraph_test` `testSumEdges` pins the Hopf chain, the paper's (i) as a
special case, the lower rule and both refusals. Certificates carry
`sum-combine` records with the edge's maps; the checker cuts the whole's
diagram at its visible spheres itself (`link_summands`), matches each
summand to its piece node under a map sending its components to the
recorded origins, checks the sites form a tree, and redoes the arithmetic
(`cascade_check_test.py`: the chain, and the record with a wrong genus,
partition or a cyclic map refused). Measured on `L8n7{0;0;1}` (yoga, 8
threads, 600 CPU-s): 7 far sides decomposed into sums of `L2a1{0}`, each
given a genus-0 profile; the target's bound unchanged (its witnesses to the
chains are not annuli), no contradiction.

## The proof graph

The cascade's graph has cycles:
- every witness is used in **both directions** (a cobordism read backwards is
  a cobordism);
- a far side can be a link met before;
- a search from a far side can find an edge back to an ancestor, which
  improves it without re-searching the ancestor.

So bounds are a **fixed point**, relaxed from a worklist whenever an edge or
a leaf arrives (`propagate()`).

**Records** are immutable and created only on a strict improvement. They are
derived from records that already exist, so every record's children have
smaller ids and every proof is a DAG. `glue()` never returns a genus below its
input's, so a cycle cannot justify a node by itself. Never store mutable
"via" pointers: once a child improves through its parent, they loop.

**Split edges.** A split link can be built from surfaces for its pieces,
placed in disjoint balls. A surface for the whole whose blocks each lie in one
piece restricts to a surface for that piece, of at most the whole's genus. A
surface mixing pieces restricts to nothing: `K ⊔ −K` bounds an annulus, and
that says nothing about K.

**Contradiction gates.** Two things halt the driver:
- a genus below a proved lower bound;
- a partition the linking numbers forbid.

They run in every run (plan divergence 6), with literature lower bounds always
loaded. A goal run judges each hop's finds in its own graph once the hop
returns, and halts with 3, even when the same find meets the goal. A search
without a goal (`verifyslicegenus`) gets a graph of its own
(`bounds/searchjudge.h`, the outside facts by `bounds/axioms.h`'s
`NodeAxioms`, shared with the cascade): the searched link with its literature
lower bound, and each new witness entered as it is kept, on a thread of its
own (`SearchRequest::judge`). A contradiction ends the search, the row's
witnesses are written, and the run halts with 2 (the FATAL banner). That graph
replaces the old in-search check (a witness's implied upper bound below the
row's literature lower bound). It reads finds on the row collared through
every layer, so a search needs `--collar-layers` equal to `--thicken-layers`
and no `--cone`, and it needs the knot and link tables.

The same graph says when a depth-0 search is constructive (plan divergence
10): once a proof of the searched link's literature lower bound rests on no
literature value, the search prints `[+] X: CONSTRUCTIVE witness found --
reaches genus g (literature [lo, hi]). Checkpointing now.` (the text is
frozen; `-- ACHIEVED` in the progress block with it) and checkpoints at once.
It reads no database, so a link constructive only through another row's
cobordisms is not announced; that is the solve's to find. The atlas solver
runs only for `--solve-only`: a search run writes each row's search record
(`exhausted_depth`, `searched_faces`, `search_outcome`) and leaves its status
and bounds as they were, and prints no verdict line and no totals.

## Assumptions and their tests

| # | assumption | test |
|---|---|---|
| A1–A2 | partition normal form; `refines` is a partial order; Bell numbers | `profile_test` `testPartition` |
| A3 | `glue()` matches an independent Euler-characteristic model on 20,000 random shapes, both sides, with closed pieces and cycles exercised | `testGlueAgainstModel` |
| A4–A10 | the paper's rules as special cases; closed pieces | `testPaperSpecialCases` |
| A11 | the linking condition is on block totals | `testLinking` |
| A12 | Pareto invariants | `testProfile` |
| A13 | `glue()` is monotone and never lowers genus (so cycles cannot self-improve, and superseded records derive nothing new) | `testGlueMonotone` |
| B1 | band to discs gives a slice knot; band to an annulus gives genus 1 | `proofgraph_test` |
| B2 | a back edge found later improves an ancestor, with a well-founded proof (John's example) | `testCycleImprovesAncestor` |
| B3–B4 | cycles without leaves derive nothing; cycles cannot self-improve | |
| B5 | component maps matter: a permuted map gives a different, correct answer | `testComponentMapsMatter` |
| B6 | split combine and restrict; no false additivity | `testSplitEdges` |
| B7 | contradiction gates fire | `testContradictionGates` |
| B8 | non-bijective maps and boundaryless components are refused | |
| B9 | on 1,500 random graphs (>300 with directed cycles): the fixed point is independent of arrival order, equals a naive closure, is saturated, and every record rechecks | `testRandomFixedPoints` |
| C1 | a relabelled diagram is found isomorphic, and the map returned realises the isomorphism | `diagramiso_test` |
| C2 | reversing ONE component is rejected whenever it changes the writhe | |
| C3 | mirror and global reversal are found only when allowed, and flagged | |
| C4 | `splitPieces()` keeps every origin | |
| C5 | a required component map is realised exactly when it is a symmetry (the Hopf link's swap is; the Whitehead diagram's is not) | `testRequiredComponentMap` |
| D1 | `simplify()` keeps origins and every pairwise linking number over 240 random r2/r3 scrambles | `nodes_test` |
| D2 | a relabelled diagram is a "diagram" hit with a correct map | |
| D3 | scrambled diagrams of Whitehead, Borromean and Conway are one node; some need the isometry; linking numbers transport through the map | |
| D4 | `L7n1{0}` ≠ `L7n1{1}` and `L4a1{0}` ≠ `L4a1{1}`: orientation variants are different nodes | |
| D5 | one shared unknot node, with its disc | |
| E1 | a real row (`10_3`) is certified, and all 9 real witnesses assemble; curve counts match the witness file; the identity far side is the row's own node; split far sides have several pieces; **`10_3` is proved slice**, constructively, and every record rechecks | `hopedges_test` |
| E2 | a 2-component row (`L11n33{1}`) is certified; its identity witness joins each component to itself | |
| E3 | an in-process hop accounts for its surfaces as `verifyslicegenus` does (the canaries' 3_1 and L2a1{0} counts at cap 3) | `hoprunner_test` |
| E4 | every surface it keeps gives the same edge read in process as read back from its pair signature (computed later from its faces), up to the row's and the far side's symmetries | |
| E5 | its keys are distinct, and a stop request ends it | |
| E6 | a kept surface rebuilt from its faces (a certificate's route) reads exactly as in process: the same far-side curves on the same surface components; faces missing a seed triangle, or any other triangle, are refused | |
| E7 | the build digest is the same for two builds of one row, and differs for another row or another layer count | `testBuildChecksum` |
| D6 | `removeNugatoryCrossings()` removes every nugatory crossing and nothing else, keeping the link (Jones polynomial), every linking number and the component order, on connected sums through a twist (knots and links, both signs) and on kinks from Regina's own type I moves; a reduced diagram comes back untouched | `nodes_test` |
| D7 | a row with a nugatory crossing cannot be certified, and its reduced diagram's row can | `testReducedRowsCertify` |
| D8 | the checker's own removal (union-find, written separately) keeps the link and linking numbers, removes kinks and twists between summands, keeps a crossing whose two loops another component joins, turns over exactly the side between the crossing's visits (pinned: deleting the crossing without the turn gives another diagram of the same link, which the invariants cannot see), and a nugatory drawing matches its relabelled reduced node | `tools/cascade_check_test.py` |

**Mutation check (2026-09-28).** Each break of `profile.cpp` is caught:

| mutation | checks failed |
|---|---|
| drop b₁ | 9,342 |
| b₁ off by one | 33,695 |
| pairwise linking | 1 |
| no dominance | 1,236 |

## Phase 0 (2026-09-28): feasibility with today's binaries

All runs on yoga, 10 threads, 100k surfaces per hop, with the campaign's
search shape plus `--exact-far-side-names`.

| row | crossings | seed faces | search + drain | CPU | peak RSS | witnesses |
|---|---|---|---|---|---|---|
| `10_3` | 10 | 120 | 124 s | 1,353 s | 542 MB | 9 (**verified slice**) |
| `13n_65` | 13 | 156 | 303 s | 2,120 s | 654 MB | 43 (bounded 1) |
| `18nh_00000601` | 18 | 216 | 617 s | 4,487 s | 891 MB | 268 |

- **Start-up is 0.4 s.** The tables and signatures load in under 0.4 s.
- **Far sides a cascade would explore.** From `13n_65`, 8 genus-0 bands reach
  2-component links with lk = 0, including the 7-crossing `L7a4{0}`. None is
  in DG's `RibbonLinks`. From `18nh_00000601`, 14 lk-0 2-component far sides,
  one of 8 crossings.
- **`10_3` had simply never been searched** by the fixed binary. Some of the
  "unverified" slice knots are rows waiting for a campaign, not hard cases.
- **Orientation round trip.** The pipeline is `farsidediagram pd=` →
  spherogram `simplify('global')` → `PD_code()` → knotbuilder, then the
  identity witness's exact name. It returns `L11n33{1}`, `L11n113{1}` and
  `L10n83{0;1}` unchanged.

  Negative controls:
  - reversing one component of `L11n33{1}` gives `L11n33{0}`;
  - reversing a linked component of `L10n83{0;1}` gives `L10n83{0;0}`;
  - reversing its unlinked component keeps the class. The class table
    predicts this: it flips both tags, and `{1;0}` ~ `{0;1}`. So choose
    controls with the class table in hand.

**Subtleties Phase 0 found** (each needs a test or a flag in the driver):
1. spherogram's `simplify` drops split unknots into `unlinked_unknot_components`;
   count them, or a split unknot vanishes.
2. Regina's `Link::simplify()` keeps component indices and orientations: r1
   and r2 null a component in place, and r3 touches none. `simplifyExhaustive()`
   may replace the diagram, so never use it without re-deriving the component
   map.
3. spherogram PD labels start at 0, Regina's at 1 (knotbuilder takes both).
   Regina's `pd()` is a formatted string; `pdData()` is the list.
4. Above 13 crossings the complement-route Pachner fallback never names
   anything and costs about 25% of a hop's CPU (97 tries, 0 named, 1,162
   thread-seconds on `18nh_00000601`). Hops pass `--no-retriangulate-on-miss`.
5. Off-table knot far sides are named by raw isomorphism signatures, which are
   not canonical, so one knot becomes many witnesses. The cascade dedupes by
   its own canonical identity before computing pairsigs.
6. The in-search names of one knot can differ (`K13n65` and `13n_65`). Node
   identity must never rest on names.

## A node's next hop carries on from its last (since 2026-09-29)

Every in-process hop records its search frontier (`../searchfrontier.h`):
exactly how far each root got. It is written to `hop_*/frontier.txt` and kept
per node.
- **A node's next hop resumes it.** That happens when the loop raises the hop
  budget and searches every useful node again deeper, or when a hub is widened.
  The hop's surface target is its breadth, so it adds only the surfaces beyond
  the frontier. Before, each doubling searched every node's prefix again.
- **Cost.** A node searched at 2,000, then 4,000, then 8,000 surfaces now
  searches 8,000 in all instead of 14,000. Measured on `3_1`, the three hops
  reached a cumulative 2,003, 4,002 and 8,003 surfaces.
- **Nodes with nothing new to search.** A node already searched past the budget
  (a hub's wide hop, say) is skipped at that budget. One whose frontier is
  complete is never chosen again.
- **When a frontier is kept.** Only if its hop's accounting balanced,
  something was examined, its drain ran to the end and its pending file is
  fsynced (plan divergence 1, every search's rule); otherwise the node's next
  hop starts afresh. The frontier records the pending file and its length;
  a later search resumes it only once `sign` has signed that far
  (`<pending>.signed`), unless the file is in its own run directory, which
  its own sign step signs.
- **The log.** Each hop prints `[+] hop k <subject>: breadth: ...; resumed
  yes|no|none`.

Frontiers live in the process: a new cascade run starts every node afresh.
Carrying them across runs needs the node's row to be the same row. The
fingerprint would check that, but runs do not yet share a node registry.

## The driver's own time (2026-09-30)

On rung 2 of `2026-09-open2` (28 runs on `a7e6a39ac`, master witnesses on),
hops were only 46% of each run's wall time. The rest ran on the driver
thread alone while the hop threads sat idle:
- **~26%: loading master rows.** Profiled over a whole run with DWARF call
  graphs of the main thread:
  - 56% of its time was in `loadMaster`;
  - 37% was reading witnesses back: `outgoingLinkFast`'s isomorphism
    enumeration from each witness's ambient onto its subject row's
    thickening;
  - 7.5% was building those rows.
- **~16%: the lower report**, one what-if `propagateLower()` per tabulated
  node, half of it rebuilding `allPartitions()` on every call.
- **~12%: signing the kept surfaces** at the end (`--pair-sig-cache`, which
  the campaigns did not pass).

Three changes, in order:

| commit | change | effect |
|---|---|---|
| `b73a49f90` | `allPartitions(n)` built once per n and shared (`std::call_once`); the lower report's what-ifs on the run's threads, written in the same order | lower report ~4x faster |
| `9c12b2112` | master rows read back on the run's threads, then assembled serially in the same order; each read-back kept across runs in `--read-back-cache DIR` (`readbackcache.h`: one file per row, keyed by the row's PD and layers, validated by the thickening's build digest; failures kept; flock-guarded appends; a torn line is cut) | master loading ~10x faster cold, the isomorphism search gone when warm |
| `41695488e` | read-back in batches of 2 x threads rows (only a batch's thickenings held at once: 5.5 GB became 1.9 GB); node choice among equal-volume nodes made reproducible (volumes compared to 1e-6, then depth, then node id) | the same graph every run |

**A/B (halcyon, 14 threads, the campaign paused, master witnesses only,
lower report on):**

| target | base `a7e6a39ac` | `b73a49f90` | `41695488e` cold | `41695488e` warm |
|---|---|---|---|---|
| `L10n112{0;1;0;1}` (12 master loads) | 90.3 s | 73.0 s | 17.1 s | 7.8 s |
| `L10a174{0;0;1;1}` (20 master loads) | 193.8 s | 99.5 s | 29.0 s | 17.5 s |
| `L11n449{1;1;1}` (1 master load) | 1.0 s | 1.1 s | 1.1 s | 0.8 s |

- **Equivalence.** `cascade.jsonl`'s master lines and `lower_report.jsonl`
  are byte-identical to the base on `L10n112` and `L11n449`, and between
  cold and warm cache everywhere. On `L10a174`, the base chose between
  `L10a123{0;1}` and `L10a123{1;0}` (one volume) by SnapPea's rounding
  noise, so it was not reproducible. The new rule loads the same 20 rows and
  assembles the same 5,854 witnesses, with the same target bounds and lower
  report; two cold runs are identical.
- **Tests.** `readbackcache_test`: the serialisation round-trips, a stale
  digest is refused, a torn line is cut, and cached equals fresh on real
  witnesses. The other cascade tests and the cascade canaries are
  unchanged.
- **Measurement.** `tools/phase_times.py ROOT [--since ISO]` breaks a
  campaign root's wall time down by phase from each hop's `cascade.jsonl`
  timers (inside hops against the run's whole wall); the split of the rest
  came from `perf record --call-graph dwarf -t <main thread>`.
- **In campaigns:** profile keys `read_back_cache` and `pair_sig_cache`,
  passed by `cascade_worker.py` and written by `cascade_ladder.sh` from
  `$READ_BACK_CACHE` and `$PAIR_SIG_CACHE`. Each master line in
  `cascade.jsonl` now has `read_back_cache_hits`.

## Where a run's time goes, and signing on every thread (2026-09-30)

**Every phase is timed now** (`ec57fa959`):
- **Hop records** carry the driver's time since the previous hop:
  `choose_s` (of which `useful_s` and `lower_slack_s`), `master_s` and
  `kept_s`.
- **Master records** carry `wall`, `readback_s` (the pool),
  `assemble_s` (serial), `name_s` and `propagate_s`.
- **The run record.** A `{"run":...}` line at the end has the run's `wall`,
  `cpu`, `cores`, `startup_s` and `loop_s`, the hops' wall and CPU, those
  driver totals, and `store_s` (split into `store_dedupe_s` and
  `store_sign_s`), `lower_report_s` and `node_bounds_s`.

**What it showed.** Campaigns before `8d2226455` ran at 9.3 of halcyon's 14
cores on close1, and yoga's knots12 rows at about 2.6 of 12:
- **close1 (eager master loads).** Half of a row's wall time fell between
  hops: `choose()`'s graph copies and eager loads, a steady 7–8 s before
  every hop.
- **knots12.** Signing took three quarters of every row, on one thread. A
  knot row keeps surfaces from one row, and `pairSigsOf` gave each row one
  thread, nearly all of it in the ambient's `isoSigDetail()` (`pairsig.h`).
  Next came the dedupe parse of the 3 GB master, also serial.

`8d2226455` (lazy loads) fixed close1: 11.2 of 14 cores, and ~12 while the
loop runs. The two commits below fixed knots12.

| commit | change |
|---|---|
| `b89858f43` | `parallelIsoSigDetail()` (`../parallelisosig.h`): start *i* of `IsoSigClassic`'s order on thread *i* mod *T*, each keeping its least (encoding, *i*), using the engine's own `fillFrom()` and `IsoSigPrintable::encode()`. The least over threads is the serial loop's first least, so it returns the same signature *and* isomorphism. `pairSigsOf` gives each row's context `threads / rows` threads. |
| `7cb7db640` | `witnessstore::witnessIdentities()`: the dedupe's identities read in byte ranges cut at line starts, on the run's threads, with `loadWitnesses()`'s own line parser (header, empty, malformed and torn lines alike) |

**A/B.** On yoga, with the campaign paused, cold pair-signature caches,
and the benchmark rows as the campaigns run them:

| run | base `ec57fa959` | `b89858f43` | `7cb7db640` |
|---|---|---|---|
| `12a_227`, 12 threads: wall / cores | 133 s / 2.6 | 57 s / 7.5 | **51 s / 9.8** |
| of which signing / dedupe | 87.5 / 11.1 s | 22.9 / 7.8 s | 24.1 / 3.1 s |
| `12a_227`, 8 threads: wall | 123 s | 57 s | **54 s** |
| `L10a174{0;0;1;1}`, 12 threads: wall / signing | 160 s / 19.3 s | **140 s / 5.5 s** | |

- **Equivalence.**
  - `parallelisosig_test`: the same signature and isomorphism at 1–16
    threads, on row thickenings, their T and Regina examples, each also
    randomly relabelled.
  - `witnessstore_test`: `witnessIdentities()` equals `loadWitnesses()`'s
    identities at every range split.
  - `--sign-only`: each base run's own `kept.csv`, re-signed by the new
    binary against the same master, gives a store and sidecar
    **byte-identical** to the base's (four runs, 108 witnesses; `L10a174`'s
    dedupe drops 8 of 19). Two searches of one row keep different surfaces
    (a multithreaded hop's first 50k), so only the same surfaces compare.
- **Scaling.** `profile_parallelisosig` on `12a_227`'s thickening (2,304
  pentachora), on yoga:

  | threads | 1 | 2 | 4 | 6 | 8 | 12 |
  |---|---|---|---|---|---|---|
  | seconds | 76.3 | 50.6 | 34.7 | 30.1 | 27.2 | 22.8 |

  yoga is a 15 W laptop part (i7-10710U: 6 cores, 12 threads), whose clock
  falls as cores load, so everything parallel scales poorly there.
  halcyon (Ryzen 7 7800X3D, 8 cores) should do better; not yet measured.
- **Not done.** Building the target's context in the background from the
  start of the run. One thread would have to finish the whole ambient
  isoSig while the loop runs, and the store would wait on it, which is
  slower than building it on every thread at the end.

## In-process hops (`hoprunner.h`, since 2026-09-28)

A hop is one row searched on a node's diagram. By default it runs in the
cascade's own process (`HopSearcher::run()`), on exactly what a
`verifyslicegenus` row runs on, from `../rowsearch.h`:
- the same thickening: the one the hop's `HopAssembler` certified, via
  `WitnessRedrawer::rowBuild()`;
- the same gates and accounting (`gateSurface()`, `RowAccounting`);
- the same namer (`farside::DiagramNamer` with exact names);
- the campaign's shape (`HopShape`: root budget 840 from c5 on).

What changes is what happens to an accepted surface:
- It is read straight off the search: `farside::orientedOutgoingLink()` over
  the boundary the drain hands out, then `HopAssembler::addRead()`. There is
  no pair signature written and read back.
- It is deduplicated by `verifyslicegenus`'s witness identity plus its
  grouping: which row components and how many far-side curves each surface
  component carries. That grouping is what profiles read, and the identity
  alone collapses it.
- It keeps its faces (`SurfaceBoundaryInfo::captureFaces`), not a pair
  signature.
- **A certificate records it by those faces**, with its row PD, its layers, and
  the thickening's digest (`WitnessRedrawer::buildChecksum()`: every gluing,
  and which pentachoron face each triangle is). Its key is its hop key
  (`hop<k>#<i>`).
  - The checker rebuilds the row's thickening from the PD and refuses a
    different digest, so a change to the construction, or to Regina's
    skeleton numbering, is caught rather than misread.
  - It then rebuilds the surface face by face with the search's own checks
    (`farsidediagram --faces`, `WitnessRedrawer::rebuild()`).
- **No pair signature is computed for a certificate.** Its cost is almost all
  the ambient's isomorphism signature: about 25 s of one thread for a
  10-crossing row on halcyon, and different for every row, since every row
  has its own thickening. It was 61% of Pool A v1's wall time (below).
  `pairSigsOf()` (each row rebuilt once, its ambient part computed once,
  rows in parallel) remains for a witness bound for the atlas, which keys
  witnesses by their signatures.

`--hop-mode child` keeps the old route: a `verifyslicegenus` child per hop,
its witnesses read back from their pair signatures.

`hoprunner_test` pins the equivalence (E3–E5 above). An in-process read and a
pair-signature read of one surface can differ in two harmless ways:
- **The far side's curves sit elsewhere in T.** The signature route carries the
  surface back by any isomorphism that sends its incoming curve onto L × {0},
  and T has automorphisms preserving L.
- **A symmetric far side is matched by either component map.** For the Hopf
  link, both maps are true statements.

The test therefore compares the edges up to the row's and the far side's
symmetries (`findDiagramIsomorphism()` with a required component map).

**A/B (2026-09-28, yoga, 10 threads, 50k surfaces per hop, the new master
search).** The same targets, one after the other; per-hop figures from
`cascade.jsonl`:

| target | mode | hop wall | hop CPU | kept / witnesses | result |
|---|---|---|---|---|---|
| `10_155` | child | 71.4 s | 182 s | 18 | slice, 1 hop |
| `10_155` | process | 25.8 s | 184 s | 18 | slice, 1 hop |
| `11n_39` | child | 104.3 s | 290 s | 27 | slice, 1 hop |
| `11n_39` | process | 33.6 s | 219 s | 26 | slice, 1 hop |

- **Hops are about three times faster in process**, for the same result.
- **Where the child's time went.** In `10_155`, the search itself took 2.9 s
  and the drain 20 s. The pair signatures took 50.6 s of one thread, almost
  all of it the ambient context.
- **In process that cost moves to the certificate.** It is paid only when
  the goal is met, once per distinct row in the proof. Both runs above signed
  one witness serially, in ~55–65 s, before `pairSigsOf()` shared the
  ambient.
- **Whole runs were therefore about equal:** 81 vs 92 s, and 116 vs 110 s.
  In-process runs gain as proofs take more hops, or as hops go without a
  success. Certificates from faces (above) have since removed the signing.

## Pool A: the cascade against the sweep (2026-09-28, halcyon, 14 threads)

The 30 rows the atlas verified only through a chain of other rows
(`depends_on`, 2–6 rows long). The cascade ran each one alone, constructive,
with the master withheld, so its only leaves were the tables and anchors. It
used 50k surfaces per hop (up to 200k) and at most 8 expansions.

**The sweep's cost, measured with the same binary and machine.** Eight rows,
each run exactly as a campaign row (`remote_run.sh`'s flags, `hosts.conf`'s
`[campaign]` and `[halcyon]` values, 1M surfaces):

| row | wall | CPU | witnesses |
|---|---|---|---|
| `6_1` | 27 s | 330 s | 12 (cap 5 exhausted at 947k) |
| `8_8` | 40 s | 460 s | 22 |
| `9_20` | 50 s | 492 s | 33 |
| `10_27` | 58 s | 517 s | 42 |
| `L10a118{0}` | 54 s | 549 s | 8 |
| `L9a27{1}` | 55 s | 546 s | 38 |
| `L10a37{0}` | 63 s | 598 s | 29 |
| `11a_213` | 71 s | 619 s | 63 |

**v1 (in-process hops, certificates signed inline).**
- 30 of 30 proved, in 91 hops (2–5 per target).
- 1,224 s of wall time, and about 3,040 CPU-s: the searches 2,160, signing 748,
  assembly 112.
- The sweep's chains for the same 30 targets span 120 rows, ≈ 61,800 CPU-s at
  the rows' ~515 CPU-s. So the cascade used about 1/20 of the CPU, or 1/29
  counting the searches alone.
- The comparison is per target. A sweep row serves many targets at once, so
  this is not a comparison of campaign totals.
- halcyon's cores were 13% busy on average. Of the wall time:
  - the certificate's signing took 61%, one thread;
  - the hops took 28%, at about 6 of 14 cores (at 50k surfaces a hop's setup
    and drain tail are a large share);
  - assembly took 9%, one thread.
- 3 nodes were refused as hop rows: their simplified diagrams kept a nugatory
  crossing, which knotbuilder's drawer cannot certify
  (`removeNugatoryCrossings()` since).

**Each change measured alone, against its parent** (same 30 targets and
settings; wall is the driver process, start to exit, summed):

| run | change | met | hops | wall | hop CPU |
|---|---|---|---|---|---|
| v1 | in-process hops, certificates signed inline | 30/30 | 91 | 1,224 s | 2,157 s |
| v2 | certificates from faces and a build digest (no pair signatures) | 30/30 | 90 | 465 s | 2,147 s |
| v3 | nugatory crossings removed from node diagrams | 30/30 | 85 | 451 s | 2,078 s |
| v4 | new nodes named on a pool; rounds and drain tail timed | 30/30 | 86 | 452 s | 2,087 s |
| v5 | the search's progress reporters and the row watchdog woken when their work ends | 30/30 | 86 | 375 s | 2,079 s |
| v6 | the exact namer's HOMFLY index built on a pool | 30/30 | 84 | 181 s | 2,085 s |
| v7 | one set of exact-naming table caches per process, shared by the node namer and every hop | 30/30 | 84 | 167 s | 2,066 s |

- v2's certificates all check (faces route).
- v3 refuses no row. The three refused nodes now lead to shorter proofs:
  `10_27` takes 2 hops instead of 4, and `L10a136{1;0}` 2 instead of 6.
- v3's checks first failed, because the checker could not match a
  non-hyperbolic piece to a reduced node. The checker now removes nugatory
  crossings itself (D8), and every check passes, v4's included.
- v4 gained nothing. Its timers show why:
  - **"Naming" is an index build.** It costs 3.8 s at the first hop of every
    target, whether that hop adds 2 new nodes or 28, and 0.0 s for the 160
    nodes named at the 56 later hops. That 3.8 s is `ExactNamer`'s HOMFLY
    index, which computes ~34,000 table polynomials on one thread the first
    time a process needs it: 110 s of v4's 452. A pool of callers cannot
    help; the index itself must be built in parallel.
  - **The search's wall is 323 s**, over 86 hops:
    - round 1: 62 s;
    - round 2: 17.5 s;
    - the drain tail: 135 s, draining 3.78M surfaces. At 50k surfaces a
      round ends in under a second, so almost every surface is drained
      after it;
    - the rest: 108 s.

    The search, the drain tail and the rest all come out near whole
    seconds. Both progress reporters sleep in 1 s steps, and each is joined
    only when it wakes. So every hop waits up to a second after its
    enumeration and again after its drain, doing nothing. That is harmless
    in a 60 s sweep row, but in a 2–6 s hop it is most of the idle cores.
- v5 wakes the reporters (and the row watchdog, which polled every 200 ms)
  the moment their work ends.
  - The search's wall falls from 323 s to 247 s: the drain tail goes
    135 → 80 s and the rest 108 → 87 s.
  - The cores busy during the search rise from 6.5 to 8.4.
  - The rest that remains is concentrated: 26 hops carry 86 of its 87 s,
    at 3.3–3.8 s each, almost all of them a target's first hop.
  - Each hop's far-side namer is a new `ExactNamer`, which builds the HOMFLY
    index again when a drain thread first needs it. The search waits for
    that build at the join, or finds it inside the drain tail (five hops
    have tails of 4.1–4.6 s against a median of 0.65 s).
  - So v5 spends ~220 of its 375 s building one index, one thread at a
    time.
- v6 builds the index on a pool: ~0.4 s instead of 3.8 s.
  - The wall halves, 375 → 181 s.
  - The rest falls from 87 s to 3.8 s (no hop over 1 s), the drain tail
    from 80 s to 61 s, and naming from 110 s to 12.8 s.
  - Serial time is 10% of hop wall. The search keeps **13.8 of 14 cores
    busy**, and the whole run 82% of them (hop CPU against 14 × wall); the
    remainder is per-process startup and the one index per process.
- v7 shares one set of table caches (`exactnaming::TableCaches`) between the
  node namer and every hop's far-side namer. So each process builds the
  index once, on hop 0's drain pool, and the node namer reuses it.
  - Naming falls from 12.8 s to 1.3 s, and serial time to 4% of hop wall.
  - The wall is 167 s for all 30 targets. About 11 s of it is outside the
    hops: startup, ~0.4 s per process, which a single long cascade pays once.
  - The sweep's chains for the same targets cost ≈ 61,800 CPU-s, and the
    cascade uses ~2,100: about 1/29. A sweep row serves many targets at
    once, so this is per target, not per campaign.

## Knots through 8 crossings (2026-09-28, halcyon, v8 build)

After c4's merge the atlas verified 33 of the 35 knots through 8 crossings.
`8_16` and `8_18` were only verified-assisted: their bounds rested on the
literature values of `10_99` and `L10a71{0}`. The cascade ran each knot alone:
constructive, master withheld, 50k surfaces per hop (up to 1M on a revisit),
at most 7,200 CPU-s.

| knot | goal | hops | wall | CPU | the proof, down to its one leaf |
|---|---|---|---|---|---|
| `8_16` | g₄ ≤ 1 | 18 | 100 s | 1,380 s | a witness to `L8a9{0}` plus one tube; `L8a9{0}` bounds an annulus: genus 0 to a 4-crossing 3-component link, genus 0 to a 2-component unlink, discs |
| `8_18` | g₄ ≤ 1 | 14 | 73 s | 1,009 s | a witness to `L7a1{0}`, whose components bound disjoint surfaces of total genus 1: genus 0 to a 4-crossing 3-component link, then a 2-component unlink, discs |
| `8_8` | slice | 1 | 2 s | 29 s | its own row's disc |
| `8_9` | slice | 1 | 2 s | 29 s | its own row's disc |

- Every certificate is CERTIFIED by `cascade_check.py`.
- Every proof's only leaf is the unknot's disc, so with the literature lower
  bounds every knot through 8 crossings now has a constructive proof of its
  g₄: 33 from the sweep, `8_16` and `8_18` from the cascade.
- The atlas's what-if chain for `8_16` (with exact far-side names) went
  through `L8a9{0}` → `9_27`. The cascade reached `L8a9{0}` and proved it by
  its own route.
- Soundness rests on composing witnesses from different triangulations (the
  plan's §1), which the paper does not yet state. Like every bound here, it
  also inherits the draft's open check on relative Wall.
- Stored in the atlas at `results/cascade/2026-09-28_knots8/`: per knot,
  `certificate.json`, `check.txt`, `cascade.jsonl`, `driver.log` and the hop
  directories.
- In `8_16`'s run one node was refused as a hop row, for a new reason: its
  triangulated row did not redraw as its own diagram (no
  orientation-preserving isomorphism). It is not nugatory; still open.

## Try other diagrams of a stuck target (2026-09-28)

**Alternative diagrams can be incredibly productive.**
- `11a_164` (slice) resisted every hop shape on its table diagram, at about
  30,000 CPU-s in all, with a best bound of 1:
  - 50k-surface hops;
  - the master's witnesses;
  - 1M-surface hops;
  - 4M-surface hops with face cap 6.
- From the first alternative diagram SnapPy offered, also 11 crossings, the
  cascade found a two-band ribbon disc in **10 s and 2 hops**, CERTIFIED.
- A hop changes the surface only near the outgoing end of one diagram's
  thickening. So a band that is short in one diagram can be out of reach in
  another, and more search on the same diagram only re-explores the same
  neighbourhood.

How:
- `tools/alt_diagrams.py <knot>` generates diagrams (spherogram's
  `many_diagrams()`, and `backtrack()` then a light simplify).
- It keeps those proved to be the knot, by an isometry of exteriors with our
  table PD's. By Gordon–Luecke that is the knot up to mirror, which g₄ does
  not see.
- Run each through `one_target.sh` at a cheap rung first. Keep the identity
  proof beside the certificate (`identity_proofs.json`).
- The cascade's own exact namer tags the target node with its table name
  independently.

Next for the scheduler: when a node stalls, revisit it with a **new
diagram** before a bigger surface target. This is the plan's "a revisit
prefers a new diagram variant", not yet implemented.

## The independent checker (`tools/cascade_check.py`)

`cascadesearch` writes `certificate.json` when a goal is met. It holds:
- the proof records;
- each witness's hop directory, row PD, cobordism shape and maps;
- each far-side piece's match (method, component map, mirror and reversal
  flags, origins);
- each node's diagram as signed Gauss data.

The checker replays the certificate without calling any cascade code:
- **Redrawing.** It rereads each witness with `farsidediagram --gauss` (the
  search pipeline's validated drawer; `--gauss` adds each diagram's signed
  Gauss data, and the output is unchanged without it). An in-process
  witness goes through `--faces`. The row's rebuilt thickening must match
  the certificate's digest, and the surface is rebuilt with the search's
  embedding, flatness and properness checks before it is read.
  `CASCADE_ATLAS` points the checker at an atlas checkout elsewhere than
  yoga's.
- **The row-to-node map.** It finds this with its own exhaustive diagram
  isomorphism.
- **Pieces.** It re-splits far sides with its own union-find, and reproduces
  splits that appear only after simplification with its own Regina
  `simplify()` runs, each an isotopy.
- **Piece identities.** It reproduces the diagram match under exactly the
  certificate's map, trying the drawn piece and 40 of its own `simplify()`
  runs, each also with its nugatory crossings removed by its own code, either
  way up (a node's diagram is reduced; since v3, non-hyperbolic pieces reach
  such nodes and have no isometry to fall back on). Failing that, it looks
  for an isometry carrying meridians with one sign (the atlas's
  `row_certificates.meridian_signs`), with cusps built in component order.
- **Literature leaves.** It re-proves the node's identity with the atlas's
  own Python pipeline: `fsid.identify_knot` for knots,
  `row_certificates.certify_hyperbolic` (uniform sign) for links.
- **Arithmetic.** It recomputes every gluing with its own Euler-characteristic
  count.

It exits 0 only if every record checks. A replay it cannot do (split
restrictions, direct witnesses, for now) counts as a failure, never a pass.

**Lower certificates** (`lower_certificate.json`, since 2026-09-30) are
recognised by their `facts`. A fact says `lower(node, partition) ≥ value`
and carries its reason; the checker replays a witness fact's edge exactly
as an upper record's (`check_witness_edge`: redraw, row map, shape, pieces
and identities), then caps a surface of the fact's partition onto it with
its own `glue()` and requires the other end's partition to be the one
recorded (refining what that fact is stored for), the addition to match,
and `value ≤ from − addition`; a literature fact's node must be proved to
be its table entry and its value at most the table's lower end (never the
target's own entry); a linking fact's partition must really be forbidden by
the node's linking matrix. Split facts are refused as "not implemented".
The upper records a split-piece fact subtracts are checked as in an upper
certificate. Mutation check (L10n74's certificate, 2026-09-30): the
addition zeroed, the value raised by one, the other end's partition
refined, the leaf's value above the table's, the leaf renamed to another
entry, the target's partition refined, and the goal raised — all seven
refused, each with its reason.

## First runs (2026-09-28, yoga, 10 threads, 50k surfaces per hop, constructive)

| target | result | hops | search CPU | checker |
|---|---|---|---|---|
| `10_22` | slice | 1 | 1,041 s | — (certificate predates replay fields) |
| `11a_28` | slice | 1 | 1,224 s | — (certificate predates replay fields) |
| `11a_35` | slice | **4** | 3,865 s | CERTIFIED |
| `11a_36` | slice | 1 | 850 s | CERTIFIED |
| `11a_58` | slice | 1 | 1,002 s | CERTIFIED |
| `11a_87` | slice | 1 | 952 s | CERTIFIED |

| `11a_96` | slice | **5** | 5,206 s | CERTIFIED |
| `L8a9{0}` (master witnesses only) | slice | **0** | 0 s | CERTIFIED |

**Master witnesses as free edges** (`--master-witnesses`, off by default and
for the benchmark).
- **What is loaded.** One read-only scan indexes the master
  `cobordisms.csv` by subject and by recorded far side. Far-side names are
  only a hint; every witness found is redrawn and identified exactly. When a
  node is a table entry the atlas searched, the driver loads that row's
  witnesses, and those of other rows whose far side names it. The latter are
  the atlas's reverse hops, which a forward search from the node cannot find.
- **The row's own identity.** Each row is certified like a hop, and must
  intern as the node. Its map, `row_node_map`, is in the certificate, and the
  checker proves it independently.
- **Scheduling.** Free loads come before any paid hop.
- **Result.** `L8a9{0}` (the atlas's chain `L8a9{0}` → `6_1` → unlink) is
  proved with no search at all: 479 witnesses assembled, 117 s.

**`11a_35` is the first multi-hop proof.** One hop on `11a_35` alone reached
only genus 1. The proof:
1. a band takes `11a_35` to a 2-component link L (node 19);
2. a cobordism keeping L's components apart takes L to the 2-component
   unlink, so L bounds disjoint discs;
3. so `11a_35` is slice.

This is exactly the partition-aware step the solvers cannot make: they would
charge `n − 1`.

The atlas has never searched any 11-crossing row, since campaigns stopped at
10. `10_3` and `10_22` fell to one hop because the fixed binary had never
searched them.

## Lower bounds (measured 2026-09-28)

Reading a witness backwards also transports lower bounds. The paper's
`g₄(L₀) ≥ g₄(L₁) − g − n₀ + 1` charges `n₀ − 1` for the worst gluing. The
exact charge, capping with a connected surface, is n₀ − c, where c is the
number of the witness's pieces. So the charge vanishes when each piece
carries one component of L₀.

On the master, this settles **0** open link entries. 40 open `[0;1]` links
have genus-0 disconnected witnesses to far sides with lower bound 1, but every
one has c = n₀ − 1: a band merging two components, where the penalty is real.
For knots, transported bounds from τ, s or signatures are redundant, since
those invariants obey the same inequality.

The payoff is concordance to knots proven non-slice by *other* means. That is
DG §7's mechanism, which settled about 2,000 knots, 1,672 through the Conway
knot.

**Why, in general (2026-09-28).** An invariant f ≤ g₄ that changes by at
most a witness's charge across it satisfies f(target) ≥ f(source) −
charge. So a lower bound resting on such an f never beats f computed at the
target itself. That covers |σ_ω|/2, |τ|, |s|/2 and ν⁺ for knots, and
Murasugi–Tristram for links: its μ − 1 term is exactly the slack splitting
bands add. The tables already record these. So a transported bound can
close an interval only from a source whose bound is NOT of this kind, and
then only across a charge-0, genus-0 witness. From KnotInfo/LinkInfo
(database_knotinfo 2026.9.1), such sources are:
- ≈2,228 knots whose lower bound exceeds max(|σ|/2, |τ|, |s|/2, Levine–
  Tristram, ν): 2,106 non-slice with every such invariant 0, and 122 with
  g₄ = 2 above a floor of at most 1;
- 653 links above the Murasugi floor ⌈(|σ| − μ + 1)/2⌉ (e.g. the
  Whitehead link, `L5a1`).

The same argument is why σ or τ leaves on intermediate nodes could never
raise a target's lower bound.

**The lower report (`--lower-report`).** At a run's end, for each tabulated
node Y, it forgets every lower bound in a copy of the graph
(`clearLowerBounds()`), seeds Y alone, relaxes, and reads lower(target):
- `carries`: seeded at Y's literature lower bound, what Y alone gives the
  target now;
- `could_carry`: seeded at `could` = min(Y's literature upper bound, its
  best proved genus), the most Y could ever give it.

Only values Y could really have are seeded. A larger one (an earlier
version seeded a huge M to read off charges) contradicts Y's own proved
surfaces, and the split rules, which read those surfaces, then pump bounds
without limit. Each line of `lower_report.jsonl` also says whether Y is a
special source (the atlas's `data/lower_bound_sources.csv`, via
`--lower-sources`). A `carries` above the target's literature lower bound
closes the entry from below. A `could_carry` above it marks a near miss
worth a deeper search.

### Implemented (`ProofGraph::propagateLower()`, `lower()`)

**The quantity.** `lower(n, Q)` bounds below the genus of EVERY surface for
n whose partition refines Q. A bound stored for Q applies to every refinement
of Q, and `lower()` reads it that way, so it is monotone by construction.

**Seeds:**
- the literature lower bound, for every Q;
- the linking condition: a forbidden Q, and so each of its refinements, gets
  `kNoSurface`.

**Across a witness, in both directions.** Cap node X's surface (partition P
refining Q) onto the witness. This gives a surface for the other end Y with a
computable partition P′ and genus at most `genus(F) + glue-with-genus-0`. So

    lower(X, Q) ≥ min over P refining Q of [lower(Y, P′(P)) − addition(P)].

**Split links:**
- a partition that never mixes pieces bounds the whole by the sum of the
  pieces' bounds; a mixing one gets nothing (`K ⊔ −K`);
- a piece gets `lower(whole, Q ∪ P_B) − g_B` for every PROVED surface for
  the other pieces. This is the paper's split-union (ii), generalised.

**Gates.** The iteration is capped at 200 passes. Sound facts cannot pump a
bound past the truth, so non-convergence is reported as inconsistent input.
Afterwards, every proved surface must satisfy every lower bound; a failure is
a contradiction.

**Driver.** It runs `propagateLower()` after every hop. Its what-if
optimism is capped at each partition's proved lower bound.

**Tests** (`proofgraph_test`):
- L1: concordance transport;
- L2: the paper's reverse inequality, band and merging band;
- L3: no penalty for disjoint annuli;
- L4: split rules;
- L5: 600 random "worlds" (the upper closure taken as the truth, true minima
  as literature) never produce a bound above an achieved surface, and never a
  contradiction;
- L6: the lower report's cleared what-if;
- L7: transport is monotone under refinement (lem:transport-monotone,
  below), on 400 random worlds with random literature bounds, every edge,
  both directions, every partition (>10,000 pairs, >100 of them strict);
- L8: `lowerIf()`, the what-if that keeps the graph's bounds;
- L9: every raised bound remembers its reason, down to a literature leaf.

**Mutation check.** Each break is caught:

| mutation | checks failed |
|---|---|
| drop the genus subtraction | 1,339 |
| split sum off by one | 3 |
| piece ignores the other genus | 1 |

**Lemma (transport is monotone under refinement; `lem:transport-monotone`,
2026-09-30).** For a witness e and its end X, let `t_e(P)` be what e
transports to X for surfaces with partition exactly P: the other end's
bound at the partition P′(P) the cap induces, minus what the cap adds. If P
refines Q then `t_e(P) ≥ t_e(Q)`. *Proof.* Cap a surface of partition P
(genus 0) onto e. Refining the cap's partition adds vertices to the gluing
graph and no edges, so its components can only split and the addition
`Σ_K b₁(K) = E − V + #components` can only fall; the induced partition P′
of the other end only refines; and `lower(Y, ·)` is monotone under
refinement (a bound stored for a partition applies to its refinements) and
the linking condition is upward-closed (a refinement of a forbidden
partition is forbidden). ∎ So the bound for surfaces refining Q is `t_e(Q)`
itself: `lowerAcross()` evaluates the one partition instead of minimising
over its refinements (Bell² → Bell per edge), and the AND over the target's
partitions that a lower proof needs is met by any single edge. L7 checks
the equality on random graphs, so a future shape or map for which it fails
is caught. (The plan first read the per-edge minimum as a gap that a
"cover rule" over edges would fix; the lemma says the rule is vacuous.)

## Lower-bound mode (2026-09-30)

Until now the cascade only *sought* upper bounds: `useful(n)` kept a node
only if its best conceivable profile could complete the target's upper
goal, and a run stopped when nothing was useful for that goal. Lower bounds
were propagated after every hop and reported at the end, never aimed at.
The one open entry closed from below, `L10n74{1;0}`, came from the
target's own ninth hop. John's idea: keep a node if its best case could
give a new upper **or** lower bound, let the search rise in complexity (a
larger far side can carry a higher lower bound back), and bound it.
`--goal-lower G` does that (plan: `~/.claude/plans/immutable-cooking-puffin.md`).

**What the data said** (open2's 238 `lower_report.jsonl`): of 37,492 nodes,
only 2 ever carried a bound above a target's literature floor (L10n74's);
203 near-misses over 122 targets all pass through 26 open `[0;1]` links of
8–11 crossings, hubs among them (`L8a21{0;0;1}` under 31 targets,
`L8n7{0;0;1}` 22, `L10a174{0;0;1;1}` 18, `L10n112{0;0;1;0}` 16). Lower
bounds flow hub → target at charge 0, so one non-slice hub closes dozens
of entries. Special sources (`data/lower_bound_sources.csv`): 3,013;
`lit_lo` 1 for 2,827, ≥ 2 for 186, ≥ 3 for 33, 4 for 7.

**The gate** (`usefulLower()`). A node n is kept for the lower goal if,
given the best lower bounds it could ever have, it would carry the goal to
the target over the edges found so far: `ProofGraph::lowerIf()` seeds n in
a copy that keeps every bound the graph has, relaxes, and reads
`lower(target, goal)`. The seeds are admissible: any transported bound is
a literature seed minus non-negative charges, so at most the largest
special source's (`L_max`, 4 in the tables; a Lipschitz-floor source
cannot beat the target's own floor across the exact charge, see "Lower
bounds"); a proved surface refining a partition caps it; the literature
upper bound, a connected surface, caps the coarsest. A seed above a proved
surface is refused, and nothing is seeded above what is known. The gate is
what makes the search terminate: a chain from a `[0;1]` target can spend
at most `L_max − 1` charge, and only charge-0 hops are unbounded in
number — bounded by `--lower-max-crossings` (16: the chain must come back
to a table entry), the expansion and CPU limits and the budget ladder.
What-ifs run for every candidate the upper gate rejects, in one parallel
batch per `choose()`, cached per graph version.

**Ordering.** Among nodes of one crossing count, the upper gate's nodes
first, then the lower gate's by slack (what they would carry beyond the
goal: the charge still affordable, which decides how many sources are in
reach); then volume, depth, id as before.

**Reasons.** Every raised lower bound remembers the fact that raised it
(`ProofGraph::lowerWhy()`): a literature or linking leaf, a witness with
the other end's partition and the cap's addition, or a split rule with the
records it subtracted. Raises are strict and no rule increases what it
reads, so reasons never cycle: a lower bound's proof is a tree. A met
lower goal prints it (`[+] LOWER GOAL MET`, then the chain) and writes
`lower_certificate.json` for the checker.

**Master rows on their recorded diagrams.** Found on the way: `loadMaster`
read every stored row on the table's PD, but a cascade hop's row is its
node's *simplified* diagram, so cascade-recorded witnesses under a table
name never read back (50 of L10n74's 52 failed). Now the store's
`.rows.csv` sidecar gives each witness its row PD, and such a row is
interned from its recorded diagram with that match's component map. Lines
stored before the sidecar existed (open2's first 6,559) still have no row;
`tools/farside/exact/cascade_row_candidates.py` recovers those from hop
logs.

**Reproduction (2026-09-30, yoga).** `L10n74{1;0}` with `--goal-genus 0
--goal-lower 1`, open2's store as master witnesses (its 20 target rows
given their recorded diagram from the packed `hop_kept_lines.csv`) and no
search: 24 witnesses assembled, `LOWER GOAL MET` in 2 s, the chain
`L10n74{0;1} {0,1,2} ≥ 1` across `cb72899e102e` (genus 0, 2 pieces) from
`L11n48{0} {0,1} ≥ 2` with the cap adding 1, leaf `literature L11n48{0}
2`; the certificate CERTIFIED by `cascade_check.py` (identity by isometry
with uniform meridian signs) and all seven mutations refused.

**Every node's bounds** go to `<work>/node_bounds.jsonl` at a run's end:
identity (table name or untabulated), diagram, linking matrix, every proved
profile entry with its record and whether its proof is constructive, and
every lower bound with its reason's kind. Composites and links beyond the
tables get bounds here that no table records (John, 2026-09-30: those are
results to keep); the atlas folds runs into a table beyond the tables.

**Master loads are lazy, and witnesses are filtered by genus before they
are read (2026-09-30, John).** Campaign close1's rows took a median
1,234 s against open2's ~250 s, a third of it hops: every table node the
graph met had its stored rows read back at once (35–45 loads, ~2,500–3,000
pair-signature read-backs per row, graphs of 1,000–2,200 nodes). Two
hereditary checks, in the enumerator's spirit: a node's rows are loaded
only when it is the one about to be expanded (the target's own rows still
first), since a proof runs through a node only when the search picks it
(`--master-loads eager` restores the old behaviour); and a stored witness
is skipped unread when its genus exceeds what any proof could spend — the
upper goal, or the largest special source's bound minus the lower goal
(`glue()` never lowers a genus; a lower chain's charge is at least the
witness's genus) — counted as `skipped_genus` per load. A/B on
`L8a21{0;0;1}` (yoga, 8 threads, 600 CPU-s, open2's store as master
witnesses): eager 1,368 s wall, 67 loads, 1,165 nodes; lazy **278 s**, 8
loads, 479 nodes; both 10 hops, ~650 s search CPU, best 1, `cpu-budget`;
0 witnesses skipped by genus (all have genus ≤ 3).

**Not built, deliberately.** A charge ladder (charge-0 hops before any
charged one, raised when nothing is useful): the slack ordering already
puts charge-0 nodes first within a budget level; measure the A/B first.
Split facts in the checker. The A*-style ordering for the upper side.

**Still to measure:** the 26 hub nodes as targets with `--goal-lower 1`
at a small rung, against the upper-only run on the same targets.
