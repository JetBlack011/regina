# surfer

Enumerating PL-embedded surfaces in triangulated 4-manifolds, and using them
to find cobordisms in S³ × I that bound the smooth 4-genus of knots and links.
The method is written up in `pl_enumeration_draft/`; the results live in the
separate `cobordism-atlas` repository.

This file is a map of the code: what each part does, what it guarantees, and
where it is tested. Each header's own `\file` comment has the details.

## Tools

| tool | what |
|---|---|
| `verifyslicegenus` | **the search driver.** One row (knot or link) at a time: thicken the row's triangulation into S³ × [0,2], seed it with the collar L × [0,2], enumerate surfaces from it, name each far side, record witnesses, and solve for bounds (`cobordismgraph`). `--solve-only` re-derives statuses from a witness file. |
| `surfer` | the general enumerator, with diagnostics (cache statistics, root coverage) that `verifyslicegenus` does not print |
| `peripheral_slopes` | drills far sides (from witnesses' pair signatures) or table links (from PD codes) keeping signed meridians, for the far-side pipeline's SnapPy half |
| `farsidediagram` | witnesses' outgoing links as oriented diagrams, via `knotbuilder/diagramdrawer` (see `knotbuilder/README.md`) |
| `knotbuilder/triangulateknot` | a PD code → the triangulation's and its complement's isomorphism signatures |
| `tools/bench_search.sh` | rough, repeatable search benchmarks, one row per run (see "Performance" below) |
| `tools/compare_surface_sets.sh` | whether two builds accept and describe exactly the same surfaces: exhaustive runs with `--surface-log`, sorted and compared |

`bogocheck.cpp` and `fillmanifold.cpp` are old scratch programs that no target
builds.

## The code

**Building the ambient manifold**

| file | what |
|---|---|
| `knotbuilder/` | a PD code → a triangulation T of S³ containing the link; the crossing block's geometry; drawing curves of T as diagrams. Has its own README. |
| `simplicialprism.{h,cpp}`, `cobordismbuilder.{h,cpp}` | T × I by prism layers (the paper's thickening); `CobordismBuilder` keeps the top prisms, so the outgoing boundary is known to be T × {2} simplex by simplex |
| `collar.{h,cpp}` | the collar L × [0, 2]: the seed every searched surface contains |

**Enumerating surfaces**

| file | what |
|---|---|
| `enumerate_cis.{h,cpp}` | connected induced subgraphs of the triangle-adjacency graph (optionally seeded, budgeted, with iterative deepening). Since 2026-09-28, relying on the prunes being anti-monotonic: a child that fails is set aside for the rest of its parent's subtree, and a seed neighbour that fails with the seed alone for good; the depth cap is the enumerator's own (`setMaxSize()`), so a node at the cap returns at once; children go back into the candidate list in place, so the list at any node is a function of the path; and a budgeted pass that runs out records a `Position`, from which the next pass carries on instead of retracing. Checked against brute force, and resumed passes against one unbudgeted pass in order, by `tests/enumerator_test` |
| `skeleton.{h,cpp}` | adjacency graphs over faces |
| `embeddedsubmanifold.{h,cpp}` | `EmbeddedSubmanifold` / `KnottedSurface`: a growing set of triangles with incremental embeddedness, local flatness, orientability; `orientedBoundaryLinks()` orients the boundary per surface component |
| `vertexlinks.{h,cpp}`, `rollbackunionfind.{h,cpp}` | the incremental local checks and their memoisation |
| `linkingnumber.{h,cpp}` | the linking number of two closed petals' traces in Lk(v), by cochains on Lk(v) itself (since 2026-09-28): push B off into the dual cells, solve δx = PD(B*) over GF(2⁶¹−1), read x(A). Checks its own answer (δβ = 0, δx = β everywhere) and declines rather than guess; `KnottedSurface::addFace()` then falls back to drilling (`linkcomplement`). `--audit-linking` runs both routes on every miss |
| `embeddingsearch.{h,cpp}` | the parallel search over roots (parallel across roots, never within one) |
| `surfacesearch.{h,cpp}` | the search as `verifyslicegenus` uses it: enumeration plus the **drain**, which describes and names each accepted surface's boundary |

**Naming boundaries**

| file | what |
|---|---|
| `linkcomplement.{h,cpp}` | drilling curves (Weeks' edge pinching) to build their complement |
| `farsidenaming.{h,cpp}` | **how the search names far sides** (since 2026-09-26): draw the curves (`knotbuilder/diagramdrawer`), simplify with Regina's `Link::simplify()`, and name by exact diagram signature against the knot and link tables. Anything a diagram cannot name falls back to `identifycomplement`. `--no-diagram-naming` restores the old route. |
| `identifycomplement.{h,cpp}` | naming a complement: unknot and handlebody recognition, the local census (`census.sqlite`), Regina's census, and a Pachner search (knots only). Now the fallback for far sides, and still how everything else is named. `BoundarySignatureCache` memoises by the curve's edge set up to automorphism, in front of either route |
| `exactnaming/` | **exact far-side names** (since 2026-09-27): an oriented far-side diagram named, with a proof, as a table entry or a split union / sum of table entries. Orientation variant and relative chirality are pinned, and every name is flagged as an identity or a description. Hyperbolic pieces are named by the SnapPea kernel's isometry test (an isometry of complements carrying meridians to meridians; the kernel's own isometry sources, which Regina ships unbuilt, are compiled into this library), the rest by Reidemeister searches. See `exactnaming/README.md` |
| `farsideredraw.{h,cpp}` | a stored witness's outgoing link recovered from its pair signature, exactly as the search saw it (the row's own thickening, an isomorphism pinned by L × {0}); shared by `farsidediagram` and `farsidename` |
| `farsidename.cpp` | the tool: pair signatures → exact names, one line per witness. Its output becomes the atlas's `data/far_side_exact.csv`, which the solvers read with `--far-side-exact` |
| `tableclasses.cpp` | the table's link classes: table names that are one oriented link up to mirror and global reversal (a shared diagram, or an isometry carrying meridians to meridians with uniform orientation signs). Writes the atlas's `data/table_link_classes.csv`, which the solvers read with `--link-classes`; fails on a class whose literature values differ |
| `linknames.h` | census names → KnotInfo / LinkInfo names |
| `peripheral.{h,cpp}` | drilling while keeping meridians, signed by the curve's orientation |

**Witnesses and the solver**

| file | what |
|---|---|
| `pairsig.{h,cpp}` | isomorphism signatures of (ambient, surface) pairs: a witness's canonical form, and `fromPairSig()` to decode one |
| `witnesskey.{h,cpp}` | `sha1(pairsig)[:12]`, the key per-witness far-side resolutions are stored under |
| `cobordismgraph.{h,cpp}` | the solver (`propagate()`: the cobordism inequalities over witnesses and literature bounds), plus the per-surface decisions: `splitBoundary()` (which boundary is the search side), `buildRowOrientation()` / `classifyRowOrientation()` (does the surface run with the row's orientation), `witnessIdentity()` (dedupe) |
| `farsidecurves.{h,cpp}` | a surface's outgoing curves carried onto knotbuilder's T through the thickening's own top prisms, oriented against the row per surface component: the far side as an **oriented** link |
| `csvwriter.{h,cpp}` | sharded CSV output |

`tools/frontier.py` in the atlas is a deliberately independent
reimplementation of the solver; `frontier.py --check` must say AGREE after any
change to either.

## Invariants the search keeps (since 2026-09-26)

An audit on 2026-09-26 found the search silently discarding qualifying
surfaces, often every surface of a row. The fixes, and the checks that now
guard against that whole class of fault:

- **The search side is identified by geometry, never by name.** In a seeded
  search it is the component holding the seed. The row setup asserts once that
  the seed's incoming edges form one closed curve per component, and that no
  searchable non-seed triangle touches the incoming boundary. Before this, one
  mismatched census name (`#17` against `#6`) discarded a whole row.
- **The row map pins L × {0}.** `buildRowOrientation()` takes the isomorphism
  that carries L onto the seed's own edges; T's automorphisms made an
  arbitrary one unreliable.
- **Orientation is judged per surface component.** Components can be oriented
  independently and tube together either way. Comparing signs across
  components rejected most disconnected surfaces of links.
- **Dedupe** includes whether a surface is resolved. Far sides that bear no
  bound are still deduped by their complement's name: keying them by edge set
  made almost every surface its own witness.
- **Only knots get the Pachner search** (`--retriangulate-links` restores it
  for links); a link's name bears nothing in the search.
- **Census writes land.** A lookup used to leave its read statement open, so
  every insert waited out the 5 s busy timeout and failed. Each
  `identification:` line now reports `census writes N ok/M failed`.
- **Far sides are named from their diagrams, never by drilling them first**
  (`farsidenaming.h`). Every name is proved:
  - no crossings left after simplifying means `Unknot` or the n-component unlink;
  - an identical table diagram (up to mirror and reversal) means that table name;
  - a knot diagram whose name the complement route proved earlier in the row gets that name again.

  The complement route remains only for:
  - a knot no table or earlier fallback knows;
  - a link whose linking numbers are all zero (it might be an unlink that
    the simplifier cannot clear, and an unlink bears a bound);
  - a degenerate drawing;
  - a drawing the drawer refuses as **not planar** (`knotbuilder::NonPlanar`,
    since 2026-09-27). That is always a drawer defect; each is counted in the
    row's `diagram naming:` line (`non-planar drawings`), and any at all prints
    a WARNING. Before the gate, the vertical-corner-edge defect
    (`knotbuilder/README.md`) let 13 c4 far sides be named `diagram:` after a
    non-planar drawing.

  A link that matches no table diagram and is provably not an unlink (a
  nonzero linking number, or a Jones polynomial other than the unlink's) is
  named `diagram:<signature>`. Like any link name, that bears no bound, but it
  tells different links apart exactly.
  - **Fallback cost:** every fallback answer is remembered against its
    diagram's signature, so there is at most one per distinct diagram per
    row, typically 0–10.
  - **Speed:** naming, which was about 94% of the drain, is now a few seconds
    of thread time per 1M-surface row. Rows run in roughly the time the
    enumeration takes: 246–326 s on halcyon for rows that took 707–1,079 s.
- **Every surface is accounted for.** Each row prints
  `accounting: accepted A, described D, recorded R, duplicate …, other-orientation …, search-side-elsewhere …, impossible I, drain complete|skipped, ok|WARNING`.
  An imbalance or an impossible state halts, after the witnesses are written.
  No exhaustion is claimed unless the row balanced. The campaign dispatcher
  fails any row without a clean accounting line.
- **An append-only witness store.** `cobordisms.csv` is appended and fsynced,
  never rewritten. Loading keeps byte offsets instead of pair signatures.
  `--solve-only` never writes. `--rewrite-witnesses` is the one explicit
  rewrite (the 12→13 column migration), and it verifies every line. A torn last
  line is ignored on load and truncated before the next append.

## Performance: measuring it, and its history (since 2026-09-28)

**The `search profile:` line.** Every row prints where its search time went,
after its `identification:` line:

```
search profile: prototype 0.0s (unknot misses 30 in 0.0s, linking misses 0 in 0.0s);
  rounds 4.1s 7.1s 0.0s; drain tail 965785 surfaces in 42.0s; nodes 121818346,
  attempts 131797610, evaluated 131797610, charged 131787694, replayed 7655;
  petal misses: unknot 4110 in 17.2s, linking 762 in 0.0s (cochains 762, fallbacks 0)
```

- **prototype:** seed commit and root filtering, single-threaded, before any
  worker starts.
- **rounds:** each IDDFS round's wall time.
- **drain tail:** what was left queued when the search ended.
- **The walk:**
  - nodes: visits;
  - attempts: `tryAdd()` calls;
  - evaluated: those that reach the embedding checks (all of them, now that
    the depth cap is the enumerator's);
  - charged: those counted against root budgets;
  - replayed: the uncharged re-adds that resume a root where its previous
    pass stopped.
- **Petal misses:** time is thread time. `cochains`/`fallbacks` say how the
  linking numbers were computed (`linkingnumber.h`), and `--audit-linking`
  adds a comparison with drilling on every miss.

In a budgeted round each root's walk is a function of the root alone, so the
walk counts are deterministic. Identical counts between two builds mean they
searched identically.

Before 1e (resumed passes), each pass re-walked its root from the start, and
"replayed" counted the charged attempts that retraced the previous pass.
Charged − replayed then equalled an unbudgeted run's attempts exactly,
checked on `6_1` and `L6a3{1}`.

**The benchmarks** (`tools/bench_search.sh`), run on halcyon with its
`hosts.conf` profile (14 threads):

- **B1:** exhaustive to 4 added faces on `10_141`, `L10a14{0}` and
  `L10a127{1;1}`, with an empty witness file.
  - Fixed work, so every correct build accepts the same surfaces.
  - `bench_search.sh ab A B` alternates the two builds (ABAB) and flags any
    difference in accepted counts or in the walk.
- **B2:** one production row (1M surfaces, caps 4 then 5, a copy of the real
  witness store).

Binaries are kept in halcyon's `~/bench-bins`, and results in
`~/bench-runs/results.tsv`.

**Equivalence.** `tools/compare_surface_sets.sh A B` runs the canary rows
exhaustively at cap 3 with both builds. The sorted surface logs (type,
triangle count and pair signature of every described surface) must be
identical.

**History.** B1 wall time per row, medians of two, in seconds. The "walk" and
"sets" columns say whether the walk counts and the surface sets matched the
build before.

| build | change | 10_141 | L10a14{0} | L10a127{1;1} | total | walk | sets |
|---|---|---|---|---|---|---|---|
| `f183dffbe` | master | 105.7 | 110.9 | 113.5 | 330.2 | — | — |
| 00-counters | the `search profile:` line | 105.6 | 110.4 | 112.8 | 328.9 (0.996×) | — | — |
| 1a | a child's neighbours join C only once it passes `tryAdd` | 74.2 | 78.9 | 81.5 | 234.7 (0.712×) | same | same |
| 3 | the seed's own report (and with it the pair-signature context) moves to the drain thread | 51.0 | 56.5 | 59.5 | 167.1 (0.715×) | same | same |
| 4 | drain embeddings keep the seed; queued entries hold only the added faces | 47.5 | 52.0 | 54.5 | 154.1 (0.933×) | same | same |
| 2 | petal linking numbers by cochains (`linkingnumber.h`) | 43.0 | 45.9 | 48.4 | 137.3 (0.903×) | same | same |
| 1b | a child that fails is set aside for its parent's subtree; a seed neighbour that fails, for good | 39.7 | 43.3 | 45.9 | 128.8 (0.939×) | differs | same |
| 1c | the depth cap is the enumerator's (`setMaxSize()`) | 38.9 | 42.0 | 45.2 | 126.1 (0.971×) | differs | same |
| 1d | candidates go back into the list in place | 38.8 | 41.8 | 44.9 | 125.5 (1.002×) | differs | same |
| 1e | budget passes resume instead of replaying | 38.3 | 41.6 | 44.5 | 124.4 (0.991×) | differs | same |
| final | the surface target stops the search exactly | 38.5 | 41.5 | 44.7 | 124.7 (0.984×) | — | same |

Each ratio is against the row above, measured in the same A/B run, so
successive rows' absolute times differ by run-to-run noise (about ±1 s).
Items 1b–1e change the order of the walk on purpose ("differs"); what must
not change is what is accepted, and it did not, including for 1e at a
budget of 2,000 attempts per pass (many suspensions and resumptions per
root), and in `enumerator_test`, where resumed passes visit exactly what one
unbudgeted pass does, in order. A budgeted and an unbudgeted run of 1e on
`10_141` visit the same 29,314,297 nodes.

B1 stops improving at 1b because it is no longer bound by the search: on
`10_141` the search round is 69 s at 00, 22 s at 2 and 4 s at 1e. What is
left is the pair-signature context (~27 s of one thread, which the first new
witness waits for) and the drain.

Item 3 helps B1, whose witness file starts empty, and any row whose seed is
new. With the real witness store the seed's witness already exists, no pair
signature is needed for it, and the prototype pass is already ~0 s.

**B2** (`10_141`, 1M surfaces, a copy of the real witness store):

| build | root budget | wall | rounds (cap 4, cap 5) | drain tail | peak RSS | attempts | linking misses (thread time) |
|---|---|---|---|---|---|---|---|
| 00-counters | 50,000 | 299.7 s | 70.7 s, 211.7 s | 409,505 in 12 s | 1,453 MB | 14.96 B | 1,057 (1,352 s) |
| 4 | 50,000 | 190.7 s | 35.4 s, 137.4 s | 571,127 in 13 s | 1,238 MB | 14.59 B | 976 (1,209 s) |
| 2 | 50,000 | 109.6 s | 22.6 s, 62.6 s | 854,384 in 19 s | 630 MB | 15.01 B | 1,022 (< 0.1 s) |
| 1e | 50,000 | 59.7 s | 4.0 s, 8.2 s | 1,043,032 in 25 s | 666 MB | 0.147 B | 956 (< 0.1 s) |
| 1e | 840 | 61.7 s | 3.9 s, 8.2 s | 1,087,392 in 25 s | 669 MB | 0.151 B | 915 (< 0.1 s) |
| final | 840 | 59.1 s | 4.1 s, 7.1 s | 965,785 in 42 s | 642 MB | 0.132 B | 762 (< 0.1 s) |

- The linking misses were about 32% of a production row's thread time;
  item 2 takes them to under a second.
- Before 1e, 55% of all attempts were replays of earlier budget passes
  (8.23 B of 14.96 B at 00).
- **The root budget.** 1e changes what a unit of it buys, so
  `root_budget_start` was recalibrated to keep c4's shape: c4's rows made a
  median 10.2 round-2 passes per root (quartiles 9.9–10.3, over all 308).
  1e at 50,000 made 4.3, and 840 gives 10.2.
  - Measure with `tools/round_passes.py <row>.err`.
  - `hosts.conf` `[campaign]` sets 840 from c5 on.
- **The surface target.** At 1e's speed, the watchdog's once-a-second stop
  overshot the target by 9% (1,090,787 accepted for 1M). The search now
  stops itself exactly (1,000,003).
- The last row's drain tail includes waiting for the pair-signature context,
  which it built in the tail rather than during the search.

**Where a production row's time goes now** (`perf` flat profile of the final
B2 row):
- allocator churn (`malloc`, `free` and friends): ~37%;
- the pair-signature context's isomorphism signature (`IsoSigData::fillFrom`,
  `encode`): ~5%, all on one thread, at the first new witness;
- the drain's per-surface work (`boundaryLinks`, `orientedBoundaryLinks`,
  `boundaryEdgeSurfaceComponent`, `singularVertices_`, and Regina building
  the boundary curves' triangulations): ~15%;
- the search itself (`addFace`, `extendFiltered`): ~9%.

At 00, the enumerator's own code alone was 60% of a B1 row. The next
targets are the pair-signature context (cache it per row, or build it off
the critical path) and the drain's allocations.

## Tests

`ctest` from the build's `utils/surfer` directory; `embeddedsubmanifold_test`
needs an idle machine. Besides each component's own unit tests:

| test | guards |
|---|---|
| `tests/enumerator_test` | the enumerator returns exactly the brute-force set of connected induced subgraphs (with seeds, budgets, depth caps, anti-monotonic filters); budgeted passes that resume one another visit exactly what one unbudgeted pass does, in the same order, down to a ration of one attempt |
| `tests/predicate_order_test` | `KnottedSurface`'s prunes (P_1, flatness, transversality) agree with a from-scratch reference on every set of triangles in small closed triangulations, whatever the order faces are added, including one-vertex ones where a triangle has several corners at a vertex |
| `tests/rowmap_test` | the row map lands exactly on L × {0}, no searchable face touches it, the bare collar classifies as matching (optionally over a whole table) |
| `tests/name_independence_test.sh` | perturbing every name (identified or drawn: it runs with diagram naming and requires that it was used) changes nothing the search accepts or records |
| `tests/interrupted_outcome_test.sh` | a search stopped by SIGINT is recorded as `interrupted`, never `exhausted` |
| `tests/census_test` | census lookups, and that an insert after a hit lands |
| `tests/cobordismgraph_test` | the solver's rules, `splitBoundary`, per-component orientation, witness identity; exact far sides (bound by their own variant, receive a bound only when exact), `m`/`r` knot marks and the slice test with them, sum pieces, and `--sum-rules` |
| `exactnaming/tests/exactnaming_test` | exact names: every table entry names itself, orientation variants pinned (L7n1{0} vs {1}), granny vs square, splits and sums, the search path forced, the table-side search and the isometry step each on its own |
| `exactnaming/tests/snappeaisometry_test` | the SnapPea isometry test: far-side drawings of 8_16, L8a1, 10_151, L11n281 match their table entries; a HOMFLY twin (10_156) and three links sharing L11n353's complement do not |
| `exactnaming/tests/isometry_validation` (not in ctest) | the isometry step on the whole table: kernel census vs KnotInfo/LinkInfo, every entry under every orientation transform scrambled and renamed, all HOMFLY- and volume-twin pairs, threads. See `exactnaming/README.md` |
| `knotbuilder/tests/diagramdrawer_test` | the drawer (see `knotbuilder/README.md`) |
| `tests/farsidenaming_test` | on real thickenings: the collar's far side is named as the row itself (a table knot, a table link's base) straight from its diagram, with no fallback; a small curve is `Unknot`; empty or missing signature tables are refused |
| `tests/linkingnumber_test` | the cochain linking number against diagrams (knotbuilder draws table links; Regina reads |lk| off the PD code: 0, 1, 2, 3), against the drilling route on random disjoint cycles, and on linked components rerouted across triangles (which cannot change lk). With a link table as argument, sweeps every 2-component row |

The atlas adds campaign-level checks: `tools/orchestrate/canaries.sh` (exact
surface accounting on fixed rows, run by `verify.sh gate` and at staging) and
`tools/orchestrate/audit_rows.py` (after every merge).
