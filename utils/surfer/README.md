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
| `enumerate_cis.{h,cpp}` | connected induced subgraphs of the triangle-adjacency graph (optionally seeded, budgeted, with iterative deepening). Checked against brute force by `tests/enumerator_test` |
| `skeleton.{h,cpp}` | adjacency graphs over faces |
| `embeddedsubmanifold.{h,cpp}` | `EmbeddedSubmanifold` / `KnottedSurface`: a growing set of triangles with incremental embeddedness, local flatness, orientability; `orientedBoundaryLinks()` orients the boundary per surface component |
| `vertexlinks.{h,cpp}`, `rollbackunionfind.{h,cpp}` | the incremental local checks and their memoisation |
| `embeddingsearch.{h,cpp}` | the parallel search over roots (parallel across roots, never within one) |
| `surfacesearch.{h,cpp}` | the search as `verifyslicegenus` uses it: enumeration plus the **drain**, which describes and names each accepted surface's boundary |

**Naming boundaries**

| file | what |
|---|---|
| `linkcomplement.{h,cpp}` | drilling curves (Weeks' edge pinching) to build their complement |
| `farsidenaming.{h,cpp}` | **how the search names far sides** (since 2026-09-26): draw the curves (`knotbuilder/diagramdrawer`), simplify with Regina's `Link::simplify()`, and name by exact diagram signature against the knot and link tables. Anything a diagram cannot name falls back to `identifycomplement`. `--no-diagram-naming` restores the old route. |
| `identifycomplement.{h,cpp}` | naming a complement: unknot and handlebody recognition, the local census (`census.sqlite`), Regina's census, and a Pachner search (knots only). Now the fallback for far sides, and still how everything else is named. `BoundarySignatureCache` memoises by the curve's edge set up to automorphism, in front of either route |
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
  - a degenerate drawing.

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

## Tests

`ctest` from the build's `utils/surfer` directory; `embeddedsubmanifold_test`
needs an idle machine. Besides each component's own unit tests:

| test | guards |
|---|---|
| `tests/enumerator_test` | the enumerator returns exactly the brute-force set of connected induced subgraphs (with seeds, budgets, iterative deepening, hereditary filters) |
| `tests/predicate_order_test` | `KnottedSurface`'s prunes (P_1, flatness, transversality) agree with a from-scratch reference on every set of triangles in small closed triangulations, whatever the order faces are added, including one-vertex ones where a triangle has several corners at a vertex |
| `tests/rowmap_test` | the row map lands exactly on L × {0}, no searchable face touches it, the bare collar classifies as matching (optionally over a whole table) |
| `tests/name_independence_test.sh` | perturbing every name (identified or drawn: it runs with diagram naming and requires that it was used) changes nothing the search accepts or records |
| `tests/interrupted_outcome_test.sh` | a search stopped by SIGINT is recorded as `interrupted`, never `exhausted` |
| `tests/census_test` | census lookups, and that an insert after a hit lands |
| `tests/cobordismgraph_test` | the solver's rules, `splitBoundary`, per-component orientation, witness identity |
| `knotbuilder/tests/diagramdrawer_test` | the drawer (see `knotbuilder/README.md`) |
| `tests/farsidenaming_test` | on real thickenings: the collar's far side is named as the row itself (a table knot, a table link's base) straight from its diagram, with no fallback; a small curve is `Unknot`; empty or missing signature tables are refused |

The atlas adds campaign-level checks: `tools/orchestrate/canaries.sh` (exact
surface accounting on fixed rows, run by `verify.sh gate` and at staging) and
`tools/orchestrate/audit_rows.py` (after every merge).
