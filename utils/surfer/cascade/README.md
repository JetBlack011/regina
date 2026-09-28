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
| `tools/cascade_check.py` | the independent checker |
| `tests/` | one test per module; `data/` holds real witnesses from the Phase 0 hops and small tables |

The first four files have no Regina dependency.

Still to come:
- the SnapPy oracle (DG `PlausibleKnots`, `RibbonLinks`, HFK) for leaf facts
  beyond the tables;
- switching hops early, depth first, when a far side could beat the demand
  (`HopSearcher::run()`'s `stop` is the hook).

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
- A certificate signs only the witnesses its proof uses (`pairSigsOf()`):
  each row rebuilt once, its ambient part computed once (`PairSigContext`),
  rows in parallel. Nearly all of a signature's cost is the ambient's, about
  50 s for a 10-crossing row. In the certificate such a witness is named by
  its signature's key, like any other, and its hop key becomes `provenance`.

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
  success.

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
  Gauss data, and the output is unchanged without it).
- **The row-to-node map.** It finds this with its own exhaustive diagram
  isomorphism.
- **Pieces.** It re-splits far sides with its own union-find, and reproduces
  splits that appear only after simplification with its own Regina
  `simplify()` runs, each an isotopy.
- **Piece identities.** It reproduces the diagram match under exactly the
  certificate's map. Failing that, it looks for an isometry carrying
  meridians with one sign (the atlas's `row_certificates.meridian_signs`),
  with cusps built in component order.
- **Literature leaves.** It re-proves the node's identity with the atlas's
  own Python pipeline: `fsid.identify_knot` for knots,
  `row_certificates.certify_hyperbolic` (uniform sign) for links.
- **Arithmetic.** It recomputes every gluing with its own Euler-characteristic
  count.

It exits 0 only if every record checks. A replay it cannot do (split
restrictions, direct witnesses, for now) counts as a failure, never a pass.

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
  contradiction.

**Mutation check.** Each break is caught:

| mutation | checks failed |
|---|---|
| drop the genus subtraction | 1,339 |
| split sum off by one | 3 |
| piece ignores the other genus | 1 |
