# surfer: PL surfaces in triangulated 4-manifolds

A library that enumerates surfaces made of the triangles of a triangulated
4-manifold: connected sets of triangles that form an embedded, locally flat
surface (or one whose only self-intersections can be resolved), optionally
containing a fixed seed and meeting the boundary as asked. It knows nothing of
knots, tables or slice genus; `cobound` builds the slice-genus search on it.
The method is the paper's §3 (`sec:enum`) and §4 (`sec:main`).

It links Regina, OpenSSL's libcrypto (SHA-1) and `linknaming_complement`
(drilling, and the unknot and unlink certificates its local checks rest on).
Only its tests and the driver's `--pd` mode also use `diagramtriangulation`.

## The code

| file | what |
|---|---|
| `submanifold/skeleton` | `Skeleton`: the adjacency graph of a triangulation's k-faces, one node per face and an edge for each pair sharing a facet. What the search walks. |
| `submanifold/submanifold` | `EmbeddedSubmanifold<dim,subdim>`: a set of faces kept up to date as faces are added and removed (embedded at every codimension, closed, proper, orientable). `KnottedSurface` (dimension 4, triangles) adds the petal checks at vertex links, resolvable self-intersections, smoothness at the boundary, and the boundary curves, oriented per surface component. |
| `submanifold/vertexlinks` | `PetalCache`: petal identities and the memoised unknot and linking answers, shared by every thread's surface. `SelfIntersectionCensus`, for measurement only. |
| `submanifold/linkingnumber` | the linking number of two petals' traces in Lk(v), from the cochains of Lk(v) itself (paper `lem:linking-cochains`). |
| `submanifold/rollbackunionfind` | `RollbackUnionFind`: union by size, no path compression, exact rollback of recent merges. |
| `enumeration/inducedsubgraphs` | `ConnectedInducedSubgraphEnumerator`: every connected induced subgraph once (after Alokshiya, Salem and Abed), seeded, pruned by anti-monotonic predicates, under a depth cap and per-root budgets. |
| `enumeration/submanifoldsearch` | `EmbeddingSearch<dim,subdim>`: the parallel search over roots, its rounds, budgets, surface target and frontier. |
| `enumeration/surfacesearch` | `SurfaceSearch`: `EmbeddingSearch<4,2>` over `KnottedSurface`s, plus the drain, which describes each accepted surface and names its boundary curves through a `BoundaryNamer`. `SurfaceSearchLimits`: the queue's and caches' sizes. |
| `enumeration/namecache` | `BoundarySignatureCache`: the drain's names of a boundary component's marked edge sets, up to that component's automorphisms. |
| `enumeration/searchstrategy` | `SearchFrontier`: how far a search got, exactly, and its file format. |
| `pairsig/pairsig` | pair signatures, isomorphism signatures of (ambient, marked subcomplex) pairs (paper `app:pairsig`), and `fromPairSig()`. `PairSigContext`: a signature's ambient part, computed once per ambient and optionally kept on disk. |
| `pairsig/parallelisosig.h` | `parallelIsoSigDetail()`: Regina's `isoSigDetail()` on several threads, giving the same signature and the same isomorphism. |
| `pairsig/sha1` | `sha1Hex()`, over OpenSSL's `SHA1()`. |
| `report/csvwriter` | `CsvWriter`, a sharded CSV writer for many threads; `csvField()` and its inverse `parseCsvLine()`. |
| `report/atomicwrite` | `atomicWrite()`: the one way a file is replaced whole. |
| `report/progress` | `RollingReport`: a block of status lines redrawn in place. Each driver formats its own block. |
| `surfer.cpp` | the proof-of-concept driver (below). |

## Enumeration

A candidate is a set of triangles connected in the triangle adjacency graph.
`ConnectedInducedSubgraphEnumerator` visits each connected induced subgraph
exactly once, by reverse search from each root in turn (paper
`sec:connected-subgraphs`). Three things shape the walk:

- **A seed.** The seed's triangles are contracted to one vertex, so every
  candidate contains the seed (paper `def:seed-contraction`,
  `alg:seeded`). A boundary component can be protected: no triangle with an
  edge on it is searchable besides the seed's, so every surface meets that
  component exactly as the seed does.
- **Predicates.** A predicate is checked as each triangle is added. Every
  predicate the search prunes on is anti-monotonic (paper
  `def:anti-monotonic`): a set that fails has no extension that passes. So
  pruning loses nothing, and a candidate that fails is set aside for the rest
  of its parent's subtree (paper `lem:pruned-candidate`).
- **Limits.** A depth cap (`setMaxSize()`) bounds the triangles added to the
  seed. A per-root budget rations each root's attempts per pass: a root that
  spends its ration is suspended at a recorded `Position`, and the next pass,
  with the ration multiplied by the growth factor, carries on from there. So a
  stopped search has covered a definite share of every root rather than a
  prefix of the root list, and the total work stays within growth/(growth − 1)
  of the last pass.

`EmbeddingSearch::search()` runs the roots on worker threads, never splitting
a root between threads. With iterative deepening it first runs capped rounds
(round i capped at `iddfs_start + (i − 1) · iddfs_step` added triangles), each
reporting only what the previous round could not reach, then a final round up
to the hard face cap. It stops on its own once every root is enumerated (a
round in which that happens gives `deepestExhaustedCap`, the search's one
exhaustive claim), at the surface target (counted by the worker that finds
each surface, so it overshoots by at most the few in flight), or on
`requestStop()`. A search handles SIGINT itself unless the caller turns that
off (`setSigintHandling(false)`), as `cobound` does.

## What a surface must satisfy

`KnottedSurface` keeps every check incremental, so adding or removing a
triangle costs work near that triangle only.

- **Pruned while adding** (each anti-monotonic):
  - P_1: every edge lies in at most two triangles (paper `def:p1`,
    `lem:p1-anti-monotonic`);
  - at an interior vertex v, a closed petal's trace in Lk(v) is unknotted
    (local flatness) and any two closed petals' traces have linking number 0
    (transversality) (paper `def:petal`, `def:trace`, `def:transverse`,
    `lem:smoothness-anti-monotonic`). Every corner of a triangle is registered
    before any petal is checked, so the answer does not depend on the order of
    addition;
  - orientability, when the search asks for it (`orientableOnly`; paper
    `lem:orientability-anti-monotonic`).
- **Checked on the finished candidate** (`isAcceptable()`):
  - embedded at every codimension (`isEmbedded()`), or, when resolvable
    surfaces are accepted, every self-intersection is at an interior vertex
    whose trace T_v(S) is certified an unlink (`isResolvable()`; paper
    `def:unlinked-self-intersection`, `thm:resolution`). The certificate is
    `complement::certifiesUnlink()`, which never says yes falsely; boundary
    vertices never qualify;
  - smooth at the boundary (`isSmoothAtBoundary()`): every petal at a boundary
    vertex, closed through the cone point of the coned-off ball, is unknotted
    (paper `def:petal-knotted`). This is a filter, not a prune: at a boundary
    vertex a knotted open petal can still close up into an unknotted one, so
    the condition is not anti-monotonic (paper `sec:pruning`);
  - the boundary condition: `all`, `closed`, `proper` (every boundary edge on
    the ambient boundary) or `connected` (proper, and at most one boundary
    curve per ambient boundary component).

Under `proper`, these make every accepted surface properly embedded and
locally flat, or perturbable near its unlinked self-intersections into such a
surface of the same topology (paper `cor:local-flat-smooth`,
`cor:search-certifies`).

Petal answers are memoised in one `PetalCache` shared by every thread, keyed by
the petal's corners. A petal linking number is computed from the cochains of
Lk(v): β = PD(B*) is solved as β = δx over GF(2⁶¹ − 1) and x(A) read off.
Every answer is checked (β a cocycle, δx = β everywhere); a failed check
returns nothing, and the caller falls back to drilling one trace out of Lk(v)
and reading the other's homology class.

## The drain

A `SurfaceSearch` that tracks boundaries (`proper` or `connected`) queues each
accepted surface's added triangles, and drain threads rebuild and describe
them while the search runs (`SurfaceBoundaryInfo`): its topology
(orientability, genus, punctures, and the genus after discarding closed
components and tubing the rest together), and per ambient boundary component
its curves' edges and names. Within the callback the caller may also take the
boundary curves oriented by the surface (each surface component oriented
independently, with a map saying which component each edge bounds), the
surface's triangles, from which its pair signature can be computed later, and
the pair signature itself.

Boundary curves are named by the caller's `BoundaryNamer`, which the search
refuses to start without. `UnlinkBoundaryNamer` names census-free
(`Unknot`, `<n>-component unlink`, else the complement's isomorphism
signature); `cobound` supplies the slice-genus search's own. A namer may say
`Unknot` or `<n>-component unlink` only with a proof. The search reads no
other name (it counts curves), so any other name answers to its reader.
`cobound solve` gates bounds by the curves' count (`cobound/README.md`, "The
atlas solver (`cobound solve`)"): one curve bears a bound whatever it is called,
so a name of one curve must denote one knot up to mirror, which a name taken
from its complement does (Gordon–Luecke; cobound's own `|` descriptions are
the exception described there); 2 or more curves bear nothing unless they are
an unlink or proved by another input, so their name need not determine
them. Names are memoised per boundary component by
`BoundarySignatureCache`, keyed by the marked edge set up to the component's
automorphisms, so a namer sees each distinct curve set once.

When more than `pendingSurfaceCap` surfaces wait, the search pauses and every
worker drains until the queue is empty. After the search, the rest of the
queue is drained on all threads; `rebuildFailures()` counts surfaces the drain
could not rebuild, which is always 0 unless something is broken.

## Frontiers

A search's traversal depends only on its graph, its roots and its schedule:
each root's walk depends on that root alone, and a budget pass carries on from
its `Position`. So where a search stopped is fully described by the round it
was in and, for each root of that round, whether it finished, how many passes
it had (its level) and where its last pass stopped. With
`setRecordFrontier(true)` the search records exactly that as a
`SearchFrontier`, with cumulative counts over every run that continued it. A
stop suspends each in-flight root where it stands
(`BudgetedPredicate::suspendOnStop()`), so a frontier never overclaims.

Handed back with `setResumeFrontier()`, a frontier resumes the search if and
only if its fingerprint is the search's; otherwise the search starts afresh
and `resumeRefusal()` says why.

- **What resuming guarantees.** A chain of runs, each resuming the last one's
  frontier, reports exactly what one uninterrupted run reports when taken to
  the end: nothing twice (the seed included), nothing missed. In the same order
  only without budgets: a root that spends its ration goes to the back of the
  root queue, while a resumed round rebuilds the queue in index order, so a
  budgeted search stopped at a surface target and resumed reaches its next
  target through other surfaces than one run would.
- **The surface target is the search's breadth.** A resumed search counts its
  frontier's surfaces, so it adds only new ones, and one already that broad
  returns at once.
- **The fingerprint** is `sha1Hex()` of a text naming everything that fixes
  the traversal and what it accepts
  (`EmbeddingSearch::frontierFingerprint_()`): `SearchFrontier::kTraversalVersion`;
  the dimension; every gluing of the triangulation; the search graph and each
  vertex's faces (as a set: a seed's faces arrive in no fixed order); whether it
  is seeded; the sorted roots; the round schedule and hard cap; the budget
  schedule; the boundary condition by name (`all`, `closed`, `proper`,
  `connected`); orientability pruning; and the subclass's context
  (`SurfaceSearch`: `resolve_unlinked 0` or `1`). These literals are frozen:
  changing any of them makes every stored frontier foreign, and nothing
  resumes again. Bump `kTraversalVersion` with any change to which candidates
  the enumeration visits or in what order; `test_traversal_pinned` fails until
  it is bumped and re-pinned.
- **The file** (`save()`, `load()`) starts `surfer-search-frontier 1`, or `2`
  when it records where the caller keeps the search's finds and that file's
  durable length (`pending`). Format 2 stores that path relative to the
  frontier's own directory, so a moved tree keeps its resume. A frontier is
  written atomically but not fsynced; `read()` refuses an empty, cut or
  malformed file by name, and a caller treats it as absent. `summary()` is the
  one-line `round r/R cap c, roots d done + s part-walked of N, max level L;
  runs n, satisfying S, found F, attempts A[, complete]`.

Skipping a frontier's prefix is sound only if the run that recorded it kept
what it found. The library records the caller's pending file and its length;
the rule for when that suffices is the caller's (`cobound/README.md`, "A
search's life").

## Pair signatures

`pairSig()` is an isomorphism invariant of an (ambient triangulation, marked
subcomplex) pair: Regina's isomorphism signature of the ambient, with the
marked faces appended, minimised over every automorphism of the canonical
ambient, so isomorphic pairs give byte-identical strings even when the ambient
has automorphisms moving the subcomplex. `fromPairSig()` decodes one. Stored
cobordisms are recorded by them, so their encoding is frozen.

Almost all of a signature's cost is the ambient's own isomorphism signature
(32.8 s of 32.9 s on a 1,728-pentachoron thickening), and every surface of a
search shares the ambient. So `PairSigContext` computes that part once,
`LazyPairSigContext` defers it to the first signature (most searches find
nothing to sign), and `PairSigContext::cached()` keeps it on disk, keyed by
`ambientKey()` (a sha1 of the ambient's gluings) and verified on loading. The
ambient part is computed on several threads by `parallelIsoSigDetail()`
(start i of `isoSigDetail()`'s order on thread i mod T, the least over threads
by encoding and start index), which returns the serial answer byte for byte.
`sig()` equals `pairSig()` byte for byte.

## The driver

`surfer` (built as `build/utils/surfer/surfer/surfer`) runs one search and
prints its progress and the surfaces' homeomorphism types:

```
surfer [-a|-c|-p|--connected] [options] <isosig>      a 4-manifold by its isomorphism signature
surfer [-a|-c|-p|--connected] [options] --pd <pd>     T x [0,k] built from a diagram
```

With `--pd`, T is `diagramtriangulation`'s triangulation of the diagram,
thickened `--thicken-layers` times (default 1) and seeded with the link's
collar through the first `--collar-layers` (default 1; 0 for no seed).
`surfer --help` lists the options (threads, rounds, face cap, budgets,
`--orientable-only`, `--resolve-unlinked`, `-o` for one CSV line per surface,
cache limits). It names boundary curves with `UnlinkBoundaryNamer`.

## Performance

**Measuring.** Each search of a `cobound run` prints a `search profile:` line:

```
[+] 10_141: search profile: prototype 0.0s (unknot misses 30 in 0.0s, linking misses 0
  in 0.0s); rounds 4.1s 7.1s 0.0s; drain tail 965785 surfaces in 42.0s; nodes 121818346,
  attempts 131797610, evaluated 131797610, charged 131787694, replayed 7655; petal
  misses: unknot 4110 in 17.2s, linking 762 in 0.0s (cochains 762, fallbacks 0)
```

(one line, wrapped here; with `audit_linking` it ends with the audit's
counts.)

- prototype: committing the seed and filtering roots, on one thread, before
  any worker starts;
- rounds: each round's wall time; drain tail: what was left to describe when
  the search ended, and how long it took;
- the walk: nodes visited, `tryAdd()` attempts, those evaluated by the
  embedding checks, those charged to root budgets, and the uncharged re-adds
  that take a suspended root back to its `Position`;
- petal misses: thread time spent on unknot and linking answers the cache did
  not have, and how the linking numbers were computed.

Within a budgeted round each root's walk is a function of the root alone, so
the walk counts are deterministic: two builds with equal counts searched
identically. `cobound/tests/bench_search.sh` runs the two benchmarks:

- **B1:** exhaustive to 4 added triangles on `10_141`, `L10a14{0}` and
  `L10a127{1;1}`, with an empty database. Fixed work, so every correct build
  accepts the same surfaces; `bench_search.sh ab A B` alternates two builds
  and flags any difference in accepted counts or in the walk.
- **B2:** one production search (`10_141`, 1M surfaces, caps 4 then 5).

`cobound/tests/compare_surface_sets.sh A B` is the equivalence check: both
builds search the canary targets exhaustively at cap 3 with a surface log, and
the sorted logs (type, triangle count and pair signature of every described
surface) must be identical.

**The root budget, 840.** A pass's ration counts attempts, and since a pass
resumes where the last one stopped rather than retracing it, a unit buys more
than it did when passes replayed. 840 keeps the breadth each root gets at a
production shape: production searches made a median 10.2 round-2 passes per
root when passes replayed (quartiles 9.9–10.3, over 308 searches); with
resuming passes, 50,000 gives 4.3 and 840 gives 10.2. A frontier records each
root's passes; the `breadth:` line prints the largest.

**mimalloc.** Every executable links mimalloc statically and whole, so it
replaces glibc's `malloc` for the engine library too; results cannot depend
on the allocator. Allocator churn was about a third of a production search's
CPU. B1 on halcyon (14 threads) took 98.7 s against glibc's 122.2 s (×0.81),
and the drain on yoga 0.11 ms per surface against 0.22. The archive must
define `malloc` as machine code: Arch's packaged one is LTO bitcode, which a
linker skips with only a warning, so CMake checks it and refuses one that does
not. Build mimalloc from source into `~/.local/lib/libmimalloc.a`, or pass
`-DSURFER_MIMALLOC=<path>`. Without it CMake warns and links glibc's.

**The pair-signature context.** Signing is dominated by the ambient's
isomorphism signature (above), which is why it is computed once per
thickening, on every thread, and kept across runs. On `12a_227`'s thickening
(2,304 pentachora, a 6-core laptop part), `profile_parallelisosig` takes
76.3 s on 1 thread, 34.7 s on 4 and 22.8 s on 12; a context loaded from the
cache takes a fraction of a second.

## Tests

`ctest` in `build/utils/surfer/surfer`. `submanifold_test` is the slow one,
deliberately exhaustive: about 20 s on an idle 8-core machine against a 300 s
timeout, which CPU contention alone can exceed, so run it on an idle machine.

| test | pins |
|---|---|
| `inducedsubgraphs_test` | on many random small graphs, the enumerator visits exactly the connected induced subgraphs a powerset sweep finds, each once: unseeded and seeded, under an anti-monotonic predicate, under the depth cap, and under budgets, whose passes together visit exactly what one unbudgeted pass does, in order, down to a ration of one attempt |
| `submanifold_test` | `isEmbedded()` against a from-scratch injectivity check over every connected subset of several triangulations' 2-skeletons, in the search's own order; boundary links' homology (H_1 = Z^n) on coned triangulations; petal linking numbers (Hopf, torus links, Whitehead); petal checks on cones over knots and links (a trefoil cone rejected, unlink cones resolvable, the Whitehead link's not, a Hopf cone still pruned); boundary-vertex self-intersections never resolvable; the boundary filter rejects a knotted boundary petal; along random walks, every accepted face can be removed and re-added (the prunes are hereditary) |
| `predicate_order_test` | P_1, flatness and transversality agree with a from-scratch reference on every set of triangles up to a size cap in small closed triangulations, whatever the order of addition, including one-vertex triangulations, where a triangle has several corners at one vertex |
| `submanifoldsearch_test` | boundary conditions; orientability pruning; the protected boundary component and the seed's exemption; the depth cap; iterative deepening and budgets against one pass; `deepestExhaustedCap`; queue back-pressure drops and doubles nothing; frontiers: the file round trip, the relative pending path, damaged files refused, stopped and resumed chains equal to one run (seeded and not, budgeted and not, with and without rounds, 1 and 2 threads), a foreign frontier refused, and the traversal pinned to `kTraversalVersion` |
| `vertexlinks_test` | `PetalCache`: petal identity independent of corner order, distinct petals distinct, misses then hits, the linking cache symmetric, the clear threshold, stale ids a clean miss |
| `rollbackunionfind_test` | rollback undoes unions exactly, and later unions behave |
| `linkingnumber_test` | the cochain linking number against diagrams (Regina's linking number of the PD code: 0 to 3) and against drilling on random disjoint cycles (a sample in ctest; `--full` compares everywhere); a link table as argument sweeps every 2-component row |
| `pairsig_test` | round trips; invariance under relabelling, including ambients with automorphisms moving the marked set; malformed input refused; a context equals the free function, decodes, and survives concurrent first use; its disk cache stores and loads, rebuilds a damaged file, and keeps a relabelled ambient apart |
| `parallelisosig_test` | `parallelIsoSigDetail()` equals `isoSigDetail()` (signature and isomorphism) at every thread count, on thickenings of table links, randomly relabelled, and on Regina's examples |
| `namecache_test` | `BoundarySignatureCache`: one entry per edge set up to the boundary component's automorphisms, checked against an independent action on edges |
| `sha1_test` | `sha1Hex()` against FIPS 180-1's vectors, length boundaries, and Python's `hashlib` values |
| `atomicwrite_test` | a replaced file has the new contents and nothing beside it; on any failure the old file is left as it was |
| `progress_test` | the block is redrawn in place, rewinding exactly the previous block's lines |

Not in ctest: `profile_driver`, `profile_boundarylinks`, `profile_pairsig`
and `profile_parallelisosig` (timing drivers), `replay_pairsigs` (stored pair
signatures re-encoded through a context must come back byte for byte) and
`review_repeated_vertices` (triangles with a repeated vertex in every table
row's thickening).
