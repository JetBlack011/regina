# exactnaming — exact names for far sides

Names an **oriented** far-side diagram — the far side drawn from knotbuilder's
T × {2} with the orientation its surface induces (`../farsidecurves.h`,
`../knotbuilder/diagramdrawer.h`) — as a table entry, or a split union or
connected sum of table entries. The name comes with a proof, with the
orientation variant and relative chirality pinned, and with no candidate sets.
Paper: §5.3 "Naming by drawing" (`sec:naming-drawing`), and Lemmas
`lem:split-union`, `lem:sum-along-components`, `cor:sum-pieces` for how the
names bound.

| file | what |
|---|---|
| `exacttables.{h,cpp}` | the knot and link tables as oriented diagrams, each entry in four versions (as written, mirrored, reversed, both), keyed by Regina's signature with neither reflection nor reversal allowed; a second, unoriented index for search results; canonical names for entries that are one link up to mirror and global reversal (the Hopf link's two variants) |
| `gaussdiagram.{h,cpp}` | an oriented diagram as signed crossings plus each component's crossing sequence (`regina::Link::fromData()`'s input). `splitPieces()` and `visibleSum()` cut it without losing orientation or chirality |
| `exactnamer.{h,cpp}` | decompose, identify each piece, compose the name |
| `snappeaisometry.{h,cpp}` | `KernelLink`: a piece's complement in Regina's SnapPea kernel, built from the diagram so its peripheral curves are the diagram's meridians; `sameLinkAs()` is the kernel's isometry test (what SnapPy's `is_isometric_to()` calls) |
| `snappeakernel_isometry{,_cusped}.cpp` | compile the kernel's own isometry sources, which Regina ships unbuilt in `engine/snappea/kernel/unused/`, against the kernel Regina does build; the engine is untouched |
| `tests/exactnaming_test.cpp`, `tests/snappeaisometry_test.cpp` | see below |

The command-line front end is `../farsidename.cpp`. It takes pair signatures
and redraws each witness exactly as the search saw it (`../farsideredraw.h`).

## How a far side is named

1. **Decompose.** `simplify()` and `simplifyExhaustive()` never reflect or
   reverse. After them, a disconnected diagram gives a split union. Two arcs
   on the same two faces whose removal disconnects the crossings give a
   connected sum along the component through both, with each side closed
   along the circle. Repeat on the pieces. Components in no piece are split
   unknots.
2. **Identify each piece**, in this order.
   - *Exactly:* after `simplify()` tries, the piece's diagram **is** one version
     of a table entry, which pins its oriented variant and its mirror and
     reversal relative to the drawing.
   - *By isometry* (hyperbolic pieces): the table entries whose HOMFLY
     polynomial, in either mirror, is the piece's and whose table diagram is
     no bigger are a shortlist (a hint only). The piece's complement is built
     by the SnapPea kernel from the diagram, so its peripheral curves are the
     diagram's meridians and longitudes. The kernel's `compute_isometries()`
     canonises both complements and tries every combinatorial isomorphism of
     the canonical triangulations; a match is proved when one takes every
     meridian to a meridian (`isometry_extends_to_link()`). That isomorphism
     is exact: floating point only steers the canonisation, so it can cause a
     miss but never a false match. Filling along the meridians extends it to
     S³, so the piece is the table link up to mirror and orientations; then
     the variant is pinned as below. Milliseconds, whatever diagram was drawn.
   - *By search* (what the isometry leaves, chiefly non-hyperbolic pieces):
     a bounded `rewrite()` at height 2 reaches a diagram of a table base; then,
     from the table's side, the flype orbit of an alternating candidate or a
     `rewrite()` outward from a candidate's table diagram (heights 1–4); then
     height 3 forward for ≤ 12 crossings. `rewrite()` may reflect and reverse,
     so this proves the base only.
   - *Pinning the variant* (isometry or search): every oriented variant V of
     the base, and mV, is ruled out if its HOMFLY polynomial or its linking
     numbers differ from the piece's, computed on the drawn orientation. Both
     are invariants up to global reversal, and the table lists every
     orientation up to global reversal, so the survivors contain the piece.
   - *Otherwise untabulated:* `diagram:<sig>`, the signature of the oriented
     diagram, with reversal disallowed and the lesser of it and its global
     reverse. A signature determines its link.
3. **Compose**, in the atlas grammar:

| shape | name | exact? |
|---|---|---|
| unknot, unlink | `Unknot`, `<n>-component unlink` | yes |
| one piece | its table name (`L7n1{1}`, `8_17`), or `A\|B` when several variants survive (so far always with one literature g₄) | yes when one name survives |
| one piece plus split unknots | `L7n1{1} u Unknot` | yes (lem:split-union (iv)) |
| composite knot | `3_1#m3_1`, `3_1#mr9_32`: each summand marked `m`/`r` where its symmetry type makes the mark matter, spelled with the fewest marks over the four global mirror/reversal choices | yes when every mark that matters is pinned |
| split union | `A u B` (Unknot last) | no |
| knots summed into a link | `3_1 #_? L2a1{0}` | no |
| links summed along components | `#{L2a1{0}[?] # L4a1{1}[?]}` | no |

**Exact** means an identity: one link up to mirror and global reversal, which
may receive a bound from the subject. Anything else is a proved
**description**: each piece is exactly the named table link, but the gluing is
not all recorded (`?` marks a component index not computed). That is enough
for the sum and split rules, but not for an identity. In `#{...}` a piece summed
at two sites is written at both. The solvers therefore count every occurrence
as a piece, which only weakens the bounds (`cor:sum-pieces`).

## Tests

`tests/exactnaming_test` (in ctest):
- every table entry, drawn as itself, mirrored or reversed, names itself;
- L7n1{0} with one component reversed is L7n1{1} (g₄ 0, not 2);
- granny (`3_1#3_1`) and square (`3_1#m3_1`) are told apart;
- a sum and its global reverse share a name (`3_1#8_17`);
- splits, split unknots, unlinks, `K #_? L` and `#{…}`;
- the search path, forced (a kinked diagram with simplification off), must
  still pin L7n1's variant;
- far sides the forward search could not reach (another minimal diagram of
  10_151, a flype of 11a_18's, 8_16 drawn two crossings above minimal) are
  named by the table-side search and by isometry, each with the other steps
  off; 8_16 must not be taken for 10_156, which shares its HOMFLY polynomial.

`tests/snappeaisometry_test` (in ctest): `sameLinkAs()` finds 8_16, L8a1,
10_151 and L11n281 from far-side drawings; it rejects 8_16's drawing against
10_156, and three far sides whose complement is L11n353's while the links are
not (an isometry exists, none carries meridians to meridians); a torus knot is
not hyperbolic and matches nothing.

## Validation of the isometry step (2026-09-27)

`tests/isometry_validation` (not in ctest: 15 minutes on 10 threads) runs on
the full tables, 17,153 entries.

- **Census.** The kernel's hyperbolicity and volume for every entry agree
  with KnotInfo (geometric type, volume) and LinkInfo (volume):
  - 16,981 hyperbolic, volumes equal to 10⁻⁶;
  - 172 not hyperbolic (torus knots, satellites, non-hyperbolic links).
- **Self-identification.** Every entry is taken under every transform: as
  written, mirrored, reversed, and each component c ≥ 1 reversed, 58,679
  cases in all. Each is scrambled by seeded random Reidemeister moves (two
  R2, up to eight R3, one R1), which keep the orientation. Every scramble is
  a diagram other than the table's, even up to reflection and reversal.
  - Named by the isometry step alone, on the raw scrambled piece: **all
    57,841 hyperbolic cases give exactly the reference name**, the full
    namer's exact diagram match on the unscrambled diagram. The 838
    non-hyperbolic cases stay untabulated, left to the searches.
  - For comparison, the search-only namer gives 17,673 identical, 4,439
    proved sets containing the reference, 36,567 untabulated (the scramble
    is beyond its reach) and 0 different.
- **Negatives.** All 3,166 pairs of distinct table links sharing a HOMFLY
  polynomial (either mirror) and all 2,784 pairs with equal volume, 4,786 in
  all. None is called the same link, including the 279 pairs whose
  complements are isometric (e.g. `L10n1` with `L7a2`, `L9n7`, `L11n16`).
- **Positives.** All 2,648 pairs of variants of one hyperbolic link are the
  same link. 526 of them are the same *oriented* link up to mirror and
  global reversal, and share one canonical name; their literature g₄ agree
  in every case.
- **Threads.** The isometry namings redone on one thread with a fresh namer:
  identical, 0 of 58,679 different.
- **Memory.** `snappeaisometry_test` and `exactnaming_test` under valgrind:
  no errors, every heap block freed.
- **Edge cases** (`snappeaisometry_test`): virtual, zero-crossing, split,
  composite and 50-crossing diagrams; orientation readings (a link to
  itself, to its global reverse, to its mirror, and a variant pair that
  differs as oriented links, `L7a6{0}` and `L7a6{1}`).

What it found, and what changed:
- **Canonical classes were per diagram.** Variants that are one oriented
  link but have different table diagrams (e.g. `L8a4{0}` and `L8a4{1}`: the
  526 pairs above) had two names. So one far side could be named either
  way, depending on where the randomised `simplify()` landed; four of c4's
  witnesses flipped between runs. Now such variants are merged by an
  isometry with uniform orientation signs (`canonicalName()`, per base on
  first use). A hyperbolic piece's variant, and a knot's mirror and
  reversal, are pinned by the isometry itself, with the HOMFLY polynomial
  and linking numbers kept as a cross-check (a disagreement throws).
- **One miss in 58,679.** Canonisation occasionally retriangulates at
  random. A hyperbolic piece that matches nothing is now retried on
  randomised triangulations.
- **A table duplicate.** LinkInfo's `L11n432` and `L11n437` are one
  non-hyperbolic link, mirror images: a height-2 Reidemeister search
  connects their diagrams, and their g₄ pair up (0, 1, 2). A far side of
  it may be named either way; no bound is wrong. The isometry cannot merge
  them (not hyperbolic).
- **The solvers need the classes too.** Merged names moved bounds from one
  row to another (`L8a4{1}`'s to `L8a4{0}`'s). `../tableclasses` writes the
  classes, 933 table names in all (408 by diagram, 525 by isometry), to the
  atlas's `data/table_link_classes.csv`. Both solvers read it opt-in
  (`--link-classes`) and treat a class as one node. With it, the what-if
  solve improves 342 rows over the baseline and worsens none (verified
  584 → 863), and the two solvers agree.

On real data:
- **Master, 249,527 witnesses:** 0 failures, so no variant the isometry
  pinned was ruled out by the HOMFLY polynomial or linking numbers. Names
  with unpinned variants fell from 2,124 to 88, all non-hyperbolic. Every
  other change from the previous run is a class merge (13,673) or a
  pinning (2,033).
- **c4, 4,416 witnesses:** every change against the pre-isometry run is a
  class merge (62) or a pinning (120); all 571 knot far sides agree with
  the search's own names.
- **In the search:** a c4-shaped row with `--exact-far-side-names`, the
  kernel called from ten drain threads: the accounting balances, with 0
  naming failures. Its 6 multi-curve far sides are named exactly as
  `farsidename` names them.
- **Tests:** all 23 `surfer` ctest targets and the canaries pass.

## Measured

On the full master (2026-09-27, 249,527 witnesses, 403 s with 10 processes).
The driver, and the numbers in detail, are in cobordism-atlas
`tools/farside/exact/`.
- **Named:** every witness, with 0 failures and 0 non-planar drawings; 4
  stay untabulated (one 12-crossing link, and three links with L11n353's
  complement that are not L11n353). Pieces: 322,544 by diagram, 23,101 by
  isometry, 1,765 by search (all non-hyperbolic). Since the isometry pins
  variants, 88 names (0.04%) have unpinned variants, all non-hyperbolic.
- **Against the SnapPy pipeline:** 211,963 of 211,984 agree; none is left
  untabulated where SnapPy named it; all 20 of its unidentified far sides
  are named; 1 disagreement, where the pipeline is wrong.
- **History of the gap.** The first run (forward search only) left 250 far
  sides untabulated that SnapPy named: right link, wrong diagram. A search
  from the table's side (flype orbits, rewrites outward from the table
  diagram) closed all 250 in 806 s; the isometry step names the same 250,
  identically, in 44 s, and is what runs first now.
- **Solvers:** a what-if solve with these names and the table's link classes
  improves 342 rows and worsens none (verified 584 → 863).

On c4's 4,416 witnesses (2026-09-27):
- 66% of names are exact identities (41% one oriented table link, 16% one
  table knot, 3.4% untabulated), and 31% are sums or splits with every piece
  pinned;
- 2.8% have unpinned variants, all with a single literature g₄;
- against an independent Python/Regina redraw, 1,671 of 1,671 table names are
  the same link, and against the search's own names, 1,255 of 1,255;
- speed: `farsidename` spends a few ms per witness, plus ~2 s per row decoding the row's ambient
  and finding its isomorphisms onto the thickening once (`farsideredraw.h`, `outgoingLinkFast()`).
  The first version took ~0.9 s per witness: 600 ms decoding each pair signature from
  scratch (the ambient, a skeleton, and a `KnottedSurface` built only to read back the
  face list), and 300 ms rebuilding a `KnottedSurface` on the thickening, whose `addFaces()`
  re-runs every embeddedness and local-flatness check with a cold cache. The isomorphism
  search itself was ~20 ms. `farsidename --reference` keeps the old path, and on a uniform
  sample of 1,000 master witnesses the two agree on 999. The other one got the same link
  pinned by the reference and left as a proved 3-variant set by the fast path (a different
  drawing).
