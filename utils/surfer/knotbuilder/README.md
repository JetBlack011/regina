# knotbuilder

Triangulating link diagrams, and drawing curves of those triangulations back
as diagrams.

| file | what |
|---|---|
| `knotbuilder.{h,cpp}` | `buildLink()`: a PD code → a triangulation T of S³ with the link L in its 1-skeleton (adapted from the Texas Experimental Geometry Lab students' construction; paper, Appendix A) |
| `blockgeometry.{h,cpp}` | the crossing block as the box [-1,1]² × [0,1], and an exact check of that embedding |
| `diagramdrawer.{h,cpp}` | `DiagramDrawer`: edge cycles of T → an oriented PD code with crossing signs and linking numbers |
| `triangulateknot.cpp` | tool: a PD code → T's and its complement's isomorphism signatures |
| `tests/knotbuilder_test.cpp` | the construction itself |
| `tests/diagramdrawer_test.cpp` | the drawer (below); given a table, sweeps it |
| `tests/triangulate_knot_table.cpp` | tool: `buildLink()` over a whole table (~13 min) |

`../farsidecurves.{h,cpp}` and `../farsidediagram.cpp` put the drawer to
work on the search's far sides; see the end of this file.

## The construction

`buildLink()` replaces each crossing by a **block** of 14 tetrahedra (a 3-ball,
isomorphism signature `oavfuGyfLaGbaahgklaggo`; tetrahedra 14k..14k+13 for
crossing k, in `Block`'s order: cores 0–5, then walls 0–7), glues the blocks
wall to wall — one wall per slot of the crossing, the two walls of each
strand label glued together — and cones off the two boundary spheres of the
resulting S² × I (`finiteToIdeal()`). The link is 3 edges per crossing: the
under-strand is one edge, the over-strand two.

`CobordismBuilder` thickens T into T × [0,2] (prism layers). A layer's top
face is a copy of the tetrahedron below it, so the search's outgoing boundary
is literally T × {2}: every far side is an edge cycle of T itself.

## The block is a box

All 13 vertices of the block lie on the surface of the box [-1,1]² × [0,1]
(integer units of `BLOCK_SCALE` = 2²⁰ in the code):

| vertices | where |
|---|---|
| 4 bottom corners | (±1, ±1, 0) |
| 4 top corners | (±1, ±1, 1), each above a bottom corner (a vertical edge joins them) |
| 4 strand points | midpoints of the bottom edges: (0,−1,0) and (0,1,0) on the under-strand, (−1,0,0) and (1,0,0) on the over-strand |
| top centre | (1/7, 1/11, 1): the over-strand's peak, nudged off the centre so that it does not project onto the under-strand |

- **Bottom face:** an octagon (the four corners, with the strand points as side midpoints). Its inner square is split by the under-strand edge.
- **Top face:** a square fanned from the top centre.
- **Walls:** each wall is a vertical rectangle fanned from its strand point. Its triangulation is symmetric left to right, so two glued blocks agree on the wall they share.
- **Strands:** the under-strand runs straight across the bottom face. The over-strand arches from one strand point up through the top centre and down to the opposite one.

`blockCoordinates()` derives each vertex's point from its combinatorial role
in a freshly built `Block` (never from Regina's vertex numbering), and
`verifyBlockModel()` proves it is an embedding, in exact integer arithmetic:

- every tetrahedron is nondegenerate and oriented consistently with the block;
- the signed volumes (1/3 or 1/6 of the unit box) sum to the box's volume;
- every boundary triangle lies in one face of the box, and each face is
  covered exactly (its triangles' areas sum to the face's).

The boundary goes onto the box's boundary with degree one and every volume
is positive, so every point of the box lies in exactly one tetrahedron: the
block tiles the box. `DiagramDrawer`'s constructor further checks, per
triangulation, that every block is oriented alike and that each vertex sits
at one point per block (true exactly when the diagram is reduced).

## Drawing a curve

Project each block vertically onto its square. The squares tile the diagram's
sphere (one square per crossing; its corners are the diagram's regions), so a
curve of T becomes a curve on S², and height says what passes over what:

- **Inside a block.** An edge through a block's interior, or in its top or
  bottom face, projects injectively into that block's square. Crossings are
  found pairwise inside each block, in that block's own coordinates.
- **Walls.** Every degenerate projection comes from an edge lying in a vertical
  wall: the shared side of two blocks, or a vertical corner line. There are no
  vertical edges inside a block. Such an edge is replaced by a small **tent**
  pushed into one adjacent block. Only the edge's middle moves, a distance of
  1/200 of the block into the box's interior, and the push grows with height,
  so that wall points differing only in height separate. A vertical corner edge
  gets a two-point tent, so that its projection does not fold back on itself.
  This is an isotopy of the curve. Near the wall the curve has nothing else:
  at a vertex it has only its own two edges, and every other edge leaves the
  wall at once.
- **Corners.** A region's bottom and top points project to the same corner of
  every square around it. If the curve passes through both (not along the
  vertical edge joining them), walk the squares around that corner
  counterclockwise. The cyclic order of the curve's four ends there decides
  whether it is a crossing (interleaved: the top pass over the bottom one) or
  two passes that merely touch.
- **Cones.** A pass through a cone apex (curve edges u → apex → w) becomes an
  arc over everything (top apex) or under everything (bottom apex). It is
  routed through the squares by a breadth-first path, crossing sides at fixed
  points that are never strand points. Its height is monotone along the arc,
  so the arc stays unknotted.

A PD code needs only the cyclic order of the four ends at each crossing. Every
crossing is found inside one block or at one corner, so the drawer never lays
the diagram out globally. All arithmetic is exact: integer coordinates, with
128-bit products for cross products and for comparing crossing parameters. A
degenerate projection (two pieces touching, collinear overlaps, equal heights)
throws `Degenerate`, and `draw()` retries with other perturbations. None has
been needed on any input so far.

**Orientation.** Component i of the result is the i-th `EdgeCycle` given,
traversed in its given direction. The PD code follows the KnotTheory
convention (the incoming under-strand first, then counterclockwise), and
`CrossingInfo::sign` is +1 for a right-handed crossing. `linkingNumber(i, j)`
is half the signed sum over crossings between i and j.

## How it is validated

| check | result |
|---|---|
| `verifyBlockModel()` | exact embedding |
| every table row (all 12,965 knots through 13 crossings, all 4,188 oriented links through 11): draw knotbuilder's own link from T | **all identical to the input diagram**: the same Regina signature with neither mirror nor reversal allowed, and the same signed linking matrix up to relabelling components. 7.6 s for the knots, 2.1 s for the links |
| negative test: one component of the Hopf link drawn reversed | caught: the signature changes and the linking number is negated |
| random closed curves in T (80 in `ctest`, 95% through a wall, corner line or cone point) | drawn unknot ⇔ the complement drilled from T is a solid torus |
| 400 real far sides from 287 rows (114 knots, 286 links, prime and knotted per the diagram pipeline), via `farsidediagram` | **400/400** named as the pipeline recorded, by its own `fsid.identify_piece` on the drawn diagram |
| 122 genus-0 witnesses whose surface is a union of annuli, pairing row components with far-side components (82 with nonzero linking numbers) | **all 122** preserve every linking number exactly, as a concordance of oriented links must |

The table sweep, the random curves and the negative test are in
`tests/diagramdrawer_test.cpp`. The witness checks are
`cobordism-atlas/tools/farside/drawfromT/test_cpp_farsides.py`.

The sweep only exercises the paths knotbuilder's own link uses: interior edges
and the bottom-face diagonal. Walls, corners and cones are covered by the
random curves and, more strongly, by the real far sides.

## Using it

```cpp
auto [tri, edges, reversed] = knotbuilder::buildLink(knotbuilder::parsePDCode(pd));
knotbuilder::DiagramDrawer drawer(tri, crossings);          // once per row
knotbuilder::Diagram d = drawer.draw(curves);               // curves: std::vector<EdgeCycle>
regina::Link link = d.link();                               // crossingless components included
long lk = d.linkingNumber(0, 1);
```

`draw()` is `const`, so one drawer can serve many threads. A far side is first
carried onto T by `farside::OutgoingMap` (`../farsidecurves.h`). That map uses
the thickening's own top prisms, with no isomorphism search, so an
orientation-reversing automorphism of T cannot mirror the far side.
`farside::orientedOutgoingLink()` orients each surface component against the
row, as the search's orientation check does. `farsidediagram '<row PD>' <
pairsigs` does all of this for witnesses decoded from their pair signatures
(`--layers 1` for the earliest runs' witnesses).

**In the search.** `verifyslicegenus` names every far side through this
drawer (`../farsidenaming.{h,cpp}`, a `BoundaryNamer` that `SurfaceSearch`
consults on each boundary-cache miss). It draws the curves, simplifies with
Regina's `Link::simplify()`, and names by exact diagram signature against the
tables. The complement route remains only as a fallback; see
`../README.md`.

**Speed.**
- About 0.5 ms per table row, including building T.
- In the search, one far side is drawn, simplified and named in about 0.15 ms.
  The complement route took roughly 20–60 ms per far side.
- On halcyon at 1M surfaces:
  - each row names 445k–720k far sides in 12–41 s of thread time, with 0–28
    fallbacks;
  - the drain now finishes seconds after the search;
  - rows went from 707–1,079 s (fixed binary, complement naming) to
    246–326 s.
- The offline `farsidediagram` spends about 0.3 s per witness, almost all of it
  in the 4-dimensional isomorphism search that puts a decoded witness back on
  the row's thickening. The search itself never needs that step.
