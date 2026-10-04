# diagramtriangulation: the triangulations a search runs in

From a link diagram (a PD code) this builds a triangulation T of S³ with the
link L in its 1-skeleton, thickens T into T × [0, k] with L's collar as the
search's seed, and draws curves of T back as oriented diagrams. Every
outgoing link a search finds is a set of curves of T (the outgoing boundary is
a copy of it), which is why drawing them back is exact and cheap. Paper:
`sec:knotbuilder` and Appendix A (`app:triangulations`, by Samantha Ward,
Angel Yuan and Jingyuan Zhang, from whose Texas Experimental Geometry Lab
construction this is adapted).

It depends only on `linknaming_complement` (the edge-cycle helpers) and
Regina. Namespace `diagramtriangulation`, except the thickening's classes.

| file | what |
|---|---|
| `pdcode` | PD codes as text: `pdLabels()`, `parsePDCode()` (the one parser T is built from), `formatPDCode()` (the one formatter) |
| `fromdiagram` | `buildLink()`: a PD code → T and L's edges, 3 per crossing, each with the direction the PD code gives it; `Block`, the crossing's 14 tetrahedra; `reduceVertices()` |
| `block` | the crossing block as the box [-1,1]² × [0,1] (`blockCoordinates()`), and `verifyBlockModel()`, an exact proof that it is an embedding |
| `todiagram` | `DiagramDrawer`: closed edge curves of T → an oriented diagram (`Diagram`: PD code, crossing signs, signed Gauss data, linking numbers, `link()`) |
| `thickening/prism` | `SimplicialPrism`: a simplex × [0,1] as simplices |
| `thickening/thickening` | `CobordismBuilder`: T × [0, k] by prism layers. `OutgoingMap`: the outgoing boundary read back as T, edge by edge. `ThickenedLink` and `buildAmbient()`: what a search runs in |
| `thickening/collar` | `CollarBuilder`: the collar L × [0, c] through the first c layers, the seed every searched surface contains |
| `triangulatediagram.cpp` | the tool `triangulateknot` (`build/utils/surfer/knotbuilder/triangulateknot`, the path it had before this part was renamed, kept for commands written against it): a PD code → the isomorphism signatures of T and of L's complement |

## PD codes

`pdLabels()` reads a PD code's integers in fours, as written, whatever the
punctuation: `[[1;5;2;4];...]`, the knot table's spaced spelling, LinkInfo's
`PD[X[4; 1; 3; 2]; ...]`. `parsePDCode()` is `pdLabels()` renumbered from 0
(a code that already holds a 0 is taken as 0-based), so T depends only on the
integers and their order, never on the spelling.

`formatPDCode()` writes the two spellings stored files use, both frozen
(`PDSpelling`): `[[a;b;c;d];...]` for every PD a goal run records (a given PD
is respelt as `formatPDCode(pdLabels(text), semicolons)`: the same integers in
the same order, so the same T), and `[[a,b,c,d],...]` for `cobound draw`'s
`pd=` field. Reading a PD code into a `regina::Link` is linknaming's
(`linkFromTablePD()`).

## The construction

`buildLink()` replaces each crossing by a **block** of 14 tetrahedra (a 3-ball;
tetrahedra 14k..14k+13 for crossing k, in `Block`'s order: cores 0–5, then
walls 0–7). It glues the blocks wall to wall, one wall per slot of the
crossing and the two walls of each strand label glued together, and cones off
the two boundary spheres of the resulting S² × I (`finiteToIdeal()`). The
link is 3 edges per crossing: the under-strand one edge, the over-strand two.
`reversed[i]` says whether edge i runs against its vertex order in the PD
code's direction, so orientation variants of one link get different
directions (`L6a3`'s variants are told apart).

`reduceVertices()` pinches away every internal edge but L's (and loops), for
tools that want a small triangulation; a search never uses it.

## The block is a box

All 13 vertices of the block lie on the surface of the box [-1,1]² × [0,1]
(integer units of `BLOCK_SCALE` = 2²⁰):

| vertices | where |
|---|---|
| 4 bottom corners | (±1, ±1, 0) |
| 4 top corners | (±1, ±1, 1), each above a bottom corner (a vertical edge joins them) |
| 4 strand points | midpoints of the bottom edges: (0,−1,0) and (0,1,0) on the under-strand, (−1,0,0) and (1,0,0) on the over-strand |
| top centre | (1/7, 1/11, 1): the over-strand's peak, nudged off the centre so that it does not project onto the under-strand |

- **Bottom face:** an octagon (the four corners, with the strand points as
  side midpoints). Its inner square is split by the under-strand edge.
- **Top face:** a square fanned from the top centre.
- **Walls:** each wall is a vertical rectangle fanned from its strand point,
  symmetric left to right, so two glued blocks agree on the wall they share.
- **Strands:** the under-strand runs straight across the bottom face; the
  over-strand arches from one strand point up through the top centre and
  down to the opposite one.

`blockCoordinates()` derives each vertex's point from its combinatorial role
in a freshly built `Block` (never from Regina's vertex numbering), and
`verifyBlockModel()` proves it an embedding in exact integer arithmetic:
- every tetrahedron is nondegenerate and oriented consistently with the block;
- the signed volumes sum to the box's volume;
- every boundary triangle lies in one face of the box, and each face is
  covered exactly (its triangles' areas sum to the face's).

The boundary goes onto the box's boundary with degree one and every volume is
positive, so every point of the box lies in exactly one tetrahedron: the block
tiles the box. `DiagramDrawer`'s constructor checks, per triangulation, that
every block is oriented alike and that each vertex sits at one point per
block, which holds exactly when the diagram is reduced; it refuses any other
triangulation.

## Drawing a curve

Project each block vertically onto its square. The squares tile the diagram's
sphere (one per crossing; their corners are the diagram's regions), so a
curve of T becomes a curve on S², and height says what passes over what:

- **Inside a block.** An edge through a block's interior, or in its top or
  bottom face, projects injectively into that block's square. Crossings are
  found pairwise inside each block, in that block's own coordinates.
- **Walls.** Every degenerate projection comes from an edge lying in a
  vertical wall: the shared side of two blocks, or a vertical corner line.
  There are no vertical edges inside a block. Such an edge is replaced by a
  small **tent** pushed into one adjacent block: only the edge's middle moves,
  by 1/200 of the block into the box's interior, and the push grows with
  height, so that wall points differing only in height separate. A vertical
  corner edge gets a two-point tent, a thin loop that leaves the corner and
  comes back to it. This is an isotopy of the curve: near the wall the curve
  has nothing else, since at a vertex it has only its own two edges and every
  other edge leaves the wall at once.
- **Corners.** A region's bottom and top points project to the same corner of
  every square around it. If the curve passes through both, walk the squares
  around that corner counterclockwise: the cyclic order of the curve's four
  ends there decides whether it is a crossing (interleaved: the top pass over
  the bottom one) or two passes that merely touch. That includes a curve
  running along the vertical edge joining them, whose tent's shadow passes
  the corner twice, on the way in at the bottom point and on the way out at
  the top one.
- **Cones.** A pass through a cone apex (curve edges u → apex → w) becomes an
  arc over everything (top apex) or under everything (bottom apex), routed
  through the squares by a breadth-first path that crosses sides at fixed
  points, never strand points. Its height is monotone along the arc, so the
  arc stays unknotted.

A PD code needs only the cyclic order of the four ends at each crossing, and
every crossing is found inside one block or at one corner, so the drawer never
lays the diagram out globally. All arithmetic is exact: integer coordinates,
with 128-bit products for cross products and for comparing crossing
parameters. A degenerate projection (two pieces touching, collinear overlaps,
equal heights, three pieces through one point) throws `Degenerate`, and
`draw()` retries with other perturbations (8 by default).

**Planarity.** Every closed curve in S³ has a planar diagram, so a drawing that
is not one (`regina::Link::isClassical()` false) is a defect of the drawer.
`draw()` throws `NonPlanar` for it and never retries, since another
perturbation could hide the defect rather than remove it. The check runs on
the raw drawing, before any simplification. A search names such curves by
their complement and counts them (`linknaming/README.md`).

**Orientation.** Component i of the result is the i-th `EdgeCycle` given,
traversed in its given direction. `Diagram::link()` builds the Regina link
from each component's signed crossing sequence (`Diagram::gauss`,
`regina::Link::fromData()`), not from the PD code: a PD code cannot fix the
orientation of a component that passes over every crossing it meets
(`regina::Link::pdAmbiguous()`). The PD code follows the KnotTheory
convention (the incoming under-strand first, then counterclockwise),
`CrossingInfo::sign` is +1 for a right-handed crossing, and
`linkingNumber(i, j)` is half the signed sum over the crossings between i and
j. Crossingless components are listed in `crossingless`: each is an unknot
split from the rest.

## The thickening

`CobordismBuilder` takes an ordered copy of T (relabelling vertices within
tetrahedra if need be; it fails rather than proceed when it cannot) and adds
prism layers: each `thicken()` adds a copy of the prism pieces over every
tetrahedron, 4 pentachora each, glued across T's gluings, so k layers
triangulate T × [0, k]. There is no coning: a search runs in S³ × I with two
boundary spheres, the incoming T × {0} and the outgoing T × {k}.

- **The outgoing boundary is T.** Over each tetrahedron σ, the top prism piece
  has σ × {k} as a facet. `OutgoingMap` uses exactly that to carry each edge of
  the outgoing boundary onto an edge of T, with no isomorphism search, so no
  automorphism of T (some reverse orientation) can relabel or mirror an
  outgoing link. `carry()` takes an oriented curve, `carryCycle()` an
  unordered closed one.
- **The seed.** `CollarBuilder` traces L's edges through the first c layers:
  the collar L × [0, c], the triangles every searched surface contains.
- **What a search runs in.** `buildAmbient(pd, layers, collarLayers, out)`
  fills a `ThickenedLink`: the PD code, T and L's edges, the builder, the
  4-dimensional ambient, the collar's faces in index order, the incoming
  boundary component and L's component count (counted by walking its
  edges). `cobound` thickens with 2 layers, collared through both.

T and its thickening are a function of the PD code's integers alone, and
must stay so byte for byte: stored certificates record surfaces by faces
against a digest of the thickening, and frontier fingerprints, pair
signatures and the read-back cache all key on it.

## How it is validated

| check | result |
|---|---|
| `verifyBlockModel()` | an exact embedding |
| every table row (12,965 knots through 13 crossings, 4,188 oriented links through 11): T's own link drawn from T | **all identical to the input diagram**: the same Regina signature with neither mirror nor reversal allowed, and the same signed linking matrix up to relabelling components; 7.6 s for the knots, 2.1 s for the links |
| one component of the Hopf link drawn reversed | caught: the signature changes and the linking number is negated |
| random closed curves in T (80 in ctest, 95% through a wall, corner line or cone point) | drawn unknot ⇔ the complement drilled from T is a solid torus |
| **every** simple closed curve of ≤ 5 edges in T for 3_1, 4_1, 5_2 and L4a1{0} (27,472 curves), each drawn with all 8 perturbations and in both directions (439,549 drawings) | **all planar**; every 25th also agrees with drilling. `SHORT_CYCLES_MAX_LEN=6` (172,749 curves, 2.76M drawings) passes too |
| pairs of disjoint curves of ≤ 5 edges, one through a region's bottom point and one through its top point (6,833 pairs at every corner of 3_1, 4_1 and L4a1{0}; 68,325 drawings) | **all planar**; every perturbation, both component orders and global reversal give the same HOMFLY polynomial and linking number, and Regina's linking number equals the drawer's |
| 400 stored outgoing links from 287 searches (114 knots, 286 links), drawn from their pair signatures | **400/400** named as SnapPy's pipeline had named them, by that pipeline's own piece identification on the drawn diagram |
| 122 genus-0 cobordisms made of annuli, pairing incoming and outgoing components (82 with nonzero linking numbers) | **all 122** preserve every linking number exactly, as a concordance of oriented links must |
| 4,416 stored cobordisms of 77 searches, redrawn | 0 non-planar drawings, 0 failures |

The table sweep, the random curves, the short curves and the corner pairs are
`tests/todiagram_test.cpp` (the sweep with a table as argument). The stored
cobordisms are checked by the atlas's
`tools/farside/drawfromT/test_cpp_farsides.py`.

The table sweep exercises only the paths T's own link uses: interior edges and
the bottom-face diagonal. Walls, corners and cones are covered by the random
curves, the stored outgoing links and, exhaustively for short curves, by the
short-curve and corner-pair tests: a short curve stays near one wall or
corner, so they cover every local configuration a curve of that length can
make, including a curve running up a vertical corner edge, which neither the
sweep nor a few random curves reach.

## Using it

```cpp
using namespace diagramtriangulation;
TriangulationWithLink t = buildLink(parsePDCode(pd));      // T, L's edges, their directions
DiagramDrawer drawer(t.tri, crossings);                     // once per T
std::vector<EdgeCycle> curves = DiagramDrawer::cyclesOf(t.edges, t.reversed); // L itself
Diagram d = drawer.draw(curves);                            // throws NonPlanar, Degenerate
regina::Link link = d.link();                               // crossingless components included
long lk = d.linkingNumber(0, 1);
```

`draw()` is `const`, so one drawer serves many threads. Drawing T's own link
and building T take about 0.5 ms per table row; drawing an outgoing link and
naming its diagram takes about 0.15 ms, against 20–60 ms for naming it by its
complement, which is why a search names its outgoing links by drawing them.

## Tests

`ctest` in `build/utils/surfer/diagramtriangulation`.

| test | pins |
|---|---|
| `pdcode_test` | `parsePDCode()` renumbers from 0, reads a code holding a 0 as 0-based, whatever the punctuation (bare integers too); `pdLabels()` keeps labels and crossing order as written; `formatPDCode()` in each stored spelling, byte for byte; respelling LinkInfo's and the spaced knot table's spellings keeps the integers and their order, so T; the round trip. The round trip over all 17,153 table PD codes is `linknaming/tests/tables_test` |
| `fromdiagram_test` | `buildLink()` gives a valid closed S³ with the link's components (trefoil, Hopf link, many named knots, a non-alternating regression); the edge directions are consistent and tell `L6a3`'s variants apart; `reduceVertices()`; the output can be ordered; T thickened for the trefoil, the Hopf link and the figure-eight |
| `todiagram_test` | the block model; T's own link drawn back exactly, orientation included, for knots and link variants; the reversed-component negative test; random curves against drilling; every short curve planar; every corner pair planar and consistent |
| `thickening_test` | prism layers give a valid product: balls, a self-glued S³ over one and two layers, doubled tetrahedra, ∂Δ⁴, every facet pair of two glued tetrahedra, a doubled triangle in dimension 2; `isOrdered()` checks every facet; the base boundary component is the bottom; `buildAmbient()`'s seed is in index order |

Not in ctest: `triangulate_knot_table <table.csv> [limit]` builds T for every
row of a table and checks it is a closed S³ with a well-formed link (minutes
for the knot table); `todiagram_test <table.csv> [max crossings]` sweeps a
table (the second row of the validation table).
