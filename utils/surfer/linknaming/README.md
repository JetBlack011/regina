# linknaming: naming a link, with proof

Names an oriented link diagram, with a proof, as a table entry or a split
union or connected sum of table entries; names a triangulated complement by
the census; and gives the census-free certificates for the unknot and unlinks
that the search's local checks rest on. It knows nothing of surfaces. Paper:
§5 (`sec:comp`), in particular `sec:naming-drawing`, `prop:drawing` and
`prop:piece-names`; how names bound is `lem:split-union`,
`lem:sum-along-components` and `cor:sum-pieces`.

Two libraries. `linknaming_complement` (`complement/`, Regina only) is the
bottom layer of the project: links as subcomplexes of a triangulated
3-manifold's 1-skeleton. `linknaming_lib` (the rest) adds SQLite for the local
census. Namespace `linknaming`, with `complement` and `census` for those
folders.

| file | what |
|---|---|
| `names` | the grammar of names and descriptions, read back: component counts, bases and tags, splits, alternatives, composites, sums, knot marks, `complement:` |
| `tables` | the knot and link tables: `readTableRows()` (the one reader), `parseTableG4()`, symmetry types, `isElementarySlice()` (the slice-composite anchors), `Tables` (every entry as four oriented diagrams, indexed), `linkFromTablePD()` |
| `linknamer` | `LinkNamer`: a diagram → a name or a description (`name()`), the route a search names its outgoing links by (`nameDrawing()`), link classes (`canonicalName()`); `NamerLimits`; `NamingStats` |
| `isometry/` | `KernelLink`: a diagram's complement in Regina's SnapPea kernel, and `sameLinkAs()`, the kernel's isometry test carrying meridians to meridians. `snappeakernel_isometry{,_cusped}.cpp` compile the kernel's own isometry sources, which Regina ships unbuilt (`engine/snappea/kernel/unused/`) |
| `diagrams/gaussdiagram` | `GaussDiagram`: an oriented diagram as signed crossings and each component's crossing sequence; `splitPieces()`, `visibleSum()` |
| `diagrams/diagramiso` | `findDiagramIsomorphism()`: diagram isomorphisms, with the component map they realise |
| `diagrams/simplification` | simplification that keeps every component in its slot and every linking number; nugatory crossings removed; components over or under everything lifted off |
| `complement/edgecycles` | walks over a set of edges: directed chaining, orienting undirected edges, counting closed curves |
| `complement/linkcomplement` | `EdgeComplement`, `Knot`, `Link`: curves of a triangulation, and their complement by drilling |
| `complement/meridians` | drilling that keeps signed meridians, and slopes in SnapPea's peripheral basis (`cobound meridians`) |
| `complement/unlinknaming` | the unknot and unlinks from their complements, census-free: `isUnknot()`, `certifiesUnlink()`, `capInCone()` |
| `complement/complementcache` | the one cache of complement answers |
| `census/censusnaming` | `census::nameComplement()`: a complement named by the local census, Regina's census and a Pachner search |
| `census/censusnames.h` | census names → table names (generated) |
| `tableclasses.cpp` | the tool `tableclasses` (`build/utils/surfer/tableclasses`): the tables' link classes |
| `gen_census_names.cpp`, `gen_census_names_header.py` | the generator pair behind `census/censusnames.h` |

## Names and descriptions

A **name** is an identity: it determines the link up to mirror and global
reversal, neither of which a slice genus sees, so it may stand as one link and
receive bounds (`LinkName::isName`). Anything weaker is a proved
**description**: what it says is true, but it does not single out one link.

| shape | written | |
|---|---|---|
| unknot, unlink | `Unknot`, `<n>-component unlink` | name |
| one piece, its variant pinned | its table name (`L7n1{1}`, `8_17`) | name |
| one piece with split unknots | `L7n1{1} u Unknot` | name |
| an untabulated piece | `diagram:<sig>`, the oriented diagram's signature (reversal disallowed, the lesser of it and its global reverse), which determines the link | name |
| a composite knot | `3_1#m3_1`, `3_1#mr9_32`: each summand marked `m`/`r` where its symmetry type makes the mark matter, spelled with the fewest marks over the four global mirror and reversal choices | name when every mark that matters is pinned |
| a piece whose variant is not pinned | `A\|B`, the variants that survive | description |
| a split union | `A u B` (Unknot last) | description |
| knots summed into a link | `3_1 #_? L2a1{0}` | description |
| links summed along components | `#{L2a1{0}[?] # L4a1{1}[?]}` | description |
| a link of 2+ components named only by its complement | `complement:<name>` | description; bears nothing |

`?` marks a component index the namer does not compute. In `#{...}` a piece
summed at two sites is written at both, so a reader counts each occurrence as
a piece, which only weakens the bounds built from the pieces
(`cor:sum-pieces`). A knot never gets `complement:`: a knot's complement names
it (Gordon–Luecke), so a knot named by its complement carries that name bare.
`names.h` reads all of this back; both solvers and the search use it.

## Naming a diagram (`LinkNamer::name()`)

1. **Decompose.** Simplify (Regina's `simplify()`, then `simplifyExhaustive()`
   up to `exhaustiveHeight`; neither reflects nor reverses). A disconnected
   diagram is a split union. Two arcs on the same two faces whose removal
   disconnects the crossings are where a circle meets the diagram: a connected
   sum along the component through both, each side closed along the circle.
   Repeat on the pieces. Components left in no piece are split unknots.
2. **Name each piece**, in this order.
   - *By its diagram:* after up to `simplifyTries` more simplifications, the
     piece's diagram **is** one version of a table entry. That pins its
     oriented variant, and its mirror and reversal relative to the drawing.
   - *By isometry* (hyperbolic pieces): the table entries whose HOMFLY
     polynomial, in either mirror, is the piece's and whose table diagram is no
     bigger are a shortlist, a hint only. The piece's complement is built by
     the SnapPea kernel from the diagram, so its peripheral curves are the
     diagram's meridians and longitudes, and a match is proved when an isometry
     carries every meridian to a meridian (below). Filling along the meridians,
     the piece is the table link up to mirror and orientations; the isometry's
     action on the oriented meridians then pins the variant (and a knot's
     mirror and reversal), with the HOMFLY polynomial and linking numbers
     checked against it (a disagreement throws). Milliseconds, whatever
     diagram was drawn.
   - *By search* (what the isometry leaves, chiefly non-hyperbolic pieces): a
     bounded `rewrite()` (`searchHeight`, `searchVisits`, then a deeper
     `deepHeight` round for small pieces) reaches a diagram of a table base;
     failing that, from the table's side, the flype orbit of an alternating
     candidate's table diagram, or a `rewrite()` outward from a candidate's
     table diagram (heights 1 to `tableSideHeight`). `rewrite()` may reflect and
     reverse, so this proves the base only. Every oriented variant V of the
     base, and mV, is then ruled out if its HOMFLY polynomial or linking
     numbers differ from the piece's (computed on the drawn orientation); the
     table lists every orientation up to global reversal, so the survivors
     contain the piece.
   - *Otherwise* the piece is untabulated: `diagram:<sig>`.
3. **Compose** the name in the grammar above.

`NamerLimits` sets every bound. `cobound name`'s defaults are its fields'
(`simplifyTries` 24, `exhaustiveHeight` 1, `searchHeight` 2 over 20,000
diagrams up to 16 crossings, `deepHeight` 3 over 200,000 up to 12, the table
side to height 4 up to 12); `namer_*` keys change them.

## The route a search names by (`LinkNamer::nameDrawing()`)

A search draws each outgoing link (`diagramtriangulation`'s drawer) and names
the drawing here, with the complement as the fallback
(`cobound/outgoing/outgoingnamer`):

- **A drawing that failed** (degenerate, or refused as not planar: a drawer
  defect) goes to the complement route: a knot is named by its complement, a
  link only described, `complement:<name>`, unless the complement proves it an
  unlink.
- **A knot:** simplified once. No crossing left: the unknot. Then its **table
  signature** (the simplified diagram is a version of a table knot,
  `Tables::exact()`); a complement answer learned earlier for that diagram;
  its pieces (cut at visible sums, each by its diagram, its isometry and, as
  the limits allow, its search). A knot some piece of which no table names is
  named by its complement, and the answer is learned against its simplified
  diagram.
- **A link:** `name()` of its pieces, oriented as drawn. One whose untabulated
  piece might be an unlink in disguise (every linking number 0 and the
  unlink's Jones polynomial) is tried on the complement route, which only ever
  replaces the name by a proved unlink.
- **Memory.** Every drawing's answer is kept against its exact signature
  (cleared at 1,000,000), so a drawing met again costs one signature and keeps
  its name.
- **Tests.** Under `SURFER_TEST_PERTURB_NAMES` every name and description gets
  a fresh suffix (`name_independence_test`: no gate may depend on a name).

The search's limits (`OutgoingNamer::limits()`) are `simplifyTries` 2, a
knot's `knotSimplifyTries` 0, no exhaustive simplification, no Reidemeister
search and no table side. Measured on a production search (`10_141`, 1M
surfaces):
- one simplification misses the table signature for 5.7% of distinct table
  knots drawn, but a retry recovered none of the ~275 misses in any run, since
  it re-simplifies an already simplified piece; the isometry names every one
  of them (no knot went to the complement in 15 runs);
- the table-side search on a link that is none of its candidates runs to the
  end, a minute or more per name; `cobound name` refines such names offline;
- naming every outgoing link costs about 1% of the search's wall time
  (20.48 s against 20.34 s without naming, peak memory 713 against 694 MB):
  4.5 s of thread time for 753,055 edge sets.

`NamingStats::summary()` is the search's `diagram naming:` line: per edge set
(the frozen `far sides drawn`), how each was named (`unknot`, `unlink`,
`table knot`, `learned knot`, `table link`, `other link`, `learned link`),
complement fallbacks (learned, non-planar drawings), oriented names per
surface (cached, by complement), the time in each route, and the slowest name.

## Tables

`Tables::load()` holds every table entry (a knot, or one oriented variant
`L7n1{1}` of a link) as the diagram its PD code draws, in four versions: as
written, mirrored, reversed and both. Each version is keyed by Regina's
signature with neither reflection nor reversal allowed (`exact()`), so a
diagram whose own such signature is a key **is** that version of that entry.
A second index, by the signature allowing both, gives a diagram's base for
the search step (`base()`). Two bases sharing an unoriented diagram are a
table error, refused on loading.

`readTableRows()` is the one reader of the table files (`Name,PD,Genus-4D`).
`parseTableG4()` reads `2` or `[0;1]`, and returns nothing for anything else,
so a malformed value never becomes a bound. `readSymmetryTable()` reads the
atlas's `knot_symmetry.csv`. `isElementarySlice()` decides the anchors: a
composite knot whose summands pair off into concordance inverses
(`K # m(Kʳ)` bounds a ribbon disc), read from each summand's marks and
symmetry type. `3_1#m3_1` and `4_1#4_1` are always anchors; any other needs
every summand's symmetry type, and a summand without one is refused, never
guessed. So `9_32#mr9_32` is an anchor, while `9_32#m9_32` and `8_17#m8_17`
(that is, `8_17 # 8_17ʳ`) are not.

**Link classes.** Some table names are one oriented link up to mirror and
global reversal: their diagrams coincide (`Tables::canonical()`), or an
isometry carries meridians to meridians with one orientation sign.
`LinkNamer::canonicalName()` names each such class by one name, computed per
base on first use, and the namer always writes the class's name. Over the
tables, 933 names are in classes (408 by diagram, 525 by isometry), and every
class's literature 4-genera agree. `tableclasses` writes them as the atlas's
`data/table_link_classes.csv` (`name,canonical,proof`), which `cobound solve`
reads with `link_classes`, and fails on a class whose members' 4-genera
differ. LinkInfo lists one non-hyperbolic link twice, `L11n432` and
`L11n437` (mirror images); no isometry can merge them, so they stay two names.

## The isometry

`KernelLink` builds a diagram's complement in the SnapPea kernel
(`regina::SnapPeaTriangulation(const Link&)`), so its peripheral curves are
the diagram's own meridians and longitudes. `sameLinkAs()` is the kernel's
`compute_isometries()`: it canonises both complements (the Epstein–Penner
cell decomposition, found numerically) and tries every combinatorial
isomorphism of the canonical triangulations, computing each one's action on
the peripheral curves in integers. A match takes every meridian to a meridian
(`isometry_extends_to_link()`). A positive answer is therefore exact: filling
every cusp along its meridian extends the isomorphism to a homeomorphism of S³
carrying one link onto the other. Floating point only steers the
canonisation; a wrong turn there makes a miss, never a false match, and a
piece that matches nothing is retried on randomised triangulations. Only
hyperbolic complements are compared. The kernel is not documented as
thread-safe, so every call into it is serialised by one mutex; each takes
milliseconds.

## Component maps

`findDiagramIsomorphism(a, b, allowMirror, allowReverse)` decides whether two
diagrams are the same up to relabelling crossings, permuting components and
choosing each component's start, optionally also up to mirror image and
reversing every component, never some components. It is exact and
exhaustive, and returns the component map the isomorphism realises; given a
required map, it considers only isomorphisms realising that map, so with
`a == b` it asks whether a permutation of components is a symmetry. Goal runs
build every map between components on it (`cobound/README.md`, "Component
maps").

`simplifyKeepingComponents()` simplifies with Reidemeister moves only, keeps
each component in its slot (a crossingless one included), removes every
nugatory crossing, and throws if the component count or any pairwise linking
number changed. `liftSplitComponents()` lifts off a component that is over
(or under) every crossing it meets: an isotopy splitting off an unknot whose
orientation a PD code could not carry.

## Complements

- **Drilling** (`linkcomplement`): a set of edges is drilled out of a
  triangulation by pinching each edge (Weeks' two-tetrahedron gadget, as He,
  Sedgwick and Spreer use it; paper `sec:drilling`), then simplified.
- **Meridians** (`meridians`): drilling that keeps each component's meridian,
  signed by the curve's orientation, and expresses it as a slope in SnapPea's
  peripheral basis. A complement does not determine a link, but the
  complement with its meridians does. The frame is common only if SnapPea
  installs its basis into our own triangulation written verbatim
  (`Triangulation<3>::snapPea()`): wrapping the drilled complement in
  `regina::SnapPeaTriangulation` removes its finite vertices and
  retriangulates it, after which its peripheral curves describe another
  triangulation. `cobound meridians` is the batch front end the atlas's SnapPy
  pipeline drives.
- **The unknot and unlinks** (`unlinknaming`): a curve is the unknot when its
  complement is a solid torus (handlebody genus, cached); several curves are
  an unlink when π₁ of the complement simplifies to a free group on as many
  generators with no relations. `certifiesUnlink()` (a petal trace T_v(S) in
  Lk(v)) also checks that the edges are disjoint closed curves off the
  boundary; it never says yes falsely. `capInCone()` closes an arc in a ball
  through the cone point, for petals at boundary vertices.
- **The cache** (`complementcache`): one map of answers per complement
  isomorphism signature (handlebody genus and census name merged), one mutex,
  one set of counters, one entry limit (`complement_cache_limit`), separate
  from the census lookup's mutex so a hit never waits on a lookup.
- **Edge walks** (`edgecycles`): every walk the project makes over a set of
  edges, one helper per kind: directed chaining (with an open-chain policy),
  orienting undirected edges, counting closed curves. Every choice is by input
  position, never by address, so the output is a function of the input alone.

## The census

`census::nameComplement()` names a complement: `Unknot` or
`<n>-component unlink` when proved (above); else the local census, Regina's
`Census::lookup()`, and, on a miss, a Pachner search over retriangulations
(`retriangulate()`, height 2, 8,000 candidates, `retriangulate_time_budget`
seconds) when `retriangulate_on_miss` is on, for knot complements only (a
link's name from its complement bears nothing); else the complement's bare
isomorphism signature. A census hit is translated to a table name where
`census/censusnames.h` knows it (`4_1 (m004 : #1)`).

- **The local census** (`census/census.sqlite`, untracked): an SQLite table
  from isomorphism signature to name, read by a connection per thread, so
  lookups run in parallel. Its path is compiled in (`SURFER_CENSUS_PATH`) and
  set by `cobound`'s `census` key. `insertCensusEntry()` only ever adds
  (`INSERT OR IGNORE`, WAL mode), and only when `census_updates` is on; a
  search inserts its own knot's complement under its table name, and every
  Pachner hit. `insertCounts()` reports writes that failed.
- **Regina's census** is not safe from several threads, so its lookups are
  serialised by one mutex.
- **Census names** are not canonical: one manifold can come back under
  different entry numbers (`#N`) and different census names on different
  runs, so nothing may gate on them; and one complement belongs to
  infinitely many links, so a link's census name bears nothing. `censusnames.h` is generated
  (`gen_knot_census_names`, then `gen_census_names_header.py`) from our own
  PD codes, so its names are KnotInfo's and LinkInfo's, verified against our
  diagrams; a census name is never turned into a table name by pattern
  (`L108019` is `5_1`'s complement, not `8_19`).

## Validation of the isometry step

`isometry_validation` runs on the full tables (17,153 entries; minutes on 10
threads, not in ctest):
- **Census.** The kernel's hyperbolicity and volume agree with KnotInfo and
  LinkInfo for every entry: 16,981 hyperbolic (volumes equal to 10⁻⁶), 172
  not (torus knots, satellites, non-hyperbolic links).
- **Self-identification.** Every entry under every transform (as written,
  mirrored, reversed, and each component c ≥ 1 reversed: 58,679 cases),
  scrambled by seeded random Reidemeister moves into another diagram of the
  same oriented link, is named by the isometry step alone exactly as the
  full namer names the unscrambled diagram: all 57,841 hyperbolic cases. The
  838 non-hyperbolic ones stay untabulated, left to the searches. The
  search-only namer, for comparison, names 17,673 exactly and leaves 36,567
  untabulated, and is never wrong.
- **Negatives.** None of the 4,786 pairs of distinct table links sharing a
  HOMFLY polynomial (either mirror) or a volume is called the same link,
  including the 279 pairs whose complements are isometric (`L10n1` with
  `L7a2`, `L9n7`, `L11n16`).
- **Positives.** All 2,648 pairs of variants of one hyperbolic link are the
  same link; 526 of them are one oriented link up to mirror and global
  reversal, and get one class name.
- **Threads.** Redone on one thread with a fresh namer: identical.

`cobound name` over the atlas's whole database (266,725 cobordisms) names
every outgoing link, with no failures.

## Tests

`ctest` in `build/utils/surfer/linknaming`. `tables_test`'s whole-table round
trip needs the atlas's tables (`SURFER_TEST_ATLAS_DATA`).

| test | pins |
|---|---|
| `names_test` | the grammar read back: component counts, bases and tags, census suffixes, splits, composites, sums, knot marks, descriptions |
| `tables_test` | the one table reader and 4-genus parser (malformed values refused), the symmetry enum and its reader (every KnotInfo spelling), `linkFromTablePD()` on both tables' spellings, the anchors (`isElementarySlice()`: inverse pairs by symmetry type and marks, `8_17#m8_17` and unknown types refused), and every table PD code's round trip through `diagramtriangulation`'s parser and formatter |
| `linknamer_test` | every table entry drawn as itself named as itself; reversing one component of `L7n1{0}` (g₄ 2) gives `L7n1{1}` (g₄ 0); granny and square told apart, a sum and its global reverse one name; splits, split unknots, unlinks and sums in the grammar, a name only when an identity; the visible-sum cut; shared table caches name alike; `nameDrawing()`: a table knot by its signature, a sum by its pieces, a link by its pieces, a drawing from memory, an untabulated knot and a failed knot drawing by the complement, a failed link drawing described `complement:` unless an unlink, names and descriptions perturbed alike |
| `isometry_test` | `sameLinkAs()` finds `8_16`, `L8a1`, `10_151` and `L11n281` from stored outgoing drawings; rejects `8_16` against `10_156` (a HOMFLY twin), and three links whose complement is `L11n353`'s (an isometry exists, none carries meridians to meridians); a torus knot is not hyperbolic and matches nothing |
| `diagramiso_test` | the component map returned is exactly the relabelling applied, over random relabellings; reversing one component is refused whenever it changes the link; mirror and global reversal only when allowed, and flagged; `splitPieces()` keeps every origin; a required map is realised exactly when it is a symmetry (the Hopf link's swap is, the Whitehead diagram's is not) |
| `simplification_test` | linking matrices; `simplify()` keeps origins and every pairwise linking number over hundreds of random scrambles; `removeNugatoryCrossings()` removes every nugatory crossing and nothing else, keeping the link, the linking numbers and the component order |
| `edgecycles_test` | each walk on closed curves, two curves, an arc, a branching, a loop edge, a repeated edge, and a real link's edges |
| `linkcomplement_test` | `EdgeComplement::edgeIndices()` and a link's split into components |
| `meridians_test` | drilling keeps a closed meridian on the right cusp; the one in-process slope round trip (a lens space's loop edge, whose complement has no finite vertex) |
| `unlinknaming_test` | `certifiesUnlink()` and `capInCone()` |
| `censusnaming_test` | the complement cache's limit and its behaviour under concurrent clears, the split-unlink fast path, the Pachner-search policy |
| `census_test` | local census lookups on a scratch database (hits translated by `censusnames.h`, misses, a missing file); an insert lands, also after a hit on the same thread; with `census_updates` off neither a direct insert nor a Pachner hit writes |

Not in ctest: `isometry_validation` (above), and `gen_knot_census_names`.
