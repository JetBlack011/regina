# surfacefinder

Research tool for enumerating embedded surfaces in triangulated 4-manifolds,
built within [Regina](https://github.com/regina-normal/regina).
Given a triangulated 4-manifold (built directly from an isomorphism
signature, or thickened up from a knot/link diagram via a PD code),
`surfacefinder` performs a DFS over a "gluing graph" of triangles to find all
embedded 2-surfaces satisfying some boundary condition (closed, properly
bounded, or unconstrained).

## Contents

### Library headers (no standalone build output)

| File | Purpose |
|---|---|
| `knotbuilder.h` / `knotbuilder.cpp` | Builds a closed triangulation of S³ from a PD code, with the knot/link embedded as a specific cycle of edges. One `Block` per crossing, glued wall-to-wall along the diagram's edges. |
| `cobordismbuilder.h` | `CobordismBuilder<dim>`: thickens a `Triangulation<dim>` into a `Triangulation<dim+1>` (`thicken()`/`thicken(n)`), and caps a thickened cobordism into a ball (`cone()`). |
| `simplicialprism.h` | `SimplicialPrism<dim>`: the "staircase" triangulation of `Δ^(dim-1) × I` that `CobordismBuilder::thicken()` builds out of, one prism per base simplex. |
| `gluinggraph.h` | `GluingGraph<dim>` and `GluingNode<dim>`: builds the graph of triangles/gluings that `surfacefinder`'s DFS runs over, and the DFS itself (`findSurfaces()`). |
| `gluing.h` | Small `Gluing<dim, subdim>` helper struct used when deferring/batching face identifications. |
| `knottedsurfaces.h` | `KnottedSurface<dim>` (the surface under construction during DFS, with its invariants), `EdgeComplement`/`Knot`/`Link` (boundary link extraction and complement recognition), and `SurfaceCondition`. |

### Programs

| Program | What it does |
|---|---|
| `surfacefinder` | The actual surface search tool. See [Programs](#programs) below. |
| `triangulateknot` | Diagnostic/demo tool: takes a PD code, prints the resulting S³ triangulation, its knot/link complement, and one thickened layer. |

### Tests

| File | Covers |
|---|---|
| `gluinggraph_test.cpp` | `GluingGraph`/`findSurfaces()` against hand-verified triangulations (single tetrahedron, two tetrahedra, minimal S³, S⁴). |
| `cobordismbuilder_test.cpp` | `CobordismBuilder::thicken()`/`cone()`: multi-layer thickening, capping, ordering, and an exhaustive sweep over `SimplicialPrism::glue()`'s facet-pairing permutations. Independent of knotbuilder. |
| `knotbuilder_test.cpp` | `knotbuilder::buildLink()` end to end: validity/closure/sphere checks, a dedicated regression test for a non-alternating diagram, a sweep over 18 named knots (3 to 13 crossings, alternating and non-alternating) cross-checked against Regina's own `Link::fromPD()`, and the full `knotbuilder` → `CobordismBuilder` → ball pipeline. |
| `triangulate_knot_table.cpp` | Standalone utility (**not** a CTest test — see [Tests](#tests)): sweeps `knotbuilder::buildLink()` over every row of a PD-code CSV file. |
| `pd_codes.csv` | PD codes for every knot up to 13 crossings (12,467 rows), from [KnotInfo](https://knotinfo.org/). Gitignored (`*.csv`); re-download from KnotInfo if missing. |

### Present but not currently built

| File | Status |
|---|---|
| `cobordismbuilder.cpp` | An older standalone CLI (`{dim} <isosig>` → thickens by 2 layers, prints the isoSig) predating the current `CobordismBuilder<dim>` header-only design. Not referenced by any `CMakeLists.txt`; kept for reference. |
| `knot_builder_simplified.py` | An earlier Python prototype of the crossing-block construction (predates `knotbuilder.h`/`.cpp`) developed during TXGL, using Regina's Python bindings. Not part of the C++ build. |
| `bogocheck.cpp` | A generic "is this the minimal triangulation?" search tool (via `Triangulation::retriangulate`), not specific to surfacefinder. Not wired into any `CMakeLists.txt`. |

## Building

All build commands run from the repository's `build/` directory.

```bash
mkdir build
cd build
cmake ..
make -j8
make surfacefinder -j8
```

Resulting binaries:

```
build/utils/surfacefinder/surfacefinder
build/utils/surfacefinder/triangulateknot
build/utils/surfacefinder/tests/gluinggraph_test
build/utils/surfacefinder/tests/cobordismbuilder_test
build/utils/surfacefinder/tests/knotbuilder_test
build/utils/surfacefinder/tests/triangulate_knot_table
```

## Programs

### `surfacefinder`

```
surfacefinder <PD code>
surfacefinder { -a, --all | -b, --boundary | -c, --closed | -l, --links } <isosig>
surfacefinder [ -v, --version | -h, --help ]
```

- **`<PD code>` form**: builds S³ from the PD code via `knotbuilder`, thickens
  it by 2 layers via `CobordismBuilder<3>` (both the surface condition and
  layer count are currently hardcoded to `--boundary` and `2` in this code
  path), then searches for surfaces.
- **`<isosig>` form**: loads a `Triangulation<4>` directly from an
  isomorphism signature, with an explicit condition flag:
  - `-a`/`--all`: any embedded surface without self-intersection.
  - `-b`/`--boundary`: proper surfaces (boundary contained in the
    4-manifold's boundary; equivalent to `--closed` when the manifold is
    closed).
  - `-c`/`--closed`: closed surfaces only.
  - `-l`/`--links`: same as `--boundary`, but additionally attempts to
    recognize the complement of each surface's boundary link.

Examples:

```bash
# A single pentachoron (B^4): small enough to finish quickly.
./build/utils/surfacefinder/surfacefinder --all baa

# Trefoil PD code: builds and thickens S³ automatically, then searches.
# The search space grows fast with crossing number, so this (and PD codes
# in general) can take a while to finish.
./build/utils/surfacefinder/surfacefinder "[[1,4,2,5],[3,6,4,1],[5,2,6,3]]"
```

### `triangulateknot`

```
triangulateknot <PD code>
```

Builds S³ from the PD code, prints its isomorphism signature, builds and
prints the isoSig of the knot/link complement (via `Link::buildComplement()`),
and prints the isoSig after one `CobordismBuilder<3>::thicken()` layer. Useful
for sanity-checking a PD code or inspecting intermediate triangulations
without running a full surface search.

## Tests

### Running via CTest

The three fast test binaries are registered with CTest under the
`surfacefinder` label, alongside the rest of Regina's test suite:

```bash
cd build
ctest -L surfacefinder --output-on-failure
```

### Running individually

```bash
./build/utils/surfacefinder/tests/gluinggraph_test
./build/utils/surfacefinder/tests/cobordismbuilder_test
./build/utils/surfacefinder/tests/knotbuilder_test
```

### Triangulating the full knot-table

`triangulate_knot_table` is a standalone utility, not part of the
routine test suite (running it against the full KnotInfo table takes roughly 10-15 minutes):

```bash
./build/utils/surfacefinder/tests/triangulate_knot_table \
    utils/surfacefinder/tests/pd_codes.csv [optional row limit]
```

It checks that every PD code in the CSV produces a valid, closed S³
triangulation with a well-formed `Link`, and reports a pass/fail summary
plus any individual failures.
