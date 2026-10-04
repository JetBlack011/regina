# utils/surfer

Enumerating PL surfaces in triangulated 4-manifolds, and using them to bound
the smooth 4-genus of knots and links by cobordisms in S³ × I. The method is
written up in the paper (`pl_enumeration_draft`); the results live in the
separate `cobordism-atlas` repository.

## The four parts

| part | what | README |
|---|---|---|
| `surfer/` | a library: surfaces made of the triangles of a triangulated 4-manifold, enumerated with local checks at vertex links; frontiers; pair signatures. Knows nothing of knots. A proof-of-concept driver, `surfer`. | `surfer/README.md` |
| `diagramtriangulation/` | a library: a link diagram's triangulation T of S³ with the link in its 1-skeleton, curves of T drawn back as oriented diagrams, and the thickening T × [0, k] a search runs in. | `diagramtriangulation/README.md` |
| `linknaming/` | a library: naming a link with proof (tables, diagram matches, the SnapPea kernel's isometry test, Reidemeister searches) or only describing it; complements, meridians and the census. Knows nothing of surfaces. | `linknaming/README.md` |
| `cobound/` | the application: one program, `cobound`, that searches for cobordisms, keeps the database of them, and bounds g₄ from them. | `cobound/README.md` |

`apidocs/` builds the Doxygen pages when Regina's engine docs are on
(`REGINA_BUILD_ENGINE_DOCS`).

## The dependency rule

Dependencies run one way, and naming never depends on the search:

```
linknaming_complement     -> Regina                     (linknaming/complement/: the bottom layer)
surfer_lib                -> Regina, OpenSSL libcrypto, linknaming_complement
diagramtriangulation_lib  -> linknaming_complement
linknaming_lib            -> Regina, SQLite, linknaming_complement
cobound_lib               -> all of the above            (the only place search meets naming)
```

Executables link what they need; surfer's tests and driver, and linknaming's
tests and census-name generator, also use `diagramtriangulation_lib`.

The build enforces the rule. Every project include is written relative to
this directory (`#include "surfer/enumeration/surfacesearch.h"`), but no
library sees this directory: each sees an include view
(`surfer_include_view()` in `CMakeLists.txt`), a directory in the build tree
holding symbolic links to its own part and nothing else. A source that
includes a header of a part its library does not link fails to compile.

## Building

The parts build with Regina, from its build tree: `utils/CMakeLists.txt`
adds this directory, and nothing outside it changes, so a rebuild compiles
project code only.

```
cmake --build build -j16
```

Besides Regina they need OpenSSL's libcrypto and SQLite 3 (development
headers), and they link a static mimalloc when CMake finds one
(`~/.local/lib/libmimalloc.a`, or `-DSURFER_MIMALLOC=<path>`; see
`surfer/README.md`, "Performance"). The programs land where scripts and
campaigns expect them:

| program | path under `build/utils/surfer/` |
|---|---|
| `cobound` | `cobound` |
| `tableclasses` | `tableclasses` |
| `triangulateknot` | `knotbuilder/triangulateknot` |
| `surfer` | `surfer/surfer` |
| `gen_knot_census_names` | `linknaming/gen_knot_census_names` |

Every program kept its name and path when the code was split into its four
parts, so commands written against the old layout still work; `surfer` alone
could not, since its part's build directory took the name. That is why
`triangulateknot` still builds into `knotbuilder/`, the directory of the part
`diagramtriangulation/` was before.

The local census (`linknaming/census/census.sqlite`) is not tracked; its
default path is compiled in (`SURFER_CENSUS_PATH`), and `cobound` takes
another with the `census` key.

## Testing

Each part's tests run under `ctest` from its own build directory (61 in
all):

| part | directory under `build/utils/surfer/` | tests |
|---|---|---|
| surfer | `surfer` | 13 |
| diagramtriangulation | `diagramtriangulation` | 4 |
| linknaming | `linknaming` | 12 |
| cobound | `cobound_part` | 32 |

Running `ctest` in `build/utils/surfer` runs all of them. Each takes under a
minute on an idle 8-core machine, but a few are exhaustive and can run past
their timeouts on a loaded one: run them idle. The atlas's full tables are found through `SURFER_TEST_ATLAS_DATA` (or
the usual checkouts); `goal_layout_test` skips without them, and
`tables_test` leaves out its whole-table round trip. Each part's README has
a table of what its tests pin.
