# The project root, in four parts:
#   surfer/                a knot-agnostic library: PL surfaces in triangulated 4-manifolds
#   diagramtriangulation/  the triangulations a search runs in, built from a diagram
#   linknaming/            naming a link, with proof
#   cobound/               bounding g4 by cobordisms
#
# Every part's implementation files are compiled into ONE static library,
# surfer_common, and every project include is written relative to ONE include
# root, this directory (e.g. #include "surfer/enumeration/surfacesearch.h").
# Each part's own CMakeLists.txt adds its executables and tests, which all just
# link surfer_common rather than relisting their own subset of implementation
# files (which silently goes stale): a traditional archive only pulls in the
# .o members a given target's unresolved symbols actually need, so e.g. a test
# that only exercises rollbackunionfind.cpp still ends up linking just that,
# with no per-target bookkeeping required.
add_library(surfer_common STATIC
  diagramtriangulation/fromdiagram.cpp diagramtriangulation/block.cpp diagramtriangulation/todiagram.cpp
  diagramtriangulation/thickening/thickening.cpp diagramtriangulation/thickening/prism.cpp
  surfer/submanifold/skeleton.cpp
  surfer/submanifold/submanifold.cpp surfer/enumeration/submanifoldsearch.cpp
  surfer/enumeration/surfacesearch.cpp
  linknaming/complement/linkcomplement.cpp surfer/submanifold/linkingnumber.cpp
  linknaming/complement/meridians.cpp linknaming/census/censusnaming.cpp
  linknaming/complement/complementcache.cpp linknaming/complement/unlinknaming.cpp
  surfer/enumeration/namecache.cpp
  surfer/enumeration/inducedsubgraphs.cpp diagramtriangulation/thickening/collar.cpp
  cobound/search/search.cpp
  surfer/submanifold/rollbackunionfind.cpp surfer/submanifold/vertexlinks.cpp
  surfer/pairsig/pairsig.cpp surfer/pairsig/sha1.cpp surfer/report/csvwriter.cpp cobound/cobordisms/cobordismkey.cpp
  cobound/cobordisms/database.cpp surfer/enumeration/searchstrategy.cpp
  cobound/solver/literature.cpp cobound/solver/solver.cpp
  cobound/cobordisms/cobordism.cpp cobound/search/incoming.cpp cobound/search/preconditions.cpp
  linknaming/names.cpp cobound/outgoing/outgoinglink.cpp
  cobound/outgoing/outgoingnamer.cpp cobound/outgoing/farsideredraw.cpp
  linknaming/tables.cpp linknaming/diagrams/gaussdiagram.cpp linknaming/linknamer.cpp
  linknaming/isometry/isometry.cpp
  linknaming/isometry/snappeakernel_isometry.cpp linknaming/isometry/snappeakernel_isometry_cusped.cpp
  cobound/bounds/partition.cpp cobound/bounds/partitiongenera.cpp cobound/bounds/cobordismgraph.cpp
  linknaming/diagrams/diagramiso.cpp cobound/bounds/nodes.cpp cobound/bounds/searchcobordisms.cpp
  cobound/bounds/axioms.cpp
  cobound/cobordisms/pairsigner.cpp cobound/cobordisms/pending.cpp cobound/outgoing/fromdatabase.cpp)
target_include_directories(surfer_common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(surfer_common PUBLIC ${ENGINE_LIBRARY})

# Where every part puts the executables that scripts and campaigns expect in
# build/utils/surfer/ (see each part's CMakeLists.txt).
set(SURFER_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}")

# The SnapPea kernel's isometry sources, which Regina ships unbuilt
# (engine/snappea/kernel/unused/), include their kernel headers by bare name.
# Only these two files see the kernel directory: its headers (homology.h,
# tables.h, ...) would shadow others anywhere else.
set_source_files_properties(
  linknaming/isometry/snappeakernel_isometry.cpp linknaming/isometry/snappeakernel_isometry_cusped.cpp
  PROPERTIES INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/engine/snappea/kernel")

# Local census SQLite lookup (censusnaming.cpp's
# census::localCensusLookup()): not used anywhere else in the engine, so
# this dependency is scoped to surfer_common alone rather than added at the
# top level.
find_package(SQLite3 REQUIRED)
target_link_libraries(surfer_common PUBLIC SQLite::SQLite3)

# SHA-1 (surfer/pairsig/sha1.h): OpenSSL's libcrypto, for the frontier
# fingerprint, the pair-signature context's ambient key and the cobordism key.
find_package(OpenSSL REQUIRED COMPONENTS Crypto)
target_link_libraries(surfer_common PUBLIC OpenSSL::Crypto)

# mimalloc in place of glibc's malloc, for every surfer binary and test
# (2026-09-29). A perf profile of a production row put about a third of all
# CPU in glibc's malloc/free, and the drain ran ~1.8x faster under mimalloc
# (utils/surfer/README.md, "Performance"). Linked statically and whole, so the
# executable itself defines malloc and interposes it for the engine library
# too, and a host without the package can use a copied archive
# (-DSURFER_MIMALLOC=<path>, or ~/.local/lib/libmimalloc.a). Results cannot
# depend on the allocator; only speed does. It is one INTERFACE target,
# surfer_mimalloc, which every executable links once: the whole-archive item
# must reach each executable exactly once, whatever libraries it links.
add_library(surfer_mimalloc INTERFACE)
find_library(SURFER_MIMALLOC NAMES libmimalloc.a
             HINTS "$ENV{HOME}/.local/lib" DOC "static mimalloc archive (optional)")
# The archive must define malloc as real machine code. Arch's packaged one is
# GCC LTO bitcode, which a linker without the matching plugin skips with only
# a warning -- leaving glibc's malloc in place. Build mimalloc from source
# instead (see README, "Performance").
if(SURFER_MIMALLOC)
  execute_process(COMMAND nm "${SURFER_MIMALLOC}" OUTPUT_VARIABLE _mi_syms
                  ERROR_QUIET)
  if(NOT _mi_syms MATCHES "[0-9a-f]+ T malloc\n")
    message(WARNING "surfer: ${SURFER_MIMALLOC} defines no malloc as machine "
                    "code (LTO bitcode?); not linking it")
    set(SURFER_MIMALLOC "")
  endif()
endif()
if(SURFER_MIMALLOC)
  message(STATUS "surfer: linking mimalloc from ${SURFER_MIMALLOC}")
  target_link_libraries(surfer_mimalloc INTERFACE
    "$<LINK_LIBRARY:WHOLE_ARCHIVE,${SURFER_MIMALLOC}>")
else()
  message(WARNING "surfer: no libmimalloc.a found; using glibc malloc, which "
                  "costs the drain about 1.8x (set SURFER_MIMALLOC)")
endif()

# Default path census::localCensusLookup() opens, overridable at runtime
# via --census-db. Untracked: extended live by verifyslicegenus's
# census::insertCensusEntry(), and moved by hand with the tree; not part of
# the default build.
target_compile_definitions(surfer_common PUBLIC
  SURFER_CENSUS_PATH="${CMAKE_CURRENT_SOURCE_DIR}/linknaming/census/census.sqlite")

add_subdirectory(surfer)
add_subdirectory(diagramtriangulation)
add_subdirectory(linknaming)
add_subdirectory(cobound)

if (REGINA_BUILD_ENGINE_DOCS)
  add_subdirectory(apidocs)
endif (REGINA_BUILD_ENGINE_DOCS)
