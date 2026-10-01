# linknaming: naming a link, with proof (tables, the namer, the SnapPea
# isometry step, complements, the census, diagrams). Two libraries (see
# ../CMakeLists.txt for the dependency rule), then its tools and tests.

# The bottom layer: links as subcomplexes of a triangulated 3-manifold's
# 1-skeleton -- drilling, meridians, the unknot and unlinks, and the one
# cache of complement answers. Regina only.
add_library(linknaming_complement STATIC
  complement/edgecycles.cpp complement/linkcomplement.cpp complement/meridians.cpp
  complement/complementcache.cpp complement/unlinknaming.cpp)
target_link_libraries(linknaming_complement PUBLIC ${ENGINE_LIBRARY})
surfer_include_view(linknaming_complement linknaming/complement)

# The rest of naming: names, tables, the namer, the isometry step, diagrams
# and the census route.
add_library(linknaming_lib STATIC
  census/censusnaming.cpp
  tables.cpp diagrams/gaussdiagram.cpp linknamer.cpp names.cpp
  isometry/isometry.cpp
  isometry/snappeakernel_isometry.cpp isometry/snappeakernel_isometry_cusped.cpp
  diagrams/diagramiso.cpp diagrams/simplification.cpp)
target_link_libraries(linknaming_lib PUBLIC ${ENGINE_LIBRARY} linknaming_complement
                      SQLite::SQLite3)
surfer_include_view(linknaming_lib linknaming)

# Default path census::localCensusLookup() opens, overridable at runtime
# via --census-db. Untracked: extended live by verifyslicegenus's
# census::insertCensusEntry(), and moved by hand with the tree; not part of
# the default build.
target_compile_definitions(linknaming_lib PUBLIC
  SURFER_CENSUS_PATH="${CMAKE_CURRENT_SOURCE_DIR}/census/census.sqlite")

# The SnapPea kernel's isometry sources, which Regina ships unbuilt
# (engine/snappea/kernel/unused/), include their kernel headers by bare name.
# Only these two files see the kernel directory: its headers (homology.h,
# tables.h, ...) would shadow others anywhere else.
set_source_files_properties(
  isometry/snappeakernel_isometry.cpp isometry/snappeakernel_isometry_cusped.cpp
  PROPERTIES INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/engine/snappea/kernel")

# The table's link classes. Keeps its build directory, build/utils/surfer/.
add_executable(tableclasses tableclasses.cpp)
target_link_libraries(tableclasses linknaming_lib surfer_mimalloc)
set_target_properties(tableclasses PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${SURFER_BUILD_DIR}")
install(TARGETS tableclasses DESTINATION "${CMAKE_INSTALL_FULL_BINDIR}" COMPONENT Runtime)

# One-off generation utility: builds the census-name <-> classical-name
# correspondence data behind knot_census_names.csv and link_census_names.csv
# (untracked, beside it), which gen_census_names_header.py turns into
# census/censusnames.h. Run it manually (see gen_census_names.cpp for usage):
#
#   ./gen_knot_census_names pd_codes.csv > knot_census_names.csv
add_executable(gen_knot_census_names gen_census_names.cpp)
target_link_libraries(gen_knot_census_names linknaming_lib diagramtriangulation_lib
                      surfer_mimalloc)

add_subdirectory(tests)
