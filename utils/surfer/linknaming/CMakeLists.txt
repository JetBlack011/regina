# linknaming: naming a link, with proof (tables, the namer, the SnapPea
# isometry step, complements, the census, diagrams). Its sources are compiled
# into surfer_common (../CMakeLists.txt) like every other part's; this file
# adds its tools and tests.

# The table's link classes. Keeps its build directory, build/utils/surfer/.
add_executable(tableclasses tableclasses.cpp)
target_link_libraries(tableclasses surfer_common surfer_mimalloc)
set_target_properties(tableclasses PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${SURFER_BUILD_DIR}")
install(TARGETS tableclasses DESTINATION "${CMAKE_INSTALL_FULL_BINDIR}" COMPONENT Runtime)

# One-off generation utility: builds the census-name <-> classical-name
# correspondence data behind knot_census_names.csv and link_census_names.csv
# (untracked, beside it), which gen_census_names_header.py turns into
# census/censusnames.h. Run it manually (see gen_census_names.cpp for usage):
#
#   ./gen_knot_census_names pd_codes.csv > knot_census_names.csv
add_executable(gen_knot_census_names gen_census_names.cpp)
target_link_libraries(gen_knot_census_names surfer_common surfer_mimalloc)

add_subdirectory(tests)
