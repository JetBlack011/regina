set(DIAGRAMTRIANGULATION_TESTS
  pdcode_test
  fromdiagram_test
  todiagram_test
  thickening_test
)

foreach(TEST ${DIAGRAMTRIANGULATION_TESTS})
  add_executable(${TEST} ${TEST}.cpp)
  target_link_libraries(${TEST} diagramtriangulation_lib surfer_mimalloc)
  add_test(NAME ${TEST} COMMAND ${TEST})
  set_tests_properties(${TEST} PROPERTIES
    LABELS "surfer"
    TIMEOUT 60)
endforeach()

# Slow, one-off validation utility: sweeps fromdiagram against every row of
# a table -- the atlas's knot table, data/4d_smooth_slice_genus_13_crossings_
# pd_codes.csv, holds every knot through 13 crossings -- checking that each
# produces a valid, closed S³ triangulation. Takes minutes, so it is *not*
# registered as a CTest test (it would blow out routine `ctest` runs); run
# it manually instead:
#
#   ./triangulate_knot_table <table.csv> [limit]
#
# todiagram_test sweeps a table the same way when given one (every row
# must be redrawn exactly; ~8 s for all 12,965 knots):
#
#   ./todiagram_test <table.csv> [max crossings]
add_executable(triangulate_knot_table triangulate_knot_table.cpp)
target_link_libraries(triangulate_knot_table diagramtriangulation_lib surfer_mimalloc)
