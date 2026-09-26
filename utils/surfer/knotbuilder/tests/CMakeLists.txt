set(KNOTBUILDER_TESTS
  knotbuilder_test
  diagramdrawer_test
)

foreach(TEST ${KNOTBUILDER_TESTS})
  add_executable(${TEST} ${TEST}.cpp)
  target_link_libraries(${TEST} surfer_common)
  add_test(NAME ${TEST} COMMAND ${TEST})
  set_tests_properties(${TEST} PROPERTIES
    LABELS "surfer"
    TIMEOUT 60)
endforeach()

# Slow, one-off validation utility: sweeps knotbuilder against every knot in
# pd_codes.csv (every knot up to 13 crossings -- 12,467 of them), checking
# that each produces a valid, closed S³ triangulation. Takes ~13 minutes, so
# it is *not* registered as a CTest test (it would blow out routine `ctest`
# runs); run it manually instead:
#
#   ./triangulate_knot_table pd_codes.csv
#
# diagramdrawer_test sweeps a table the same way when given one (every row
# must be redrawn exactly; ~8 s for all 12,965 knots):
#
#   ./diagramdrawer_test pd_codes.csv [max crossings]
add_executable(triangulate_knot_table triangulate_knot_table.cpp)
target_link_libraries(triangulate_knot_table surfer_common)
