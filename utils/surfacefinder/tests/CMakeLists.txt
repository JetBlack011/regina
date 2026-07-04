# Fast, always-run correctness tests for surfacefinder's components. Each
# test binary covers one program/module:
#   gluinggraph_test        - GluingGraph (surfacefinder's core search)
#   cobordismbuilder_test   - CobordismBuilder (thicken()/cone())
#   knotbuilder_test        - knotbuilder (PD code -> triangulated S³) and
#                             its integration with CobordismBuilder
#
# Registered with CTest, so `ctest` (or `ctest --output-on-failure`) from
# the build directory runs and tracks all of them in one place, alongside
# the rest of Regina's test suite.
set(TESTS
  gluinggraph_test
  cobordismbuilder_test
  knotbuilder_test
)
set(knotbuilder_test_EXTRA_SOURCES ../knotbuilder.cpp)

foreach(TEST ${TESTS})
  add_executable(${TEST} ${TEST}.cpp ${${TEST}_EXTRA_SOURCES})
  target_link_libraries(${TEST} ${ENGINE_LIBRARY})
  add_test(NAME ${TEST} COMMAND ${TEST})
  set_tests_properties(${TEST} PROPERTIES
    LABELS "surfacefinder"
    TIMEOUT 60)
endforeach()

# Slow, one-off validation utility: sweeps knotbuilder against every knot in
# pd_codes.csv (every knot up to 13 crossings -- 12,467 of them), checking
# that each produces a valid, closed S³ triangulation. Takes ~13 minutes, so
# it is *not* registered as a CTest test (it would blow out routine `ctest`
# runs); run it manually instead:
#
#   ./triangulate_knot_table pd_codes.csv
add_executable(triangulate_knot_table triangulate_knot_table.cpp ../knotbuilder.cpp)
target_link_libraries(triangulate_knot_table ${ENGINE_LIBRARY})
