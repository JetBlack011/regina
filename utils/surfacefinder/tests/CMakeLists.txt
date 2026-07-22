set(TESTS
  embeddingsearch_test
  cobordismbuilder_test
  knotbuilder_test
  embeddedsubmanifold_test
)
set(knotbuilder_test_EXTRA_SOURCES ../knotbuilder.cpp ../cobordismbuilder.cpp ../simplicialprism.cpp ../linkcomplement.cpp)
set(cobordismbuilder_test_EXTRA_SOURCES ../cobordismbuilder.cpp ../simplicialprism.cpp)
set(embeddingsearch_test_EXTRA_SOURCES ../cobordismbuilder.cpp ../simplicialprism.cpp ../skeleton.cpp ../embeddedsubmanifold.cpp ../embeddingsearch.cpp ../linkcomplement.cpp ../enumerate_cis.cpp)
set(embeddedsubmanifold_test_EXTRA_SOURCES ../cobordismbuilder.cpp ../simplicialprism.cpp ../skeleton.cpp ../embeddedsubmanifold.cpp ../linkcomplement.cpp ../enumerate_cis.cpp ../knotbuilder.cpp)

foreach(TEST ${TESTS})
  add_executable(${TEST} ${TEST}.cpp ${${TEST}_EXTRA_SOURCES})
  target_link_libraries(${TEST} ${ENGINE_LIBRARY})
  add_test(NAME ${TEST} COMMAND ${TEST})
  set_tests_properties(${TEST} PROPERTIES
    LABELS "surfacefinder"
    TIMEOUT 60)
endforeach()

# Deliberately exhaustive -- exercises every connected induced subgraph of
# several triangulations' 2-skeletons (see embeddedsubmanifold_test.cpp)
# rather than a narrow set of hand-picked cases, so it can take on the order
# of a couple of minutes; correctness coverage is the priority here.
set_tests_properties(embeddedsubmanifold_test PROPERTIES TIMEOUT 300)

# Slow, one-off validation utility: sweeps knotbuilder against every knot in
# pd_codes.csv (every knot up to 13 crossings -- 12,467 of them), checking
# that each produces a valid, closed S³ triangulation. Takes ~13 minutes, so
# it is *not* registered as a CTest test (it would blow out routine `ctest`
# runs); run it manually instead:
#
#   ./triangulate_knot_table pd_codes.csv
add_executable(triangulate_knot_table triangulate_knot_table.cpp ../knotbuilder.cpp ../linkcomplement.cpp)
target_link_libraries(triangulate_knot_table ${ENGINE_LIBRARY})
