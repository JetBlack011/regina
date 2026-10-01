set(LINKNAMING_TESTS
  names_test
  tables_test
  edgecycles_test
  linkcomplement_test
  censusnaming_test
  unlinknaming_test
  census_test
  meridians_test
  linknamer_test
  simplification_test
)

foreach(TEST ${LINKNAMING_TESTS})
  add_executable(${TEST} ${TEST}.cpp)
  target_link_libraries(${TEST} linknaming_lib diagramtriangulation_lib surfer_mimalloc)
  add_test(NAME ${TEST} COMMAND ${TEST})
  set_tests_properties(${TEST} PROPERTIES LABELS "surfer" TIMEOUT 60)
endforeach()
# The round trip over every table PD code (when the atlas tables are there).
set_tests_properties(tables_test PROPERTIES TIMEOUT 300)

add_executable(isometry_test isometry_test.cpp)
target_link_libraries(isometry_test linknaming_lib diagramtriangulation_lib surfer_mimalloc)
add_test(NAME isometry_test COMMAND isometry_test)
set_tests_properties(isometry_test PROPERTIES LABELS "surfer" TIMEOUT 120)

add_executable(diagramiso_test diagramiso_test.cpp)
target_link_libraries(diagramiso_test linknaming_lib diagramtriangulation_lib surfer_mimalloc)
add_test(NAME diagramiso_test COMMAND diagramiso_test)
set_tests_properties(diagramiso_test PROPERTIES LABELS "surfer" TIMEOUT 300)

# Whole-table validation of the isometry step; minutes on the full tables, so
# not a ctest target. See the file's header for what it checks.
add_executable(isometry_validation isometry_validation.cpp)
target_link_libraries(isometry_validation linknaming_lib diagramtriangulation_lib surfer_mimalloc)
