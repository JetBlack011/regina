add_executable(exactnaming_test exactnaming_test.cpp)
target_link_libraries(exactnaming_test surfer_common)
add_test(NAME exactnaming_test COMMAND exactnaming_test)
set_tests_properties(exactnaming_test PROPERTIES LABELS "surfer" TIMEOUT 60)

add_executable(snappeaisometry_test snappeaisometry_test.cpp)
target_link_libraries(snappeaisometry_test surfer_common)
add_test(NAME snappeaisometry_test COMMAND snappeaisometry_test)
set_tests_properties(snappeaisometry_test PROPERTIES LABELS "surfer" TIMEOUT 120)

# Whole-table validation of the isometry step; minutes on the full tables, so
# not a ctest target. See the file's header for what it checks.
add_executable(isometry_validation isometry_validation.cpp)
target_link_libraries(isometry_validation surfer_common)
