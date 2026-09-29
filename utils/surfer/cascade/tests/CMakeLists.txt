# Seconds each. See ../README.md, "Assumptions and their tests", for what
# each pins.
foreach(TEST profile_test proofgraph_test diagramiso_test nodes_test hopedges_test
             hoprunner_test leaves_test keptstore_test)
  add_executable(cascade_${TEST} ${TEST}.cpp)
  target_link_libraries(cascade_${TEST} surfer_common)
  add_test(NAME cascade_${TEST} COMMAND cascade_${TEST})
  set_tests_properties(cascade_${TEST} PROPERTIES LABELS "surfer" TIMEOUT 300)
endforeach()
# Real witnesses from the Phase 0 hops, and small tables to name far sides by
# (data/), for the end-to-end tests.
foreach(TEST hopedges_test hoprunner_test keptstore_test)
  target_compile_definitions(cascade_${TEST} PRIVATE
    CASCADE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data")
endforeach()
