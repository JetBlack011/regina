set(TESTS
  solver_test
  cobordismkey_test
  database_test
  rowmap_test
  preconditions_test
  outgoingnamer_test
)

# Every test links the cobound library (../CMakeLists.txt), and through it
# every part's, rather than relisting its own subset of implementation files
# (which silently goes stale): a static archive only pulls in the members a
# target's unresolved symbols need.
foreach(TEST ${TESTS})
  add_executable(${TEST} ${TEST}.cpp)
  target_link_libraries(${TEST} cobound_lib surfer_mimalloc)
  add_test(NAME ${TEST} COMMAND ${TEST})
  set_tests_properties(${TEST} PROPERTIES
    LABELS "surfer"
    TIMEOUT 60)
endforeach()

# Two exhaustive cap-3 searches (the canaries' 3_1 and L2a1{0}) and a few
# pair signatures; ~25 s on an idle machine.
set_tests_properties(preconditions_test PROPERTIES TIMEOUT 120)

# The cascade's tests: seconds each. See ../README.md, "Assumptions and their
# tests", for what each pins.
foreach(TEST partitiongenera_test cobordismgraph_test links_test searchcobordisms_test
             search_test axioms_test pending_test fromdatabase_test thickening_rows_test
             json_test parallelfor_test)
  add_executable(${TEST} ${TEST}.cpp)
  target_link_libraries(${TEST} cobound_lib surfer_mimalloc)
  add_test(NAME ${TEST} COMMAND ${TEST})
  set_tests_properties(${TEST} PROPERTIES LABELS "surfer" TIMEOUT 300)
endforeach()
# Real witnesses from the Phase 0 hops, and small tables to name far sides by
# (data/), for the end-to-end tests.
foreach(TEST searchcobordisms_test search_test pending_test fromdatabase_test
             thickening_rows_test)
  target_compile_definitions(${TEST} PRIVATE
    CASCADE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data")
endforeach()

# No gate may depend on an identified name: the same exhaustive search, with
# every name perturbed, must accept exactly the same surfaces. See the
# script's header.
add_test(NAME name_independence_test
         COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/name_independence_test.sh
                 $<TARGET_FILE:verifyslicegenus>)
set_tests_properties(name_independence_test PROPERTIES
  LABELS "surfer" TIMEOUT 600)

# A search stopped by SIGINT must not be recorded as having run out of
# candidates. See the script's header.
add_test(NAME interrupted_outcome_test
         COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/interrupted_outcome_test.sh
                 $<TARGET_FILE:verifyslicegenus>)
set_tests_properties(interrupted_outcome_test PROPERTIES
  LABELS "surfer" TIMEOUT 300 SKIP_RETURN_CODE 77)

# A goal run must explore the same links whatever the heap layout (phase
# 3.0). Needs the atlas's tables; skipped without them. See the script's
# header.
add_test(NAME goal_layout_test
         COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/goal_layout_test.sh
                 $<TARGET_FILE:cascadesearch>)
set_tests_properties(goal_layout_test PROPERTIES
  LABELS "surfer" TIMEOUT 600 SKIP_RETURN_CODE 77)

# peripheral_slopes answers each input line on its own: the same records in
# the reverse order give the same output (phase 3.0). See the script's header.
add_test(NAME meridians_order_test
         COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/meridians_order_test.sh
                 $<TARGET_FILE:peripheral_slopes> ${CMAKE_CURRENT_SOURCE_DIR}/data)
set_tests_properties(meridians_order_test PROPERTIES
  LABELS "surfer" TIMEOUT 300)

# --surface-log over two rows on one thread must not reuse the previous
# row's (destroyed) CsvWriter shard. See the script's header.
add_test(NAME surface_log_rows_test
         COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/surface_log_rows_test.sh
                 $<TARGET_FILE:verifyslicegenus>)
set_tests_properties(surface_log_rows_test PROPERTIES
  LABELS "surfer" TIMEOUT 300)
