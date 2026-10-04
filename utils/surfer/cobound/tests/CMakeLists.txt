set(TESTS
  solver_test
  cobordismkey_test
  database_test
  incomingmap_test
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

# The cobordism graph's and goal runs' tests: seconds each. See ../README.md, "Assumptions and their
# tests", for what each pins.
foreach(TEST partitiongenera_test cobordismgraph_test links_test searchcobordisms_test
             search_test axioms_test pending_test fromdatabase_test thickening_rows_test
             json_test parallelfor_test appendonly_test timers_test config_test)
  add_executable(${TEST} ${TEST}.cpp)
  target_link_libraries(${TEST} cobound_lib surfer_mimalloc)
  add_test(NAME ${TEST} COMMAND ${TEST})
  set_tests_properties(${TEST} PROPERTIES LABELS "surfer" TIMEOUT 300)
endforeach()
# Real cobordisms from the Phase 0 searches, and small tables to name outgoing
# links by (data/), for the end-to-end tests.
foreach(TEST searchcobordisms_test search_test pending_test fromdatabase_test
             thickening_rows_test)
  target_compile_definitions(${TEST} PRIVATE
    COBOUND_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data")
endforeach()

# No gate may depend on a name: the same exhaustive search, with
# every name perturbed, must accept exactly the same surfaces. See the
# script's header.
add_test(NAME name_independence_test
         COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/name_independence_test.sh
                 $<TARGET_FILE:cobound>)
set_tests_properties(name_independence_test PROPERTIES
  LABELS "surfer" TIMEOUT 600)

# SIGINT and SIGTERM end the running search cleanly and the run searches no
# further, at depth 0 and with a goal; an interrupted search is never
# recorded as having run out of candidates (plan divergence 8). See the
# script's header.
add_test(NAME interrupted_outcome_test
         COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/interrupted_outcome_test.sh
                 $<TARGET_FILE:cobound> ${CMAKE_CURRENT_SOURCE_DIR}/data)
set_tests_properties(interrupted_outcome_test PROPERTIES
  LABELS "surfer" TIMEOUT 600 SKIP_RETURN_CODE 77)

# A goal run must explore the same links whatever the heap layout (phase
# 3.0). Needs the atlas's tables; skipped without them. See the script's
# header.
add_test(NAME goal_layout_test
         COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/goal_layout_test.sh
                 $<TARGET_FILE:cobound>)
set_tests_properties(goal_layout_test PROPERTIES
  LABELS "surfer" TIMEOUT 600 SKIP_RETURN_CODE 77)

# cobound meridians answers each input line on its own: the same records in
# the reverse order give the same output (phase 3.0). See the script's header.
add_test(NAME meridians_order_test
         COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/meridians_order_test.sh
                 $<TARGET_FILE:cobound> ${CMAKE_CURRENT_SOURCE_DIR}/data)
set_tests_properties(meridians_order_test PROPERTIES
  LABELS "surfer" TIMEOUT 300)

# One set of search defaults (2 layers, proper), no default for
# resolve_unlinked or work, the retired options refused, and a run's written
# configuration running the same search (plan divergence 3, phase 5). See the
# script's header.
add_test(NAME search_defaults_test
         COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/search_defaults_test.sh
                 $<TARGET_FILE:cobound> ${CMAKE_CURRENT_SOURCE_DIR}/data)
set_tests_properties(search_defaults_test PROPERTIES
  LABELS "surfer" TIMEOUT 300)

# An imbalanced search ends only itself: outcome unaccounted, the run goes
# on, and the exit code reports what the run was for (plan divergence 2).
# See the script's header.
add_test(NAME unaccounted_search_test
         COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/unaccounted_search_test.sh
                 $<TARGET_FILE:cobound> ${CMAKE_CURRENT_SOURCE_DIR}/data)
set_tests_properties(unaccounted_search_test PROPERTIES
  LABELS "surfer" TIMEOUT 300)

# Targets are searched on their given diagrams, certified in-binary: a table
# row that does not certify is refused, an untabulated goal target falls back
# to its simplified diagram, logged (plan, phase 7.2). See the script's header.
add_test(NAME given_diagram_test
         COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/given_diagram_test.sh
                 $<TARGET_FILE:cobound> ${CMAKE_CURRENT_SOURCE_DIR}/data)
set_tests_properties(given_diagram_test PROPERTIES
  LABELS "surfer" TIMEOUT 300)

# A failed output write ends its search as an I/O error, never the process:
# no frontier, no exhaustion claim, what it kept signed; exit 2 (plan, phase
# 7.3). See the script's header.
add_test(NAME io_failure_test
         COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/io_failure_test.sh
                 $<TARGET_FILE:cobound> ${CMAKE_CURRENT_SOURCE_DIR}/data)
set_tests_properties(io_failure_test PROPERTIES
  LABELS "surfer" TIMEOUT 300)

# A contradiction halts every run, with its code: 2 at depth 0, 3 with a goal
# (plan divergence 6). See the script's header.
add_test(NAME contradiction_halt_test
         COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/contradiction_halt_test.sh
                 $<TARGET_FILE:cobound> ${CMAKE_CURRENT_SOURCE_DIR}/data)
set_tests_properties(contradiction_halt_test PROPERTIES
  LABELS "surfer" TIMEOUT 300)

# A frontier records its pending file and that file's fsynced length; a
# resume refuses, with a named reason, unless `sign` signed that far or the
# file is the run's own (plan divergence 1); a moved work tree keeps its
# resume (phase 5). See the script's header.
add_test(NAME frontier_pending_test
         COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/frontier_pending_test.sh
                 $<TARGET_FILE:cobound> ${CMAKE_CURRENT_SOURCE_DIR}/data)
set_tests_properties(frontier_pending_test PROPERTIES
  LABELS "surfer" TIMEOUT 300)

# The atlas's canaries (canaries.sh, cascade_canaries.sh) as a cobound test:
# phase 0's accounting numbers on nine exhaustive searches and two goal runs. See
# the script's header.
add_test(NAME canaries_test
         COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/canaries_test.sh
                 $<TARGET_FILE:cobound> ${CMAKE_CURRENT_SOURCE_DIR}/data)
set_tests_properties(canaries_test PROPERTIES
  LABELS "surfer" TIMEOUT 600)

# surface_log over two searches on one thread must not reuse the previous
# search's (destroyed) CsvWriter shard. See the script's header.
add_test(NAME surface_log_searches_test
         COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/surface_log_searches_test.sh
                 $<TARGET_FILE:cobound>)
set_tests_properties(surface_log_searches_test PROPERTIES
  LABELS "surfer" TIMEOUT 300)
