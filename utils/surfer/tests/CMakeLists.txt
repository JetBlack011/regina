set(TESTS
  embeddingsearch_test
  cobordismbuilder_test
  embeddedsubmanifold_test
  rollbackunionfind_test
  vertexlinks_test
  pairsig_test
  linkcomplement_test
  linkingnumber_test
  identifycomplement_test
  census_test
  cobordismgraph_test
  peripheral_test
  witnesskey_test
  rowmap_test
  enumerator_test
  farsidenaming_test
  predicate_order_test
)

# Every test just links surfer_common (defined in ../CMakeLists.txt) rather
# than relisting its own subset of implementation files -- see that file's
# comment on surfer_common for why.
foreach(TEST ${TESTS})
  add_executable(${TEST} ${TEST}.cpp)
  target_link_libraries(${TEST} surfer_common)
  add_test(NAME ${TEST} COMMAND ${TEST})
  set_tests_properties(${TEST} PROPERTIES
    LABELS "surfer"
    TIMEOUT 60)
endforeach()

# Deliberately exhaustive -- exercises every connected induced subgraph of
# several triangulations' 2-skeletons (see embeddedsubmanifold_test.cpp)
# rather than a narrow set of hand-picked cases, so it can take on the order
# of a couple of minutes; correctness coverage is the priority here.
set_tests_properties(embeddedsubmanifold_test PROPERTIES TIMEOUT 300)
# Every set of triangles up to a size cap in several small closed
# triangulations, each checked against a from-scratch reference.
set_tests_properties(predicate_order_test PROPERTIES TIMEOUT 300)
# Compares with the old (drilling) route, ~1 s a call, on a sample; about a
# minute. See linkingnumber_test.cpp's --full for the whole comparison.
set_tests_properties(linkingnumber_test PROPERTIES TIMEOUT 300)


# One-off generation utility -- not a CTest test: builds the Rolfsen <->
# Regina-census-name correspondence data behind ../knot_census_names.csv.
# Run it manually (see gen_knot_census_names.cpp for usage):
#
#   ./gen_knot_census_names pd_codes.csv > ../knot_census_names.csv
add_executable(gen_knot_census_names gen_knot_census_names.cpp)
target_link_libraries(gen_knot_census_names surfer_common)

# Standalone profiling driver (see profile_driver.cpp) -- not a CTest test.
add_executable(profile_driver profile_driver.cpp)
target_link_libraries(profile_driver surfer_common)
target_compile_options(profile_driver PRIVATE -g)

# Measures where the boundary-identification drain spends its time; see
# profile_boundarylinks.cpp. Not a CTest test.
add_executable(profile_boundarylinks profile_boundarylinks.cpp)
target_link_libraries(profile_boundarylinks surfer_common)
target_compile_options(profile_boundarylinks PRIVATE -g)

# Decomposes pairSig<4,2>() into its ambient-only and surface-dependent
# parts, which is what justifies PairSigContext; see profile_pairsig.cpp.
# Not a CTest test.
add_executable(profile_pairsig profile_pairsig.cpp)
target_link_libraries(profile_pairsig surfer_common)
target_compile_options(profile_pairsig PRIVATE -g)

# Replays the pair signatures recorded in results/cobordisms.csv (all written
# by the pre-context code path) through PairSigContext, requiring
# byte-identical output. Needs a corpus file as input, so not a CTest test:
#
#   replay_pairsigs group0.txt group1.txt ...
add_executable(replay_pairsigs replay_pairsigs.cpp)
target_link_libraries(replay_pairsigs surfer_common)

# Review utility: counts triangles with a repeated vertex in each table row's
# thickening, built as verifyslicegenus builds it. Not a CTest test:
#
#   review_repeated_vertices knots.csv links.csv
add_executable(review_repeated_vertices review_repeated_vertices.cpp)
target_link_libraries(review_repeated_vertices surfer_common)

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

# --surface-log over two rows on one thread must not reuse the previous
# row's (destroyed) CsvWriter shard. See the script's header.
add_test(NAME surface_log_rows_test
         COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/surface_log_rows_test.sh
                 $<TARGET_FILE:verifyslicegenus>)
set_tests_properties(surface_log_rows_test PROPERTIES
  LABELS "surfer" TIMEOUT 300)
