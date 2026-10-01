set(TESTS
  submanifoldsearch_test
  submanifold_test
  rollbackunionfind_test
  vertexlinks_test
  pairsig_test
  parallelisosig_test
  linkingnumber_test
  inducedsubgraphs_test
  predicate_order_test
  sha1_test
  namecache_test
)

# Every test just links surfer_common (defined in ../../CMakeLists.txt) rather
# than relisting its own subset of implementation files -- see that file's
# comment on surfer_common for why.
foreach(TEST ${TESTS})
  add_executable(${TEST} ${TEST}.cpp)
  target_link_libraries(${TEST} surfer_common surfer_mimalloc)
  add_test(NAME ${TEST} COMMAND ${TEST})
  set_tests_properties(${TEST} PROPERTIES
    LABELS "surfer"
    TIMEOUT 60)
endforeach()

# Deliberately exhaustive -- exercises every connected induced subgraph of
# several triangulations' 2-skeletons (see submanifold_test.cpp)
# rather than a narrow set of hand-picked cases, so it can take on the order
# of a couple of minutes; correctness coverage is the priority here.
set_tests_properties(submanifold_test PROPERTIES TIMEOUT 300)
# Every set of triangles up to a size cap in several small closed
# triangulations, each checked against a from-scratch reference.
set_tests_properties(predicate_order_test PROPERTIES TIMEOUT 300)
# Compares with the old (drilling) route, ~1 s a call, on a sample; about a
# minute. See linkingnumber_test.cpp's --full for the whole comparison.
set_tests_properties(linkingnumber_test PROPERTIES TIMEOUT 300)

# Standalone profiling driver (see profile_driver.cpp) -- not a CTest test.
add_executable(profile_driver profile_driver.cpp)
target_link_libraries(profile_driver surfer_common surfer_mimalloc)
target_compile_options(profile_driver PRIVATE -g)

# Measures where the boundary-identification drain spends its time; see
# profile_boundarylinks.cpp. Not a CTest test.
add_executable(profile_boundarylinks profile_boundarylinks.cpp)
target_link_libraries(profile_boundarylinks surfer_common surfer_mimalloc)
target_compile_options(profile_boundarylinks PRIVATE -g)

# Decomposes pairSig<4,2>() into its ambient-only and surface-dependent
# parts, which is what justifies PairSigContext; see profile_pairsig.cpp.
# Not a CTest test.
add_executable(profile_pairsig profile_pairsig.cpp)
target_link_libraries(profile_pairsig surfer_common surfer_mimalloc)
target_compile_options(profile_pairsig PRIVATE -g)

# How parallelIsoSigDetail() scales with threads on a row's thickening.
# Not a CTest test.
add_executable(profile_parallelisosig profile_parallelisosig.cpp)
target_link_libraries(profile_parallelisosig surfer_common surfer_mimalloc)

# Replays the pair signatures recorded in results/cobordisms.csv (all written
# by the pre-context code path) through PairSigContext, requiring
# byte-identical output. Needs a corpus file as input, so not a CTest test:
#
#   replay_pairsigs group0.txt group1.txt ...
add_executable(replay_pairsigs replay_pairsigs.cpp)
target_link_libraries(replay_pairsigs surfer_common surfer_mimalloc)

# Review utility: counts triangles with a repeated vertex in each table row's
# thickening, built as verifyslicegenus builds it. Not a CTest test:
#
#   review_repeated_vertices knots.csv links.csv
add_executable(review_repeated_vertices review_repeated_vertices.cpp)
target_link_libraries(review_repeated_vertices surfer_common surfer_mimalloc)
