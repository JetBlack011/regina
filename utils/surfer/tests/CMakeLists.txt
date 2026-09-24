set(TESTS
  embeddingsearch_test
  cobordismbuilder_test
  knotbuilder_test
  embeddedsubmanifold_test
  rollbackunionfind_test
  vertexlinks_test
  pairsig_test
  linkcomplement_test
  identifycomplement_test
  census_test
  cobordismgraph_test
  peripheral_test
  witnesskey_test
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

# Slow, one-off validation utility: sweeps knotbuilder against every knot in
# pd_codes.csv (every knot up to 13 crossings -- 12,467 of them), checking
# that each produces a valid, closed S³ triangulation. Takes ~13 minutes, so
# it is *not* registered as a CTest test (it would blow out routine `ctest`
# runs); run it manually instead:
#
#   ./triangulate_knot_table pd_codes.csv
add_executable(triangulate_knot_table triangulate_knot_table.cpp)
target_link_libraries(triangulate_knot_table surfer_common)

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
