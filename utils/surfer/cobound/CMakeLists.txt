# cobound: bounding g4 by cobordisms. Its library, which links every other
# part's (see ../CMakeLists.txt for the dependency rule), then the one
# program and the tests.
add_library(cobound_lib STATIC
  search/search.cpp search/searchreport.cpp
  cobordisms/cobordismkey.cpp cobordisms/database.cpp
  solver/literature.cpp solver/solver.cpp
  cobordisms/cobordism.cpp search/incoming.cpp search/preconditions.cpp
  outgoing/outgoinglink.cpp outgoing/outgoingnamer.cpp
  bounds/partition.cpp bounds/partitiongenera.cpp bounds/cobordismgraph.cpp
  bounds/links.cpp bounds/searchcobordisms.cpp bounds/axioms.cpp bounds/searchjudge.cpp
  cobordisms/pairsigner.cpp cobordisms/pending.cpp outgoing/fromdatabase.cpp
  cobordisms/appendonly.cpp driver/signals.cpp
  solver/verdicts.cpp solver/solverinputs.cpp
  bounds/certificate.cpp bounds/databasecobordisms.cpp
  driver/config.cpp driver/setup.cpp driver/fatal.cpp driver/targets.cpp
  driver/runrecords.cpp driver/scheduler.cpp
  driver/run.cpp driver/solve.cpp driver/sign.cpp driver/draw.cpp driver/name.cpp
  driver/meridians.cpp)
target_link_libraries(cobound_lib PUBLIC ${ENGINE_LIBRARY} surfer_lib
                      diagramtriangulation_lib linknaming_lib linknaming_complement)
surfer_include_view(cobound_lib cobound)

# The one program: cobound run | solve | sign | draw | name | meridians. It
# goes to build/utils/surfer/, where scripts and campaigns expect the
# project's programs.
add_executable(cobound main.cpp)
target_link_libraries(cobound cobound_lib surfer_mimalloc)
set_target_properties(cobound PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${SURFER_BUILD_DIR}")
install(TARGETS cobound DESTINATION "${CMAKE_INSTALL_FULL_BINDIR}" COMPONENT Runtime)

add_subdirectory(tests)
