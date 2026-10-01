# cobound: bounding g4 by cobordisms. Its library, which links every other
# part's (see ../CMakeLists.txt for the dependency rule), then the drivers
# and the tests.
add_library(cobound_lib STATIC
  search/search.cpp
  cobordisms/cobordismkey.cpp cobordisms/database.cpp
  solver/literature.cpp solver/solver.cpp
  cobordisms/cobordism.cpp search/incoming.cpp search/preconditions.cpp
  outgoing/outgoinglink.cpp outgoing/outgoingnamer.cpp
  bounds/partition.cpp bounds/partitiongenera.cpp bounds/cobordismgraph.cpp
  bounds/links.cpp bounds/searchcobordisms.cpp bounds/axioms.cpp
  cobordisms/pairsigner.cpp cobordisms/pending.cpp outgoing/fromdatabase.cpp
  cobordisms/appendonly.cpp)
target_link_libraries(cobound_lib PUBLIC ${ENGINE_LIBRARY} surfer_lib
                      diagramtriangulation_lib linknaming_lib linknaming_complement)
surfer_include_view(cobound_lib cobound)
#
# Each driver keeps its build directory: verifyslicegenus and the tools in
# build/utils/surfer/ (hosts.conf's surfer=, remote_run.sh, stage.sh,
# verify.sh and the gates expect them there), cascadesearch in
# build/utils/surfer/cascade/.

set(COBOUND_TOOLS
  verifyslicegenus  driver/verifyslicegenus.cpp
  farsidediagram    driver/draw.cpp
  farsidename       driver/name.cpp
  peripheral_slopes driver/meridians.cpp
)
while(COBOUND_TOOLS)
  list(POP_FRONT COBOUND_TOOLS EXE SRC)
  add_executable(${EXE} ${SRC})
  target_link_libraries(${EXE} cobound_lib surfer_mimalloc)
  set_target_properties(${EXE} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${SURFER_BUILD_DIR}")
  install(TARGETS ${EXE} DESTINATION "${CMAKE_INSTALL_FULL_BINDIR}" COMPONENT Runtime)
endwhile()

add_executable(cascadesearch driver/cascadesearch.cpp)
target_link_libraries(cascadesearch cobound_lib surfer_mimalloc)
set_target_properties(cascadesearch PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${SURFER_BUILD_DIR}/cascade")

add_subdirectory(tests)
