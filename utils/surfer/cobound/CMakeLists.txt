# cobound: bounding g4 by cobordisms. Its sources are compiled into
# surfer_common (../CMakeLists.txt) like every other part's; this file adds
# the drivers and the tests.
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
  target_link_libraries(${EXE} surfer_common surfer_mimalloc)
  set_target_properties(${EXE} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${SURFER_BUILD_DIR}")
  install(TARGETS ${EXE} DESTINATION "${CMAKE_INSTALL_FULL_BINDIR}" COMPONENT Runtime)
endwhile()

add_executable(cascadesearch driver/cascadesearch.cpp)
target_link_libraries(cascadesearch surfer_common surfer_mimalloc)
set_target_properties(cascadesearch PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${SURFER_BUILD_DIR}/cascade")

add_subdirectory(tests)
