# diagramtriangulation: triangulating link diagrams (fromdiagram.cpp), the
# geometry of its crossing block (block.cpp), drawing curves of its
# triangulation as diagrams (todiagram.cpp) and the thickening (thickening/).
# Those are compiled into surfer_common (../CMakeLists.txt) like every other
# part's sources; this file adds its own tool and tests.
#
# triangulateknot keeps its build directory, build/utils/surfer/knotbuilder/.

add_executable(triangulateknot triangulatediagram.cpp)
target_link_libraries(triangulateknot surfer_common surfer_mimalloc)
set_target_properties(triangulateknot PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${SURFER_BUILD_DIR}/knotbuilder")
install(TARGETS triangulateknot DESTINATION "${CMAKE_INSTALL_FULL_BINDIR}" COMPONENT Runtime)

add_subdirectory(tests)
