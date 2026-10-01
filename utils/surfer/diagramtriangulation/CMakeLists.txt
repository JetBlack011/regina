# diagramtriangulation: triangulating link diagrams (fromdiagram.cpp), the
# geometry of its crossing block (block.cpp), drawing curves of its
# triangulation as diagrams (todiagram.cpp) and the thickening (thickening/).
# Its library (see ../CMakeLists.txt for the dependency rule), then its own
# tool and tests.
add_library(diagramtriangulation_lib STATIC
  pdcode.cpp fromdiagram.cpp block.cpp todiagram.cpp
  thickening/thickening.cpp thickening/prism.cpp thickening/collar.cpp)
target_link_libraries(diagramtriangulation_lib PUBLIC ${ENGINE_LIBRARY}
                      linknaming_complement)
surfer_include_view(diagramtriangulation_lib diagramtriangulation)

# triangulateknot keeps its build directory, build/utils/surfer/knotbuilder/.
add_executable(triangulateknot triangulatediagram.cpp)
target_link_libraries(triangulateknot diagramtriangulation_lib surfer_mimalloc)
set_target_properties(triangulateknot PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${SURFER_BUILD_DIR}/knotbuilder")
install(TARGETS triangulateknot DESTINATION "${CMAKE_INSTALL_FULL_BINDIR}" COMPONENT Runtime)

add_subdirectory(tests)
