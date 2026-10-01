# knotbuilder: triangulating link diagrams (knotbuilder.cpp), the geometry of
# its crossing block (blockgeometry.cpp) and drawing curves of its
# triangulation as diagrams (diagramdrawer.cpp). Those three are compiled
# into surfer_common (../CMakeLists.txt) like every other surfer source; this
# file adds knotbuilder's own tool and tests.

add_executable(triangulateknot triangulateknot.cpp)
target_link_libraries(triangulateknot surfer_common)
install(TARGETS triangulateknot DESTINATION "${CMAKE_INSTALL_FULL_BINDIR}" COMPONENT Runtime)

add_subdirectory(tests)
