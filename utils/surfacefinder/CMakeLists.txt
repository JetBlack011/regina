set(UTILS
  surfacefinder
  triangulateknot
)

set(TESTS
  gluinggraph_test
  cobordismbuilder_test 
)

set(surfacefinder_EXTRA_SOURCES        knotbuilder.cpp)
set(triangulateknot_EXTRA_SOURCES      knotbuilder.cpp)
set(cobordismbuilder_test_EXTRA_SOURCES knotbuilder.cpp)

foreach(EXE ${UTILS} ${TESTS})
  add_executable(${EXE} ${EXE}.cpp ${${EXE}_EXTRA_SOURCES})
  target_link_libraries(${EXE} ${ENGINE_LIBRARY})
endforeach()

install(TARGETS ${UTILS} DESTINATION "${CMAKE_INSTALL_FULL_BINDIR}" COMPONENT Runtime)
