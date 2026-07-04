set(UTILS
  #bogocheck
  #buildknot
  #surfacefinder
  fillmanifold
  )

add_executable(buildknot buildknot.cpp knotbuilder.cpp)
target_link_libraries(buildknot ${ENGINE_LIBRARY})

add_executable(surfacefinder surfacefinder.cpp knotbuilder.cpp)
target_link_libraries(surfacefinder ${ENGINE_LIBRARY})

add_executable(gluinggraph_test gluinggraph_test.cpp)
target_link_libraries(gluinggraph_test ${ENGINE_LIBRARY})

add_executable(knotbuilder_test knotbuilder_test.cpp knotbuilder.cpp)
target_link_libraries(knotbuilder_test ${ENGINE_LIBRARY})

foreach(UTIL ${UTILS})
  add_executable(${UTIL} ${UTIL}.cpp)
  target_link_libraries(${UTIL} ${ENGINE_LIBRARY})
  #foreach(SF_SOURCE ${SF_SOURCES})
  #  message(Adding_"${SF_SOURCE}"_as_a_source_to_"${UTIL}")
  #  target_sources(${UTIL} PRIVATE ${SF_SOURCE}.cpp)
  #endforeach(SF_SOURCE)
endforeach(UTIL)

install(TARGETS ${UTILS} DESTINATION "${CMAKE_INSTALL_FULL_BINDIR}" COMPONENT Runtime)
