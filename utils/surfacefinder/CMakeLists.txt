set(UTILS
  bogocheck
  surfacefinder
  fillmanifold
  knotbuilder
  )

foreach(UTIL ${UTILS})
  add_executable(${UTIL} ${UTIL}.cpp)
  target_link_libraries(${UTIL} ${ENGINE_LIBRARY})
  #foreach(SF_SOURCE ${SF_SOURCES})
  #  message(Adding_"${SF_SOURCE}"_as_a_source_to_"${UTIL}")
  #  target_sources(${UTIL} PRIVATE ${SF_SOURCE}.cpp)
  #endforeach(SF_SOURCE)
endforeach(UTIL)

install(TARGETS ${UTILS} DESTINATION "${CMAKE_INSTALL_FULL_BINDIR}" COMPONENT Runtime)
