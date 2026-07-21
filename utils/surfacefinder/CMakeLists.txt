set(UTILS
  triangulateknot
  searchforembedded
)

set(triangulateknot_EXTRA_SOURCES knotbuilder.cpp cobordismbuilder.cpp simplicialprism.cpp linkcomplement.cpp)

set(searchforembedded_EXTRA_SOURCES
  knotbuilder.cpp cobordismbuilder.cpp simplicialprism.cpp skeleton.cpp
  embeddedsubmanifold.cpp embeddingsearch.cpp linkcomplement.cpp
  enumerate_cis.cpp)

foreach(EXE ${UTILS})
  add_executable(${EXE} ${EXE}.cpp ${${EXE}_EXTRA_SOURCES})
  target_link_libraries(${EXE} ${ENGINE_LIBRARY})
endforeach()

install(TARGETS ${UTILS} DESTINATION "${CMAKE_INSTALL_FULL_BINDIR}" COMPONENT Runtime)

# Let files under tests/ resolve bare quoted includes (e.g.
# #include "cobordismbuilder.h") against this directory.
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(tests)
