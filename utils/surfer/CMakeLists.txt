# Every utils/surfer executable (and test/utility under tests/) links
# against some subset of these implementation files. Rather than each
# target relisting its own subset by hand (which silently goes stale --
# see surfacesearch.cpp needing to be added in three places at once), they
# all just link this one static library: a traditional archive only pulls
# in the .o members a given target's unresolved symbols actually need, so
# e.g. a test that only exercises rollbackunionfind.cpp still ends up
# linking just that, with no per-target bookkeeping required.
add_library(surfer_common STATIC
  knotbuilder.cpp cobordismbuilder.cpp simplicialprism.cpp skeleton.cpp
  embeddedsubmanifold.cpp embeddingsearch.cpp surfacesearch.cpp
  linkcomplement.cpp enumerate_cis.cpp collar.cpp rollbackunionfind.cpp
  vertexlinks.cpp pairsig.cpp)
target_link_libraries(surfer_common PUBLIC ${ENGINE_LIBRARY})

set(UTILS
  triangulateknot
  surfer
)

foreach(EXE ${UTILS})
  add_executable(${EXE} ${EXE}.cpp)
  target_link_libraries(${EXE} surfer_common)
endforeach()

install(TARGETS ${UTILS} DESTINATION "${CMAKE_INSTALL_FULL_BINDIR}" COMPONENT Runtime)

# Let files under tests/ resolve bare quoted includes (e.g.
# #include "cobordismbuilder.h") against this directory.
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory(tests)

if (REGINA_BUILD_ENGINE_DOCS)
  add_subdirectory(apidocs)
endif (REGINA_BUILD_ENGINE_DOCS)
