# surfer: PL surfaces in triangulated 4-manifolds. Its library (see
# ../CMakeLists.txt for the dependency rule), then its proof-of-concept
# driver and its tests.
add_library(surfer_lib STATIC
  submanifold/skeleton.cpp submanifold/submanifold.cpp
  enumeration/submanifoldsearch.cpp enumeration/surfacesearch.cpp
  submanifold/linkingnumber.cpp enumeration/namecache.cpp
  enumeration/inducedsubgraphs.cpp
  submanifold/rollbackunionfind.cpp submanifold/vertexlinks.cpp
  pairsig/pairsig.cpp pairsig/sha1.cpp report/csvwriter.cpp report/atomicwrite.cpp report/progress.cpp
  enumeration/searchstrategy.cpp)
target_link_libraries(surfer_lib PUBLIC ${ENGINE_LIBRARY} linknaming_complement
                      PRIVATE OpenSSL::Crypto)
surfer_include_view(surfer_lib surfer)

# The driver lands in this part's own build directory,
# build/utils/surfer/surfer/surfer: build/utils/surfer/surfer is this
# directory's build directory, so it cannot also be the executable. Its
# --pd mode builds the ambient through diagramtriangulation.
add_executable(surfer surfer.cpp)
target_link_libraries(surfer surfer_lib diagramtriangulation_lib surfer_mimalloc)
install(TARGETS surfer DESTINATION "${CMAKE_INSTALL_FULL_BINDIR}" COMPONENT Runtime)

add_subdirectory(tests)
