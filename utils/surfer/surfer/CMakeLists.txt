# surfer: PL surfaces in triangulated 4-manifolds. Its sources are compiled
# into surfer_common (../CMakeLists.txt) like every other part's; this file
# adds its proof-of-concept driver and its tests.
#
# The driver lands in this part's own build directory,
# build/utils/surfer/surfer/surfer: build/utils/surfer/surfer is this
# directory's build directory, so it cannot also be the executable.

add_executable(surfer surfer.cpp)
target_link_libraries(surfer surfer_common surfer_mimalloc)
install(TARGETS surfer DESTINATION "${CMAKE_INSTALL_FULL_BINDIR}" COMPONENT Runtime)

add_subdirectory(tests)
