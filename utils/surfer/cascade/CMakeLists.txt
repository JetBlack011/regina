# cascade: goal-directed chained searches (README.md). Its library sources
# are compiled into surfer_common (../CMakeLists.txt) like every other surfer
# source; this file adds the driver and the tests.

add_executable(cascadesearch cascadesearch.cpp)
target_link_libraries(cascadesearch surfer_common)

add_subdirectory(tests)
