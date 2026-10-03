# The project root, in four parts, each its own static library:
#   surfer/                surfer_lib: a knot-agnostic library, PL surfaces in
#                          triangulated 4-manifolds
#   diagramtriangulation/  diagramtriangulation_lib: the triangulations a search
#                          runs in, built from a diagram
#   linknaming/            linknaming_lib: naming a link, with proof; its
#                          complement/ folder is the bottom layer,
#                          linknaming_complement
#   cobound/               cobound_lib: bounding g4 by cobordisms
#
# The dependency rule, enforced by the build: each library links only what it
# may depend on, and sees only those parts' headers.
#   linknaming_complement     -> Regina
#   surfer_lib                -> Regina, OpenSSL's libcrypto, linknaming_complement
#   diagramtriangulation_lib  -> linknaming_complement
#   linknaming_lib            -> Regina, SQLite, linknaming_complement
#   cobound_lib               -> all of the above
# Executables link the libraries they need; surfer's tests and driver may also
# use diagramtriangulation_lib, and linknaming's tests and census-name
# generator too (none of them is a library another part links).
#
# Every project include is written relative to this directory
# (#include "surfer/enumeration/surfacesearch.h"). A library does not see this
# directory itself: it sees an include view (surfer_include_view() below), a
# directory in the build tree holding a symbolic link to its own part only.
# So a source including a header from a part its library does not link fails
# to compile. .scratch/gates/g_deps.py checks the same rule include by
# include, for executables too.

# Regina turns this on for the whole tree; here it would put each part's
# directory (and this one) on every target's include path, and a part-local
# census/ or link/ header would shadow Regina's.
set(CMAKE_INCLUDE_CURRENT_DIR OFF)

# Where every part puts the executables that scripts and campaigns expect in
# build/utils/surfer/ (see each part's CMakeLists.txt).
set(SURFER_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}")

# surfer_include_view(<target> <part path>...): gives <target> (PUBLIC) an
# include directory holding a symbolic link to each <part path>, a path
# relative to this directory such as "surfer" or "linknaming/complement", at
# that same relative path.
set(SURFER_SOURCE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")
set(SURFER_INCLUDE_VIEWS "${CMAKE_CURRENT_BINARY_DIR}/include-views")
function(surfer_include_view target)
  set(view "${SURFER_INCLUDE_VIEWS}/${target}")
  foreach(part IN LISTS ARGN)
    get_filename_component(parent "${view}/${part}" DIRECTORY)
    file(MAKE_DIRECTORY "${parent}")
    if(NOT IS_SYMLINK "${view}/${part}")
      file(CREATE_LINK "${SURFER_SOURCE_ROOT}/${part}" "${view}/${part}" SYMBOLIC)
    endif()
  endforeach()
  target_include_directories(${target} PUBLIC "${view}")
endfunction()

# Local census SQLite lookup (censusnaming.cpp's
# census::localCensusLookup()): not used anywhere else in the engine, so
# this dependency is scoped to linknaming_lib alone rather than added at the
# top level.
find_package(SQLite3 REQUIRED)

# SHA-1 (surfer/pairsig/sha1.h): OpenSSL's libcrypto, for the frontier
# fingerprint, the pair-signature context's ambient key and the cobordism key.
find_package(OpenSSL REQUIRED COMPONENTS Crypto)

# mimalloc in place of glibc's malloc, for every surfer binary and test
# (2026-09-29). A perf profile of a production row put about a third of all
# CPU in glibc's malloc/free, and the drain ran ~1.8x faster under mimalloc
# (utils/surfer/README.md, "Performance"). Linked statically and whole, so the
# executable itself defines malloc and interposes it for the engine library
# too, and a host without the package can use a copied archive
# (-DSURFER_MIMALLOC=<path>, or ~/.local/lib/libmimalloc.a). Results cannot
# depend on the allocator; only speed does. It is one INTERFACE target,
# surfer_mimalloc, which every executable links once: the whole-archive item
# must reach each executable exactly once, whatever libraries it links.
add_library(surfer_mimalloc INTERFACE)
find_library(SURFER_MIMALLOC NAMES libmimalloc.a
             HINTS "$ENV{HOME}/.local/lib" DOC "static mimalloc archive (optional)")
# The archive must define malloc as real machine code. Arch's packaged one is
# GCC LTO bitcode, which a linker without the matching plugin skips with only
# a warning -- leaving glibc's malloc in place. Build mimalloc from source
# instead (see README, "Performance").
if(SURFER_MIMALLOC)
  execute_process(COMMAND nm "${SURFER_MIMALLOC}" OUTPUT_VARIABLE _mi_syms
                  ERROR_QUIET)
  if(NOT _mi_syms MATCHES "[0-9a-f]+ T malloc\n")
    message(WARNING "surfer: ${SURFER_MIMALLOC} defines no malloc as machine "
                    "code (LTO bitcode?); not linking it")
    set(SURFER_MIMALLOC "")
  endif()
endif()
if(SURFER_MIMALLOC)
  message(STATUS "surfer: linking mimalloc from ${SURFER_MIMALLOC}")
  target_link_libraries(surfer_mimalloc INTERFACE
    "$<LINK_LIBRARY:WHOLE_ARCHIVE,${SURFER_MIMALLOC}>")
else()
  message(WARNING "surfer: no libmimalloc.a found; using glibc malloc, which "
                  "costs the drain about 1.8x (set SURFER_MIMALLOC)")
endif()

# Bottom layer first; each part adds its library, then its executables and
# tests (which may name libraries defined later).
add_subdirectory(linknaming)
add_subdirectory(surfer)
add_subdirectory(diagramtriangulation)
# The cobound part builds in cobound_part/: its program, cobound, goes to
# build/utils/surfer/cobound, where scripts and campaigns expect the
# project's programs (a build directory of that name would be in its way).
add_subdirectory(cobound cobound_part)

if (REGINA_BUILD_ENGINE_DOCS)
  add_subdirectory(apidocs)
endif (REGINA_BUILD_ENGINE_DOCS)
