# Generates standalone HTML API docs for utils/surfer, separate from
# Regina's own engine docs (see apidocs/CMakeLists.txt). Not part of the
# default build or install -- run explicitly via the surfer_docs target.

CONFIGURE_FILE (
  Doxyfile.in
  "${PROJECT_BINARY_DIR}/utils/surfer/apidocs/Doxyfile"
  ESCAPE_QUOTES @ONLY
)

ADD_CUSTOM_COMMAND(OUTPUT index.html VERBATIM
  COMMAND "${DOXYGEN_EXECUTABLE}" "${PROJECT_BINARY_DIR}/utils/surfer/apidocs/Doxyfile"
  )

ADD_CUSTOM_TARGET(surfer_docs DEPENDS index.html)
