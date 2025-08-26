# This config file tells DuckDB how to build the raster extension when using the
# in-tree build system or extension-ci-tools.  The extension is located in
# this repository (the same directory as this file) and uses its own
# CMakeLists.txt to define the targets.  DONT_LINK ensures the extension is
# built as a loadable module rather than linked into the DuckDB binary.

duckdb_extension_load(raster
    DONT_LINK
    SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}"
)