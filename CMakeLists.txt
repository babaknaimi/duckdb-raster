cmake_minimum_required(VERSION 3.16)
project(raster)

# DuckDB macros must exist (we're being built from duckdb root via EXTENSION_CONFIGS)
if(NOT COMMAND build_loadable_extension)
  message(FATAL_ERROR "DuckDB build macros not loaded. Configure DuckDB with EXTENSION_CONFIGS/DUCKDB_EXTENSION_CONFIGS.")
endif()

set(EXT_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")

set(EXTENSION_SOURCES
  "${EXT_ROOT}/src/raster.cpp"
  "${EXT_ROOT}/src/raster_functions.cpp"
  # Uncomment if you also build the static-in-CLI flavor:
  # "${EXT_ROOT}/src/raster_extension.cpp"
)

# Make headers visible to both targets without referencing targets that may be skipped
include_directories("${EXT_ROOT}/src/include")

# Always create the targets (DuckDB will only build what it needs)
build_loadable_extension(raster "" ${EXTENSION_SOURCES})
if(EXTENSION_STATIC_BUILD)
  build_static_extension(raster ${EXTENSION_SOURCES})
endif()

# Try GDAL softly (no REQUIRED): if missing (e.g. tidy job), we just skip linking.
find_package(GDAL CONFIG QUIET)
if(GDAL_FOUND)
  if(TARGET raster_loadable_extension)
    target_link_libraries(raster_loadable_extension GDAL::GDAL)
  endif()
  if(TARGET raster_extension)
    target_link_libraries(raster_extension GDAL::GDAL)
  endif()

  # RPATH so the loadable looks next to itself for libs
  if(APPLE AND TARGET raster_loadable_extension)
    set_target_properties(raster_loadable_extension PROPERTIES
      BUILD_RPATH "@loader_path" INSTALL_RPATH "@loader_path")
  elseif(UNIX AND TARGET raster_loadable_extension)
    set_target_properties(raster_loadable_extension PROPERTIES
      BUILD_RPATH "\$ORIGIN" INSTALL_RPATH "\$ORIGIN")
  endif()
else()
  message(STATUS "GDAL not found - skipping link (tidy/format mode).")
endif()