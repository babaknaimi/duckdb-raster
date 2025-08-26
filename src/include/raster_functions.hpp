#pragma once

#include "duckdb.hpp"

namespace duckdb {

/// Register all helper functions and table functions for working with raster
/// datasets.  In addition to metadata and simple statistics (width, height,
/// band count, min/max/mean), this also registers functions for pixel
/// retrieval, conversion between cell identifiers and 2D coordinates,
/// temporal aggregation over lists, summation across entire bands, and
/// reading multiple rasters as a long-format data cube.  The implementation
/// lives in raster_functions.cpp.  See RegisterRaster() in raster.cpp for
/// usage.
void RegisterRasterFunctions(DatabaseInstance &db);

} // namespace duckdb