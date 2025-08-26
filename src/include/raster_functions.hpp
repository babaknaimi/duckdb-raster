#pragma once

#include "duckdb.hpp"

namespace duckdb {

/// Register all scalar helper functions for working with raster datasets.
/// These functions provide metadata and simple statistics for raster files,
/// including width, height, band count, and basic statistical summaries
/// (minimum, maximum and mean) over a given band.  The implementation lives
/// in raster_functions.cpp.  See RegisterRaster() in raster.cpp for usage.
void RegisterRasterFunctions(DatabaseInstance &db);

} // namespace duckdb