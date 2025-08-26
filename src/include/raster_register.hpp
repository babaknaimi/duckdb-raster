#pragma once
#include "duckdb.hpp"

namespace duckdb {

/// Register the core raster functions and table functions.  This function is
/// called by the extension loader to register all functionality provided by
/// this extension.  It registers the `read_raster` table function as well as
/// additional scalar functions exposed via RegisterRasterFunctions.
void RegisterRaster(DatabaseInstance &db);

/// Forward declaration of the helper that registers additional scalar
/// functions.  Implementation lives in raster_functions.cpp.
void RegisterRasterFunctions(DatabaseInstance &db);

} // namespace duckdb
