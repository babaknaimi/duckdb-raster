#include "duckdb/main/extension.hpp"
#include "raster_register.hpp"

namespace duckdb {

class RasterExtension final : public Extension {
public:
    void Load(DuckDB &db) override {
        RegisterRaster(*db.instance);
    }
    std::string Name() override {
        // The name of the extension as exposed to DuckDB.  Changing this from
        // "geotiff" to "raster" reflects that the extension now supports a
        // variety of raster formats and functions beyond GeoTIFF.
        return "raster";
    }
    std::string Version() const override {
        return DuckDB::LibraryVersion();
    }
};

} // namespace duckdb