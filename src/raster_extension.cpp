#include "duckdb/main/extension.hpp"
#include "raster_register.hpp"

namespace duckdb {
class RasterExtension final : public Extension {
public:
	void Load(DuckDB &db) override {
		RegisterRaster(*db.instance);
	}
	std::string Name() override {
		return "geotiff";
	}
	std::string Version() const override {
		return DuckDB::LibraryVersion();
	}
};
} // namespace duckdb
