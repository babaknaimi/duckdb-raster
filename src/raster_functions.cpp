#include "duckdb.hpp"
#include "duckdb/function/scalar/scalar_function.hpp"
#include "duckdb/main/extension_util.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/common/types/validity_mask.hpp"

#include "raster_register.hpp"
#include "raster_functions.hpp"

#include "gdal_priv.h"

using namespace duckdb;

namespace {

// Helper to open a GDAL dataset given a path.  Returns a unique_ptr
// managing the dataset, or throws a IOException if opening fails.
static std::unique_ptr<GDALDataset> OpenDataset(const std::string &path) {
    GDALAllRegister();
    GDALDataset *raw = static_cast<GDALDataset *>(GDALOpen(path.c_str(), GA_ReadOnly));
    if (!raw) {
        throw IOException("GDALOpen failed for '%s'", path.c_str());
    }
    return std::unique_ptr<GDALDataset>(raw);
}

// raster_width(path) -> BIGINT
static void RasterWidthFunction(DataChunk &input, ExpressionState &state, Vector &result) {
    auto count = input.size();
    result.SetVectorType(VectorType::FLAT_VECTOR);
    for (idx_t row = 0; row < count; row++) {
        auto &vec = input.data[0];
        Value path_val = vec.GetValue(row);
        if (path_val.IsNull()) {
            result.SetValue(row, Value());
            continue;
        }
        std::string path = StringValue::Get(path_val);
        try {
            auto ds = OpenDataset(path);
            int width = ds->GetRasterXSize();
            result.SetValue(row, Value::BIGINT(width));
        } catch (...) {
            // If unable to open, return NULL
            result.SetValue(row, Value());
        }
    }
}

// raster_height(path) -> BIGINT
static void RasterHeightFunction(DataChunk &input, ExpressionState &state, Vector &result) {
    auto count = input.size();
    result.SetVectorType(VectorType::FLAT_VECTOR);
    for (idx_t row = 0; row < count; row++) {
        auto &vec = input.data[0];
        Value path_val = vec.GetValue(row);
        if (path_val.IsNull()) {
            result.SetValue(row, Value());
            continue;
        }
        std::string path = StringValue::Get(path_val);
        try {
            auto ds = OpenDataset(path);
            int height = ds->GetRasterYSize();
            result.SetValue(row, Value::BIGINT(height));
        } catch (...) {
            result.SetValue(row, Value());
        }
    }
}

// raster_band_count(path) -> BIGINT
static void RasterBandCountFunction(DataChunk &input, ExpressionState &state, Vector &result) {
    auto count = input.size();
    result.SetVectorType(VectorType::FLAT_VECTOR);
    for (idx_t row = 0; row < count; row++) {
        Value path_val = input.data[0].GetValue(row);
        if (path_val.IsNull()) {
            result.SetValue(row, Value());
            continue;
        }
        std::string path = StringValue::Get(path_val);
        try {
            auto ds = OpenDataset(path);
            int bands = ds->GetRasterCount();
            result.SetValue(row, Value::BIGINT(bands));
        } catch (...) {
            result.SetValue(row, Value());
        }
    }
}

// Helper for statistics: returns min/max/mean for a given band.  Throws on error.
static void GetBandStatistics(GDALDataset *ds, int band_index, double &min, double &max, double &mean) {
    GDALRasterBand *band = ds->GetRasterBand(band_index);
    if (!band) {
        throw IOException("Invalid band index");
    }
    double stddev = 0.0;
    // approximate=FALSE, force=TRUE => compute statistics on demand
    CPLErr err = band->GetStatistics(FALSE, TRUE, &min, &max, &mean, &stddev);
    if (err != CE_None) {
        throw IOException("GetStatistics failed");
    }
}

// raster_band_min(path, band) -> DOUBLE
static void RasterMinFunction(DataChunk &input, ExpressionState &state, Vector &result) {
    auto count = input.size();
    result.SetVectorType(VectorType::FLAT_VECTOR);
    for (idx_t row = 0; row < count; row++) {
        Value path_val = input.data[0].GetValue(row);
        Value band_val = input.data[1].GetValue(row);
        if (path_val.IsNull() || band_val.IsNull()) {
            result.SetValue(row, Value());
            continue;
        }
        std::string path = StringValue::Get(path_val);
        int band_index = band_val.GetValue<int32_t>();
        try {
            auto ds = OpenDataset(path);
            double min = 0.0, max = 0.0, mean = 0.0;
            GetBandStatistics(ds.get(), band_index, min, max, mean);
            result.SetValue(row, Value::DOUBLE(min));
        } catch (...) {
            result.SetValue(row, Value());
        }
    }
}

// raster_band_max(path, band) -> DOUBLE
static void RasterMaxFunction(DataChunk &input, ExpressionState &state, Vector &result) {
    auto count = input.size();
    result.SetVectorType(VectorType::FLAT_VECTOR);
    for (idx_t row = 0; row < count; row++) {
        Value path_val = input.data[0].GetValue(row);
        Value band_val = input.data[1].GetValue(row);
        if (path_val.IsNull() || band_val.IsNull()) {
            result.SetValue(row, Value());
            continue;
        }
        std::string path = StringValue::Get(path_val);
        int band_index = band_val.GetValue<int32_t>();
        try {
            auto ds = OpenDataset(path);
            double min = 0.0, max = 0.0, mean = 0.0;
            GetBandStatistics(ds.get(), band_index, min, max, mean);
            result.SetValue(row, Value::DOUBLE(max));
        } catch (...) {
            result.SetValue(row, Value());
        }
    }
}

// raster_band_mean(path, band) -> DOUBLE
static void RasterMeanFunction(DataChunk &input, ExpressionState &state, Vector &result) {
    auto count = input.size();
    result.SetVectorType(VectorType::FLAT_VECTOR);
    for (idx_t row = 0; row < count; row++) {
        Value path_val = input.data[0].GetValue(row);
        Value band_val = input.data[1].GetValue(row);
        if (path_val.IsNull() || band_val.IsNull()) {
            result.SetValue(row, Value());
            continue;
        }
        std::string path = StringValue::Get(path_val);
        int band_index = band_val.GetValue<int32_t>();
        try {
            auto ds = OpenDataset(path);
            double min = 0.0, max = 0.0, mean = 0.0;
            GetBandStatistics(ds.get(), band_index, min, max, mean);
            result.SetValue(row, Value::DOUBLE(mean));
        } catch (...) {
            result.SetValue(row, Value());
        }
    }
}

// raster_srs(path) -> VARCHAR
// Return the spatial reference system (projection) of a raster dataset as a WKT
static void RasterSRSFunction(DataChunk &input, ExpressionState &state, Vector &result) {
    auto count = input.size();
    result.SetVectorType(VectorType::FLAT_VECTOR);
    for (idx_t row = 0; row < count; row++) {
        Value path_val = input.data[0].GetValue(row);
        if (path_val.IsNull()) {
            result.SetValue(row, Value());
            continue;
        }
        std::string path = StringValue::Get(path_val);
        try {
            auto ds = OpenDataset(path);
            const char *proj = ds->GetProjectionRef();
            if (proj && proj[0] != '\0') {
                result.SetValue(row, Value(proj));
            } else {
                result.SetValue(row, Value());
            }
        } catch (...) {
            result.SetValue(row, Value());
        }
    }
}

} // anonymous namespace

namespace duckdb {

void RegisterRasterFunctions(DatabaseInstance &db) {
    // width
    ScalarFunction width_fn("raster_width", {LogicalType::VARCHAR}, LogicalType::BIGINT, RasterWidthFunction);
    ExtensionUtil::RegisterFunction(db, width_fn);

    // height
    ScalarFunction height_fn("raster_height", {LogicalType::VARCHAR}, LogicalType::BIGINT, RasterHeightFunction);
    ExtensionUtil::RegisterFunction(db, height_fn);

    // band count
    ScalarFunction band_count_fn("raster_band_count", {LogicalType::VARCHAR}, LogicalType::BIGINT, RasterBandCountFunction);
    ExtensionUtil::RegisterFunction(db, band_count_fn);

    // min
    ScalarFunction min_fn("raster_band_min", {LogicalType::VARCHAR, LogicalType::INTEGER}, LogicalType::DOUBLE, RasterMinFunction);
    ExtensionUtil::RegisterFunction(db, min_fn);

    // max
    ScalarFunction max_fn("raster_band_max", {LogicalType::VARCHAR, LogicalType::INTEGER}, LogicalType::DOUBLE, RasterMaxFunction);
    ExtensionUtil::RegisterFunction(db, max_fn);

    // mean
    ScalarFunction mean_fn("raster_band_mean", {LogicalType::VARCHAR, LogicalType::INTEGER}, LogicalType::DOUBLE, RasterMeanFunction);
    ExtensionUtil::RegisterFunction(db, mean_fn);

    // spatial reference system (projection) as WKT
    ScalarFunction srs_fn("raster_srs", {LogicalType::VARCHAR}, LogicalType::VARCHAR, RasterSRSFunction);
    ExtensionUtil::RegisterFunction(db, srs_fn);
}

} // namespace duckdb