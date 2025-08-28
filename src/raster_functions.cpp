#include "duckdb.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension_util.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/common/types/validity_mask.hpp"

#include "raster_register.hpp"
#include "raster_functions.hpp"
#include "raster_gdal.hpp"

#include "gdal_priv.h"

#include <cmath> // for std::isfinite

using namespace duckdb;

namespace {



// raster_width(path) -> BIGINT
static void RasterWidthFunction(DataChunk &input, ExpressionState &state,
                                Vector &result) {
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
static void RasterHeightFunction(DataChunk &input, ExpressionState &state,
                                 Vector &result) {
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
static void RasterBandCountFunction(DataChunk &input, ExpressionState &state,
                                    Vector &result) {
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

// Helper for statistics: returns min/max/mean for a given band.  Throws on
// error.
static void GetBandStatistics(GDALDataset *ds, int band_index, double &min,
                              double &max, double &mean) {
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
static void RasterMinFunction(DataChunk &input, ExpressionState &state,
                              Vector &result) {
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
static void RasterMaxFunction(DataChunk &input, ExpressionState &state,
                              Vector &result) {
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
static void RasterMeanFunction(DataChunk &input, ExpressionState &state,
                               Vector &result) {
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
static void RasterSRSFunction(DataChunk &input, ExpressionState &state,
                              Vector &result) {
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

// =============================================================================
//  Scalar Function: raster_pixel_value(path, band, row, col)
//
//  Returns the value of a single pixel from the specified raster dataset.  The
//  row and column are zero-based indices.  If any input is NULL or the
//  coordinates are out of bounds, NULL is returned.  Only the first band of
//  the dataset is used if band is 1.  This function is useful for retrieving
//  values from a raster table given row/col coordinates computed from a
//  cell_id.
static void RasterPixelValueFunction(DataChunk &input, ExpressionState &,
                                     Vector &result) {
  auto count = input.size();
  result.SetVectorType(VectorType::FLAT_VECTOR);
  for (idx_t row_idx = 0; row_idx < count; row_idx++) {
    Value path_val = input.data[0].GetValue(row_idx);
    Value band_val = input.data[1].GetValue(row_idx);
    Value row_val = input.data[2].GetValue(row_idx);
    Value col_val = input.data[3].GetValue(row_idx);
    if (path_val.IsNull() || band_val.IsNull() || row_val.IsNull() ||
        col_val.IsNull()) {
      result.SetValue(row_idx, Value());
      continue;
    }
    std::string path = StringValue::Get(path_val);
    int band_idx = band_val.GetValue<int32_t>();
    int64_t r = row_val.GetValue<int64_t>();
    int64_t c = col_val.GetValue<int64_t>();
    try {
      auto ds = OpenDataset(path);
      int width = ds->GetRasterXSize();
      int height = ds->GetRasterYSize();
      if (r < 0 || c < 0 || r >= height || c >= width) {
        result.SetValue(row_idx, Value());
        continue;
      }
      // Clamp band index to [1, band_count]
      int band_count = ds->GetRasterCount();
      if (band_idx < 1 || band_idx > band_count) {
        result.SetValue(row_idx, Value());
        continue;
      }
      double val = 0.0;
      GDALRasterBand *band = ds->GetRasterBand(band_idx);
      if (!band) {
        result.SetValue(row_idx, Value());
        continue;
      }
      CPLErr err =
          band->RasterIO(GF_Read, c, r, 1, 1, &val, 1, 1, GDT_Float64, 0, 0);
      if (err != CE_None) {
        result.SetValue(row_idx, Value());
      } else {
        result.SetValue(row_idx, Value::DOUBLE(val));
      }
    } catch (...) {
      result.SetValue(row_idx, Value());
    }
  }
}

// =============================================================================
//  Scalar Functions: raster_cell_row and raster_cell_col
//
//  Given a cell_id and the raster width, compute the zero-based row and
//  column indexes.  These helpers are useful to convert the output of
//  read_raster (which yields a linear cell id) into 2D indices that can
//  be passed to raster_pixel_value.
static void RasterCellRowFunction(DataChunk &input, ExpressionState &,
                                  Vector &result) {
  auto count = input.size();
  result.SetVectorType(VectorType::FLAT_VECTOR);
  for (idx_t row_idx = 0; row_idx < count; row_idx++) {
    Value cell_val = input.data[0].GetValue(row_idx);
    Value width_val = input.data[1].GetValue(row_idx);
    if (cell_val.IsNull() || width_val.IsNull()) {
      result.SetValue(row_idx, Value());
      continue;
    }
    int64_t cell_id = cell_val.GetValue<int64_t>();
    int64_t width = width_val.GetValue<int64_t>();
    if (width <= 0) {
      result.SetValue(row_idx, Value());
      continue;
    }
    int64_t row = cell_id / width;
    result.SetValue(row_idx, Value::BIGINT(row));
  }
}

static void RasterCellColFunction(DataChunk &input, ExpressionState &,
                                  Vector &result) {
  auto count = input.size();
  result.SetVectorType(VectorType::FLAT_VECTOR);
  for (idx_t row_idx = 0; row_idx < count; row_idx++) {
    Value cell_val = input.data[0].GetValue(row_idx);
    Value width_val = input.data[1].GetValue(row_idx);
    if (cell_val.IsNull() || width_val.IsNull()) {
      result.SetValue(row_idx, Value());
      continue;
    }
    int64_t cell_id = cell_val.GetValue<int64_t>();
    int64_t width = width_val.GetValue<int64_t>();
    if (width <= 0) {
      result.SetValue(row_idx, Value());
      continue;
    }
    int64_t col = cell_id % width;
    result.SetValue(row_idx, Value::BIGINT(col));
  }
}

// =============================================================================
//  Scalar Function: raster_temporal_mean(values)
//
//  Computes the mean of a list of numeric values, ignoring NULLs.  Returns
//  NULL if no non-null values are present.
static void RasterTemporalMeanFunction(DataChunk &input, ExpressionState &,
                                       Vector &result) {
  auto count = input.size();
  result.SetVectorType(VectorType::FLAT_VECTOR);
  for (idx_t row_idx = 0; row_idx < count; row_idx++) {
    Value list_val = input.data[0].GetValue(row_idx);
    if (list_val.IsNull()) {
      result.SetValue(row_idx, Value());
      continue;
    }
    auto &children = ListValue::GetChildren(list_val);
    double sum = 0.0;
    idx_t n = 0;
    for (auto &child : children) {
      if (child.IsNull())
        continue;
      sum += child.GetValue<double>();
      n++;
    }
    if (n == 0) {
      result.SetValue(row_idx, Value());
    } else {
      result.SetValue(row_idx, Value::DOUBLE(sum / (double)n));
    }
  }
}

// =============================================================================
//  Scalar Function: raster_band_sum(path, band)
//
//  Computes the sum of all pixel values in a given band of the raster.  If the
//  dataset cannot be opened or the band is invalid, returns NULL.  This
//  function can be used for simple spatial aggregation across the entire
//  raster.  Note: reading the entire raster into memory may be expensive for
//  large datasets.
static void RasterBandSumFunction(DataChunk &input, ExpressionState &,
                                  Vector &result) {
  auto count = input.size();
  result.SetVectorType(VectorType::FLAT_VECTOR);
  for (idx_t row_idx = 0; row_idx < count; row_idx++) {
    Value path_val = input.data[0].GetValue(row_idx);
    Value band_val = input.data[1].GetValue(row_idx);
    if (path_val.IsNull() || band_val.IsNull()) {
      result.SetValue(row_idx, Value());
      continue;
    }
    std::string path = StringValue::Get(path_val);
    int band_idx = band_val.GetValue<int32_t>();
    try {
      auto ds = OpenDataset(path);
      int band_count = ds->GetRasterCount();
      if (band_idx < 1 || band_idx > band_count) {
        result.SetValue(row_idx, Value());
        continue;
      }
      GDALRasterBand *band = ds->GetRasterBand(band_idx);
      if (!band) {
        result.SetValue(row_idx, Value());
        continue;
      }
      int width = ds->GetRasterXSize();
      int height = ds->GetRasterYSize();
      // Read the entire band row by row to avoid allocating huge buffers
      double sum = 0.0;
      std::vector<double> buffer;
      buffer.resize(width);
      for (int y = 0; y < height; y++) {
        CPLErr err = band->RasterIO(GF_Read, 0, y, width, 1, buffer.data(),
                                    width, 1, GDT_Float64, 0, 0);
        if (err != CE_None) {
          throw IOException("RasterIO failed when summing band");
        }
        for (int x = 0; x < width; x++) {
          // Treat NaN as NULL
          double v = buffer[x];
          if (std::isfinite(v)) {
            sum += v;
          }
        }
      }
      result.SetValue(row_idx, Value::DOUBLE(sum));
    } catch (...) {
      result.SetValue(row_idx, Value());
    }
  }
}

// =============================================================================
//  Table Function: read_raster_cube(paths, vars, times)
//
//  Reads a collection of raster datasets (paths) that represent different
//  variables and/or time steps but share the same spatial extent.  It
//  produces a long-format table with one row per pixel per dataset, with
//  columns: cell_id, var (VARCHAR), time (VARCHAR) and value (DOUBLE).  The
//  user must supply three lists of equal length: paths, variable names and
//  time identifiers.  All datasets must have identical dimensions; an error
//  is thrown otherwise.  This function can be used to build a wide or long
//  data cube by pivoting the resulting table.
struct CubeBindData : public FunctionData {
  std::vector<std::string> paths;
  std::vector<std::string> vars;
  std::vector<std::string> times;
  idx_t width = 0;
  idx_t height = 0;
  unique_ptr<FunctionData> Copy() const override {
    return make_uniq<CubeBindData>(*this);
  }
  bool Equals(const FunctionData &other) const override {
    auto &o = other.Cast<CubeBindData>();
    return paths == o.paths && vars == o.vars && times == o.times &&
           width == o.width && height == o.height;
  }
};

struct CubeGlobalState : public GlobalTableFunctionState {
  std::vector<GDALDatasetPtr> ds;
  idx_t width = 0;
  idx_t height = 0;
  idx_t dataset_index = 0; // which dataset we are currently scanning
  int64_t pixel_index = 0; // 0..(width*height - 1)
  idx_t MaxThreads() const override { return 1; }
};

// Bind for read_raster_cube: parse three lists and ensure all files have same
// dimensions
static unique_ptr<FunctionData> CubeBind(ClientContext &,
                                         TableFunctionBindInput &input,
                                         vector<LogicalType> &types,
                                         vector<string> &names) {
  if (input.inputs.size() < 3) {
    throw BinderException(
        "read_raster_cube(paths, vars, times) requires three arguments");
  }
  std::vector<std::string> paths;
  std::vector<std::string> vars;
  std::vector<std::string> times;
  // paths: accept LIST<VARCHAR> or single VARCHAR
  auto parse_list = [&](const Value &v, std::vector<std::string> &out) {
    if (v.IsNull())
      return;
    if (v.type().id() == LogicalTypeId::LIST) {
      const auto &children = ListValue::GetChildren(v);
      out.reserve(children.size());
      for (auto &c : children) {
        out.push_back(StringValue::Get(c));
      }
    } else {
      // treat scalar as single-element list
      out.push_back(StringValue::Get(v));
    }
  };
  parse_list(input.inputs[0], paths);
  parse_list(input.inputs[1], vars);
  parse_list(input.inputs[2], times);
  if (paths.empty()) {
    throw BinderException("paths list must not be empty");
  }
  // If vars or times lists are empty, create default names equal to index
  if (vars.empty()) {
    vars.reserve(paths.size());
    for (size_t i = 0; i < paths.size(); i++) {
      vars.push_back("v" + std::to_string(i));
    }
  }
  if (times.empty()) {
    times.reserve(paths.size());
    for (size_t i = 0; i < paths.size(); i++) {
      times.push_back("t" + std::to_string(i));
    }
  }
  if (paths.size() != vars.size() || paths.size() != times.size()) {
    throw BinderException(
        "paths, vars and times lists must have the same length");
  }
  // Determine dimensions from first dataset and ensure all others match
  idx_t width = 0;
  idx_t height = 0;
  for (size_t i = 0; i < paths.size(); i++) {
    auto ds = OpenDataset(paths[i]);
    int w = ds->GetRasterXSize();
    int h = ds->GetRasterYSize();
    if (i == 0) {
      width = w;
      height = h;
    } else {
      if ((idx_t)w != width || (idx_t)h != height) {
        throw BinderException(
            "All rasters in read_raster_cube must have identical dimensions");
      }
    }
  }
  types = {LogicalType::BIGINT, LogicalType::VARCHAR, LogicalType::VARCHAR,
           LogicalType::DOUBLE};
  names = {"cell_id", "var", "time", "value"};
  auto bd = make_uniq<CubeBindData>();
  bd->paths = std::move(paths);
  bd->vars = std::move(vars);
  bd->times = std::move(times);
  bd->width = width;
  bd->height = height;
  return std::move(bd);
}

// Init: open all datasets and set state
static unique_ptr<GlobalTableFunctionState>
CubeInit(ClientContext &, TableFunctionInitInput &in) {
  auto &bd = in.bind_data->Cast<CubeBindData>();
  auto st = make_uniq<CubeGlobalState>();
  st->width = bd.width;
  st->height = bd.height;
  st->dataset_index = 0;
  st->pixel_index = 0;
  st->ds.reserve(bd.paths.size());
  try {
    for (auto &p : bd.paths) {
      st->ds.push_back(OpenDataset(p));
    }
  } catch (...) {
    // If any dataset fails to open, leave st->ds empty so scan returns no rows
    st->ds.clear();
  }
  return std::move(st);
}

// Scan: produce rows of (cell_id, var, time, value) across all datasets
static void CubeScan(ClientContext &, TableFunctionInput &in, DataChunk &out) {
  auto &st = in.global_state->Cast<CubeGlobalState>();
  if (st.ds.empty()) {
    out.SetCardinality(0);
    return;
  }
  auto &bd = in.bind_data->Cast<CubeBindData>();
  const idx_t total_pixels = st.width * st.height;
  idx_t out_idx = 0;
  auto cell_data = FlatVector::GetData<int64_t>(out.data[0]);
  auto var_data = out.data[1];
  auto time_data = out.data[2];
  // We'll set the value vector using SetValue to properly handle NULLs
  while (out_idx < STANDARD_VECTOR_SIZE && st.dataset_index < st.ds.size()) {
    if (st.pixel_index >= (int64_t)total_pixels) {
      // move to next dataset
      st.dataset_index++;
      st.pixel_index = 0;
      continue;
    }
    // compute row/col for current pixel
    int64_t cell_id = st.pixel_index;
    int64_t r = st.pixel_index / (int64_t)st.width;
    int64_t c = st.pixel_index % (int64_t)st.width;
    // read pixel value from current dataset
    double val = 0.0;
    bool ok = false;
    if (st.dataset_index < st.ds.size()) {
      GDALDataset *ds = st.ds[st.dataset_index].get();
      GDALRasterBand *band = ds->GetRasterBand(1);
      if (band) {
        CPLErr err = band->RasterIO(GF_Read, (int)c, (int)r, 1, 1, &val, 1, 1,
                                    GDT_Float64, 0, 0);
        if (err == CE_None) {
          ok = true;
        }
      }
    }
    // set output values
    cell_data[out_idx] = cell_id;
    // set var and time using SetValue to handle string allocation
    var_data.SetValue(out_idx, Value(bd.vars[st.dataset_index]));
    time_data.SetValue(out_idx, Value(bd.times[st.dataset_index]));
    if (ok && std::isfinite(val)) {
      out.data[3].SetValue(out_idx, Value::DOUBLE(val));
    } else {
      out.data[3].SetValue(out_idx, Value());
    }
    out_idx++;
    st.pixel_index++;
  }
  out.SetCardinality(out_idx);
}

// =============================================================================
//  Table Function: raster_metadata(path)
//
//  Returns a single-row table containing basic metadata about a raster.  This
//  includes its dimensions (width and height), the number of bands, the
//  top-left corner (origin) of the dataset, the pixel size in the X/Y
//  directions, and the spatial reference system in WKT.  This table function
//  can be used to populate a metadata table when loading rasters into DuckDB.

struct MetaState : public GlobalTableFunctionState {
  bool emitted = false;
  int64_t width = 0;
  int64_t height = 0;
  int64_t bands = 0;
  double origin_x = 0.0;
  double origin_y = 0.0;
  double pixel_x = 0.0;
  double pixel_y = 0.0;
  std::string srs;
  idx_t MaxThreads() const override { return 1; }
};

// The metadata table function requires its own bind data because
// the BindData struct used by read_raster (defined in raster.cpp) is
// not visible here.  MetaBindData simply stores the path to the
// raster file provided by the user.  During binding, we also set
// the output schema: eight columns describing width/height/bands,
// origin, pixel size and spatial reference system.
struct MetaBindData : public FunctionData {
  std::string path;
  unique_ptr<FunctionData> Copy() const override {
    return make_uniq<MetaBindData>(*this);
  }
  bool Equals(const FunctionData &other) const override {
    auto &o = other.Cast<MetaBindData>();
    return path == o.path;
  }
};

// Bind: define output schema and record the raster path
static unique_ptr<FunctionData> MetaBind(ClientContext &,
                                         TableFunctionBindInput &input,
                                         vector<LogicalType> &types,
                                         vector<string> &names) {
  if (input.inputs.empty()) {
    throw BinderException("raster_metadata(path) requires a file path");
  }
  // Define the column types and names for the metadata table
  types = {
      LogicalType::BIGINT, // width
      LogicalType::BIGINT, // height
      LogicalType::BIGINT, // bands
      LogicalType::DOUBLE, // origin_x
      LogicalType::DOUBLE, // origin_y
      LogicalType::DOUBLE, // pixel_x
      LogicalType::DOUBLE, // pixel_y
      LogicalType::VARCHAR // srs
  };
  names = {"width",    "height",  "bands",   "origin_x",
           "origin_y", "pixel_x", "pixel_y", "srs"};
  // Save the path into our custom bind data
  auto bd = make_uniq<MetaBindData>();
  bd->path = StringValue::Get(input.inputs[0]);
  return std::move(bd);
}

// Init: open dataset and record metadata
static unique_ptr<GlobalTableFunctionState>
MetaInit(ClientContext &, TableFunctionInitInput &in) {
  // Retrieve our custom bind data to get the path
  auto &bd = in.bind_data->Cast<MetaBindData>();
  auto st = make_uniq<MetaState>();
  try {
    auto ds = OpenDataset(bd.path);
    st->width = ds->GetRasterXSize();
    st->height = ds->GetRasterYSize();
    st->bands = ds->GetRasterCount();
    double gt[6] = {0};
    if (ds->GetGeoTransform(gt) == CE_None) {
      st->origin_x = gt[0];
      st->origin_y = gt[3];
      st->pixel_x = gt[1];
      st->pixel_y = gt[5];
    }
    const char *proj = ds->GetProjectionRef();
    if (proj && proj[0] != '\0') {
      st->srs = std::string(proj);
    }
  } catch (...) {
    // If anything goes wrong opening the dataset, leave all values at defaults.
  }
  return std::move(st);
}

// Scan: emit one row then finish
static void MetaScan(ClientContext &, TableFunctionInput &in, DataChunk &out) {
  auto &st = in.global_state->Cast<MetaState>();
  if (st.emitted) {
    out.SetCardinality(0);
    return;
  }
  auto r_width = FlatVector::GetData<int64_t>(out.data[0]);
  auto r_height = FlatVector::GetData<int64_t>(out.data[1]);
  auto r_bands = FlatVector::GetData<int64_t>(out.data[2]);
  auto r_origin_x = FlatVector::GetData<double>(out.data[3]);
  auto r_origin_y = FlatVector::GetData<double>(out.data[4]);
  auto r_pixel_x = FlatVector::GetData<double>(out.data[5]);
  auto r_pixel_y = FlatVector::GetData<double>(out.data[6]);
  auto &r_srs = out.data[7];
  r_width[0] = st.width;
  r_height[0] = st.height;
  r_bands[0] = st.bands;
  r_origin_x[0] = st.origin_x;
  r_origin_y[0] = st.origin_y;
  r_pixel_x[0] = st.pixel_x;
  r_pixel_y[0] = st.pixel_y;
  r_srs.SetValue(0, st.srs.empty() ? Value() : Value(st.srs));
  out.SetCardinality(1);
  st.emitted = true;
}

// =============================================================================
//  Scalar Function: ndvi(nir, red)
//
//  Computes the normalized difference vegetation index from two bands.
static void NdviFunction(DataChunk &input, ExpressionState &, Vector &result) {
  auto count = input.size();
  result.SetVectorType(VectorType::FLAT_VECTOR);
  for (idx_t row = 0; row < count; row++) {
    Value nir_val = input.data[0].GetValue(row);
    Value red_val = input.data[1].GetValue(row);
    if (nir_val.IsNull() || red_val.IsNull()) {
      result.SetValue(row, Value());
      continue;
    }
    double nir = nir_val.GetValue<double>();
    double red = red_val.GetValue<double>();
    double denom = nir + red;
    if (denom == 0.0) {
      result.SetValue(row, Value());
    } else {
      double ndvi = (nir - red) / denom;
      result.SetValue(row, Value::DOUBLE(ndvi));
    }
  }
}

// =============================================================================
//  Scalar Function: trend(values)
//
//  Computes the slope of a simple linear regression through a list of values.
//  The list is interpreted as a time series at regular intervals; the x values
//  are 0..n-1.  Nulls are skipped.  Returns NULL if fewer than two
//  non-null observations.
static void TrendFunction(DataChunk &input, ExpressionState &, Vector &result) {
  auto count = input.size();
  result.SetVectorType(VectorType::FLAT_VECTOR);
  for (idx_t row = 0; row < count; row++) {
    Value list_val = input.data[0].GetValue(row);
    if (list_val.IsNull()) {
      result.SetValue(row, Value());
      continue;
    }
    auto &children = ListValue::GetChildren(list_val);
    size_t n = children.size();
    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_x2 = 0.0;
    idx_t valid_count = 0;
    for (idx_t i = 0; i < n; i++) {
      auto &c = children[i];
      if (c.IsNull())
        continue;
      double y = c.GetValue<double>();
      double x = (double)valid_count; // treat sequential index for valid values
      sum_x += x;
      sum_y += y;
      sum_xy += x * y;
      sum_x2 += x * x;
      valid_count++;
    }
    if (valid_count < 2) {
      result.SetValue(row, Value());
      continue;
    }
    double n_d = (double)valid_count;
    double num = n_d * sum_xy - sum_x * sum_y;
    double den = n_d * sum_x2 - sum_x * sum_x;
    if (den == 0.0) {
      result.SetValue(row, Value());
    } else {
      double slope = num / den;
      result.SetValue(row, Value::DOUBLE(slope));
    }
  }
}

} // anonymous namespace

namespace duckdb {

void RegisterRasterFunctions(DatabaseInstance &db) {
  // width
  ScalarFunction width_fn("raster_width", {LogicalType::VARCHAR},
                          LogicalType::BIGINT, RasterWidthFunction);
  ExtensionUtil::RegisterFunction(db, width_fn);

  // height
  ScalarFunction height_fn("raster_height", {LogicalType::VARCHAR},
                           LogicalType::BIGINT, RasterHeightFunction);
  ExtensionUtil::RegisterFunction(db, height_fn);

  // band count
  ScalarFunction band_count_fn("raster_band_count", {LogicalType::VARCHAR},
                               LogicalType::BIGINT, RasterBandCountFunction);
  ExtensionUtil::RegisterFunction(db, band_count_fn);

  // min
  ScalarFunction min_fn("raster_band_min",
                        {LogicalType::VARCHAR, LogicalType::INTEGER},
                        LogicalType::DOUBLE, RasterMinFunction);
  ExtensionUtil::RegisterFunction(db, min_fn);

  // max
  ScalarFunction max_fn("raster_band_max",
                        {LogicalType::VARCHAR, LogicalType::INTEGER},
                        LogicalType::DOUBLE, RasterMaxFunction);
  ExtensionUtil::RegisterFunction(db, max_fn);

  // mean
  ScalarFunction mean_fn("raster_band_mean",
                         {LogicalType::VARCHAR, LogicalType::INTEGER},
                         LogicalType::DOUBLE, RasterMeanFunction);
  ExtensionUtil::RegisterFunction(db, mean_fn);

  // spatial reference system (projection) as WKT
  ScalarFunction srs_fn("raster_srs", {LogicalType::VARCHAR},
                        LogicalType::VARCHAR, RasterSRSFunction);
  ExtensionUtil::RegisterFunction(db, srs_fn);

  // metadata table function
  TableFunction meta_tf("raster_metadata", {LogicalType::VARCHAR}, MetaScan,
                        MetaBind, MetaInit);
  ExtensionUtil::RegisterFunction(db, meta_tf);

  // NDVI
  ScalarFunction ndvi_fn("ndvi", {LogicalType::DOUBLE, LogicalType::DOUBLE},
                         LogicalType::DOUBLE, NdviFunction);
  ExtensionUtil::RegisterFunction(db, ndvi_fn);

  // trend over a list of doubles
  ScalarFunction trend_fn("trend", {LogicalType::LIST(LogicalType::DOUBLE)},
                          LogicalType::DOUBLE, TrendFunction);
  ExtensionUtil::RegisterFunction(db, trend_fn);

  // pixel value retrieval: raster_pixel_value(path, band, row, col) -> DOUBLE
  ScalarFunction pixel_fn("raster_pixel_value",
                          {LogicalType::VARCHAR, LogicalType::INTEGER,
                           LogicalType::BIGINT, LogicalType::BIGINT},
                          LogicalType::DOUBLE, RasterPixelValueFunction);
  ExtensionUtil::RegisterFunction(db, pixel_fn);

  // cell_id conversion helpers
  ScalarFunction cell_row_fn("raster_cell_row",
                             {LogicalType::BIGINT, LogicalType::BIGINT},
                             LogicalType::BIGINT, RasterCellRowFunction);
  ExtensionUtil::RegisterFunction(db, cell_row_fn);
  ScalarFunction cell_col_fn("raster_cell_col",
                             {LogicalType::BIGINT, LogicalType::BIGINT},
                             LogicalType::BIGINT, RasterCellColFunction);
  ExtensionUtil::RegisterFunction(db, cell_col_fn);

  // temporal mean: mean of a list of doubles
  ScalarFunction tmean_fn("raster_temporal_mean",
                          {LogicalType::LIST(LogicalType::DOUBLE)},
                          LogicalType::DOUBLE, RasterTemporalMeanFunction);
  ExtensionUtil::RegisterFunction(db, tmean_fn);

  // band sum: sum of all values in a band
  ScalarFunction bsum_fn("raster_band_sum",
                         {LogicalType::VARCHAR, LogicalType::INTEGER},
                         LogicalType::DOUBLE, RasterBandSumFunction);
  ExtensionUtil::RegisterFunction(db, bsum_fn);

  // read_raster_cube table function
  TableFunction cube_tf("read_raster_cube",
                        {LogicalType::ANY, LogicalType::ANY, LogicalType::ANY},
                        CubeScan, CubeBind, CubeInit);
  ExtensionUtil::RegisterFunction(db, cube_tf);
}

} // namespace duckdb
