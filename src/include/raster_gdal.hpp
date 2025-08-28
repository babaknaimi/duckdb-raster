#pragma once
#include <memory>
#include <string>

#include "duckdb/common/exception.hpp"
#include "gdal_priv.h"

struct GDALDatasetDeleter {
  void operator()(GDALDataset *p) const {
    if (p)
      GDALClose(p);
  }
};

using GDALDatasetPtr = std::unique_ptr<GDALDataset, GDALDatasetDeleter>;

inline GDALDatasetPtr OpenDataset(const std::string &path) {
  GDALAllRegister();
  GDALDataset *raw =
      static_cast<GDALDataset *>(GDALOpen(path.c_str(), GA_ReadOnly));
  if (!raw) {
    throw duckdb::IOException("GDALOpen failed for '%s'", path.c_str());
  }
  return GDALDatasetPtr(raw);
}
