#pragma once
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "base/raster.h"
#include "cpl_conv.h"
#include "gdal_priv.h"

/// Register GDAL drivers exactly once per process.
inline void EnsureGdalRegistered() {
  static const bool done = [] {
    GDALAllRegister();
    return true;
  }();
  (void)done;
}

/// Read band 1 of a GeoTIFF into a Raster<T>.
/// Values are read as Float64 then cast to T, so T may be float or double.
/// Throws std::runtime_error on any GDAL failure.
template <typename T>
Raster<T> ReadBand(const std::string& path) {
  EnsureGdalRegistered();

  GDALDataset* ds = static_cast<GDALDataset*>(GDALOpen(path.c_str(), GA_ReadOnly));
  if (ds == nullptr) {
    throw std::runtime_error("GDAL: cannot open " + path);
  }

  Raster<T> r;
  r.width = ds->GetRasterXSize();
  r.height = ds->GetRasterYSize();

  double gt[6];
  if (ds->GetGeoTransform(gt) == CE_None) {
    // gt[1] is pixel width; magnitude is the cell size in CRS units
    r.cell_size = std::abs(gt[1]);
  }

  GDALRasterBand* band = ds->GetRasterBand(1);
  int has_nd = 0;
  const double nd = band->GetNoDataValue(&has_nd);
  r.has_nodata = has_nd != 0;
  r.nodata = nd;

  const std::size_t n = static_cast<std::size_t>(r.width) * static_cast<std::size_t>(r.height);
  std::vector<double> buf(n);
  const CPLErr err = band->RasterIO(GF_Read, 0, 0, r.width, r.height, buf.data(), r.width, r.height, GDT_Float64, 0, 0);
  GDALClose(ds);

  if (err != CE_None) {
    throw std::runtime_error("GDAL: read failed for " + path);
  }

  r.data.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    r.data[i] = static_cast<T>(buf[i]);
  }
  return r;
}
