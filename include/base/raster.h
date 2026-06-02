#pragma once
#include <cstddef>
#include <memory>
#include <vector>

#include "base/core.h"
namespace lms {
namespace raster {
/// A single 2D raster band loaded into memory (row-major).
/// Decoupled from GDAL: any reader can populate this.
template <typename T>
struct Raster {
  int width = 0;
  int height = 0;
  // pixel size in CRS units (m), taken from the geotransform
  double cell_size = 0.0;
  double nodata = 0.0;
  bool has_nodata = false;
  // row-major, size == width * height
  std::vector<T> data;

  std::size_t size() const { return data.size(); }

  T& at(int row, int col) { return data[static_cast<std::size_t>(row) * width + col]; }
  const T& at(int row, int col) const { return data[static_cast<std::size_t>(row) * width + col]; }

  bool is_nodata_at(std::size_t i) const { return has_nodata && static_cast<double>(data[i]) == nodata; }
};

/// The const-parameter grid consumed by the model: per-cell ConstParam plus
/// grid metadata and an active mask (cells inside the basin / not NoData).
template <typename T>
struct ConstRaster {
  ModelMeta meta {};
  // row-major, size == meta.width_ * meta.heigh_
  std::unique_ptr<std::vector<ConstParam<T>>> cells;
  // 1 if the cell is inside the basin (label is not NoData), else 0
  std::vector<char> active;

  const ConstParam<T>& operator[](std::size_t idx) const { return (*cells)[idx]; }
};
template <typename T>
struct StateRaster {
  ModelMeta meta_;
  // row-major, size == meta.width_ * meta.heigh_
  std::unique_ptr<std::vector<StateParam<T>>> cells;
  // 1 if the cell is inside the basin (label is not NoData), else 0
  std::shared_ptr<std::vector<char>> active_;

  StateRaster(ModelMeta meta, std::shared_ptr<std::vector<char>> active)
      : meta_(meta), active_(active), cells(std::nullptr_t) {
    std::size_t n_cells = meta_.width_ * meta_.heigh_;
    cells.reset(new std::vector<StateParam<T>>(n_cells));
  };
  StateRaster(const StateRaster& stateraster) = delete;
  StateRaster& operator=(const StateRaster& stateraster) = delete;

  StateParam<T>& operator[](std::size_t idx) const { return (*cells)[idx]; }
};
}  // namespace raster
}  // namespace lms