#pragma once
#include <cstddef>
#include <memory>
#include <vector>

#include "base/core.h"
namespace lms {
namespace raster {

using lms::core::ConstParam;
using lms::core::Label;
using lms::core::ModelMeta;
using lms::core::StateParam;
using lms::direct::Direct8;

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
  ModelMeta<T> meta {};
  // row-major, size == meta.width_ * meta.heigh_
  std::vector<ConstParam<T>> cells;
  // 1 if the cell is inside the basin (label is not NoData), else 0
  std::vector<char> active;

  ConstRaster() = default;
  ConstRaster(const ConstRaster&) = default;
  ConstRaster& operator=(const ConstRaster&) = default;
  ConstRaster(ConstRaster&&) = default;
  ConstRaster& operator=(ConstRaster&&) = default;

  ConstParam<T>& operator[](std::size_t idx) { return cells[idx]; }
  const ConstParam<T>& operator[](std::size_t idx) const { return cells[idx]; }
  size_t size() const { return cells.size(); }
};
template <typename T>
struct StateRaster {
  ModelMeta<T> meta_;
  // row-major, size == meta.width_ * meta.heigh_
  std::vector<StateParam<T>> cells;
  // 1 if the cell is inside the basin (label is not NoData), else 0
  std::vector<char> active_;

  StateRaster(ModelMeta<T> meta, std::vector<char> active) : meta_(meta), active_(active) {
    std::size_t n_cells = meta_.width_ * meta_.heigh_;
    cells = std::vector<StateParam<T>>(n_cells);
  };
  StateRaster(const StateRaster& stateraster) = default;
  StateRaster& operator=(const StateRaster& stateraster) = default;
  StateRaster(StateRaster&&) = default;
  StateRaster& operator=(StateRaster&&) = default;

  StateParam<T>& operator[](std::size_t idx) { return cells[idx]; }
  const StateParam<T>& operator[](std::size_t idx) const { return cells[idx]; }

  size_t size() const { return cells.size(); }
};

// ---------------------------------------------------------------------------
// Raster-value -> enum decoders. These encode domain assumptions about how the
// Meizhou tiles are coded; verify them against the histograms printed by main.
// ---------------------------------------------------------------------------

/// label.tif encoding (assumed): 1=Soil, 2=Channel, 3=Reservoir,
/// 4=Channel(outlet). NOTE: label.tif here ranges 1..4 while Label has 3
/// categories — value 4 is provisionally folded into Channel. Confirm via the
/// raw-value histogram.
inline Label LabelFromRaster(double v) {
  switch (static_cast<int>(v)) {
    case 1:
      return Label::Soil;
    case 2:
      return Label::Channel;
    case 3:
      return Label::Reservoir;
    case 4:
      return Label::Channel;
    default:
      return Label::Soil;
  }
}

/// ESRI D8 encoding: powers of two clockwise from East.
/// 1=E, 2=SE, 4=S, 8=SW, 16=W, 32=NW, 64=N, 128=NE — which maps to Direct8 by
/// log2.
inline Direct8 D8FromRaster(double v) {
  switch (static_cast<int>(v)) {
    case 1:
      return Direct8::Right;
    case 2:
      return Direct8::DownRight;
    case 4:
      return Direct8::Down;
    case 8:
      return Direct8::DownLeft;
    case 16:
      return Direct8::Left;
    case 32:
      return Direct8::UpLeft;
    case 64:
      return Direct8::Up;
    case 128:
      return Direct8::UpRight;
    default:
      // NoData, sink, or outlet with no defined direction
      return Direct8::Right;
  }
}
}  // namespace raster
}  // namespace lms
