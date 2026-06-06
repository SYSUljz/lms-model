#pragma once
// Factory that reads a directory of GeoTIFF tiles and assembles the per-cell
// ConstParam grid (ConstRaster<T>) consumed by the model.
//
// This I/O layer is deliberately decoupled from the (still-WIP) model/algorithm
// files: it depends only on core.h data structures + GDAL, so it compiles and
// runs on its own and you can inspect what was read before the solver exists.

#include <array>
#include <filesystem>
#include <queue>
#include <stack>
#include <stdexcept>
#include <string>

#include "base/core.h"
#include "base/direct.h"
#include "base/raster.h"
#include "io/gdal_reader.h"

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

template <typename T>
class ModelBuilder {
  std::filesystem::path dir_;

 public:
  ModelBuilder& from_directory(const std::string& dir) {
    dir_ = dir;
    return *this;
  }

  /// Read every required tile and build the per-cell ConstParam grid.
  ConstRaster<T> build_const_raster() const {
    if (dir_.empty()) {
      throw std::runtime_error("ModelBuilder: model directory not set");
    }

    auto load = [&](const char* name) { return ReadBand<T>((dir_ / name).string()); };

    // topology / classification
    Raster<T> label = load("label.tif");
    Raster<T> d8 = load("d8.tif");
    Raster<T> slope = load("slope.tif");
    // soil parameters
    Raster<T> sat = load("sat.tif");
    Raster<T> fc = load("fc.tif");
    Raster<T> wl = load("wl.tif");
    Raster<T> zs = load("zs.tif");
    Raster<T> ks = load("ks.tif");
    Raster<T> b = load("b.tif");
    Raster<T> n = load("n.tif");
    Raster<T> v = load("v.tif");
    // channel parameters
    Raster<T> bs = load("bs.tif");
    Raster<T> bw = load("bw.tif");
    Raster<T> manning = load("manning.tif");

    const int W = label.width;
    const int H = label.height;

    const std::array<const Raster<T>*, 14> all = {&label, &d8, &slope, &sat, &fc, &wl, &zs,
                                                  &ks,    &b,  &n,     &v,   &bs, &bw, &manning};
    for (const Raster<T>* r : all) {
      if (r->width != W || r->height != H) {
        throw std::runtime_error("ModelBuilder: raster dimension mismatch");
      }
    }

    ConstRaster<T> out;
    out.meta.width_ = static_cast<std::size_t>(W);
    out.meta.heigh_ = static_cast<std::size_t>(H);
    out.meta.cell_size_ = static_cast<std::size_t>(label.cell_size);

    const std::size_t n_cells = static_cast<std::size_t>(W) * static_cast<std::size_t>(H);
    out.cells.reset(new std::vector<ConstParam<T>>(n_cells));
    out.active.assign(n_cells, 0);

    for (std::size_t i = 0; i < n_cells; ++i) {
      // A cell is part of the basin iff its label is not NoData.
      if (label.is_nodata_at(i)) {
        continue;
      }
      out.active[i] = 1;

      ConstParam<T>& c = (*out.cells)[i];
      c.label = LabelFromRaster(static_cast<double>(label.data[i]));
      c.d8 = D8FromRaster(static_cast<double>(d8.data[i]));
      c.slop = slope.data[i];
      c.sat = sat.data[i];
      c.fc = fc.data[i];
      c.wl = wl.data[i];
      c.zs = zs.data[i];
      c.ks = ks.data[i];
      c.b = b.data[i];
      c.n = n.data[i];
      c.v = v.data[i];
      c.bs = bs.data[i];
      c.bw = bw.data[i];
      c.manning = manning.data[i];
      c.ep = T {0};  // no potential-evapotranspiration tile yet
    }

    return out;
  }
  [[nodiscard]] std::vector<int> BuildOrder() {
    auto load = [&](const char* name) { return ReadBand<T>((dir_ / name).string()); };
    Raster<int> d8 = load("d8.tif");
    auto width = d8.width;
    auto height = d8.height;

    std::vector<int> outlets;

    for (size_t i = 0; i < d8.data.size(); ++i) {
      int dir = d8.data[i];

      if (dir != 1 && dir != 2 && dir != 4 && dir != 8 && dir != 16 && dir != 32 && dir != 64 && dir != 128) {
        outlets.push_back(i);
        continue;
      }

      int x = i % width;
      int y = i / width;
      int nx = x, ny = y;

      switch (dir) {
        case 1:
          nx++;
          break;  // E
        case 2:
          nx++;
          ny++;
          break;  // SE
        case 4:
          ny++;
          break;  // S
        case 8:
          nx--;
          ny++;
          break;  // SW
        case 16:
          nx--;
          break;  // W
        case 32:
          nx--;
          ny--;
          break;  // NW
        case 64:
          ny--;
          break;  // N
        case 128:
          nx++;
          ny--;
          break;  // NE
      }

      if (nx < 0 || nx >= width || ny < 0 || ny >= height || d8.is_nodata_at(nx + width * ny)) {
        outlets.push_back(i);
      }
    }

    if (outlets.size() > 1) {
      throw std::runtime_error("BuildOrder: more than one outlet found");
    }

    std::vector<int> order;
    if (outlets.empty()) {
      return order;
    }

    // Map a neighbour offset (dx, dy) to the D8 value that would cause that
    // neighbour to flow INTO the current cell.  Index: [dy+1][dx+1].
    // e.g. neighbour at (dx=-1, dy=-1) — upper-left — flows into us iff its
    // D8 == 2 (SE); neighbour at (dx=+1, dy=0) — right — flows in iff D8 == 16 (W).
    static const int up_d8[3][3] = {
        {2,  4,  8},     // dy=-1  (upper-left, above, upper-right)
        {1,  0,  16},    // dy= 0  (left,      self,  right)
        {128, 64, 32}    // dy=+1  (lower-left, below, lower-right)
    };

    // BFS upstream from the single outlet.
    // Push every visited cell onto a stack; popping yields upstream→downstream
    // order so that all inflows to a cell are computed before the cell itself.
    std::queue<int> q;
    std::vector<bool> visited(d8.data.size(), false);
    std::stack<int> stk;

    q.push(outlets[0]);
    visited[outlets[0]] = true;

    while (!q.empty()) {
      int front = q.front();
      stk.push(front);
      q.pop();

      int x = front % width;
      int y = front / width;

      // Examine all 8 neighbours; if a neighbour's D8 direction points INTO
      // (x, y), that neighbour is immediately upstream and gets enqueued.
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
          if (dx == 0 && dy == 0) continue;

          int nx = x + dx;
          int ny = y + dy;
          if (nx < 0 || nx >= width || ny < 0 || ny >= height) continue;

          int nidx = nx + width * ny;
          if (visited[nidx]) continue;
          if (d8.is_nodata_at(nidx)) continue;

          int ndir = static_cast<int>(d8.data[nidx]);
          if (ndir != up_d8[dy + 1][dx + 1]) continue;

          visited[nidx] = true;
          q.push(nidx);
        }
      }
    }

    // Pop from stack: most-upstream cells first, outlet last.
    order.reserve(stk.size());
    while (!stk.empty()) {
      order.push_back(stk.top());
      stk.pop();
    }

    return order;
  }
  model<T> BuildModel() {}

 private:
  unique_ptr<> BuildStateRaster() {}
  void BuildStreamOrder() {}
};
