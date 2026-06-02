#include <filesystem>
#include <stdexcept>

#include "base/rain.h"
template <typename T>
class RainBuilder {
  std::filesystem::path dir_;

 public:
  RainBuilder& from_directory(const std::string& dir) {
    dir_ = dir;
    return *this;
  }

  /// Read a .
  ConstRaster<T> build_rain_csv() const {
    if (dir_.empty()) {
      throw std::runtime_error("RainBuilder: model directory not set");
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
        throw std::runtime_error("RainBuilder: raster dimension mismatch");
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
  model<T> BuildRain() {}

 private:
  unique_ptr<> BuildStateRaster() {}
  void BuildStreamOrder() {}
};
