// Driver to exercise the GDAL I/O layer against real data and print a summary,
// so raster input can be inspected without the (WIP) solver.
//
//   ./lms_io_test [model_dir]   (default: data/meizhou/model)

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <string>

#include "base/model.h"
#include "base/rain.h"
#include "builder/raster_builder.h"
#include "io/gdal_reader.h"

using T = double;
using lms::core::Label;
using lms::core::ConstParam;
using lms::raster::Raster;
using lms::raster::ConstRaster;
static const char* LabelName(Label l) {
  switch (l) {
    case Label::Soil:
      return "Soil";
    case Label::Channel:
      return "Channel";
    case Label::Reservoir:
      return "Reservoir";
  }
  return "?";
}

// Stats over active cells, optionally restricted to a single Label (used for
// channel-only parameters, which are NoData on soil cells).
template <typename Acc>
static void PrintFieldStats(const char* name, const ConstRaster<T>& cr, Acc acc,
                            std::optional<Label> only = std::nullopt) {
  double mn = std::numeric_limits<double>::max();
  double mx = std::numeric_limits<double>::lowest();
  double sum = 0.0;
  std::size_t cnt = 0;
  for (std::size_t i = 0; i < (*cr.cells).size(); ++i) {
    if (!cr.active[i]) continue;
    if (only && (*cr.cells)[i].label != *only) continue;
    const double val = static_cast<double>(acc((*cr.cells)[i]));
    mn = std::min(mn, val);
    mx = std::max(mx, val);
    sum += val;
    ++cnt;
  }
  if (cnt == 0) {
    std::printf("  %-10s : (no matching cells)\n", name);
    return;
  }
  std::printf("  %-10s : min=%-12.5g max=%-12.5g mean=%-12.5g\n", name, mn, mx, sum / cnt);
}

int main(int argc, char** argv) {
  const std::string dir = (argc > 1) ? argv[1] : "data/meizhou/model";

  try {
    // Raw label histogram: verify the label.tif -> Label mapping assumption.
    Raster<double> raw_label = ReadBand<double>(dir + "/label.tif");
    std::map<int, std::size_t> raw_hist;
    for (std::size_t i = 0; i < raw_label.data.size(); ++i) {
      if (raw_label.is_nodata_at(i)) continue;
      raw_hist[static_cast<int>(raw_label.data[i])]++;
    }

    ModelBuilder<T> builder;
    const ConstRaster<T> cr = builder.from_directory(dir).build_const_raster();

    std::size_t active = 0;
    std::map<int, std::size_t> label_counts;
    for (std::size_t i = 0; i < (*cr.cells).size(); ++i) {
      if (!cr.active[i]) continue;
      ++active;
      label_counts[static_cast<int>((*cr.cells)[i].label)]++;
    }

    std::printf("== Model I/O summary (%s) ==\n", dir.c_str());
    std::printf("grid        : %zu x %zu  (cell_size = %zu m)\n", cr.meta.width_, cr.meta.heigh_, cr.meta.cell_size_);
    std::printf("cells       : %zu total, %zu active (in basin)\n", (*cr.cells).size(), active);

    std::printf("raw label.tif values:\n");
    for (const auto& [val, c] : raw_hist) std::printf("  value %d : %zu cells\n", val, c);

    std::printf("decoded labels:\n");
    for (const auto& [val, c] : label_counts)
      std::printf("  %-10s : %zu cells\n", LabelName(static_cast<Label>(val)), c);

    std::printf("field stats over active cells:\n");
    PrintFieldStats("slop", cr, [](const ConstParam<T>& c) { return c.slop; });
    PrintFieldStats("ks", cr, [](const ConstParam<T>& c) { return c.ks; });
    PrintFieldStats("sat", cr, [](const ConstParam<T>& c) { return c.sat; });
    PrintFieldStats("fc", cr, [](const ConstParam<T>& c) { return c.fc; });
    PrintFieldStats("wl", cr, [](const ConstParam<T>& c) { return c.wl; });
    PrintFieldStats("zs", cr, [](const ConstParam<T>& c) { return c.zs; });
    PrintFieldStats("b", cr, [](const ConstParam<T>& c) { return c.b; });
    PrintFieldStats("n", cr, [](const ConstParam<T>& c) { return c.n; });
    PrintFieldStats("v", cr, [](const ConstParam<T>& c) { return c.v; });

    std::printf("channel-only params (over Channel cells):\n");
    PrintFieldStats("bs", cr, [](const ConstParam<T>& c) { return c.bs; }, Label::Channel);
    PrintFieldStats("bw", cr, [](const ConstParam<T>& c) { return c.bw; }, Label::Channel);
    PrintFieldStats("manning", cr, [](const ConstParam<T>& c) { return c.manning; }, Label::Channel);

    // Dump one active cell so the assembled ConstParam can be eyeballed.
    for (std::size_t i = 0; i < (*cr.cells).size(); ++i) {
      if (!cr.active[i]) continue;
      const ConstParam<T>& c = (*cr.cells)[i];
      std::printf(
          "sample active cell #%zu: label=%s d8=%d slop=%g ks=%g sat=%g fc=%g "
          "bw=%g manning=%g\n",
          i, LabelName(c.label), static_cast<int>(c.d8), static_cast<double>(c.slop), static_cast<double>(c.ks),
          static_cast<double>(c.sat), static_cast<double>(c.fc), static_cast<double>(c.bw),
          static_cast<double>(c.manning));
      break;
    }

    // =======================================================================
    // Section 2: BuildStateRaster test
    // =======================================================================
    std::printf("\n== BuildStateRaster test ==\n");

    auto sr = builder.BuildStateRaster(cr);
    assert(sr.meta_.width_ == cr.meta.width_);
    assert(sr.meta_.heigh_ == cr.meta.heigh_);
    assert(sr.meta_.cell_size_ == cr.meta.cell_size_);
    std::printf("  meta match: %zu x %zu (cell_size=%zu m)  PASSED\n",
                sr.meta_.width_, sr.meta_.heigh_, sr.meta_.cell_size_);

    // Verify active mask is copied correctly
    std::size_t sr_active_cnt = 0;
    for (std::size_t i = 0; i < (*sr.cells).size(); ++i) {
      if ((*sr.active_)[i]) ++sr_active_cnt;
    }
    assert(sr_active_cnt == active);
    std::printf("  active mask: %zu active cells  PASSED\n", sr_active_cnt);

    // Verify StateParam default-initialization: all fields should be zero
    {
      bool all_zero = true;
      for (std::size_t i = 0; i < (*sr.cells).size(); ++i) {
        if (!(*sr.active_)[i]) continue;
        const auto& sp = (*sr.cells)[i];
        if (sp.soil_moisture != T{0} || sp.actual_evaporate != T{0} || sp.runoff != T{0} ||
            sp.lat_mm != T{0} || sp.per_mm != T{0} || sp.groundwater_mm != T{0} ||
            sp.prev_t_flow != T{0} || sp.upstream_in_flow != T{0} || sp.water_level != T{0} ||
            sp.temp != T{0} || sp.lateral_in_flow_mm != T{0}) {
          all_zero = false;
          break;
        }
      }
      assert(all_zero);
      std::printf("  StateParam default-init: all fields zero  PASSED\n");
    }

    // =======================================================================
    // Section 3: BuildOrder test
    // =======================================================================
    std::printf("\n== BuildOrder test ==\n");

    std::vector<int> order = builder.BuildOrder();
    assert(order.size() == active);
    std::printf("  order size == active cells (%zu)  PASSED\n", order.size());

    // Every cell in order must be active
    {
      bool all_active = true;
      for (int idx : order) {
        if (!cr.active[static_cast<std::size_t>(idx)]) {
          all_active = false;
          break;
        }
      }
      assert(all_active);
      std::printf("  all ordered cells are active  PASSED\n");
    }

    // No duplicates in order (topological sort guarantees this)
    {
      std::vector<bool> seen((*cr.cells).size(), false);
      bool no_dup = true;
      for (int idx : order) {
        if (seen[static_cast<std::size_t>(idx)]) { no_dup = false; break; }
        seen[static_cast<std::size_t>(idx)] = true;
      }
      assert(no_dup);
      std::printf("  no duplicates in order  PASSED\n");
    }

    // Print first/last few cells in the order for inspection
    std::printf("  order[0..4]:   ");
    for (std::size_t i = 0; i < std::min(std::size_t{5}, order.size()); ++i)
      std::printf("%d ", order[i]);
    std::printf("\n  order[last-4..]: ");
    for (std::size_t i = (order.size() > 5 ? order.size() - 5 : 0); i < order.size(); ++i)
      std::printf("%d ", order[i]);
    std::printf("\n");

    // =======================================================================
    // Section 4: BuildTarget test
    // =======================================================================
    std::printf("\n== BuildTarget test ==\n");

    auto target = builder.BuildTarget(order);
    assert(target.size() == (*cr.cells).size());
    std::printf("  target vector size == n_cells (%zu)  PASSED\n", target.size());

    // Outlet cells should have no target (nullopt)
    // Non-outlet cells in the order should have a valid target
    {
      std::size_t with_target = 0, no_target = 0;
      for (int idx : order) {
        if (target[static_cast<std::size_t>(idx)].has_value()) {
          ++with_target;
          auto tgt = *target[static_cast<std::size_t>(idx)];
          assert(tgt >= 0);
          assert(static_cast<std::size_t>(tgt) < (*cr.cells).size());
        } else {
          ++no_target;
        }
      }
      assert(no_target >= 1);  // at least one outlet (may be more)
      std::printf("  cells with target: %zu, without target (outlet): %zu  PASSED\n",
                  with_target, no_target);
    }

    // =======================================================================
    // Section 5: Station building test
    // =======================================================================
    std::printf("\n== Station test ==\n");

    using lms::rain::Station;
    using lms::core::ModelMeta;

    // Set up meta with known geo bounds for deterministic raster position test
    ModelMeta<T> station_meta{};
    station_meta.width_ = 100;
    station_meta.heigh_ = 200;
    station_meta.cell_size_ = 30;
    station_meta.raster_min_lat_ = static_cast<T>(20.0);
    station_meta.raster_max_lat_ = static_cast<T>(30.0);
    station_meta.raster_min_long_ = static_cast<T>(110.0);
    station_meta.raster_max_long_ = static_cast<T>(120.0);

    // Create stations at known locations
    Station<T> s_center(1, static_cast<T>(25.0), static_cast<T>(115.0), "center", station_meta);
    Station<T> s_corner(2, static_cast<T>(20.0), static_cast<T>(110.0), "sw_corner", station_meta);
    Station<T> s_northeast(3, static_cast<T>(30.0), static_cast<T>(120.0), "ne_corner", station_meta);

    s_center.GeoPos2RasterPos(station_meta);
    s_corner.GeoPos2RasterPos(station_meta);
    s_northeast.GeoPos2RasterPos(station_meta);

    std::printf("  Station '%s': geo(%.4f, %.4f) -> raster(%zu, %zu)\n",
                s_center.name_.c_str(), static_cast<double>(s_center.geo_lat_),
                static_cast<double>(s_center.geo_long_), s_center.raster_x_, s_center.raster_y_);
    std::printf("  Station '%s': geo(%.4f, %.4f) -> raster(%zu, %zu)\n",
                s_corner.name_.c_str(), static_cast<double>(s_corner.geo_lat_),
                static_cast<double>(s_corner.geo_long_), s_corner.raster_x_, s_corner.raster_y_);
    std::printf("  Station '%s': geo(%.4f, %.4f) -> raster(%zu, %zu)\n",
                s_northeast.name_.c_str(), static_cast<double>(s_northeast.geo_lat_),
                static_cast<double>(s_northeast.geo_long_), s_northeast.raster_x_, s_northeast.raster_y_);

    // Verify raster positions
    // center: (115-110)/(120-110)*100 = 50, (25-20)/(30-20)*200 = 100
    assert(s_center.raster_x_ == 50);
    assert(s_center.raster_y_ == 100);
    // sw_corner: (110-110)/(120-110)*100 = 0, (20-20)/(30-20)*200 = 0
    assert(s_corner.raster_x_ == 0);
    assert(s_corner.raster_y_ == 0);
    // ne_corner: (120-110)/(120-110)*100 = 100, (30-20)/(30-20)*200 = 200
    assert(s_northeast.raster_x_ == 100);
    assert(s_northeast.raster_y_ == 200);

    std::printf("  GeoPos2RasterPos coordinate conversion  PASSED\n");

    // Test Station field access
    assert(s_center.id_ == 1);
    assert(s_center.name_ == "center");
    assert(s_center.idx_ == 0);  // default-initialized
    std::printf("  Station field access  PASSED\n");

    // Test Rainfall / RainfallEvent data structures (compile-time test)
    {
      constexpr std::size_t kStations = 3;
      lms::rain::Rainfall<T, kStations> rf;
      rf.rainfall_row_.fill(T{0});
      assert(rf.rainfall_row_[0] == T{0});

      lms::rain::RainfallEvent<T, kStations> event;
      event.duration_ = 0;
      event.time_interval_s_ = 3600;
      assert(event.rainfall_event_.empty());
      std::printf("  Rainfall / RainfallEvent data structures  PASSED\n");
    }

    std::printf("\n===== All tests PASSED =====\n");

    return 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
}
