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

#include "base/factor.h"
#include "base/model.h"
#include "base/rain.h"
#include "builder/rain_builder.h"
#include "builder/raster_builder.h"
#include "io/gdal_reader.h"

using T = double;
using lms::core::ConstParam;
using lms::core::Label;
using lms::factor::Factor;
using lms::raster::ConstRaster;
using lms::raster::Raster;
int main(int argc, char** argv) {
  const std::string base_dir = "data/meizhou";
  const std::string model_dir = (argc > 1) ? argv[1] : base_dir + "/model";

  try {
    std::printf("== Full Model Construction Start ==\n");

    // 1. Build Raster Data
    ModelBuilder<T> builder;
    builder.from_directory(model_dir);
    ConstRaster<T> cr = builder.build_const_raster();

    // Set up Meta Data
    lms::core::ModelMeta<T> meta {};
    meta.width_ = cr.meta.width_;
    meta.heigh_ = cr.meta.heigh_;
    meta.cell_size_ = cr.meta.cell_size_;
    meta.time_interval_s_ = 3600;  // 1 hour
    meta.confluence_steps_ = 6;
    meta.runoff_dt_s_ = 600;
    // Bounds for coordinate conversion (placeholder values)
    meta.raster_min_lat_ = 2580000.0;
    meta.raster_max_lat_ = 2610000.0;
    meta.raster_min_long_ = 630000.0;
    meta.raster_max_long_ = 660000.0;

    auto order = builder.BuildOrder();
    auto targets = builder.BuildTarget(order);
    auto sr = builder.BuildStateRaster(cr);

    std::printf("Raster data built: %zu active cells, %zu order length\n",
                std::count(cr.active.begin(), cr.active.end(), 1), order.size());

    // 2. Build Rainfall Data
    const std::string rain_file = "2007060608.csv";
    const std::string station_file = "st_old.CSV";
    constexpr std::size_t kStationCnt = 20;

    lms::rain::RainBuilder<T, kStationCnt> rain_builder;
    rain_builder.from_directory(base_dir);
    rain_builder.BuildAll(station_file, rain_file, meta);

    auto rainfall_matrix = rain_builder.GetRainfallMatrix();
    std::printf("Rainfall data loaded: %zu time steps, %zu stations\n", rainfall_matrix.size(),
                rain_builder.stations().size());

    // 3. Assemble Model
    lms::core::GlobalParam<T> global_param {};
    global_param.soil_alpha_ = 4.0;
    global_param.baseflow_coff = 0.998;
    global_param.v = 0.7;
    global_param.manning = 0.025;
    global_param.ss = 60.0;
    global_param.init_soil_water = 0.3;

    lms::model::Model<T> model(std::move(sr), std::move(cr), meta, global_param, std::move(order), std::move(targets),
                               std::move(rainfall_matrix), rain_builder.stations());

    std::printf("Model assembled successfully.\n");

    model.BuildStationID();
    // 4. Run Simulation
    std::printf("Starting simulation...\n");
    Factor<T> factor;
    auto factor_model = model.BuildWithFactor(factor);
    model.SimulateAll();
    std::printf("Simulation completed.\n");

    // 5. Verify results (simple check)
    const auto& results = model.GetResult();
    std::printf("Simulation Results (outlet flow at each time step):\n");
    for (std::size_t i = 0; i < results.size(); ++i) {
      std::printf("  Time step %zu: %g\n", i + 1, results[i] / 3600);
    }

    auto& final_state = model.state_param();
    double total_runoff = 0;
    for (std::size_t i = 0; i < final_state.meta_.width_ * final_state.meta_.heigh_; ++i) {
      if (final_state.active_[i]) {
        total_runoff += final_state[i].runoff;
      }
    }
    std::printf("Final Total Runoff (sum across all cells): %g\n", total_runoff);

    std::printf("== Full Model Construction End ==\n");

    return 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 1;
  }
}
