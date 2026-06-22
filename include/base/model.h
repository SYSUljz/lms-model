#pragma once

#include <limits>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "base/core.h"
#include "base/direct.h"
#include "base/factor.h"
#include "base/rain.h"
#include "base/raster.h"
// Forward declarations: free template functions defined in other translation
// units. These should eventually move into dedicated headers.

template <typename T>
void FlowGeneration(lms::core::StateParam<T>& state_param_loc, const lms::core::ConstParam<T>& const_param_loc,
                    lms::core::StateParam<T>& target_state, T rainfall, const lms::core::ModelMeta<T>& meta_data,
                    const lms::core::GlobalParam<T>& global_param);

template <typename T>
void FlowConfluenceStepOnce(lms::core::StateParam<T>& state_param_loc, const lms::core::ConstParam<T>& const_param_loc,
                            lms::core::StateParam<T>& target_state, T rainfall,
                            const lms::core::ModelMeta<T>& meta_data, const lms::core::GlobalParam<T>& global_param);

namespace lms {
namespace model {

template <typename U>
using Station = lms::rain::Station<U>;
template <typename U>
using Factor = lms::factor::Factor<U>;
template <typename T>
class Model {
 public:
  Model(lms::raster::StateRaster<T> state_param, lms::raster::ConstRaster<T> const_param,
        lms::core::ModelMeta<T> model_meta, lms::core::GlobalParam<T> global_param, std::vector<int> iter_order,
        std::vector<std::optional<int>> target_idx, std::vector<std::vector<T>> rainfall,
        std::vector<Station<T>> stations)
      : state_param_(std::move(state_param)),
        const_param_(std::move(const_param)),
        model_meta_(model_meta),
        global_param_(global_param),
        iter_order_(std::move(iter_order)),
        target_idx_(std::move(target_idx)),
        rainfall_(std::move(rainfall)),
        rainfall_data_length_(static_cast<int>(rainfall_.size())),
        // copy stations
        stations_(stations) {
    InitializeState();
  }

  // Thiessen polygon (Voronoi) interpolation: assign each cell to its nearest
  // station so that every cell's rainfall is driven by the closest station.
  Model& BuildStationID() {
    station_id_.resize(iter_order_.size());
    int width = static_cast<int>(model_meta_.width_);

    for (std::size_t i = 0; i < iter_order_.size(); ++i) {
      int item = iter_order_[i];
      int x = item % width;
      int y = item / width;

      int nearest_station_id = -1;
      double min_dist_sq = std::numeric_limits<double>::max();

      for (const auto& station : stations_) {
        double dx = static_cast<double>(x) - static_cast<double>(station.raster_x_);
        double dy = static_cast<double>(y) - static_cast<double>(station.raster_y_);
        double dist_sq = dx * dx + dy * dy;
        if (dist_sq < min_dist_sq) {
          min_dist_sq = dist_sq;
          nearest_station_id = static_cast<int>(station.id_);
        }
      }

      station_id_[i] = nearest_station_id;
    }

    return *this;
  }

  // flowGeneration doesn't require child time step iteration
  void SimulateOneStep(std::vector<T> station_rain) {
    for (int item : iter_order_) {
      state_param_[item].runoff = 0;
      state_param_[item].lateral_in_flow_mm = 0;
    }

    for (std::size_t idx = 0; idx < iter_order_.size(); ++idx) {
      auto item = iter_order_[idx];
      auto target_opt = target_idx_[idx];
      if (!target_opt.has_value()) continue;

      int target_idx = *target_opt;
      std::size_t rain_idx = station_id_[idx];

      FlowGeneration(state_param_[item], const_param_[item], state_param_[target_idx],
                     station_rain.empty() ? static_cast<T>(0) : station_rain[rain_idx], model_meta_, global_param_);
    }
    results_.push_back(FlowConfluenceMultiStep());
  }

  auto SimulateAll() {
    for (auto& item : rainfall_) {
      SimulateOneStep(item);
    }
  }
  // return a copy of this model
  Model BuildWithFactor(Factor<T> factor_) {
    lms::raster::StateRaster<T> state_param = state_param_;
    lms::core::GlobalParam<T> global_param = global_param_;
    lms::raster::ConstRaster<T> const_param = const_param_;
    ApplyGlobalFactor(global_param, factor_);
    for (size_t i = 0; i < state_param.size(); i++) {
      ApplyFactor(const_param[i], factor_);
    }
    Model<T> factormodel(state_param, const_param, model_meta_, global_param, iter_order_, target_idx_, rainfall_,
                         stations_);
    return factormodel;
  }

  // Particle Swarm Optimization
  auto PsoStep() {}

  bool BuildOrder() { return false; }

  bool BuildTargetOrder() { return false; }

  bool BuildRain() { return false; }

  int GetTargetIdx(int this_idx) { return 0; }

  lms::raster::StateRaster<T>& state_param() { return state_param_; }

  const std::vector<T>& GetResult() const { return results_; }

 private:
  void InitializeState() {
    for (int item : iter_order_) {
      auto& state = state_param_[item];
      const auto& constant = const_param_[item];

      state.soil_moisture = constant.sat * global_param_.init_soil_water;

      state.lat_mm = 0.0;
      state.per_mm = 0.0;
      state.lateral_in_flow_mm = 0.0;
      state.groundwater_mm = 0.0;
    }
  }
  lms::raster::StateRaster<T> state_param_;
  lms::raster::ConstRaster<T> const_param_;
  lms::core::ModelMeta<T> model_meta_;
  lms::core::GlobalParam<T> global_param_;
  // During the simulation, a mock channel cell is added after the outlet to
  // ensure the calculation accuracy of the outlet channel element.
  // iter_order_[a]=b means the a'th step iter we process cell with b idx in
  // param raster (both state and const)
  std::vector<int> iter_order_;
  // iter_order_[a]=c means the ath step we are processing b cell the flow
  // direction of b cell is c cell
  std::vector<std::optional<int>> target_idx_;
  std::vector<int> station_id_;
  // rainfall_[time][station]
  std::vector<std::vector<T>> rainfall_;
  int rainfall_data_length_;
  // station meta_data collection
  std::vector<Station<T>> stations_;
  // store simulate result
  std::vector<T> results_;

  // return water flow of outlet cell
  T FlowConfluenceMultiStep() {
    lms::core::StateParam<T> exit_sink;

    for (int i = 0; i < model_meta_.confluence_steps_; ++i) {
      std::printf("  Confluence step %d/%zu\n", i + 1, model_meta_.confluence_steps_);
      std::fflush(stdout);

      // Clear upstream_in_flow for all cells at the start of each confluence sub-timestep.
      for (int item : iter_order_) {
        state_param_[item].upstream_in_flow = 0;
      }
      exit_sink.upstream_in_flow = 0;

      for (std::size_t idx = 0; idx < iter_order_.size(); ++idx) {
        auto item = iter_order_[idx];
        auto target_opt = target_idx_[idx];

        if (target_opt.has_value()) {
          FlowConfluenceStepOnce(state_param_[item], const_param_[item], state_param_[*target_opt], static_cast<T>(0),
                                 model_meta_, global_param_);
        } else {
          // Route the outlet cell to the exit_sink
          FlowConfluenceStepOnce(state_param_[item], const_param_[item], exit_sink, static_cast<T>(0), model_meta_,
                                 global_param_);
        }
      }

      // Clear runoff and lateral inflow after they have been routed in the first sub-timestep
      for (int item : iter_order_) {
        state_param_[item].runoff = 0;
        state_param_[item].lateral_in_flow_mm = 0;
      }
    }

    if (!iter_order_.empty()) {
      return state_param_[iter_order_.back()].prev_t_flow;
    }
    return 0;
  }

  void ApplyFactor(lms::raster::ConstParam<T>& const_param, lms::factor::Factor<T>& factor_) {
    const_param.sat *= factor_.sat;
    const_param.fc *= factor_.fc;
    const_param.wl *= factor_.wl;
    const_param.ks *= factor_.ks;
    const_param.zs *= factor_.zs;
    const_param.b *= factor_.b;
    const_param.n *= factor_.n;
    const_param.v *= factor_.v;
    const_param.bs *= factor_.bs;
    const_param.bw *= factor_.bw;
    const_param.ep *= factor_.ep;
  }

  void ApplyGlobalFactor(lms::core::GlobalParam<T>& global_param, lms::factor::Factor<T>& factor_) {
    global_param.manning *= factor_.manning;
    global_param.soil_alpha_ *= factor_.soil_alpha;
    global_param.init_soil_water *= factor_.init_soil_water;
    global_param.ss *= factor_.ss;
  }

  void CompressRaster() {};
};

}  // namespace model
}  // namespace lms
