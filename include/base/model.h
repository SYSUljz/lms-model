#pragma once

#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "algorithm/object_function.hpp"
#include "base/core.h"
#include "base/direct.h"
#include "base/factor.h"
#include "base/flow.h"
#include "base/rain.h"
#include "base/raster.h"
// Forward declarations: free template functions defined in other translation
// units. These should eventually move into dedicated headers.

template <typename T>
void FlowGeneration(lms::core::StateParam<T>& state_param_loc, const lms::core::ConstParam<T>& const_param_loc,
                    lms::core::StateParam<T>& target_state, T rainfall, const lms::core::ModelMeta<T>& meta_data,
                    const lms::core::GlobalParam<T>& global_param, const lms::factor::Factor<T>& factor);

template <typename T>
void FlowConfluenceStepOnce(lms::core::StateParam<T>& state_param_loc, const lms::core::ConstParam<T>& const_param_loc,
                            lms::core::StateParam<T>& target_state, T rainfall,
                            const lms::core::ModelMeta<T>& meta_data, const lms::core::GlobalParam<T>& global_param,
                            const lms::factor::Factor<T>& factor);

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
        std::vector<Station<T>> stations, lms::flow::Flow<T> flow)
      : state_param_(std::move(state_param)),
        const_param_(std::move(const_param)),
        model_meta_(model_meta),
        global_param_(global_param),
        iter_order_(std::move(iter_order)),
        target_idx_(std::move(target_idx)),
        rainfall_(std::move(rainfall)),
        rainfall_data_length_(static_cast<int>(rainfall_.size())),
        // copy stations
        stations_(stations),
        flow_(std::move(flow)) {
    InitializeState();
    PrecomputeHelpers();
    global_param_.soil_alpha_exp_minus_one_inv_ = 
        (std::abs(global_param_.soil_alpha_) > static_cast<T>(1e-6))
        ? (static_cast<T>(1) / (std::exp(global_param_.soil_alpha_) - static_cast<T>(1)))
        : static_cast<T>(1);
    objective_func_ = [](const std::vector<T>& sim, const std::vector<T>& obs) {
      return T {1} - lms::objfunc::CalculateNSE(sim, obs);
    };
  }

  // Delete copy constructor and copy assignment operator
  Model(const Model&) = delete;
  Model& operator=(const Model&) = delete;

  // Implement move constructor and move assignment operator
  Model(Model&&) noexcept = default;
  Model& operator=(Model&&) noexcept = default;

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

  // Run simulation with optional factor. Returns results without mutating Model.
  std::vector<T> Simulate(const Factor<T>& factor = {}) const {
    // 1. Scale global parameters
    lms::core::GlobalParam<T> scaled_global_param = global_param_;
    scaled_global_param.manning *= factor.manning;
    scaled_global_param.soil_alpha_ *= factor.soil_alpha;
    scaled_global_param.init_soil_water *= factor.init_soil_water;
    scaled_global_param.ss *= factor.ss;

    scaled_global_param.soil_alpha_exp_minus_one_inv_ = 
        (std::abs(scaled_global_param.soil_alpha_) > static_cast<T>(1e-6))
        ? (static_cast<T>(1) / (std::exp(scaled_global_param.soil_alpha_) - static_cast<T>(1)))
        : static_cast<T>(1);

    // 2. Initialize running state_param by copying initial state_param_
    lms::raster::StateRaster<T> state_param = state_param_;
    // Re-initialize soil moisture with scaled global init_soil_water and scaled sat
    T scaled_init_soil_water = global_param_.init_soil_water * factor.init_soil_water;
    for (int item : iter_order_) {
      auto& state = state_param[item];
      const auto& constant = const_param_[item];
      state.soil_moisture = (constant.sat * factor.sat) * scaled_init_soil_water;
    }

    std::vector<T> results;
    results.reserve(rainfall_.size());

    // 3. Run step-by-step
    for (const auto& rain_step : rainfall_) {
      SimulateOneStep(state_param, scaled_global_param, factor, rain_step, results);
    }

    return results;
  }

  // Run simulation and evaluate against observed flow. Returns objective value.
  T SimulateEval(const Factor<T>& factor = {}) const {
    return objective_func_(Simulate(factor), flow_.data);
  }

  bool BuildOrder() { return false; }

  bool BuildTargetOrder() { return false; }

  bool BuildRain() { return false; }

  int GetTargetIdx(int this_idx) { return 0; }

  const lms::raster::StateRaster<T>& state_param() const { return state_param_; }

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

  void PrecomputeHelpers() {
    T cell_size = static_cast<T>(model_meta_.cell_size_);
    for (int item : iter_order_) {
      auto& constant = const_param_[item];
      T manning = constant.n;
      T slope = constant.slop;
      if (std::isnan(manning) || manning <= 0) manning = static_cast<T>(0.03);
      if (std::isnan(cell_size) || cell_size <= 0) cell_size = static_cast<T>(30);
      if (std::isnan(slope) || slope <= 0) [[unlikely]] {
        slope = static_cast<T>(1e-6);
      }
      T alpha = std::pow(manning * std::pow(cell_size, static_cast<T>(0.666667)) * std::pow(slope, static_cast<T>(-0.5)) /
                             static_cast<T>(3600),
                          static_cast<T>(0.6));
      if (std::isnan(alpha)) alpha = static_cast<T>(1);
      constant.cached_soil_alpha = alpha;
      constant.cached_dx = GetDirectFactor<T>(constant.d8) * cell_size;
    }
  }

  void SimulateOneStep(lms::raster::StateRaster<T>& state_param,
                       const lms::core::GlobalParam<T>& scaled_global_param,
                       const Factor<T>& factor,
                       const std::vector<T>& station_rain,
                       std::vector<T>& results) const {
    for (int item : iter_order_) {
      state_param[item].runoff = 0;
      state_param[item].lateral_in_flow_mm = 0;
    }

    for (std::size_t idx = 0; idx < iter_order_.size(); ++idx) {
      auto item = iter_order_[idx];
      auto target_opt = target_idx_[idx];
      if (!target_opt.has_value()) continue;

      int target_idx = *target_opt;
      std::size_t rain_idx = station_id_[idx];

      FlowGeneration(state_param[item], const_param_[item], state_param[target_idx],
                     station_rain.empty() ? static_cast<T>(0) : station_rain[rain_idx], model_meta_,
                     scaled_global_param, factor);
    }
    results.push_back(FlowConfluenceMultiStep(state_param, scaled_global_param, factor));
  }

  T FlowConfluenceMultiStep(lms::raster::StateRaster<T>& state_param,
                            const lms::core::GlobalParam<T>& scaled_global_param,
                            const Factor<T>& factor) const {
    lms::core::StateParam<T> exit_sink;

    for (int i = 0; i < model_meta_.confluence_steps_; ++i) {
      // Clear upstream_in_flow for all cells at the start of each confluence sub-timestep.
      for (int item : iter_order_) {
        state_param[item].upstream_in_flow = 0;
      }
      exit_sink.upstream_in_flow = 0;

      for (std::size_t idx = 0; idx < iter_order_.size(); ++idx) {
        auto item = iter_order_[idx];
        auto target_opt = target_idx_[idx];

        if (target_opt.has_value()) {
          FlowConfluenceStepOnce(state_param[item], const_param_[item], state_param[*target_opt], static_cast<T>(0),
                                 model_meta_, scaled_global_param, factor);
        } else {
          // Route the outlet cell to the exit_sink
          FlowConfluenceStepOnce(state_param[item], const_param_[item], exit_sink, static_cast<T>(0), model_meta_,
                                 scaled_global_param, factor);
        }
      }

      // Clear runoff and lateral inflow after they have been routed in the first sub-timestep
      for (int item : iter_order_) {
        state_param[item].runoff = 0;
        state_param[item].lateral_in_flow_mm = 0;
      }
    }

    if (!iter_order_.empty()) {
      return state_param[iter_order_.back()].prev_t_flow;
    }
    return 0;
  }

  lms::raster::StateRaster<T> state_param_;
  lms::raster::ConstRaster<T> const_param_;
  lms::core::ModelMeta<T> model_meta_;
  lms::core::GlobalParam<T> global_param_;
  std::vector<int> iter_order_;
  std::vector<std::optional<int>> target_idx_;
  std::vector<int> station_id_;
  std::vector<std::vector<T>> rainfall_;
  int rainfall_data_length_;
  std::vector<Station<T>> stations_;
  lms::flow::Flow<T> flow_;
  std::function<T(const std::vector<T>&, const std::vector<T>&)> objective_func_;
};

}  // namespace model
}  // namespace lms
