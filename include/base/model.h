#pragma once

#include <limits>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "base/core.h"
#include "base/direct.h"
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
        stations_(stations) {}

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
  auto SimulateOneStep(std::vector<T> station_rain) {
    for (std::size_t idx = 0; idx < iter_order_.size(); ++idx) {
      auto item = iter_order_[idx];
      auto target_opt = target_idx_[idx];
      if (!target_opt.has_value()) continue;

      int target_idx = *target_opt;
      std::size_t rain_idx = station_id_[idx];
      // TODO: map cell to station for correct rainfall index
      FlowGeneration(state_param_[item], const_param_[item], state_param_[target_idx],
                     station_rain.empty() ? static_cast<T>(0) : station_rain[rain_idx], model_meta_, global_param_);
    }

    return FlowConfluenceMultiStep();
  }

  auto SimulateAll() {
    for (auto& item : rainfall_) {
      SimulateOneStep(item);
    }
  }

  // Particle Swarm Optimization
  auto PsoStep() {}

  bool BuildOrder() { return false; }

  bool BuildTargetOrder() { return false; }

  bool BuildRain() { return false; }

  int GetTargetIdx(int this_idx) { return 0; }

  const lms::raster::StateRaster<T>& state_param() const { return state_param_; }

 private:
  lms::raster::StateRaster<T> state_param_;
  const lms::raster::ConstRaster<T> const_param_;
  const lms::core::ModelMeta<T> model_meta_;
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

  auto FlowConfluenceMultiStep() {
    for (int i = 0; i < model_meta_.confluence_steps_; ++i) {
      for (std::size_t idx = 0; idx < iter_order_.size(); ++idx) {
        auto item = iter_order_[idx];
        if (target_idx_[idx].has_value()) [[likely]] {
          FlowConfluenceStepOnce(state_param_[item], const_param_[item], state_param_[*target_idx_[idx]],
                                 // TODO: pass actual rainfall for this cell / time step
                                 static_cast<T>(0), model_meta_, global_param_);
        } else {
          // todo: find a way to process pourpoint branch
          return state_param_[idx];
        }
      }
    }
    // TODO: determine proper return value for the normal path
    return state_param_[0];
  }
};

}  // namespace model
}  // namespace lms
