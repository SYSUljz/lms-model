#pragma once

#include <memory>
#include <optional>
#include <ranges>
#include <vector>

#include "base/core.h"
#include "base/direct.h"
#include "base/raster.h"

// Forward declarations: free template functions defined in other translation
// units. These should eventually move into dedicated headers.

template <typename T>
void FlowGeneration(lms::core::StateParam<T>& state_param_loc,
                    const lms::core::ConstParam<T>& const_param_loc,
                    lms::core::StateParam<T>& target_state, T rainfall,
                    const lms::core::ModelMeta<T>& meta_data,
                    const lms::core::GlobalParam<T>& global_param);

template <typename T>
void FlowConfluenceStepOnce(lms::core::StateParam<T>& state_param_loc,
                            const lms::core::ConstParam<T>& const_param_loc,
                            lms::core::StateParam<T>& target_state, T rainfall,
                            const lms::core::ModelMeta<T>& meta_data,
                            const lms::core::GlobalParam<T>& global_param);

namespace lms {
namespace model {

template <typename T>
class Model {
 public:
  Model() {};

  // flowGeneration doesn't require child time step iteration
  auto SimulateOneStep(std::vector<T> station_rain) {
    for (auto& [idx, item] : std::views::enumerate(iter_order_)) {
      int target_idx = target_idx_[idx];
      // TODO: map cell to station for correct rainfall index
      FlowGeneration(state_param_[item], const_param_[item],
                     state_param_[target_idx],
                     station_rain.empty() ? static_cast<T>(0) : station_rain[0],
                     model_meta_, global_param_);
    }

    return FlowConfluenceMultiStep();
  }

  auto Simulate() {
    // NOTE: StateRaster copy is deleted; using reference as a workaround.
    // A proper deep-copy / clone method is needed on StateRaster.
    auto& temp_state_param_ = state_param_;
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

  auto FlowConfluenceMultiStep() {
    for (int i = 0; i < model_meta_.confluence_steps_; ++i) {
      for (auto& [idx, item] : std::views::enumerate(iter_order_)) {
        if (target_idx_[idx].has_value()) [[likely]] {
          FlowConfluenceStepOnce(
              state_param_[item], const_param_[item],
              state_param_[*target_idx_[idx]],
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
