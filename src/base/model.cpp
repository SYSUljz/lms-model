#include "base/model.h"

#include <memory>

#include "./utils.cpp"
namespace lxh {
namespace model {
template <typename T>

// flowGeneration doesn't require child time step iteration
template <typename T>
Task<T> Model::SimulateOneStep(std::vector<T> station_rain) {
  for (auto& [idx, item] : std::views::enumerate(iter_order_)) {
    int target_idx = target_idx_[idx];
    FlowGeneration(int idx, int target_idx, StateRaster<T>& state_param_, const ConstRaster<T>& const_param_,
                   time_interval);
  }

  auto result = FlowConfluenceMultiStep();
};

auto Simulate() {
  StateRaster<T> temp_state_param_ = state_param_;
  for (auto& item : rainfall_) {
    SimulateOneStep();
  }
};

// Particle Swarm Optimization
Task<T> PSO() {};
bool BuildOrder() {};
bool BuildTargetOrder() {};
bool BuildRain() {};

int GetTargetIdx(int this_idx) {};

template <typename T>
auto Model::FlowConfluenceMultiStep() {
  for (int i = 0; i < model_meta_.steps; ++i) {
    for (auto& [idx, item] : std::views::enumerate(iter_order_)) {
      if (target_idx_[idx].has_value()) [[likely]] {
        FlowConfluenceStepOnce(&state_param_[item], const&(*const_param_)[item], &state_param_[target_idx_[idx]],
                               &rainfall_, const& meta_data_, const& global_param_);
      } else {
        // todo: find a way to process pourpoint branch
        return this.state_param_[idx];
      }
    }
  }
}
}  // namespace model
}  // namespace lxh
