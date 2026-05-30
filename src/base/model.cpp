#include "./utils.cpp"
#include "base/direct.h"
#include "base/core.h"
#include "base/raster.h"
#include <memory>
template <typename T>
class model
{

public:
  // flowGeneration doesn't require child time step iteration
  Task<T> SimulateOneStep(std::vector<T> station_rain)
  {

    for (auto &[idx, item] : std::views::enumerate(iter_order_))
    {
      int target_idx = target_idx_[idx];
      FlowGeneration(int idx, int target_idx, StateRaster<T> &state_param_, const ConstRaster<T> &const_param_, time_interval);
    }

    auto result = FlowConfluenceMultiStep();
  };

  auto Simulate()
  {
    StateRaster<T> temp_state_param_ = state_param_;
    for (auto &item : rainfall_)
    {
      SimulateOneStep();
    }
  };

  // Particle Swarm Optimization
  Task<T> PSO() {};
  bool BuildOrder() {};
  bool BuildTargetOrder() {};
  bool BuildRain() {};

  int GetTargetIdx(int this_idx) {};

private:
  StateRaster state_param_;
  const ConstRaster const_param_;
  const ModelMeta<T> model_meta_;
  GlobalParam<T> global_param_;
  // During the simulation, a mock channel cell is added after the outlet to ensure the calculation accuracy of the outlet channel element.
  // iter_order_[a]=b means the a'th step iter we process cell with b idx in param raster (both state and const)
  std::vector<int> iter_order_;
  // iter_order_[a]=c means the ath step we are processing b cell the flow direction of b cell is c cell
  std::vector<std::optional<int>> target_idx_;
  std::vector<int> station_id_;
  // rainfall_[time][station]
  std::vector<std::vector<T>> rainfall_;
  int rainfall_data_length_;

  template <typename T>
  auto FlowConfluenceMultiStep()
  {
    for (int i = 0; i < model_meta_.steps; ++i)
    {
      for (auto &[idx, item] : std::views::enumerate(iter_order_))
      {
        if (target_idx_[idx].has_value()) [[likely]]
        {
          FlowConfluenceStepOnce(&state_param_[item], const &(*const_param_)[item], &state_param_[target_idx_[idx]], &rainfall_, const &meta_data_, const &global_param_);
        }
        else
        {
          // todo: find a way to process pourpoint branch
          return this.state_param_[idx];
        }
      }
    }
  }
};
