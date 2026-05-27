#include "./utils.cpp"
#include "base/direct.h"
#include "base/core.h"
template <typename T>
class model
{

public:
  Task<T> SimulateOneStep(StateRaster<T> &temp_state_param_, const ConstRaster<T> &const_param_, std::vector<T> station_rain)
  {

    for (auto &[idx, item] : std::views::enumerate(iter_order_))
    {
      int target_idx = target_idx_[idx];
      FlowGeneration(int idx, int target_idx, StateRaster<T> &state_param_, const ConstRaster<T> &const_param_, time_interval);
    }
  };
  Simulate()
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
  std::vector<ConstParam<T>> state_param_;
  const std::vector<StateParam<T>> const_param_;
  const ModelMeta<T> model_meta_;
  GlobalParam<T> global_param_;
  std::vector<int> iter_order_;
  std::vector<int> target_idx_;
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
        FlowConfluenceStepOnce(&state_param_[i], const &const_param_[i], &state_param_[target_idx_[i]], &rainfall_, const &meta_data_, const &global_param_);
      }
    }
  }
};

template <typename T>
auto ProcessCell(StateParam<T> &loc_state, ConstParam<T> &const_param, StateParam<T> &target_state, int time_interval)
{
  // copy state_param raster here
  auto local_state_param = state_param;
  auto target_idx = GetTargetIdx(idx);

  FlowConfluenceStepOnce(int idx, int target_idx, StateRaster<T> &state_param_, const ConstRaster<T> &const_param_, time_interval);
};
