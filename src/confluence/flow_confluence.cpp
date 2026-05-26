#include "src/base/core.h"
#include "utils.cpp"
#include "direct.h"
template <typename T>
void FlowConfluence(StateParam<T> &state_param_loc, const ConstParam<T> &const_param_loc, StateParam<T> &target_state, T rainfall, const ModelMeta &meta_data, const GlobalParam<T> &global_param)
{
  auto cell_size = meta_data.cell_size_;
  auto flow = state_param_loc.current_flow;
  // Slope Kinematic Wave Routing for Soil cell
  if (const_param_loc.label == Label::Soil)
  {
    auto alpha = GetSoilAlpha<T>(const_param_loc.n, cell_size, const_param_loc.slop);
    auto beta = static_cast<T>(0.6);
    auto dx = GetDirectFactor<T>(const_param_loc.label) * cell_size;
    auto upstream_flow_in = state_param_loc.upstream_in_flow;
    // todo：create  a GetFlow function calculate cell flow in and Distribute evenly to each time step
  }
};