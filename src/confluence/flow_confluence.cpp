#include "base/core.h"
#include "base/direct.h"
#include "../utils.hpp"
template <typename T>
void FlowConfluenceStepOnce(StateParam<T>& state_param_loc, const ConstParam<T>& const_param_loc,
                            StateParam<T>& target_state, T rainfall, const ModelMeta& meta_data,
                            const GlobalParam<T>& global_param) {
  auto cell_size = meta_data.cell_size_;
  auto steps = meta_data.confluence_steps_;
  auto dt = meta_data.runoff_dt_s_;
  auto alpha = GetSoilAlpha<T>(const_param_loc.n, cell_size, const_param_loc.slop);
  auto beta = static_cast<T>(0.6);
  auto dx = GetDirectFactor<T>(const_param_loc.label) * cell_size;

  auto upstream_flow_in = state_param_loc.upstream_in_flow;
  auto prev_t_flow = state_param_loc.prev_t_flow;
  auto lateralrunoff = state_param_loc.lateral_in_flow_mm;
  auto underground_flow_mm = state_param_loc.groundwater_mm;
  auto total_runoff_mm = state_param_loc.runoff + lateralrunoff + underground_flow_mm;
  // q is Lateral inflow intensity per unit length (m^3/h/m)
  // (mm * 0.001 * ArRea_m2) / (dX_m * dT_h) = m^3 / h / m
  auto q = total_runoff_mm * L * L * 0.001 / steps / dx / dt;
  // Slope Kinematic Wave Routing for Soil cell

  if (const_param_loc.label == Label::Soil) {
    auto ground_flow = state_param_loc.groundwater_mm;
    q += ground_flow * cell_size * cell_size * 0.001 / dx / steps / dt;
    state_param_loc.prev_t_flow = SolveSaintVenant(prev_t_flow, q, alpha, beta, dt, dx, upstream_flow_in, prev_t_flow);
    target_state.upstream_in_flow += state_param_loc.prev_t_flow;
  }

  else if (const_param_loc.label == Label::Channel) {
    auto bw = const_param_loc.bw;
    auto ss = global_param.ss;
    auto bs = const_param_loc.bs;
    auto manning = global_param.manning;
    auto water_depth = state_param_loc.water_level;
    auto channel_x = GetChannelX(water_depth, bw, ss);
    auto next_h = target_state.water_level;
    auto Sf = bs - (next_h - water_depth) / dx;

    auto alpha = GetChannelAlpha(manning, channel_x, Sf);
    T beta = static_cast<T>(0.6);
    // We no longer use variables like current_flow as a temp_container once we
    // calculate flow out , in will be stored in prev_t_flow.
    state_param_loc.prev_t_flow = SolveSaintVenant(prev_t_flow, q, alpha, beta, dt, dx, upstream_flow_in, prev_t_flow);
    state_param_loc.water_level = solveQtoH(state_param_loc.prev_t_flow / 3600, manning, bs, bw, ss);
    target_state.upstream_in_flow += state_param_loc.prev_t_flow;
  }
};