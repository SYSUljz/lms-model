#include "../utils.hpp"
#include "base/core.h"
#include "base/direct.h"
#include "base/factor.h"

using lms::core::ConstParam;
using lms::core::GlobalParam;
using lms::core::Label;
using lms::core::ModelMeta;
using lms::core::StateParam;
using lms::direct::GetDirectFactor;

template <typename T>
void FlowConfluenceStepOnce(StateParam<T>& state_param_loc, const ConstParam<T>& const_param_loc,
                            StateParam<T>& target_state, T rainfall, const ModelMeta<T>& meta_data,
                            const GlobalParam<T>& global_param, const lms::factor::Factor<T>& factor) {
  if (std::isnan(state_param_loc.runoff)) state_param_loc.runoff = 0;
  if (std::isnan(state_param_loc.lateral_in_flow_mm)) state_param_loc.lateral_in_flow_mm = 0;
  if (std::isnan(state_param_loc.groundwater_mm)) state_param_loc.groundwater_mm = 0;
  if (std::isnan(state_param_loc.prev_t_flow)) state_param_loc.prev_t_flow = 0;
  if (std::isnan(state_param_loc.upstream_in_flow)) state_param_loc.upstream_in_flow = 0;
  if (std::isnan(state_param_loc.water_level)) state_param_loc.water_level = 0;
  if (std::isnan(target_state.upstream_in_flow)) target_state.upstream_in_flow = 0;
  if (std::isnan(target_state.water_level)) target_state.water_level = 0;

  auto cell_size = meta_data.cell_size_;
  auto steps = meta_data.confluence_steps_;
  // runoff_dt_s_ is the sub-step time in seconds (e.g. 600s = 10 min)
  auto dt_sub = static_cast<T>(meta_data.runoff_dt_s_) / static_cast<T>(3600.0);
  // total confluence time = sub-step × steps (e.g. 10 min × 6 = 1 h)
  auto dt_full = dt_sub * static_cast<T>(steps);
  auto alpha = const_param_loc.cached_soil_alpha * std::pow(factor.n, static_cast<T>(0.6));
  auto beta = static_cast<T>(0.6);
  auto dx = const_param_loc.cached_dx;

  auto upstream_flow_in = state_param_loc.upstream_in_flow;
  auto prev_t_flow = state_param_loc.prev_t_flow;
  auto underground_flow_mm = state_param_loc.groundwater_mm;

  // Calculate lateral runoff contribution from upstream cells' lateral inflow
  T lateral_runoff = 0;
  if (const_param_loc.label == Label::Soil) {
    auto sat = const_param_loc.sat * factor.sat;
    auto cur = state_param_loc.soil_moisture;
    auto ks = const_param_loc.ks * factor.ks * static_cast<T>(meta_data.time_interval_s_) / static_cast<T>(3600.0);
    auto k = GetK(cur, sat, const_param_loc.b * factor.b, ks);
    auto zs = const_param_loc.zs * factor.zs;
    auto depth = state_param_loc.lateral_in_flow_mm;

    if (depth > k) {
      auto soil_alpha = global_param.soil_alpha_;
      if (std::abs(soil_alpha) > static_cast<T>(1e-6)) {
        auto a_val =
            (std::exp(soil_alpha * cur / sat) - static_cast<T>(1)) * global_param.soil_alpha_exp_minus_one_inv_;
        lateral_runoff += (depth - k) * a_val;
        depth = k + (depth - k) * (static_cast<T>(1) - a_val);
      }
    }

    state_param_loc.soil_moisture += depth / zs;
    if (state_param_loc.soil_moisture > sat) {
      lateral_runoff += zs * (state_param_loc.soil_moisture - sat);
      state_param_loc.soil_moisture = sat;
    }
  } else {
    lateral_runoff = state_param_loc.lateral_in_flow_mm;
  }

  // q is Lateral inflow intensity per unit length (m^3/h/m)
  // Runoff enters concentrated in the sub-step: volume / dt_sub / dx
  auto runoff_q = (state_param_loc.runoff + lateral_runoff) * cell_size * cell_size * static_cast<T>(0.001) / dx / dt_sub;
  // Groundwater is steady baseflow spread over the full timestep: volume / dt_full / dx
  auto groundwater_q = underground_flow_mm * cell_size * cell_size * static_cast<T>(0.001) / dx / dt_full;
  auto q = runoff_q + groundwater_q;

  // Slope Kinematic Wave Routing for Soil cell
  if (const_param_loc.label == Label::Soil) {
    state_param_loc.prev_t_flow =
        SolveSaintVenant(prev_t_flow, q, alpha, beta, dt_sub, dx, upstream_flow_in, prev_t_flow);
    target_state.upstream_in_flow += state_param_loc.prev_t_flow;
  }

  else if (const_param_loc.label == Label::Channel) {
    auto bw = const_param_loc.bw * factor.bw;
    auto ss = global_param.ss;
    auto bs = const_param_loc.bs * factor.bs;
    auto manning = global_param.manning;
    auto water_depth = state_param_loc.water_level;
    auto channel_x = GetChannelX(water_depth, bw, ss);
    auto next_h = target_state.water_level;
    auto Sf = bs - (next_h - water_depth) / dx;

    auto channel_alpha = GetChannelAlpha(manning, channel_x, Sf);
    T channel_beta = static_cast<T>(0.6);
    state_param_loc.prev_t_flow =
        SolveSaintVenant(prev_t_flow, q, channel_alpha, channel_beta, dt_sub, dx, upstream_flow_in, prev_t_flow);
    state_param_loc.water_level = solveQtoH(state_param_loc.prev_t_flow / 3600, manning, bs, bw, ss);
    target_state.upstream_in_flow += state_param_loc.prev_t_flow;
  }
}

template void FlowConfluenceStepOnce<double>(StateParam<double>& state_param_loc,
                                             const ConstParam<double>& const_param_loc,
                                             StateParam<double>& target_state, double rainfall,
                                             const ModelMeta<double>& meta_data,
                                             const GlobalParam<double>& global_param,
                                             const lms::factor::Factor<double>& factor);

template void FlowConfluenceStepOnce<float>(StateParam<float>& state_param_loc,
                                            const ConstParam<float>& const_param_loc,
                                            StateParam<float>& target_state, float rainfall,
                                            const ModelMeta<float>& meta_data,
                                            const GlobalParam<float>& global_param,
                                            const lms::factor::Factor<float>& factor);
