#include <cmath>

#include "base/core.h"
#include "base/direct.h"
#include "../utils.hpp"

using lms::core::StateParam;
using lms::core::ConstParam;
using lms::core::ModelMeta;
using lms::core::GlobalParam;
using lms::core::Label;
using lms::direct::GetDirectFactor;

template <typename T>
T GetStepInFlow(StateParam<T>& state_param_loc, const ConstParam<T>& const_param_loc, T rainfall,
                const ModelMeta<T>& meta_data, const GlobalParam<T>& global_param) { return T{}; }

template <typename T>
void FlowGeneration(StateParam<T>& state_param_loc, const ConstParam<T>& const_param_loc, StateParam<T>& target_state,
                    T rainfall, const ModelMeta<T>& meta_data, const GlobalParam<T>& global_param) {
  // calculate flow generate in soil cell
  if (const_param_loc.label == Label::Soil) {
    auto sat = const_param_loc.sat;
    auto fc = const_param_loc.fc;
    auto soil_moisture = const_param_loc.soil_moisture;
    auto slop = const_param_loc.slop;
    auto ks = const_param_loc.ks * meta_data.time_interval_s_;
    auto b = const_param_loc.b;
    auto k = GetK(soil_moisture, sat, b, ks);
    auto zs = const_param_loc.zs;
    auto lateral_in_q = state_param_loc.lateral_in_flow_mm;
    auto direct_factor = GetDirectFactor<T>(const_param_loc.d8);
    auto ep = const_param_loc.ep;
    auto v = const_param_loc.v;
    auto cell_size = meta_data.cell_size_;
    auto soil_alpha = global_param.soil_alpha_;
    // calculate evaporation
    if (soil_moisture > fc) {
      state_param_loc.actual_evaporate = ep * v;
    } else if (soil_moisture > const_param_loc.wl) {
      state_param_loc.actual_evaporate = (1 - v) * ep * (soil_moisture - const_param_loc.wl) / (fc - const_param_loc.wl);
    } else {
      state_param_loc.actual_evaporate = 0;
    }
    // calculate produce runoff
    T percolate_out_q = 0;
    T lateral_out_q = 0;
    if (soil_moisture > fc) {
      // Percolate out quantity
      percolate_out_q = 0.001 * (k + state_param_loc.per_mm) / 2 * cell_size * cell_size;
      lateral_out_q =
          0.001 * (k * slop + state_param_loc.lateral_in_flow_mm) / 2 * direct_factor * cell_size * zs * 0.001;
      auto excess_q = 0.001 * zs * (soil_moisture - fc) * cell_size * cell_size;
      if (percolate_out_q + lateral_out_q > excess_q) {
        percolate_out_q = excess_q * percolate_out_q / (percolate_out_q + lateral_out_q);
        lateral_out_q = excess_q - percolate_out_q;
      }
    }
    auto depth = rainfall - state_param_loc.actual_evaporate +
                 (lateral_in_q - percolate_out_q - lateral_out_q) / cell_size / cell_size * 1000;
    if (depth > k) {
      auto a = (std::exp(soil_alpha * soil_moisture / sat) - 1) / (std::exp(soil_alpha) - 1);
      state_param_loc.runoff += (depth - k) * a;
      depth = k + (depth - k) * (1 - a);
    }

    // update downstream grid's lateral_in_q_
    target_state.lateral_in_flow_mm += lateral_out_q / cell_size / cell_size * 1000;
    state_param_loc.per_mm = percolate_out_q / cell_size / cell_size * 1000;
    state_param_loc.soil_moisture += depth / zs;

    if (state_param_loc.soil_moisture > sat) {
      state_param_loc.runoff += zs * (state_param_loc.soil_moisture - sat);
      state_param_loc.soil_moisture = sat;
    }

    state_param_loc.groundwater_mm = global_param.baseflow_coff * state_param_loc.groundwater_mm +
                                     (1 - global_param.baseflow_coff) * state_param_loc.per_mm;
    // state_param_loc.runoff += state_param_loc.groundwater_mm;
  }

  else {
    auto depth = rainfall - state_param_loc.actual_evaporate;
    if (depth > 0) {
      state_param_loc.runoff += depth;
    } else {
      state_param_loc.actual_evaporate = rainfall;
    }
  }
}
