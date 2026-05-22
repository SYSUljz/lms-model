#include "model.cpp"
#include "utils.cpp"
template <typename T>
void FlowGeneration(StateParam<T> &loc_state_param, ConstParam<T> &loc_const_param, StateParam<T> &target_state, int time_interval, size_t cell_size)
{
  // calculate flow generate in soil cell
  if (loc_const_param.label == Label::Soil)
  {
    auto sat = loc_const_param.sat;
    auto fc = loc_const_param.fc;
    auto soil_moisture = loc_const_param.soil_moisture;
    auto slop = loc_const_param.slop;
    auto ks = loc_const_param.ks * time_interval;
    auto b = loc_const_param.b;
    auto k = GetK(soil_moisture, sat, b, ks);
    auto zs = loc_const_param.zs;
    auto lateral_in_q = loc_state_param.lateral_in_flow;
    auto direct_factor = GetDirectFactor(loc_const_param.d8);
    auto ep = loc_const_param.ep;
    auto v = loc_const_param.v;
    // calculate evaporation
    if (soil_moisture > fc)
    {
      loc_state_param.Ea = ep * v;
    }
    else if (soil_moisture > loc_const_param.wl)
    {
      loc_state_param.Ea = (1 - v) * Ep * (curr - loc_const_param.wl) / (fc - loc_const_param.wl);
    }
    else
    {
      loc_state_param.Ea = 0;
    }
    if (soil_moisture > fc)
    {
      // Percolate out quantity
      auto percolate_out_q = 0.001 * (k + loc_state_param.per_mm) / 2 * cell_size * cell_size;
      auto lateral_out_q = 0.001 * (k * slop + loc_state_param.lat_mm) / 2 * direct_factor * cell_size * zs * 0.001;
      auto excess_q = 0.001 * zs * (cur - fc) * cell_size * cell_size;
    }
  }
}
