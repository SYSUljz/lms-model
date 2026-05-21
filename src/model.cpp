#include "./utils.cpp"
enum class Label
{
  Soil,
  Channel,
  Reservoir

};
struct StationAttribute
{
  size_t station_id;
};
template <typename T>
struct ConstParam
{
  // Topographic factors
  Label label;
  T d8;
  T slop;
  /* soil factors
  sat      — saturated water content
  fc       — field capacity
  wl       — wilting coefficient
  zs       — soil layer thickness
  ks       — saturated hydraulic conductivity
  b        — pore size distribution index
  n        — vegetation coverage
  v        — vegetation interception coefficient
  ss       — slope roughness
  bs       — channel slope
  bw       — channel width
  manning  — Manning's roughness coefficient
  Ep       — potential evapotranspiration (optional) */
  T sat;
  T fc;
  T wl;
  T zs;
  T ks;
  T b;
  T n;
  T v;
  T ss;
  T bs;
  T bw;
  T manning;
  T Ep;
};
/* state factors:
  cur        — current soil moisture
Ea         — actual evapotranspiration
runoff     — surface runoff
lat_mm     — lateral flow
per_mm     — percolation
Qg         — groundwater flow
QPrevT     — flow at previous time step
QPrevX     — flow at previous spatial step
Qcurr      — current flow
waterLevel — water level
temp       — intermediate computation variable
*/
template <typename T>
struct StateParam
{
  T cur;
  T Ea;
  T runoff;
  T lat_mm;
  T per_mm;
  T Qg;
  T QPrevT;
  T Qcurr;
  T waterLevel;
  T temp;
};

template <typename T>
class ConstRaster
{
  ConstRaster(const ConstRaster &raster) = delete;
  ConstRaster &operator=(const ConstRaster &raster) = delete;

private:
  std::vector<ConstParam> raster_;
  size_t weith_;
  size_t heigh_;
};

template <typename T>
class StateRaster
{
private:
  std::vector<StateParam> raster_;
  size_t weith_;
  size_t heigh_;
};
template <typename T>
class model
{

public:
  Task<T> SimulateOneStep(StateRaster<T> &temp_state_param_, const ConstRaster<T> &const_param_, std::vector<T> station_rain)
  {

    for (auto &[idx, item] : std::views::enumerate(iter_order_))
    {
      int target_idx = target_idx_[idx];
      ProcessCell(&temp_state_param_.raster_[idx], &const_param_.raster_[idx], &temp_state_param_.raster_[target_idx], time_interval_s_);
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
  auto FlowGeneration(int idx, int target_idx, StateRaster<T> &state_param_, const ConstRaster<T> &const_param_) {

  };
  auto FlowConfluence(int idx, int target_idx, StateRaster<T> &state_param_, const ConstRaster<T> &const_param_) {};

private:
  StateRaster<T> state_param_;
  const ConstRaster<T> const_param_;
  std::vector<int> iter_order_;
  std::vector<int> target_idx_;
  std::vector<int> station_id_;
  // rainfall_[time][station]
  std::vector<std::vector<T>> rainfall_;
  int time_interval_s_;
  int rainfall_data_length_;
};

template <typename T>
auto ProcessCell(StateParam<T> &loc_state, ConstParam<T> &const_param, StateParam<T> &target_state, int time_interval)
{
  // copy state_param raster here
  auto local_state_param = state_param;
  auto target_idx = GetTargetIdx(idx);
  FlowGeneration(int idx, int target_idx, StateRaster<T> &state_param_, const ConstRaster<T> &const_param_, time_interval);
  FlowConfluence(int idx, int target_idx, StateRaster<T> &state_param_, const ConstRaster<T> &const_param_, time_interval);
};