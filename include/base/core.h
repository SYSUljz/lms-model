#pragma once
#include <cstddef>

#include "direct.h"

namespace lms {
namespace core {
enum class Label { Soil, Channel, Reservoir };

template <typename T>
struct ModelMeta {
  size_t width_;
  size_t heigh_;
  // raster resolution (m)
  size_t cell_size_;
  size_t time_interval_s_;
  size_t runoff_dt_s_;
  size_t confluence_steps_;
  T raster_min_lat_;
  T raster_max_lat_;
  T raster_min_long_;
  T raster_max_long_;
};

template <typename T>
struct GlobalParam {
  T soil_alpha_;
  // groundwater recession coefficient
  T baseflow_coff;
  T v;
  T manning;
  T ss;
  T b;
};

template <typename T>
struct ConstParam {
  // Topographic factors
  Label label;
  Direct8 d8;
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
  ep       — potential evapotranspiration (optional) */
  T sat;
  T fc;
  T wl;
  T zs;
  T ks;
  // Clapp-Hornberger pore size distribution index
  T b;
  T n;
  T v;
  T bs;
  // channel bottom width (m)
  T bw;
  T manning;
  T ep;
};
/* state factors:
  soil_moisture        — current soil moisture
actual_evaporate         — actual evapotranspiration
runoff     — surface runoff
lat_mm     — lateral flow
per_mm     — percolation
groundwater_mm         — groundwater flow
QPrevT     — flow at previous time step
QPrevX     — flow at previous spatial step
Qcurr      — current flow
waterLevel — water level
temp       — intermediate computation variable
*/

template <typename T>
struct StateParam {
  T soil_moisture {0};
  T actual_evaporate {0};
  T runoff {0};
  T lat_mm {0};
  T per_mm {0};
  T groundwater_mm {0};
  T prev_t_flow {0};
  T upstream_in_flow {0};
  T water_level {0};
  T temp {0};
  T lateral_in_flow_mm {0};
};
}  // namespace core
}  // namespace lms