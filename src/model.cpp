template <typename T>
struct lms_cell
{
  // Topographic factors
  T label;
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
  // Ep       — potential evapotranspiration (optional) */
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

struct model
{
  using std::vector;

public:
  Task<T> simulate() {};
  T get_result() {};

private:
  vector<lms_cell> raster;
};