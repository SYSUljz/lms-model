#include <cmath>
/// @brief Clapp-Hornberger unsaturated hydraulic conductivity formula
/// @param cur Current volumetric soil moisture θ
/// @param sat Saturated soil moisture θ_sat
/// @param b Clapp-Hornberger pore-size distribution index (soil texture parameter）
/// @param ks Saturated hydraulic conductivity (mm/step)
/// @return Hydraulic conductivity K at current moisture (mm/step)
template <typename T>
T inline GetK(T cur, T sat, T b, T ks)
{
  return ks * std::pow(cur / sat, 2 * b + 3);
};

template <typename T>
T GetSoilAlpha(T manning, T cell_size, T slope)
{
  if (slope == 0)
      [[unlikely]]
  {
    slope = static_cast<T>(1e-6);
  }
  T alpha = std::pow(manning * std::pow(cell_size, static_cast<T>(0.666667)) * std::pow(slope, static_cast<T>(-0.5)) / static_cast<T>(3600), static_cast<T>(0.6));
  return alpha;
}