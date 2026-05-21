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