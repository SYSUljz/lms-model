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

/// Solving the simplified form of Saint-Venant's equations using Newton's iteration method (evolution of moving/spreading waves)
/// @param iQ Initial guess/Target flow at the current time [m^3/h]
/// @param q Lateral inflow per unit length [m^3/h / m]
/// @param alpha Evolution equation parameter alpha
/// @param beta Evolution equation parameter beta (usually 0.6 for slopes)
/// @param dT Time step [h]
/// @param dX Space step/grid spacing [m]
/// @param iQPrevX Flow at upstream node at current time [m^3/h]
/// @param iQPrevT Flow at current node at previous time [m^3/h]
/// @return The current grid flow rate is [m^3/h]
///
template <typename T>
T SolveSaintVenant(T iQ, T q, T alpha, T beta, T dT, T dX, T iQPrevX, T iQPrevT)
{
  int cnt = 0;
  while (cnt < 20)
  {
    cnt++;
    auto diff = iQ * dT / dX + alpha * std::pow(iQ, beta) - iQPrevX * dT / dX - alpha * std::pow(iQPrevT, beta) - q * dT;
    duu = dT / dX + alpha * beta * std::pow(iQ, beta - 1);

    dff /= duu;

    if (std::abs(dff) < static_cast<T>(1e-6))
    {
      return iQ;
    }

    iQ -= dff;

    if (iQ < 0) [[unlikely]]
    {
      iQ = std::abs(iQ);
    }
  }
}

/// calculate Wetted Perimeter
//
/// @param iWaterLevel  (Water depth) [m]
/// @param ibw         (Bottom width) [m]
/// @param iss        (Side slope) [radians]
/// @return            Wetted Perimeter [m]

template <typename T>
T GetChannelX(T iWaterLevel, T ibw, T iss)
{
  T X = ibw + 2 * iWaterLevel / Math.sin(iss);

  return X;
}
