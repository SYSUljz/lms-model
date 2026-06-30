#include <cmath>
/// @brief Clapp-Hornberger unsaturated hydraulic conductivity formula
/// @param cur Current volumetric soil moisture θ
/// @param sat Saturated soil moisture θ_sat
/// @param b Clapp-Hornberger pore-size distribution index (soil texture
/// parameter）
/// @param ks Saturated hydraulic conductivity (mm/step)
/// @return Hydraulic conductivity K at current moisture (mm/step)
template <typename T>
T inline GetK(T cur, T sat, T b, T ks) {
  return ks * std::pow(cur / sat, 2 * b + 3);
};

template <typename T>
T GetSoilAlpha(T manning, T cell_size, T slope) {
  if (std::isnan(manning) || manning <= 0) manning = static_cast<T>(0.03);
  if (std::isnan(cell_size) || cell_size <= 0) cell_size = static_cast<T>(30);
  if (std::isnan(slope) || slope <= 0) [[unlikely]] {
    slope = static_cast<T>(1e-6);
  }
  T alpha = std::pow(manning * std::pow(cell_size, static_cast<T>(0.666667)) * std::pow(slope, static_cast<T>(-0.5)) /
                         static_cast<T>(3600),
                      static_cast<T>(0.6));
  if (std::isnan(alpha)) alpha = static_cast<T>(1);
  return alpha;
}

template <typename T>
T GetChannelAlpha(T manning, T cell_size, T slope) {
  if (std::isnan(manning) || manning <= 0) manning = static_cast<T>(0.03);
  if (std::isnan(cell_size) || cell_size <= 0) cell_size = static_cast<T>(30);
  if (std::isnan(slope) || slope <= 0) [[unlikely]] {
    slope = static_cast<T>(1e-6);
  }
  T alpha = std::pow(manning * std::pow(cell_size, static_cast<T>(0.666667)) * std::pow(slope, static_cast<T>(-0.5)) /
                         static_cast<T>(3600),
                      static_cast<T>(0.6));
  if (std::isnan(alpha)) alpha = static_cast<T>(1);
  return alpha;
}

/// Solving the simplified form of Saint-Venant's equations using Newton's
/// iteration method (evolution of moving/spreading waves)
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
T SolveSaintVenant(T iQ, T q, T alpha, T beta, T dT, T dX, T iQPrevX, T iQPrevT) {
  if (std::isnan(iQ) || iQ < 0) iQ = 0;
  if (std::isnan(q) || q < 0) q = 0;
  if (std::isnan(alpha) || alpha <= 0) alpha = 1.0;
  if (std::isnan(beta) || beta <= 0) beta = 0.6;
  if (std::isnan(dT) || dT <= 0) dT = 1.0;
  if (std::isnan(dX) || dX <= 0) dX = 30.0;
  if (std::isnan(iQPrevX) || iQPrevX < 0) iQPrevX = 0;
  if (std::isnan(iQPrevT) || iQPrevT < 0) iQPrevT = 0;

  T dff = iQ;
  T duu;
  int cnt = 0;
  T iQPrevT_safe = std::max(iQPrevT, static_cast<T>(1e-10));
  T const_term = alpha * std::pow(iQPrevT_safe, beta);

  while (cnt < 20) {
    cnt++;
    T iQ_safe = std::max(iQ, static_cast<T>(1e-10));
    T pow_val = std::pow(iQ_safe, beta);
    auto diff =
        iQ * dT / dX + alpha * pow_val - iQPrevX * dT / dX - const_term - q * dT;
    duu = dT / dX + alpha * beta * (pow_val / iQ_safe);

    dff = diff;
    dff /= duu;

    if (std::isnan(dff)) {
      break;
    }

    if (std::abs(dff) < static_cast<T>(1e-6)) {
      return iQ;
    }

    iQ -= dff;

    if (std::isnan(iQ) || iQ < 0) [[unlikely]] {
      iQ = std::isnan(iQ) ? 0 : std::abs(iQ);
    }
  }
  return std::isnan(iQ) ? static_cast<T>(0) : iQ;
}

/// calculate Wetted Perimeter
//
/// @param iWaterLevel  (Water depth) [m]
/// @param ibw         (Bottom width) [m]
/// @param iss        (Side slope) [radians]
/// @return            Wetted Perimeter [m]

template <typename T>
T GetChannelX(T iWaterLevel, T ibw, T iss) {
  if (std::isnan(iss) || iss <= 0) {
    iss = static_cast<T>(0.5); // safe default
  }
  T X = ibw + 2 * iWaterLevel / std::sin(iss);

  return X;
}

///
///
/// @param iQ       河道流量 (Flow rate) [m^3/s]
/// @param iManning 曼宁粗糙系数 (Manning's n)
/// @param Sf       能坡 (Friction slope) [ratio]
/// @param ibw      河底宽度 (Bottom width) [m]
/// @param iss      岸坡坡度 (Side slope) [radians]
/// @return         计算得到的水位/水深 [m]
///
template <typename T>
T solveQtoH(T iQ, T iManning, T Sf, T ibw, T iss) {
  if (std::isnan(iQ) || iQ < 0) iQ = 0;
  if (std::isnan(iManning) || iManning <= 0) iManning = 0.03;
  if (std::isnan(Sf) || Sf <= 0) Sf = static_cast<T>(1e-6);
  if (std::isnan(ibw) || ibw <= 0) ibw = 1.0;
  if (std::isnan(iss) || iss <= 0) iss = 0.5;

  T dff;
  T duu;
  T h = static_cast<T>(0.5);
  int MAXITER = 20;
  for (int j = 0; j <= MAXITER; j++) {
    T h_safe = std::max(h, static_cast<T>(1e-10));
    double a = ibw * h_safe + h_safe * h_safe / std::tan(iss);
    if (std::isnan(a) || a < 0) a = 1e-10;

    dff = std::pow(a, static_cast<T>(1.5)) - iQ * iManning * std::pow(Sf, static_cast<T>(-0.5)) *
                                                  std::pow(ibw + 2.0 * h_safe / std::sin(iss), static_cast<T>(0.66667));
    duu = static_cast<T>(1.5) * std::pow(a, static_cast<T>(0.5)) * (ibw + static_cast<T>(2.0) * h_safe / std::tan(iss)) -
          iQ * iManning * std::pow(Sf, static_cast<T>(-0.5)) *
              std::pow(ibw + static_cast<T>(2) * h_safe / std::sin(iss), static_cast<T>(-0.3333)) * static_cast<T>(1.3333) /
              std::sin(iss);
    
    if (std::isnan(duu) || duu == 0) {
      break;
    }
    
    dff /= duu;
    
    if (std::isnan(dff)) {
      break;
    }

    h -= dff;
    if (std::abs(dff) < static_cast<T>(0.000001)) {
      break;
    }
    if (std::isnan(h) || h < 0) {
      h = static_cast<T>(0.5);
    }
  }
  return std::isnan(h) ? static_cast<T>(0.5) : h;
}
