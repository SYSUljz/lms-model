#pragma once
#include <vector>
#include <numeric>
#include <cmath>
namespace lms {
namespace objfunc {

template <typename T>
T CalculateNSE(const std::vector<T>& sim, const std::vector<T>& obs) {
  T mean_obs = std::accumulate(obs.begin(), obs.end(), T {0}) / obs.size();
  T num = 0.0;
  T den = 0.0;
  for (size_t i = 0; i < sim.size(); ++i) {
    num += std::pow(sim[i] - obs[i], 2);
    den += std::pow(obs[i] - mean_obs, 2);
  }
  return T {1} - (num / den);
}
}  // namespace objfunc
}  // namespace lms
