#pragma once

#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "base/core.h"

namespace lms {
namespace rain {

using lms::core::ModelMeta;

template <typename T>
class Station {
 public:
  std::size_t id_;
  T geo_lat_;
  T geo_long_;
  std::size_t raster_x_ {0};
  std::size_t raster_y_ {0};
  std::string name_;

  Station(size_t id, T geo_lat, T geo_long, std::string name, const ModelMeta<T>& meta)
      : id_(id), geo_lat_(geo_lat), geo_long_(geo_long), name_(std::move(name)) {}
  void GeoPos2RasterPos(const ModelMeta<T>& meta) {
    double x_val = std::round((static_cast<double>(geo_long_) - static_cast<double>(meta.raster_min_long_)) /
                              (static_cast<double>(meta.raster_max_long_) - static_cast<double>(meta.raster_min_long_)) *
                              static_cast<double>(meta.width_));
    double y_val = std::round((static_cast<double>(geo_lat_) - static_cast<double>(meta.raster_min_lat_)) /
                              (static_cast<double>(meta.raster_max_lat_) - static_cast<double>(meta.raster_min_lat_)) *
                              static_cast<double>(meta.heigh_));
    long long rx = static_cast<long long>(x_val);
    long long ry = static_cast<long long>(meta.heigh_) - 1 - static_cast<long long>(y_val);

    // Clamp to valid range [0, width - 1] and [0, heigh - 1]
    if (rx < 0) rx = 0;
    if (rx >= static_cast<long long>(meta.width_)) rx = static_cast<long long>(meta.width_) - 1;
    if (ry < 0) ry = 0;
    if (ry >= static_cast<long long>(meta.heigh_)) ry = static_cast<long long>(meta.heigh_) - 1;

    raster_x_ = static_cast<size_t>(rx);
    raster_y_ = static_cast<size_t>(ry);
  }
};

template <typename T, std::size_t station_cnt>
class Rainfall {
 public:
  std::array<T, station_cnt> rainfall_row_ {};
};

template <typename T, std::size_t station_cnt>
class RainfallEvent {
 public:
  size_t duration_;
  std::vector<Rainfall<T, station_cnt>> rainfall_event_ {};
  size_t time_interval_s_;
};

}  // namespace rain
}  // namespace lms
