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
    raster_x_ = static_cast<size_t>(
        std::round((static_cast<double>(geo_long_) - static_cast<double>(meta.raster_min_long_)) /
                   (static_cast<double>(meta.raster_max_long_) - static_cast<double>(meta.raster_min_long_)) *
                   static_cast<double>(meta.width_)));
    raster_y_ = static_cast<size_t>(
        std::round((static_cast<double>(geo_lat_) - static_cast<double>(meta.raster_min_lat_)) /
                   (static_cast<double>(meta.raster_max_lat_) - static_cast<double>(meta.raster_min_lat_)) *
                   static_cast<double>(meta.heigh_)));
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
