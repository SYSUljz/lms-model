#include <array>
#include <string>

#include "base/core.h"
namespace lms {
namespace rain {
template <typename T>
class Station {
  std::size_t id_;
  std::size_t idx_;
  T geo_lat_;
  T geo_long_;
  std::size_t raster_x_ {0};
  std::size_t raster_y_ {0};
  std::string name_;

  Station(size_t id, size_t idx, T geo_lat, T geo_long, std::string& name, const ModelMeta<T>& meta)
      : id_(id), idx_(idx), geo_lat_(geo_lat), geo_long_(geo_long), name_(name) {}
  void GeoPos2RasterPos(const ModelMeta<T>& meta);
};
template <typename T, std::size_t station_cnt>
class Rainfall {
  std::array<T, station_cnt> rainfall_row_ {};
};
template <typename T, std::size_t station_cnt>
class RainfallEvent {
  size_t duration_;
  std::vector<Rainfall<T, station_cnt>> rainfall_event_ {};
  size_t time_interval_s_;
};

}  // namespace rain
}  // namespace lms