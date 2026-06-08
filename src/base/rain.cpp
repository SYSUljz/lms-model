#include <cmath>

#include "base/rain.h"

namespace lms {
namespace rain {

using lms::core::ModelMeta;

template <typename T>
void Station<T>::GeoPos2RasterPos(const ModelMeta<T>& meta) {
  raster_x_ = static_cast<size_t>(
      std::round((geo_long_ - meta.raster_min_long_) / (meta.raster_max_long_ - meta.raster_min_long_) * meta.width_));
  raster_y_ = static_cast<size_t>(
      std::round((geo_lat_ - meta.raster_min_lat_) / (meta.raster_max_lat_ - meta.raster_min_lat_) * meta.heigh_));
}

}  // namespace rain
}  // namespace lms
