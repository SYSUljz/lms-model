#include <base/rain.h>
namespace lxh {
namespace rain {
using lxh::core::ModelMeta;
template <typename T>
void Station::GeoPis2RasterPos(const ModelMeta<T>& meta) {
  raster_x_ = static_cast<size_t>(
      std::round((geo_long_ - meta.raster_min_long_) / (meta.raster_max_long_ - meta.raster_min_long_) * meta.width_));
  raster_y_ = static_cast<size_t>(
      std::round((geo_lat_ - meta.raster_min_lat_) / (meta.raster_max_lat_ - meta.raster_min_lat_) * meta.heigh_));
}
}  // namespace rain
}  // namespace lxh