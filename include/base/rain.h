#include <array>
#include <string>
namespace lms {
namespace rain {
class Station {
  std::size_t id_;
  std::size_t idx_;
  std::string name_;
};
template <typename T, std::size_t station_cnt>
class Rainfall {
  std::array<T, station_cnt> rainfall_row_ {};
};
template <typename T, std::size_t station_cnt>
class RainfallEvent {
  size_t duration_;
  std::vector<Rainfall<T, station_cnt>> rainfall_row_ {};
};

}  // namespace rain
}  // namespace lms