#pragma once
namespace lms {
namespace flow {
template <typename T>
struct Flow {
  size_t duration;
  std::vector<T> data {};
  size_t time_interval_s;
};
}  // namespace flow
}  // namespace lms
