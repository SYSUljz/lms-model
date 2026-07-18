#pragma once
#include "core.h"

namespace lms {
namespace factor {
using lms::core::ConstParam;
using lms::core::GlobalParam;

template <typename T>
struct Factor {
  T sat {1};
  T fc {1};
  T wl {1};
  T ks {1};
  T zs {1};
  T b {1};
  T n {1};
  T v {1};
  T bs {1};
  T bw {1};
  T ep {1};
  T manning {1};
  T soil_alpha {1};
  T init_soil_water {1};
  T ss {1};
};

}  // namespace factor
}  // namespace lms
