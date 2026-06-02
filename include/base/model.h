#include <memory>

#include "base/core.h"
#include "base/direct.h"
#include "base/raster.h"
namespace lxh {
namespace model {
template <typename T>
class Model {
 public:
  Model() {};
  // flowGeneration doesn't require child time step iteration
  auto SimulateOneStep(std::vector<T> station_rain);

  auto Simulate();

  // Particle Swarm Optimization
  auto PsoStep();
  bool BuildOrder();
  bool BuildTargetOrder();
  bool BuildRain();

  int GetTargetIdx(int this_idx);

 private:
  StateRaster state_param_;
  const ConstRaster const_param_;
  const ModelMeta<T> model_meta_;
  GlobalParam<T> global_param_;
  // During the simulation, a mock channel cell is added after the outlet to
  // ensure the calculation accuracy of the outlet channel element.
  // iter_order_[a]=b means the a'th step iter we process cell with b idx in
  // param raster (both state and const)
  std::vector<int> iter_order_;
  // iter_order_[a]=c means the ath step we are processing b cell the flow
  // direction of b cell is c cell
  std::vector<std::optional<int>> target_idx_;
  std::vector<int> station_id_;
  // rainfall_[time][station]
  std::vector<std::vector<T>> rainfall_;
  int rainfall_data_length_;

  template <typename T>
  auto FlowConfluenceMultiStep();
};
}  // namespace model
}  // namespace lxh
