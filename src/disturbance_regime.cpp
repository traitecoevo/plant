#include <plant/disturbance_regime.h>

namespace plant {

std::vector<double> DisturbanceRegime::r_density(std::vector<double> time) const {
  std::vector<double> ret;
  ret.reserve(time.size());
  for (const auto& t : time) {
    ret.push_back(density(t));
  }
  return ret;
}

}
