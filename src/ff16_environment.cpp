#include <plant/models/ff16_environment.h>

namespace plant {

const double get_k_I(const FF16_Environment environment) {
  return environment.canopy.k_I;
}

}
