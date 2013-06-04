#include "patch.h"

namespace model {

template <>
void Patch<CohortTop>::r_step() {
  r_step_deterministic();
}

}
