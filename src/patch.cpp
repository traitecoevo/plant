#include "patch.h"

namespace model {

template <>
void Patch<CohortTop>::step() {
  step_deterministic();
}

}
