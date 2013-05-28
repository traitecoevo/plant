#include "ebt.h"

namespace model {

EBT::EBT(Parameters p)
  : patch(p),
    // ode_solver(*patch),
    schedule(p.size()) {
}

EBT::EBT(Parameters *p)
  : patch(p),
    // ode_solver(*patch),
    schedule(p->size()) {
}

Patch<CohortTop> EBT::r_patch() const {
  return patch;
}

}
