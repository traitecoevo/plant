#include "ebt.h"

namespace model {

EBT::EBT(Parameters p)
  : patch(p),
    ode_solver(&patch),
    schedule(p.size()),
    time(0.0) {
}

EBT::EBT(Parameters *p)
  : patch(p),
    ode_solver(&patch),
    schedule(p->size()),
    time(0.0) {
}

void EBT::step() {
  std::vector<double> y(patch.ode_size());
  patch.ode_values(y.begin());
  ode_solver.set_state(y, time);
  ode_solver.step();
  // or get this from environment?
  time = ode_solver.get_time();
}

Patch<CohortTop> EBT::r_patch() const {
  return patch;
}

}
