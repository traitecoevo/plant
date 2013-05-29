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

void EBT::run_next() {
  const CohortSchedule::Event e = schedule.next_event();
  if (e.time == R_PosInf)
    ::Rf_error("Already reached end of schedule");
  ode_solver.advance(e.time);
  patch.add_seedling(e.cohort);
}

Patch<CohortTop> EBT::r_patch() const {
  return patch;
}

CohortSchedule EBT::r_cohort_schedule() const {
  return schedule;
}

void EBT::r_set_cohort_schedule(CohortSchedule x) {
  schedule = x;
}

}
