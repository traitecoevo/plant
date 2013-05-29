#include "ebt.h"

namespace model {

EBT::EBT(Parameters p)
  : patch(p),
    ode_solver(&patch),
    schedule(patch.size()) {
}

EBT::EBT(Parameters *p)
  : patch(p),
    ode_solver(&patch),
    schedule(patch.size()) {
}

void EBT::step() {
  std::vector<double> y(patch.ode_size());
  patch.ode_values(y.begin());
  ode_solver.set_state(y, get_time());
  ode_solver.step();
}

void EBT::run_next() {
  if (schedule.remaining() == 0)
    ::Rf_error("Already reached end of schedule");
  const CohortSchedule::Event e = schedule.next_event();
  ode_solver.advance(e.time);
  patch.add_seedling(e.cohort);
  schedule.pop();
}

double EBT::get_time() const {
  return ode_solver.get_time();
}

Patch<CohortTop> EBT::r_patch() const {
  return patch;
}

CohortSchedule EBT::r_cohort_schedule() const {
  return schedule;
}

void EBT::r_set_cohort_schedule(CohortSchedule x) {
  util::check_length(x.types(), patch.size());
  schedule = x;
}

}
