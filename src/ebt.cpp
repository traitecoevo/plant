#include "ebt.h"

namespace model {

EBT::EBT(Parameters p)
  : patch(p),
    ode_solver(&patch, p.control.ode_control),
    schedule(patch.size()) {
}

EBT::EBT(Parameters *p)
  : patch(p),
    ode_solver(&patch, p->control.ode_control),
    schedule(patch.size()) {
}

void EBT::run_next() {
  const bool events_remaining = schedule.remaining() > 0;
  const CohortSchedule::Event e = schedule.next_event();
  const double next_time = e.time_introduction();
  const bool done = !events_remaining &&
    (!util::is_finite(next_time) || get_time() >= next_time);
  if (done)
    ::Rf_error("Already reached end of schedule");

  advance(next_time);
  if (events_remaining) {
    add_seedling(e.cohort);
    schedule.pop();
  }
}

double EBT::get_time() const {
  return ode_solver.get_time();
}

// NOTE: ode_solver.reset() will set time within the solver to zero.
// However, there is no other current way of setting the time within
// the solver.  It might be better to add a set_time method within
// ode::Solver, and then here do explicitly ode_solver.set_time(0)?
void EBT::reset() {
  patch.reset();
  schedule.reset();
  ode_solver.reset();
}

Patch<CohortTop> EBT::r_patch() const {
  return patch;
}

CohortSchedule EBT::r_cohort_schedule() const {
  return schedule;
}

void EBT::r_set_cohort_schedule(CohortSchedule x) {
  util::check_length(x.get_n_species(), patch.size());
  schedule = x;
}

void EBT::add_seedling(size_t species_index) {
  patch.add_seedling(species_index);
}

void EBT::advance(double time) {
  ode_solver.set_state_from_problem();
  ode_solver.advance(time);
}

}
