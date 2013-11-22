#include "ebt.h"

namespace model {

EBT::EBT(Parameters p)
  : parameters(p),
    patch(*parameters.get()),
    ode_solver(&patch, parameters->control.ode_control),
    schedule(patch.size()) {
}

EBT::EBT(Parameters *p)
  : parameters(p),
    patch(*parameters.get()),
    ode_solver(&patch, parameters->control.ode_control),
    schedule(patch.size()) {
}

void EBT::run() {
  reset();
  while (schedule.remaining() > 0) {
    run_next();
  }
}

void EBT::run_next() {
  const CohortSchedule::Event e = schedule.next_event();
  if (!util::identical(get_time(), e.time_introduction()))
    ::Rf_error("Start time not what was expected");
  patch.add_seedling(e.species_index);
  ode_solver.set_state_from_problem();
  if (schedule.fixed_times())
    ode_solver.advance_fixed(e.times);
  else
    ode_solver.advance(e.time_end());
  schedule.pop(); // or do at next_event()?  Only matters on error.
}

std::vector<double> EBT::fitness() const {
  std::vector<double> w = fitness_raw();
  const size_t n_spp = patch.size();
  const double Pi_0 = parameters->Pi_0;
  for (size_t i = 0; i < n_spp; ++i) {
    w[i] *= Pi_0 * parameters->seed_rain[i];
  }
  return w;
}

std::vector<double> EBT::fitness_raw() const {
  const size_t n_spp = patch.size();
  std::vector<double> w;
  const Disturbance& d = patch.get_disturbance_regime();
  for (size_t i = 0; i < n_spp; ++i) {
    w.push_back(model::fitness(patch.at(i), schedule.times(i), d));
  }
  return w;
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

// * Ode interface.
size_t EBT::ode_size() const {
  return patch.ode_size();
}

ode::iterator_const EBT::set_ode_values(double time,
					ode::iterator_const it) {
  return patch.set_ode_values(time, it);
}

ode::iterator EBT::ode_values(ode::iterator it) const {
  return patch.ode_values(it);
}

ode::iterator EBT::ode_rates(ode::iterator it) const {
  return patch.ode_rates(it);
}

Patch<CohortTop> EBT::r_patch() const {
  return patch;
}

CohortSchedule EBT::r_cohort_schedule() const {
  return schedule;
}

void EBT::r_set_cohort_schedule(CohortSchedule x) {
  if (patch.ode_size() > 0)
    ::Rf_error("Cannot set schedule without resetting first");
  util::check_length(x.get_n_species(), patch.size());
  schedule = x;
}

std::vector<double> EBT::r_ode_times() const {
  return ode_solver.get_times();
}

Parameters EBT::r_parameters() const {
  return *parameters.get();
}

std::vector<double> EBT::r_times(size_t species_index) const {
  return schedule.r_times(species_index);
}
void EBT::r_set_times(std::vector<double> times, size_t species_index) {
  if (patch.ode_size() > 0)
    ::Rf_error("Cannot set schedule without resetting first");
  schedule.r_set_times(times, species_index);
}

Rcpp::List EBT::r_get_state() const {
  return Rcpp::List::create(Rcpp::_["patch"]    = patch.r_get_state(),
			    Rcpp::_["schedule"] = schedule.r_get_state());
}

// Note that we use force here because when using this we are
// generally happy about changing the dimensions.
void EBT::r_set_state(Rcpp::List x) {
  Rcpp::List patch_state = x["patch"],
    schedule_state = x["schedule"];
  patch.r_force_state(patch_state);
  schedule.r_set_state(schedule_state);
}

}
