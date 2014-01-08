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
    Rcpp::stop("Start time not what was expected");
  patch.add_seedling(e.species_index);
  ode_solver.set_state_from_problem();
  if (schedule.fixed_times())
    ode_solver.advance_fixed(e.times);
  else
    ode_solver.advance(e.time_end());
  schedule.pop(); // or do at next_event()?  Only matters on error.
}

// * Fitness calculations
std::vector<double> EBT::r_fitness_cohort(size_t species_index) const {
  return fitness_cohort(util::check_bounds_r(species_index, patch.size()));
}
std::vector<double> EBT::fitness_cohort(size_t species_index) const {
  const std::vector<double> times = schedule.times(species_index);
  const Disturbance& disturbance_regime = patch.get_disturbance_regime();
  const double scale =
    parameters->Pi_0 * parameters->seed_rain[species_index];
  std::vector<double> seeds = patch.at(species_index).seeds();
  for (size_t i = 0; i < seeds.size(); ++i)
    seeds[i] *= disturbance_regime.density(times[i]) * scale;
  return seeds;
}

double EBT::fitness(size_t species_index) const {
  return util::trapezium(schedule.times(species_index),
			 fitness_cohort(species_index));
}

std::vector<double> EBT::leaf_area_error(size_t species_index) const {
  return patch.at(species_index).leaf_area_error();
}

std::vector<double> EBT::fitness_error(size_t species_index) const {
  return util::local_error_integration(schedule.times(species_index),
				       fitness_cohort(species_index));
}

std::vector<double> EBT::fitnesses() const {
  std::vector<double> w;
  for (size_t i = 0; i < patch.size(); ++i)
    w.push_back(fitness(i));
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

// * R interface

double EBT::r_fitness(size_t species_index) const {
  return fitness(util::check_bounds_r(species_index, patch.size()));
}

std::vector<double> EBT::r_fitness_error(size_t species_index) const {
  return fitness_error(util::check_bounds_r(species_index, patch.size()));
}

std::vector<double> EBT::r_leaf_area_error(size_t species_index) const {
  return leaf_area_error(util::check_bounds_r(species_index, patch.size()));
}

Patch<CohortTop> EBT::r_patch() const {
  return patch;
}

CohortSchedule EBT::r_cohort_schedule() const {
  return schedule;
}

void EBT::r_set_cohort_schedule(CohortSchedule x) {
  if (patch.ode_size() > 0)
    Rcpp::stop("Cannot set schedule without resetting first");
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
    Rcpp::stop("Cannot set schedule without resetting first");
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
