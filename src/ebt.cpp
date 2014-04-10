#include "ebt.h"

namespace model {

EBT::EBT(Parameters p)
  : parameters(p),
    patch(*parameters.get()),
    ode_solver(this, parameters->control.ode_control),
    schedule(patch.size()) {
  if (!util::identical(parameters->patch_area, 1.0))
    Rcpp::stop("Patch area must be exactly 1 for the EBT");
}

EBT::EBT(Parameters *p)
  : parameters(p),
    patch(*parameters.get()),
    ode_solver(this, parameters->control.ode_control),
    schedule(patch.size()) {
  if (!util::identical(parameters->patch_area, 1.0))
    Rcpp::stop("Patch area must be exactly 1 for the EBT");
}

void EBT::run() {
  reset();
  while (!complete()) {
    run_next();
  }
}

// TODO: The return value here is a bit of a hack for
// build.schedule().  It would be nice to fix this.  One alternative
// would be to update the CohortSchedule to naturally deal with the
// case of multiple introductions per time.  The other would be to say
// how many events we skipped over perhaps.  This will do for now
// though.
std::vector<int> EBT::run_next() {
  std::vector<int> ret;
  const double t0 = get_time();
  CohortSchedule::Event e = schedule.next_event();
  while (true) {
    if (!util::identical(t0, e.time_introduction()))
      Rcpp::stop("Start time not what was expected");
    ret.push_back(e.r_species_index());
    // TODO: This causes the light environment to be computed multiple
    // times when multiple residents are added.  The environment is
    // not recomputed when mutants are added though.
    //
    // The other option would be to create a new vector and pass that
    // through to EBT::add_seedlings().
    patch.add_seedling(e.species_index);
    schedule.pop();
    if (e.time_end() > t0 || complete()) {
      break;
    } else {
      e = schedule.next_event();
    }
  }

  ode_solver.set_state_from_problem();
  if (schedule.fixed_times())
    ode_solver.advance_fixed(e.times);
  else
    ode_solver.advance(e.time_end());

  return ret;
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
  // NOTE: Possibly wasteful, but at the moment we only return these
  // one at a time, so maybe not.
  const double tot_leaf_area  = patch.leaf_area_above(0.0);
  return patch.at(species_index).leaf_area_error(tot_leaf_area);
}

std::vector<double> EBT::fitness_error(size_t species_index) const {
  // NOTE: Possibly wasteful, possibly nbd.  If we grab all fitness at
  // once, then we should only compute this scaling factor once.
  double tot_seed_out = 0.0;
  for (size_t i = 0; i < patch.size(); ++i)
    tot_seed_out += fitness(i);
  return util::local_error_integration(schedule.times(species_index),
				       fitness_cohort(species_index),
				       tot_seed_out);
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

bool EBT::complete() const {
  return schedule.remaining() == 0;
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
