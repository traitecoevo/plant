// -*-c++-*-
#ifndef PLANT_PLANT_SCM_H_
#define PLANT_PLANT_SCM_H_

#include <plant/patch.h>
#include <plant/cohort_schedule.h>
#include <plant/ode_solver.h>
#include <plant/scm_utils.h>

using namespace Rcpp;

namespace plant {

template <typename T, typename E>
class SCM {
public:
  typedef T             strategy_type;
  typedef E             environment_type;
  typedef Individual<T,E>    individual_type;
  typedef Cohort<T,E>   cohort_type;
  typedef Species<T,E>  species_type;
  typedef Patch<T,E>    patch_type;
  typedef Parameters<T,E> parameters_type;


  SCM(parameters_type p);

  void run();
  std::vector<size_t> run_next();

  double time() const;
  void reset();
  bool complete() const;

  // * Output total seed rain calculation (not per capita)
  double seed_rain(size_t species_index) const;
  std::vector<double> seed_rains() const;

  // * R interface
  std::vector<util::index> r_run_next();
  parameters_type r_parameters() const {return parameters;}
  const patch_type& r_patch() const {return patch;}
  // TODO: These are liable to change to return all species at once by
  // default.  The pluralisation difference between
  // SCM::r_competition_effect_error and Species::r_competition_effects_error will get
  // dealt with then.
  double              r_seed_rain(util::index species_index) const;
  std::vector<double> r_seed_rain_cohort(util::index species_index) const;
  std::vector<double> r_seed_rain_error(util::index species_index) const;
  std::vector<std::vector<double> > r_seed_rain_error() const;
  std::vector<double> r_competition_effect_error(util::index species_index) const;
  std::vector<double> r_ode_times() const;
  bool r_use_ode_times() const;
  void r_set_use_ode_times(bool x);

  CohortSchedule r_cohort_schedule() const {return cohort_schedule;}
  void r_set_cohort_schedule(CohortSchedule x);
  void r_set_cohort_schedule_times(std::vector<std::vector<double> > x);

private:
  double seed_rain_total() const;
  std::vector<double> seed_rain_cohort(size_t species_index) const;

  parameters_type parameters;
  patch_type patch;
  CohortSchedule cohort_schedule;
  ode::Solver<patch_type> solver;
};

template <typename T, typename E>
SCM<T,E>::SCM(parameters_type p)
  : parameters(p),
    patch(parameters),
    cohort_schedule(make_cohort_schedule(parameters)),
    solver(patch, make_ode_control(p.control)) {
  parameters.validate();
  if (!util::identical(parameters.patch_area, 1.0)) {
    util::stop("Patch area must be exactly 1 for the SCM");
  }
}

template <typename T, typename E>
void SCM<T,E>::run() {
  reset();
  while (!complete()) {
    run_next();
  }
}

template <typename T, typename E>
std::vector<size_t> SCM<T,E>::run_next() {
  std::vector<size_t> ret;
  const double t0 = time();

  CohortSchedule::Event e = cohort_schedule.next_event();
  while (true) {
    if (!util::identical(t0, e.time_introduction())) {
      util::stop("Start time not what was expected");
    }
    ret.push_back(e.species_index);
    cohort_schedule.pop();
    if (e.time_end() > t0 || complete()) {
      break;
    } else {
      e = cohort_schedule.next_event();
    }
  }
  patch.add_seeds(ret);

  const bool use_ode_times = cohort_schedule.using_ode_times();
  solver.set_state_from_system(patch);
  if (use_ode_times) {
    solver.advance_fixed(patch, e.times);
  } else {
    solver.advance(patch, e.time_end());
  }

  return ret;
}

template <typename T, typename E>
double SCM<T,E>::time() const {
  return patch.time();
}

// NOTE: solver.reset() will set time within the solver to zero.
// However, there is no other current way of setting the time within
// the solver.  It might be better to add a set_time method within
// ode::Solver, and then here do explicitly ode_solver.set_time(0)?
template <typename T, typename E>
void SCM<T,E>::reset() {
  patch.reset();
  cohort_schedule.reset();
  solver.reset(patch);
}

template <typename T, typename E>
bool SCM<T,E>::complete() const {
  return cohort_schedule.remaining() == 0;
}

template <typename T, typename E>
double SCM<T,E>::seed_rain(size_t species_index) const {
  return util::trapezium(cohort_schedule.times(species_index),
                         seed_rain_cohort(species_index));
}

template <typename T, typename E>
std::vector<double> SCM<T,E>::seed_rains() const {
  std::vector<double> ret;
  for (size_t i = 0; i < patch.size(); ++i) {
    ret.push_back(seed_rain(i));
  }
  return ret;
}

template <typename T, typename E>
std::vector<util::index> SCM<T,E>::r_run_next() {
  return util::index_vector(run_next());
}

template <typename T, typename E>
double SCM<T,E>::r_seed_rain(util::index species_index) const {
  return seed_rain(species_index.check_bounds(patch.size()));
}

template <typename T, typename E>
std::vector<double>
SCM<T,E>::r_seed_rain_cohort(util::index species_index) const {
  return seed_rain_cohort(species_index.check_bounds(patch.size()));
}

template <typename T, typename E>
std::vector<double> SCM<T,E>::r_seed_rain_error(util::index species_index) const {
  // TODO: This causes this to happen too often, given we usually get
  // all the errors I think? (see TODO in class definition)
  double tot_seed_out = seed_rain_total();
  const size_t idx = species_index.check_bounds(patch.size());
  return util::local_error_integration(cohort_schedule.times(idx),
                                       seed_rain_cohort(idx),
                                       tot_seed_out);
}

template <typename T, typename E>
std::vector<std::vector<double> > SCM<T,E>::r_seed_rain_error() const {
  std::vector<std::vector<double> > ret;
  double tot_seed_out = seed_rain_total();
  for (size_t i = 0; i < patch.size(); ++i) {
    ret.push_back(util::local_error_integration(cohort_schedule.times(i),
                                                seed_rain_cohort(i),
                                                tot_seed_out));
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> SCM<T,E>::r_competition_effect_error(util::index species_index) const {
  // TODO: I think we need to scale this by total area; that should be
  // computed for everything so will get passed in as an argument.
  // const double tot_competition_effect  = patch.compute_competition(0.0);
  const size_t idx = species_index.check_bounds(patch.size());
  return patch.r_competition_effect_error(idx);
}

template <typename T, typename E>
std::vector<double> SCM<T,E>::r_ode_times() const {
  return solver.get_times();
}

template <typename T, typename E>
bool SCM<T,E>::r_use_ode_times() const {
  return cohort_schedule.using_ode_times();
}

template <typename T, typename E>
void SCM<T,E>::r_set_use_ode_times(bool x) {
  cohort_schedule.r_set_use_ode_times(x);
}

template <typename T, typename E>
void SCM<T,E>::r_set_cohort_schedule(CohortSchedule x) {
  if (patch.ode_size() > 0) {
    util::stop("Cannot set schedule without resetting first");
  }
  util::check_length(x.get_n_species(), patch.size());
  cohort_schedule = x;

  // Update these here so that extracting Parameters would give the
  // new schedule, this making Parameters sufficient.
  parameters.cohort_schedule_max_time = cohort_schedule.get_max_time();
  parameters.cohort_schedule_times = cohort_schedule.get_times();
}

template <typename T, typename E>
void SCM<T,E>::r_set_cohort_schedule_times(std::vector<std::vector<double> > x) {
  if (patch.ode_size() > 0) {
    util::stop("Cannot set schedule without resetting first");
  }
  cohort_schedule.set_times(x);
  parameters.cohort_schedule_times = x;
}

template <typename T, typename E>
double SCM<T,E>::seed_rain_total() const {
  double tot = 0.0;
  for (size_t i = 0; i < patch.size(); ++i) {
    tot += seed_rain(i);
  }
  return tot;
}

template <typename T, typename E>
std::vector<double> SCM<T,E>::seed_rain_cohort(size_t species_index) const {
  const std::vector<double> times = cohort_schedule.times(species_index);
  const Disturbance& disturbance_regime = patch.disturbance_regime();
  const double S_D = parameters.strategies[species_index].S_D;
  const double scal = S_D * parameters.seed_rain[species_index];
  std::vector<double> seeds = patch.at(species_index).seeds();
  for (size_t i = 0; i < seeds.size(); ++i) {
    seeds[i] *= disturbance_regime.density(times[i]) * scal;
  }
  return seeds;
}

}

#endif
