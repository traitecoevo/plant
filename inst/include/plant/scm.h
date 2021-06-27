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

  // * Output total offspring calculation (not per capita)
  std::vector<double> net_reproduction_ratio_by_cohort_weighted(size_t species_index) const;
  double net_reproduction_ratio_for_species(size_t species_index) const;
  std::vector<double> net_reproduction_ratios() const;
  std::vector<double> offspring_production() const;

  // * R interface
  std::vector<util::index> r_run_next();
  parameters_type r_parameters() const {return parameters;}
  const patch_type& r_patch() const {return patch;}

  void r_set_state(double time,
                   const std::vector<double>& state,
                   const std::vector<size_t>& n)
                   {patch.r_set_state(time, state, n, std::vector<double>());};

  // TODO: These are liable to change to return all species at once by
  // default.  The pluralisation difference between
  // SCM::r_competition_effect_error and Species::r_competition_effects_error will get
  // dealt with then.
  double r_net_reproduction_ratio_for_species(util::index species_index) const;
  std::vector<std::vector<double> > r_net_reproduction_ratio_errors() const;
  std::vector<double> r_competition_effect_error(util::index species_index) const;
  std::vector<double> r_ode_times() const;
  bool r_use_ode_times() const;
  void r_set_use_ode_times(bool x);

  CohortSchedule r_cohort_schedule() const {return cohort_schedule;}
  void r_set_cohort_schedule(CohortSchedule x);
  void r_set_cohort_schedule_times(std::vector<std::vector<double> > x);

private:
  double total_offspring_production() const;

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
  patch.introduce_new_cohorts(ret);

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
std::vector<util::index> SCM<T,E>::r_run_next() {
  return util::index_vector(run_next());
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


// Offspring production, equal to overall fitness scaled by the birth rate
template <typename T, typename E>
std::vector<double> SCM<T,E>::offspring_production() const {
  std::vector<double> ret;
  for (size_t i = 0; i < patch.size(); ++i) {
    ret.push_back(net_reproduction_ratio_for_species(i) *
                  parameters.birth_rate[i]);
  }
  return ret;
}

// Overall fitness
template <typename T, typename E>
std::vector<double> SCM<T,E>::net_reproduction_ratios() const {
  std::vector<double> ret;
  for (size_t i = 0; i < patch.size(); ++i) {
    ret.push_back(net_reproduction_ratio_for_species(i));
  }
  return ret;
}

// Integrate over lifetime fitness of individual cohorts
template <typename T, typename E>
double SCM<T,E>::net_reproduction_ratio_for_species(size_t species_index) const {
  return util::trapezium(cohort_schedule.times(species_index),
                         net_reproduction_ratio_by_cohort_weighted(species_index));
}

// R interface method
template <typename T, typename E>
double SCM<T,E>::r_net_reproduction_ratio_for_species(util::index species_index) const {
  return net_reproduction_ratio_for_species(species_index.check_bounds(patch.size()));
}

// Cohort fitness within a meta-population of patches
template <typename T, typename E>
std::vector<double> SCM<T,E>::net_reproduction_ratio_by_cohort_weighted(size_t species_index) const {
  // cohort introduction times
  const std::vector<double> times = cohort_schedule.times(species_index);

  // retrieve lifetime fitness for each cohort
  std::vector<double> net_reproduction_ratio_by_cohort_weighted =
                        patch.at(species_index).net_reproduction_ratio_by_cohort();

  // weight by probabilty of reproduction
  for (size_t i = 0; i < net_reproduction_ratio_by_cohort_weighted.size(); ++i) {
    net_reproduction_ratio_by_cohort_weighted[i] *=
      patch.survival_weighting->density(times[i]) * // probability of landing in patch of a given age
      parameters.strategies[species_index].S_D; // probability of survival during dispersal (assumed constant)
  }

  return net_reproduction_ratio_by_cohort_weighted;
}

// Sum up all offspring produced
template <typename T, typename E>
double SCM<T,E>::total_offspring_production() const {
  double total = 0.0;
  std::vector<double> offspring = offspring_production();
  for (size_t i = 0; i < patch.size(); ++i) {
    total += offspring[i];
  }
  return total;
}

// Check integration errors
template <typename T, typename E>
std::vector<std::vector<double> > SCM<T,E>::r_net_reproduction_ratio_errors() const {
  std::vector<std::vector<double> > ret;
  double total_offspring = total_offspring_production();
  for (size_t i = 0; i < patch.size(); ++i) {
    ret.push_back(util::local_error_integration(cohort_schedule.times(i),
                                                net_reproduction_ratio_by_cohort_weighted(i),
                                                total_offspring));
  }
  return ret;
}

}

#endif
