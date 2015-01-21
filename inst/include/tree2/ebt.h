// -*-c++-*-
#ifndef TREE_EBT_H_
#define TREE_EBT_H_

#include <tree2/patch.h>
#include <tree2/cohort_schedule.h> // may move to Parameters?
#include <tree2/ode_solver.h>

namespace tree2 {

template <typename T>
class EBT {
public:
  typedef T plant_type;
  typedef Cohort<plant_type> cohort_type;
  typedef Patch<cohort_type> patch_type;
  EBT(Parameters p);

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
  Parameters r_parameters() const {return parameters;}
  patch_type r_patch()      const {return patch;}
  double              r_seed_rain(util::index species_index) const;
  std::vector<double> r_seed_rain_cohort(util::index species_index) const;
  std::vector<double> r_seed_rain_error(util::index species_index) const;
  std::vector<double> r_leaf_area_error(util::index species_index) const;
  std::vector<double> r_ode_times() const;

  CohortSchedule r_cohort_schedule() const {return cohort_schedule;}
  void r_set_cohort_schedule(CohortSchedule x);

private:
  double seed_rain_total() const;
  std::vector<double> seed_rain_cohort(size_t species_index) const;

  Parameters parameters;
  patch_type patch;
  CohortSchedule cohort_schedule;
  ode::Solver<patch_type> solver;
};

template <typename T>
EBT<T>::EBT(Parameters p)
  : parameters(p),
    patch(parameters),
    cohort_schedule(patch.size()),
    solver(patch, make_ode_control(p.control)) {
  parameters.validate();
  if (!util::identical(parameters.patch_area, 1.0)) {
    util::stop("Patch area must be exactly 1 for the EBT");
  }
}

template <typename T>
void EBT<T>::run() {
  reset();
  while (!complete()) {
    run_next();
  }
}

template <typename T>
std::vector<size_t> EBT<T>::run_next() {
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

  // TODO: Using r_ function -- fix at some point...
  const bool use_ode_times = cohort_schedule.r_use_ode_times();
  solver.set_state_from_system(patch);
  if (use_ode_times) {
    solver.advance_fixed(patch, e.times);
  } else {
    solver.advance(patch, e.time_end());
  }

  return ret;
}

template <typename T>
double EBT<T>::time() const {
  return solver.get_time();// TODO: ?patch.environment.time;
}

// NOTE: solver.reset() will set time within the solver to zero.
// However, there is no other current way of setting the time within
// the solver.  It might be better to add a set_time method within
// ode::Solver, and then here do explicitly ode_solver.set_time(0)?
template <typename T>
void EBT<T>::reset() {
  patch.reset();
  cohort_schedule.reset();
  solver.reset(patch);
}

template <typename T>
bool EBT<T>::complete() const {
  return cohort_schedule.remaining() == 0;
}

template <typename T>
double EBT<T>::seed_rain(size_t species_index) const {
  return util::trapezium(cohort_schedule.times(species_index),
                         seed_rain_cohort(species_index));
}

// TODO: Use ranged iterator?
template <typename T>
std::vector<double> EBT<T>::seed_rains() const {
  std::vector<double> ret;
  for (size_t i = 0; i < patch.size(); ++i) {
    ret.push_back(seed_rain(i));
  }
  return ret;
}

template <typename T>
std::vector<util::index> EBT<T>::r_run_next() {
  return util::index_vector(run_next());
}

template <typename T>
double EBT<T>::r_seed_rain(util::index species_index) const {
  return seed_rain(species_index.check_bounds(patch.size()));
}

template <typename T>
std::vector<double>
EBT<T>::r_seed_rain_cohort(util::index species_index) const {
  return seed_rain_cohort(species_index.check_bounds(patch.size()));
}

template <typename T>
std::vector<double> EBT<T>::r_seed_rain_error(util::index species_index) const {
  // TODO: This causes this to happen too often, given we usually get
  // all the errors I think?
  double tot_seed_out = seed_rain_total();
  const size_t idx = species_index.check_bounds(patch.size());
  return util::local_error_integration(cohort_schedule.times(idx),
                                       seed_rain_cohort(idx),
                                       tot_seed_out);
}

template <typename T>
std::vector<double> EBT<T>::r_leaf_area_error(util::index species_index) const {
  // const double tot_leaf_area  = patch.leaf_area_above(0.0);
  const size_t idx = species_index.check_bounds(patch.size());
  return patch.r_leaf_area_error(idx);
}

template <typename T>
std::vector<double> EBT<T>::r_ode_times() const {
  return solver.get_times();
}

template <typename T>
void EBT<T>::r_set_cohort_schedule(CohortSchedule x) {
  if (patch.ode_size() > 0) {
    util::stop("Cannot set schedule without resetting first");
  }
  util::check_length(x.get_n_species(), patch.size());
  cohort_schedule = x;
}

template <typename T>
double EBT<T>::seed_rain_total() const {
  double tot = 0.0;
  for (size_t i = 0; i < patch.size(); ++i) {
    tot += seed_rain(i);
  }
  return tot;
}

template <typename T>
std::vector<double> EBT<T>::seed_rain_cohort(size_t species_index) const {
  const std::vector<double> times = cohort_schedule.times(species_index);
  const Disturbance& disturbance_regime = patch.disturbance_regime();
  const double scal = parameters.Pi_0 * parameters.seed_rain[species_index];
  std::vector<double> seeds = patch.at(species_index).seeds();
  for (size_t i = 0; i < seeds.size(); ++i) {
    seeds[i] *= disturbance_regime.density(times[i]) * scal;
  }
  return seeds;
}


}

#endif
