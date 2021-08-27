// -*-c++-*-
#ifndef PLANT_PLANT_STOCHASTIC_PATCH_RUNNER_H_
#define PLANT_PLANT_STOCHASTIC_PATCH_RUNNER_H_

#include <plant/stochastic_patch.h>
#include <plant/stochastic_utils.h>

namespace plant {

// The name here is likely to change.
//
// The interface mirrors that of SCM; we'll use the same basic
// approach for a schedule too perhaps.
//
// One thing that is different though is that we need to take care in
// tracking who is who; deaths are going to make it hard to plot sizes
// vs time without some care.
//
// One option is to make "StochasticCohort<T,E>" that would include an
// ID.  Another option is to track required bits of data that within
// the patch somehow?
template <typename T, typename E>
class StochasticPatchRunner {
public:
  typedef T                      strategy_type;
  typedef E                 environment_type;
  typedef Individual<T,E>             individual_type;
  typedef StochasticSpecies<T,E> species_type;
  typedef StochasticPatch<T,E>   patch_type;
  typedef SpeciesParameters<T>        species_params_type;

  StochasticPatchRunner(species_params_type s, environment_type e, Control c);

  void run();
  size_t run_next();
  void advance(double time_);

  double time() const {return patch.time();}
  void reset();
  bool complete() const;

  // * R interface
  util::index r_run_next();
  species_params_type r_species_parameters() const {return species_parameters;}
  const patch_type& r_patch() const {return patch;}

  CohortSchedule r_schedule() const {return cohort_schedule;}
  void r_set_schedule(CohortSchedule x);
  void r_set_schedule_times(std::vector<std::vector<double> > x);

private:
  bool deaths();

  species_params_type species_parameters;
  patch_type patch;
  CohortSchedule cohort_schedule;
  ode::Solver<patch_type> solver;
};

template <typename T, typename E>
StochasticPatchRunner<T,E>::StochasticPatchRunner(species_params_type s, environment_type e, Control c)
  : species_parameters(s),
    patch(s, e, c),
    cohort_schedule(make_empty_stochastic_schedule(species_parameters, c)),
    solver(patch, make_ode_control(c)) {

  species_parameters.validate();
}

template <typename T, typename E>
void StochasticPatchRunner<T,E>::run() {
  reset();
  while (!complete()) {
    run_next();
  }
}

template <typename T, typename E>
size_t StochasticPatchRunner<T,E>::run_next() {
  const double t0 = time();

  // NOTE: Unlike SCM::run_next(), this assumes that there is only a
  // single event at a given time.  That's not all bad -- multiple
  // events could occur at a single time but the time-saving trick of
  // not computing the light environment would not work.
  CohortSchedule::Event e = cohort_schedule.next_event();
  if (!util::identical(t0, e.time_introduction())) {
    util::stop("Start time not what was expected");
  }
  const size_t idx = e.species_index;
  cohort_schedule.pop();

  if (patch.introduce_new_cohort(idx)) {
    solver.set_state_from_system(patch);
  }
  advance(e.time_end());

  return idx;
}

template <typename T, typename E>
void StochasticPatchRunner<T,E>::advance(double time_) {
  // Clones some of Solver<T,E>::advance()
  solver.set_time_max(time_);
  while (solver.get_time() < time_) {
    solver.step(patch);
    if (deaths()) {
      solver.set_state_from_system(patch);
    }
  }
}

template <typename T, typename E>
bool StochasticPatchRunner<T,E>::deaths() {
  const auto ret = patch.deaths();
  return std::any_of(ret.begin(), ret.end(), [](size_t i) {return i > 0;});
}

// NOTE: solver.reset() will set time within the solver to zero.
// However, there is no other current way of setting the time within
// the solver.  It might be better to add a set_time method within
// ode::Solver, and then here do explicitly ode_solver.set_time(0)?
template <typename T, typename E>
void StochasticPatchRunner<T,E>::reset() {
  patch.reset();
  cohort_schedule.reset();
  solver.reset(patch);
  if (cohort_schedule.size() > 0) {
    const double t = cohort_schedule.next_event().time_introduction();
    if (t >= 0.0) {
      solver.step_to(patch, t);
    }
  }
}

template <typename T, typename E>
bool StochasticPatchRunner<T,E>::complete() const {
  return cohort_schedule.remaining() == 0;
}

template <typename T, typename E>
util::index StochasticPatchRunner<T,E>::r_run_next() {
  return util::index(run_next());
}

template <typename T, typename E>
void StochasticPatchRunner<T,E>::r_set_schedule(CohortSchedule x) {
  if (patch.ode_size() > 0) {
    util::stop("Cannot set schedule without resetting first");
  }
  util::check_length(x.get_n_species(), patch.size());
  cohort_schedule = x;

  // Update these here so that extracting SpeciesParameters would give the
  // new schedule, this making SpeciesParameters sufficient.
  species_parameters.cohort_schedule_times = cohort_schedule.get_times();
  reset();
}

template <typename T, typename E>
void StochasticPatchRunner<T,E>::r_set_schedule_times(std::vector<std::vector<double> > x) {
  if (patch.ode_size() > 0) {
    util::stop("Cannot set schedule without resetting first");
  }
  cohort_schedule.set_times(x);
  species_parameters.cohort_schedule_times = x;
  reset();
}




}

#endif
