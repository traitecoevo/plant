// -*-c++-*-
#ifndef PLANT_PLANT_STOCHASTIC_PATCH_RUNNER_H_
#define PLANT_PLANT_STOCHASTIC_PATCH_RUNNER_H_

#include <plant/runner.h>
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
// One option is to make "StochasticNode<T,E>" that would include an
// ID.  Another option is to track required bits of data that within
// the patch somehow?
template <typename T, typename E>
class StochasticPatchRunner
    : public ScheduleDrivenRunner<StochasticPatchRunner<T, E>> {
public:
  typedef T                       strategy_type;
  typedef E                       environment_type;
  typedef Individual<T, E>        individual_type;
  typedef StochasticSpecies<T, E> species_type;
  typedef StochasticPatch<T, E>   patch_type;
  typedef Parameters<T, E>        parameters_type;

  StochasticPatchRunner(parameters_type p, environment_type e, Control c);

  // run() (reset, then step to completion) and complete() are inherited from
  // ScheduleDrivenRunner -- the stochastic runner needs no per-step bookkeeping,
  // so the base defaults apply unchanged. node_schedule also lives in the base;
  // the using-declaration lets the unqualified references below resolve.
  using ScheduleDrivenRunner<StochasticPatchRunner<T, E>>::node_schedule;

  size_t run_next();
  void advance(double time_);

  double time() const { return patch.time(); }
  void reset();

  // * R interface
  util::index r_run_next();
  parameters_type r_parameters() const { return parameters; }
  const patch_type &r_patch() const { return patch; }

  // TODO(#479): consider renaming NodeSchedule -> Schedule
  NodeSchedule r_schedule() const { return node_schedule; }
  void r_set_schedule(NodeSchedule x);
  void r_set_schedule_times(std::vector<std::vector<double>> x);
  Rcpp::List r_get_state() const { return patch.r_get_state(); };

private:
  bool deaths();

  parameters_type parameters;
  patch_type patch;
  odelia::ode::Solver<patch_type> solver;
};

template <typename T, typename E>
StochasticPatchRunner<T, E>::StochasticPatchRunner(parameters_type p,
                                                   environment_type e,
                                                   Control c)
    : ScheduleDrivenRunner<StochasticPatchRunner<T, E>>(
          make_empty_stochastic_schedule(p)),
      parameters(p), patch(parameters, e, c),
      solver(patch, make_ode_control(c)) {
  parameters.validate();
  solver.set_collect(false);
}

template <typename T, typename E>
size_t StochasticPatchRunner<T, E>::run_next() {
  const double t0 = time();
  auto& patch_solver = solver.get_system_ref();

  // NOTE: Unlike SCM::run_next(), this assumes that there is only a
  // single event at a given time.  That's not all bad -- multiple
  // events could occur at a single time but the time-saving trick of
  // not computing the light environment would not work.
  NodeSchedule::Event e = node_schedule.next_event();
  if (!util::identical(t0, e.time_introduction())) {
    util::stop("Start time not what was expected");
  }
  const size_t idx = e.species_index;
  node_schedule.pop();

  if (patch_solver.introduce_new_node(idx)) {
    solver.set_state_from_system();
  }
  advance(e.time_end());
  patch = solver.get_system_ref();

  return idx;
}

template <typename T, typename E>
void StochasticPatchRunner<T, E>::advance(double time_) {
  solver.advance_adaptive({solver.time(), time_});
  patch = solver.get_system_ref();
  if (deaths()) {
    solver.get_system_ref() = patch;
    solver.set_state_from_system();
  }
}

template <typename T, typename E> bool StochasticPatchRunner<T, E>::deaths() {
  const auto ret = patch.deaths();
  return std::any_of(ret.begin(), ret.end(), [](size_t i) { return i > 0; });
}

// complete() is inherited from ScheduleDrivenRunner (node_schedule.remaining()
// == 0); the stochastic runner adds no completion condition of its own.

// NOTE: solver.reset() will set time within the solver to zero.
// However, there is no other current way of setting the time within
// the solver.  It might be better to add a set_time method within
// odelia::ode::Solver, and then here do explicitly ode_solver.set_time(0)?
template <typename T, typename E> void StochasticPatchRunner<T, E>::reset() {
  patch.reset();
  node_schedule.reset();
  solver.get_system_ref() = patch;
  solver.reset();
  if (node_schedule.size() > 0) {
    const double t = node_schedule.next_event().time_introduction();
    if (t >= 0.0) {
      solver.advance_fixed({solver.time(), t});
      patch = solver.get_system();
    }
  }
}

template <typename T, typename E>
util::index StochasticPatchRunner<T, E>::r_run_next() {
  return util::index(run_next());
}

template <typename T, typename E>
void StochasticPatchRunner<T, E>::r_set_schedule(NodeSchedule x) {
  if (patch.ode_size() > 0) {
    util::stop("Cannot set schedule without resetting first");
  }
  util::check_length(x.get_n_species(), patch.size());
  node_schedule = x;

  // Update these here so that extracting Parameters would give the
  // new schedule, this making Parameters sufficient.
  parameters.node_schedule_times = node_schedule.get_times();
  reset();
}

template <typename T, typename E>
void StochasticPatchRunner<T, E>::r_set_schedule_times(
    std::vector<std::vector<double>> x) {
  if (patch.ode_size() > 0) {
    util::stop("Cannot set schedule without resetting first");
  }
  node_schedule.set_times(x);
  parameters.node_schedule_times = x;
  reset();
}

} // namespace plant

#endif
