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
// One option is to make "StochasticNode<T,E>" that would include an
// ID.  Another option is to track required bits of data that within
// the patch somehow?
template <typename T, typename E> class StochasticPatchRunner {
public:
  typedef T                       strategy_type;
  typedef E                       environment_type;
  typedef Individual<T, E>        individual_type;
  typedef StochasticSpecies<T, E> species_type;
  typedef StochasticPatch<T, E>   patch_type;
  typedef Parameters<T, E>        parameters_type;

  StochasticPatchRunner(parameters_type p, environment_type e, Control c);

  void run();
  size_t run_next();
  void advance(double time_);

  double time() const { return patch.time(); }
  void reset();
  bool complete() const;

  // * R interface
  util::index r_run_next();
  parameters_type r_parameters() const { return parameters; }
  const patch_type &r_patch() const { return patch; }

  // R-facing names mirror SCM (node_schedule / set_node_schedule[_times]).
  NodeSchedule r_node_schedule() const { return node_schedule; }
  void r_set_node_schedule(NodeSchedule x);
  void r_set_node_schedule_times(std::vector<std::vector<double>> x);
  Rcpp::List r_get_state() const { return patch.r_get_state(); };

private:
  bool deaths();

  parameters_type parameters;
  patch_type patch;
  NodeSchedule node_schedule;
  odelia::ode::Solver<patch_type> solver;
};

template <typename T, typename E>
StochasticPatchRunner<T, E>::StochasticPatchRunner(parameters_type p,
                                                   environment_type e,
                                                   Control c)
    : parameters(p), patch(parameters, e, c),
      node_schedule(make_empty_stochastic_schedule(parameters)),
      solver(patch, make_ode_control(c)) {
  parameters.validate();
  solver.set_collect(false);
}

template <typename T, typename E> void StochasticPatchRunner<T, E>::run() {
  reset();
  while (!complete()) {
    run_next();
  }
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
  // The stochastic tower only knows about node introductions; a schedule
  // carrying any other event type belongs to the deterministic solver until
  // this runner is migrated onto the shared queue too (#601).
  if (!e.is_node_introduction()) {
    util::stop("The stochastic solver does not yet apply scheduled events "
               "other than node introductions");
  }
  const size_t idx = e.target_index;
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
  // deaths() acts on the solver's live system in place (see below). If any
  // individual died the ODE system shrank, so re-pull its now-smaller state
  // into the solver. The `patch` member snapshot is refreshed once, by the
  // caller (run_next), rather than copied in and back out around the deaths.
  if (deaths()) {
    solver.set_state_from_system();
  }
}

// Apply stochastic deaths to the solver's owned system in place. Mirrors the way
// SCM operates on solver.get_system_ref() directly instead of round-tripping
// through the `patch` member -- avoids two full Patch copies per step.
template <typename T, typename E> bool StochasticPatchRunner<T, E>::deaths() {
  const auto ret = solver.get_system_ref().deaths();
  return std::any_of(ret.begin(), ret.end(), [](size_t i) { return i > 0; });
}

template <typename T, typename E>
bool StochasticPatchRunner<T, E>::complete() const {
  return node_schedule.remaining() == 0;
}

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
      // One step of this length would be tens of years for a late first arrival,
      // and the environment's own states are integrated over it, so this leg
      // needs error control like every other advance in the runner.
      //
      // Only when there is something to integrate, though. On an empty patch
      // whose environment carries no state the system is zero-width, and
      // advance_adaptive still walks the step-size controller: with no elements
      // to measure, every step is accepted and step_size_last ratchets up by 5x
      // until it saturates at ode_step_size_max. That survives into the first
      // real step, because only SolverInternal::step() writes step_size_last and
      // step_to() -- which advance_fixed drives -- does not, and nothing between
      // here and the first advance() restores it (set_state_from_system() resizes
      // and re-reads state and rates, but not the step size). The first step
      // after the first arrival would then begin from a rejected step of
      // ode_step_size_max rather than ode_step_size_initial. Error control
      // recovers, but the realised step sequence differs: FF16's collected
      // trajectory moves by up to 1e-5 relative, inside ode_tol_rel but not
      // identical. Taking the fixed step when the system is empty is what keeps
      // FF16 and K93 bit-identical.
      if (patch.ode_size() > 0) {
        solver.advance_adaptive({solver.time(), t});
      } else {
        solver.advance_fixed({solver.time(), t});
      }
      patch = solver.get_system();
    }
  }
}

template <typename T, typename E>
util::index StochasticPatchRunner<T, E>::r_run_next() {
  return util::index(run_next());
}

template <typename T, typename E>
void StochasticPatchRunner<T, E>::r_set_node_schedule(NodeSchedule x) {
  if (patch.node_ode_size() > 0) {
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
void StochasticPatchRunner<T, E>::r_set_node_schedule_times(
    std::vector<std::vector<double>> x) {
  if (patch.node_ode_size() > 0) {
    util::stop("Cannot set schedule without resetting first");
  }
  node_schedule.set_times(x);
  parameters.node_schedule_times = x;
  reset();
}

} // namespace plant

#endif
