// -*-c++-*-
#ifndef PLANT_PLANT_STOCHASTIC_PATCH_RUNNER_H_
#define PLANT_PLANT_STOCHASTIC_PATCH_RUNNER_H_

#include <cmath>
#include <odelia/ode_solver.hpp>
#include <plant/patch_runner_base.h>
#include <plant/stochastic_patch.h>
#include <plant/stochastic_utils.h>

namespace plant {

// The name here is likely to change.
//
// The interface mirrors that of SCM, and the two share their state, their
// schedule handling and their solver reset through PatchRunnerBase.
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
  : public PatchRunnerBase<StochasticPatchRunner<T, E>, StochasticPatch<T, E>> {
  typedef PatchRunnerBase<StochasticPatchRunner<T, E>, StochasticPatch<T, E>>
    base_type;
  friend base_type;
public:
  typedef T                       strategy_type;
  typedef E                       environment_type;
  typedef Individual<T, E>        individual_type;
  typedef StochasticSpecies<T, E> species_type;
  using typename base_type::patch_type;
  using typename base_type::parameters_type;

  StochasticPatchRunner(parameters_type p, environment_type e, Control c);

  void run();
  size_t run_next();
  void advance(double time_);

  // time(), complete(), r_parameters(), r_patch(), r_get_state(), r_schedule()
  // and the schedule setters are PatchRunnerBase's, as SCM's are.
  using base_type::time;
  using base_type::complete;
  void reset();

  // * R interface
  util::index r_run_next();

private:
  using base_type::parameters;
  using base_type::patch;
  using base_type::schedule;
  using base_type::solver;

  // PatchRunnerBase's hooks: a schedule set from R is taken as it is, and
  // takes effect at once.
  void prepare_schedule(Schedule&) const {}
  void schedule_changed() { reset(); }

  void land(double t);
  void draw_deaths();

  double death_interval;
  // Draws taken so far on the death-check grid; the next is at
  // (death_checks + 1) * death_interval.
  size_t death_checks;
};

template <typename T, typename E>
StochasticPatchRunner<T, E>::StochasticPatchRunner(parameters_type p,
                                                   environment_type e,
                                                   Control c)
    : base_type(p, e, c,
                [](const parameters_type& q) {
                  return make_empty_stochastic_schedule(q);
                }),
      death_interval(c.stochastic_death_interval), death_checks(0) {
  parameters.validate();
  if (!(death_interval > 0.0 && std::isfinite(death_interval))) {
    util::stop("stochastic_death_interval must be positive and finite");
  }
  solver.set_collect(false);
  // Its patch declares no tolerance scales, so the solver never reads them.
  if (c.ode_tol_scale_soil != 1.0) {
    util::warning("The stochastic runner ignores ode_tol_scale_soil");
  }
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

  double t_end;
  const schedule_entry& intro = this->take_entry(t0, t_end);
  // The stochastic tower only knows about node introductions; an instant
  // carrying any other event type belongs to the deterministic solver until
  // this runner is migrated onto the shared queue too (#601).
  if (!intro.actions.empty()) {
    util::stop("The stochastic solver does not yet apply scheduled events "
               "other than node introductions");
  }

  // Arrival times are drawn from a continuous distribution, so an introduction
  // here names one species; the loop is what makes that an observation rather
  // than an assumption. Only the first reaches R, which is why this stays a
  // walk over the set rather than the set itself: the return is one arrival.
  bool widened = false;
  for (const size_t idx : intro.species) {
    if (patch_solver.introduce_new_node(idx)) {
      widened = true;
    }
  }
  if (widened) {
    solver.set_state_from_system();
  }
  advance(t_end);
  patch = solver.get_system_ref();

  return intro.species.front();
}

// Integrates to `time_`, drawing deaths at every point of the death-check grid
// on the way and at the end of the patch's lifetime. When deaths are drawn
// does not depend on when seeds arrive. The `patch` member snapshot is
// refreshed once, by the caller (run_next).
template <typename T, typename E>
void StochasticPatchRunner<T, E>::advance(double time_) {
  double t = solver.time();
  // Skip the grid points the run has already passed (reset() integrates to
  // the first arrival without drawing, as nothing is alive before it), and any
  // that rounds onto t.
  const double passed = std::floor(t / death_interval);
  if (passed > static_cast<double>(death_checks)) {
    death_checks = static_cast<size_t>(passed);
  }
  bool drawn_at_end = false;
  while (t < time_) {
    const double t_check = (death_checks + 1) * death_interval;
    if (!(t_check > t)) {
      ++death_checks;
      continue;
    }
    const bool check = t_check <= time_;
    const double t_next = check ? t_check : time_;
    solver.advance_adaptive({t, t_next});
    land(t_next);
    if (check) {
      ++death_checks;
      draw_deaths();
    }
    drawn_at_end = check && t_next == time_;
    t = t_next;
  }
  if (complete() && !drawn_at_end) {
    draw_deaths();
  }
}

template <typename T, typename E>
void StochasticPatchRunner<T, E>::land(double t) {
  solver.get_system_ref().set_time(t);
}

// Acts on the solver's own system, as SCM does, rather than round-tripping
// through the `patch` member. Every draw re-syncs the solver, not only one in
// which somebody died: the draw also resets each survivor's accumulated hazard
// in the system, and without the re-sync the solver would carry on from its
// own copy, so the next draw would count this interval's hazard again.
template <typename T, typename E>
void StochasticPatchRunner<T, E>::draw_deaths() {
  solver.get_system_ref().deaths();
  solver.set_state_from_system();
}

template <typename T, typename E> void StochasticPatchRunner<T, E>::reset() {
  this->restart();
  death_checks = 0;
  if (schedule.size() > 0) {
    const double t = schedule.next().time;
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
      // recovers, but the realised step sequence differs, so the collected
      // trajectory moves within the ODE tolerance rather than staying
      // identical. Taking the fixed step when the system is empty is what keeps
      // FF16 and K93 bit-identical.
      if (patch.ode_size() > 0) {
        solver.advance_adaptive({solver.time(), t});
      } else {
        solver.advance_fixed({solver.time(), t});
      }
      land(t);
      patch = solver.get_system();
    }
  }
}

template <typename T, typename E>
util::index StochasticPatchRunner<T, E>::r_run_next() {
  return util::index(run_next());
}

} // namespace plant

#endif
