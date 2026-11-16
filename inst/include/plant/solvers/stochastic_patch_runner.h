// -*-c++-*-
#ifndef PLANT_PLANT_STOCHASTIC_PATCH_RUNNER_H_
#define PLANT_PLANT_STOCHASTIC_PATCH_RUNNER_H_

#include <algorithm>
#include <cmath>
#include <odelia/ode_solver.hpp>
#include <plant/solvers/patch_runner.h>
#include <plant/community/stochastic_patch.h>
#include <plant/solvers/stochastic_utils.h>

namespace plant {

// The name here is likely to change.
//
// The interface mirrors that of DeterministicPatchRunner, and the two share
// their state, their
// schedule handling and their solver reset through PatchRunner.
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
  : public PatchRunner<StochasticPatchRunner<T, E>, StochasticPatch<T, E>> {
  typedef PatchRunner<StochasticPatchRunner<T, E>, StochasticPatch<T, E>>
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

  // time(), complete(), r_parameters(), r_patch(), r_get_state(), r_schedule()
  // and the schedule setters are PatchRunner's, as DeterministicPatchRunner's are.
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

  // PatchRunner's hooks: a schedule set from R gets this runner's death checks
  // in place of any it carries, and takes effect at once. Its last check is at
  // its max_time, so one without a finite max_time is refused before it is
  // taken.
  void prepare_schedule(Schedule& x) const {
    if (x.size() > 0 && !std::isfinite(x.get_max_time())) {
      util::stop("The stochastic runner needs a schedule with a finite "
                 "max_time, where its last death check falls");
    }
  }
  void schedule_changed() {
    place_death_checks();
    reset();
  }

  // The death-check grid, on the schedule: a check at every multiple of
  // death_interval strictly after the first introduction and before the end
  // of the patch's lifetime, and one at the end. Nothing is alive before the
  // first introduction, and a check there would still split the leg that
  // integrates the environment up to it.
  void place_death_checks();
  // Takes the next entry, applies it and integrates to the one after.
  const schedule_entry& step();
  void land(double t);
  void draw_deaths();

  double death_interval;
};

template <typename T, typename E>
StochasticPatchRunner<T, E>::StochasticPatchRunner(parameters_type p,
                                                   environment_type e,
                                                   Control c)
    : base_type(p, e, c,
                [](const parameters_type& q) {
                  return make_empty_stochastic_schedule(q);
                }),
      death_interval(c.stochastic_death_interval) {
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

// One arrival to the next: the entry at the solver's time, which introduces,
// and then every entry up to the next introduction (or to the end), which
// draws deaths. Stepping arrival to arrival is what R collects a row for.
template <typename T, typename E>
size_t StochasticPatchRunner<T, E>::run_next() {
  // Its death checks all follow the first introduction, so a walk that has
  // begun stands at one.
  if (schedule.next().species.empty()) {
    util::stop("The stochastic runner stands at no introduction");
  }
  const schedule_entry& first = step();
  while (!complete() && schedule.next().species.empty()) {
    step();
  }
  patch = solver.get_system_ref();
  return first.species.front();
}

// Acts on the solver's own system, as DeterministicPatchRunner does, rather
// than round-tripping through the `patch` member. The schedule holds the
// entry's death checks ahead of its introductions (event_type_rank), so the
// individuals a check removes died before anyone arriving here.
template <typename T, typename E>
const schedule_entry& StochasticPatchRunner<T, E>::step() {
  const double t0 = time();
  auto& patch_solver = solver.get_system_ref();

  // The stochastic tower draws deaths and introduces nodes; an instant
  // carrying any other event type belongs to the deterministic solver until
  // this runner applies those too (#601). Refused before the entry is taken,
  // so the refusal leaves the schedule where it was.
  const schedule_entry& next = schedule.next();
  if (util::identical(t0, next.time)) {
    for (const auto& a : next.actions) {
      if (a.type != EventType::DeathCheck) {
        util::stop("The stochastic solver does not yet apply scheduled events "
                   "other than node introductions and its death checks");
      }
    }
  }
  double t_end;
  const schedule_entry& entry = this->take_entry(t0, t_end);

  for (const auto& a : entry.actions) {
    if (a.type == EventType::DeathCheck) {
      draw_deaths();
    }
  }
  // Arrival times are drawn from a continuous distribution, so an introduction
  // here names one species; the loop is what makes that an observation rather
  // than an assumption. Only the first reaches R, which is why this stays a
  // walk over the set rather than the set itself: the return is one arrival.
  bool widened = false;
  for (const size_t idx : entry.species) {
    if (patch_solver.introduce_new_node(idx)) {
      widened = true;
    }
  }
  if (widened) {
    solver.set_state_from_system();
  }
  const double t = solver.time();
  if (t < t_end) {
    solver.advance_adaptive({t, t_end});
    land(t_end);
  }
  return entry;
}

template <typename T, typename E>
void StochasticPatchRunner<T, E>::place_death_checks() {
  std::vector<ScheduleEvent> checks;
  const std::vector<ScheduleEvent> events = schedule.get_events();
  const auto first = std::find_if(
    events.begin(), events.end(),
    [](const ScheduleEvent& e) { return e.is_node_introduction(); });
  if (first != events.end()) {
    const double t_first = first->time_introduction();
    const double t_max = schedule.get_max_time();
    if (!std::isfinite(t_max)) {
      util::stop("The stochastic runner needs a schedule with a finite "
                 "max_time, where its last death check falls");
    }
    // Multiples computed as k * interval, each from k, so a long grid does not
    // accumulate rounding and every point is the same double however it is
    // reached.
    size_t k = static_cast<size_t>(std::floor(t_first / death_interval));
    for (;; ++k) {
      const double t = static_cast<double>(k) * death_interval;
      if (!(t > t_first)) {
        continue;
      }
      if (!(t < t_max)) {
        break;
      }
      checks.emplace_back(t, 0, EventType::DeathCheck, EventTarget::Patch);
    }
    checks.emplace_back(t_max, 0, EventType::DeathCheck, EventTarget::Patch);
  }
  schedule.replace_actions(EventType::DeathCheck, checks);
}

template <typename T, typename E>
void StochasticPatchRunner<T, E>::land(double t) {
  solver.get_system_ref().set_time(t);
}

// Every draw re-syncs the solver, not only one in which somebody died: the
// draw also resets each survivor's accumulated hazard in the system, and
// without the re-sync the solver would carry on from its own copy, so the next
// draw would count this interval's hazard again.
template <typename T, typename E>
void StochasticPatchRunner<T, E>::draw_deaths() {
  solver.get_system_ref().deaths();
  solver.set_state_from_system();
}

template <typename T, typename E> void StochasticPatchRunner<T, E>::reset() {
  this->restart();
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
      // here and the first leg restores it (set_state_from_system() resizes
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
