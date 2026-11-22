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
// Its arrivals are drawn by draw_arrivals(), from each species' birth rate and
// the patch's area; a schedule set from R gives them instead, one per time.
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
  void draw_arrivals();

  // time(), complete(), r_parameters(), r_patch(), r_get_state(), r_schedule()
  // and the schedule setters are PatchRunner's, as DeterministicPatchRunner's are.
  using base_type::time;
  using base_type::complete;
  void reset();

  // * R interface
  util::index r_run_next();
  // What the run's scheduled events did; its death checks are not logged.
  EventLog r_event_log() const { return event_log_from_records(event_log); }

private:
  using base_type::parameters;
  using base_type::patch;
  using base_type::schedule;
  using base_type::solver;

  // PatchRunner's hooks: a schedule set from R is checked and gets this
  // runner's death checks in place of any it carries, on the copy and before
  // it is taken, so a schedule refused here leaves the runner as it was. It
  // then takes effect at once.
  void prepare_schedule(Schedule& x) const {
    check_events(x);
    place_death_checks(x);
  }
  void schedule_changed() { reset(); }

  // The death-check grid, on `x`: a check at every multiple of death_interval
  // strictly after the first introduction and before the schedule's max_time,
  // and one at max_time. Nothing is alive before the first introduction, and a
  // check there would still split the leg that integrates the environment up
  // to it.
  void place_death_checks(Schedule& x) const;
  // What would otherwise be refused only when its entry is reached, after the
  // entry is taken: a schedule with no end, or a pulse into a pool this
  // environment does not have.
  void check_events(const Schedule& x) const;
  // Takes the next entry, applies it and integrates to the one after.
  const schedule_entry& step();
  // An entry's actions, in their rank order: death checks, then the rest.
  void apply_actions(const schedule_entry& entry);
  // A leg before the first arrival; see reset().
  void lead_in(double t);
  void land(double t);
  void draw_deaths();

  double death_interval;
  // Lives on the runner rather than the patch, as DeterministicPatchRunner's
  // does: the patch is copied out after every run_next().
  std::vector<EventRecord> event_log;
  // Ten million checks is an interval of 1e-5 yr (about five minutes) over a
  // 100-year patch, and a few gigabytes of schedule.
  static constexpr double max_death_checks = 1e7;
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
  if (!(parameters.max_patch_lifetime / death_interval <= max_death_checks)) {
    util::stop("stochastic_death_interval is too small for this patch's "
               "lifetime: it would place more than " +
               util::to_string(max_death_checks) + " death checks");
  }
  solver.set_collect(false);
  // Its patch declares no tolerance scales, so the solver never reads them.
  if (c.ode_tol_scale_soil != 1.0) {
    util::warning("The stochastic runner ignores ode_tol_scale_soil");
  }
}

// Each species' arrivals over the schedule's span, drawn from its birth rate
// and the patch's area, in place of its introductions; the schedule's other
// events stay. Set as a schedule from R is, so its parameters record the times.
template <typename T, typename E>
void StochasticPatchRunner<T, E>::draw_arrivals() {
  Schedule drawn = schedule;
  if (!std::isfinite(drawn.get_max_time())) {
    util::stop("Arrivals are drawn up to the schedule's max_time, which must "
               "be finite");
  }
  for (size_t i = 0; i < patch.size(); ++i) {
    drawn.set_times(patch.at_species(i).draw_arrivals(
                      0.0, drawn.get_max_time(), patch.get_area()),
                    i);
  }
  this->r_set_schedule(drawn);
}

template <typename T, typename E> void StochasticPatchRunner<T, E>::run() {
  reset();
  while (!complete()) {
    run_next();
  }
}

// One arrival to the next: the entry at the solver's time, which introduces,
// and then every entry up to the next introduction (or to the end), which
// draws deaths and applies the run's other events. Stepping arrival to arrival
// is what R collects a row for.
template <typename T, typename E>
size_t StochasticPatchRunner<T, E>::run_next() {
  // reset() walks every entry before the first arrival, and each call walks
  // to the next, so a run that has begun stands at one.
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

  double t_end;
  const schedule_entry& entry = this->take_entry(t0, t_end);

  apply_actions(entry);
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

// Every action re-syncs the solver, since each moves the system: a draw
// removes individuals and resets the survivors' hazards, a removal shrinks the
// system, and a pulse moves the environment's pools.
template <typename T, typename E>
void StochasticPatchRunner<T, E>::apply_actions(const schedule_entry& entry) {
  bool applied = false;
  for (const auto& a : entry.actions) {
    if (a.type == EventType::DeathCheck) {
      draw_deaths();
    } else {
      event_log.push_back(solver.get_system_ref().apply_event(a));
      applied = true;
    }
  }
  if (applied) {
    solver.set_state_from_system();
  }
}

template <typename T, typename E>
void StochasticPatchRunner<T, E>::check_events(const Schedule& x) const {
  if (x.size() > 0 && !std::isfinite(x.get_max_time())) {
    util::stop("The stochastic runner needs a schedule with a finite "
               "max_time, where its last death check falls");
  }
  for (const ScheduleEvent& e : x.get_events()) {
    if (e.type == EventType::ResourcePulse &&
        e.target_index >= patch.n_resources()) {
      util::stop("The resource pulse at t = " +
                 util::to_string(e.time_introduction()) + " is into resource " +
                 util::to_string(e.target_index + 1) +
                 ", which does not exist: this environment has " +
                 util::to_string(patch.n_resources()) + " resources");
    }
  }
}

template <typename T, typename E>
void StochasticPatchRunner<T, E>::place_death_checks(Schedule& x) const {
  std::vector<ScheduleEvent> checks;
  const double t_max = x.get_max_time();
  std::vector<double> arrivals;
  for (const ScheduleEvent& e : x.get_events()) {
    if (e.is_node_introduction()) {
      arrivals.push_back(e.time_introduction());
    }
  }
  if (!arrivals.empty()) {
    // check_events() has refused an endless schedule. The runner's own setter takes times unchecked, and a time outside the
    // patch's life would put a check before the first arrival or none after
    // the last.
    for (const double t : arrivals) {
      if (!(t >= 0.0 && t <= t_max)) {
        util::stop("Arrival time " + util::to_string(t) + " is outside [0, " +
                   util::to_string(t_max) + "], the patch's lifetime");
      }
    }
    // Every check is held on the schedule, so an interval far below the run's
    // length is refused rather than allocated. Bounding the grid from zero
    // rather than from the first arrival also bounds the index k below, which
    // starts near t_first / death_interval.
    const double n_checks = t_max / death_interval;
    if (!(n_checks <= max_death_checks)) {
      util::stop("stochastic_death_interval " + util::to_string(death_interval) +
                 " would place " + util::to_string(n_checks) +
                 " death checks; at most " + util::to_string(max_death_checks) +
                 " are allowed");
    }
    const double t_first = arrivals.front();
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
  x.replace_actions(EventType::DeathCheck, checks);
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

// The leg up to the first arrival, or to an event before it.
//
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
template <typename T, typename E>
void StochasticPatchRunner<T, E>::lead_in(double t) {
  if (solver.get_system_ref().ode_size() > 0) {
    solver.advance_adaptive({solver.time(), t});
  } else {
    solver.advance_fixed({solver.time(), t});
  }
  land(t);
}

template <typename T, typename E> void StochasticPatchRunner<T, E>::reset() {
  this->restart();
  event_log.clear();
  if (schedule.size() > 0 && schedule.next().time >= 0.0) {
    lead_in(schedule.next().time);
    // Events before the first arrival, such as rain on bare ground, apply on
    // the way to it. No death check comes before it.
    while (schedule.next().species.empty()) {
      double t_end;
      const schedule_entry& entry = this->take_entry(time(), t_end);
      apply_actions(entry);
      lead_in(t_end);
      if (complete()) {
        break;
      }
    }
    patch = solver.get_system();
  }
}

template <typename T, typename E>
util::index StochasticPatchRunner<T, E>::r_run_next() {
  return util::index(run_next());
}

} // namespace plant

#endif
