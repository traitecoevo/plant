// -*-c++-*-
#ifndef PLANT_PLANT_PATCH_RUNNER_H_
#define PLANT_PLANT_PATCH_RUNNER_H_

#include <vector>
#include <odelia/ode_solver.hpp>
#include <plant/control.h>
#include <plant/solvers/schedule.h>
#include <plant/util.h>

namespace plant {

// What the two runners share:
//   * DeterministicPatchRunner<T,E>                   -- runs a deterministic Patch; and
//   * StochasticPatchRunner<T,E> -- runs a StochasticPatch.
//
// Both walk a Schedule of node introductions with an ODE solver that owns the
// live patch: take the entry the solver stands at, act on the patch, then
// integrate to the next entry. This CRTP base holds the parts of that both use
// (the parameters, the patch, the schedule and the solver) and what they do
// identically with them: the time, whether the schedule is spent, the schedule
// setters, putting the solver back at the start, and taking the next entry.
// What each does at an entry and how it integrates between entries stay in the
// derived class (DeterministicPatchRunner's events, pinned schedules and
// recordings; the stochastic
// runner's draws of establishment and death).
//
// The hooks each derived class declares, which the base calls on the derived
// type. There are no defaults, so a hook misspelt or left out is a compile
// error rather than a silent no-op:
//   * prepare_schedule(x)  -- adjusts a schedule set from R before it is taken.
//   * schedule_changed()   -- what follows setting the schedule from R.
template <typename Derived, typename P>
class PatchRunner {
public:
  typedef P                                patch_type;
  typedef typename P::parameters_type      parameters_type;
  typedef typename P::environment_type     environment_type;

  // The solver's time. The patch's own is wherever the solver last evaluated
  // it, which on a leg's final stage can round one ulp away from where the leg
  // landed, and an entry's time is compared to this exactly.
  double time() const { return solver.time(); }
  // True once every scheduled entry has been taken.
  bool complete() const { return schedule.remaining() == 0; }

  // * R interface
  parameters_type r_parameters() const { return parameters; }
  const patch_type& r_patch() const { return patch; }
  Rcpp::List r_get_state() const { return patch.r_get_state(); }
  Schedule r_schedule() const { return schedule; }

  void r_set_schedule(Schedule x) {
    refuse_if_started();
    util::check_length(x.get_n_species(), patch.size());
    d().prepare_schedule(x);
    schedule = x;
    // Update here so that extracting Parameters reflects the new schedule,
    // keeping Parameters self-sufficient.
    parameters.node_schedule_times = schedule.get_times();
    d().schedule_changed();
  }
  void r_set_node_schedule_times(std::vector<std::vector<double>> x) {
    refuse_if_started();
    schedule.set_times(x);
    parameters.node_schedule_times = x;
    d().schedule_changed();
  }

protected:
  // The schedule is made from the runner's own copy of the parameters, by
  // `schedule_from`, once that copy and the patch exist.
  template <typename ScheduleFrom>
  PatchRunner(parameters_type p, environment_type e, Control c,
                  ScheduleFrom schedule_from)
    : parameters(p), patch(parameters, e, c),
      schedule(schedule_from(parameters)),
      solver(patch, make_ode_control(c)) {}

  // The patch and the schedule back at their start, and the solver's system a
  // copy of that patch with its time and step state reset.
  //
  // NOTE: solver.reset() sets the solver's internal time to zero. There is
  // currently no other way to set that time; it might be cleaner to add an
  // odelia::ode::Solver::set_time and call set_time(0) explicitly here.
  void restart() {
    patch.reset();
    schedule.reset();
    solver.get_system_ref() = patch;
    solver.reset();
  }

  // The schedule's next entry, which must stand at `t0`, taken off the
  // schedule. `t_end` is set to the time the run integrates to from it: the
  // next entry's, or the end of the patch's lifetime. The entry stays valid,
  // since taking it only moves the schedule's cursor.
  const schedule_entry& take_entry(double t0, double& t_end) {
    const schedule_entry& entry = schedule.next();
    if (!util::identical(t0, entry.time)) {
      util::stop("Start time not what was expected");
    }
    t_end = schedule.time_end();
    schedule.pop();
    return entry;
  }

  parameters_type parameters;
  patch_type patch;
  Schedule schedule;
  odelia::ode::Solver<patch_type> solver;

private:
  // A schedule is set on a patch holding no nodes: one reset to empty, not one
  // started or initialised from an exported state.
  void refuse_if_started() const {
    if (patch.has_nodes()) {
      util::stop("Cannot set schedule without resetting first");
    }
  }

  Derived&       d()       { return static_cast<Derived&>(*this); }
  const Derived& d() const { return static_cast<const Derived&>(*this); }
};

} // namespace plant

#endif
