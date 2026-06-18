// -*-c++-*-
#ifndef PLANT_PLANT_SCM_H_
#define PLANT_PLANT_SCM_H_

#include <plant/node_schedule.h>
#include <plant/ode_solver/ode_solver.h>
#include <plant/patch.h>
#include <plant/scm_utils.h>

#include <algorithm>
#include <limits>

using namespace Rcpp;

namespace plant {

// SCM: the "Solver for Characteristics Method" driver.
//
// Owns a Patch (the population being integrated), a NodeSchedule (when each
// species' nodes are introduced), and an ODE Solver. It steps the patch
// forward by repeatedly: introducing all nodes due at the current time, then
// integrating the patch ODE system up to the next introduction time.
//
// The patch owns all of the ecology (fitness, offspring, competition, error
// computations); the SCM is the time-stepping/scheduling layer on top of it.
// Most r_* members are thin facades that expose the C++ API to R via RcppR6.
template <typename T, typename E> class SCM {
public:
  // ---- Type aliases ------------------------------------------------------
  typedef T                strategy_type;
  typedef E                environment_type;
  typedef Individual<T, E> individual_type;
  typedef Node<T, E>       node_type;
  typedef Species<T, E>    species_type;
  typedef Patch<T, E>      patch_type;
  typedef Parameters<T, E> parameters_type;

  // ---- Construction ------------------------------------------------------
  SCM(parameters_type p, environment_type e, plant::Control c);

  // ---- Simulation lifecycle ----------------------------------------------

  // Run the whole schedule from t = 0 to completion.
  void run();

  // Advance one step: introduce every node due at the current time, then
  // integrate the patch forward to the next introduction (or over the fixed
  // ode times). Returns the species indices introduced this step.
  std::vector<size_t> run_next();

  // Replay the resident's saved environment for a mutant strategy: swap in
  // mutant parameters, reuse the cached environment/ode times, and run.
  void run_mutant(parameters_type p);

  // Adaptively refine the node-introduction schedule entirely in C++:
  // repeatedly run, flag nodes whose combined error exceeds schedule_eps,
  // and bisect the interval below each flagged node (upwind scheme), up to
  // schedule_nsteps times. Replaces the R build_schedule loop.
  void refine_schedule();

  // Return patch, schedule and solver to their t = 0 state; clear history.
  void reset();

  // True once every scheduled node introduction has been consumed.
  bool complete() const;

  // Current patch time.
  double time() const;

  // ---- Outputs -----------------------------------------------------------
  // Total (not per-capita) offspring. These delegate to the patch, which owns
  // the fitness/offspring computations.
  std::vector<double> net_reproduction_ratios() const { return patch.net_reproduction_ratios(); }
  std::vector<double> offspring_production() const { return patch.offspring_production(); }

  // ---- R interface -------------------------------------------------------

  // Run / parameters / state access
  parameters_type r_parameters() const { return parameters; }
  const patch_type &r_patch() const { return patch; }
  const std::vector<patch_type> &r_history() const { return history; }
  Rcpp::List r_get_state() const { return patch.r_get_state(); };

  // Fitness / reproduction
  double r_net_reproduction_ratio_for_species(util::index species_index) const;
  std::vector<std::vector<double>> r_net_reproduction_ratio_errors() const;

  // Schedule-refinement error signals.
  // Per-node refinement error: element-wise max of the competition error
  // (sampled during the run) and the reproduction error (computed at the end).
  // This is the signal that drives refine_schedule().
  std::vector<std::vector<double>> refinement_error_by_node() const;
  std::vector<double>
  r_compute_competition_effect_error_by_node_for_species_i(util::index species_index) const;

  // Node schedule access
  NodeSchedule r_node_schedule() const { return node_schedule; }
  void r_set_node_schedule(NodeSchedule x);
  void r_set_node_schedule_times(std::vector<std::vector<double>> x);

  // ODE times: the step times the solver actually used on the last run.
  // (Whether to *pin* integration to a fixed set of times is controlled on the
  // NodeSchedule via its own use_ode_times flag.)
  std::vector<double> r_ode_times() const;

  // ---- Public state ------------------------------------------------------
  // The two toggles are exposed to R directly (access: field), so they need
  // no getter/setter wrappers.
  bool collect;                    // record a patch snapshot after each step
  bool collect_refinement_errors;  // accumulate competition errors during run
  std::vector<patch_type> history; // per-step patch snapshots when collect

private:
  // Upwind bisection: insert the midpoint of the interval below each flagged
  // node. Mirrors split_times() in build_schedule.R.
  static std::vector<double> bisect_flagged_intervals(const std::vector<double>& times,
                                                      const std::vector<bool>& split);

  parameters_type parameters;
  Control control;
  patch_type patch;
  NodeSchedule node_schedule;
  ode::Solver<patch_type> solver;
};

// ---- Construction --------------------------------------------------------

template <typename T, typename E>
SCM<T, E>::SCM(parameters_type p, environment_type e, Control c)
    : parameters(p), control(c), patch(parameters, e, c),
      node_schedule(make_node_schedule(parameters)),
      solver(patch, make_ode_control(c)) {

  parameters.validate();

  collect = false;
  collect_refinement_errors = false;

  if (!util::identical(parameters.patch_area, 1.0)) {
    util::warning("We recommened keeping patch_area = 1 for the SCM, as need to check units for all other sizes");
  }
}

// ---- Simulation lifecycle ------------------------------------------------

template <typename T, typename E> void SCM<T, E>::run() {
  reset();
  if (collect) {
    history.push_back(patch);
  }

  while (!complete()) {
    std::vector<size_t> added = run_next();
    if (collect_refinement_errors) {
      patch.collect_competition_errors(added);
    }
    if (collect) {
      history.push_back(patch);
    }
  }
}

template <typename T, typename E> std::vector<size_t> SCM<T, E>::run_next() {
  std::vector<size_t> ret;
  const double t0 = time();

  // Consume every event scheduled at the current time t0: each contributes a
  // species to introduce. Stop once the next event ends later than t0 (i.e. it
  // belongs to a later introduction) or the schedule is exhausted.
  NodeSchedule::Event e = node_schedule.next_event();
  while (true) {
    if (!util::identical(t0, e.time_introduction())) {
      util::stop("Start time not what was expected");
    }
    ret.push_back(e.species_index);
    node_schedule.pop();
    if (e.time_end() > t0 || complete()) {
      break;
    } else {
      e = node_schedule.next_event();
    }
  }

  patch.introduce_new_nodes(ret);
  solver.set_state_from_system(patch);

  // Some schedules pin the integration points (e.g. when replaying a resident
  // run for a mutant); otherwise integrate adaptively to the next event time.
  if (node_schedule.using_ode_times()) {
    solver.advance_fixed(patch, e.times);
  } else {
    solver.advance_adaptive(patch, e.time_end());
  }

  return ret;
}

template <typename T, typename E>
void SCM<T, E>::run_mutant(parameters_type p) {

  // Switch the patch to its cached (resident) environment.
  patch.set_mutant();

  // Destructive: overwrite the resident parameters with the mutant's.
  parameters = p;

  // Swap in the mutant strategies.
  patch.overwrite_strategies(parameters.strategies);

  // Rebuild the schedule for the new parameters, then pin its integration
  // points to the resident's step history so the mutant sees the same
  // environment trajectory.
  node_schedule = make_node_schedule(parameters);
  node_schedule.r_set_ode_times(patch.step_history);
  node_schedule.r_set_use_ode_times(true);
  node_schedule.reset();

  // Re-initialise solver/patch and run.
  reset();
  run();
}

// Upwind bisection of flagged intervals. For each flagged node j (j >= 1; the
// first and last nodes are never flagged), insert the midpoint of the interval
// (t[j-1], t[j]). Equivalent to sort(c(times, times[i] - dt[i-1]/2)) in R.
template <typename T, typename E>
std::vector<double> SCM<T, E>::bisect_flagged_intervals(const std::vector<double>& times,
                                                        const std::vector<bool>& split) {
  std::vector<double> ret = times;
  for (size_t j = 1; j < split.size(); ++j) {
    if (split[j]) {
      ret.push_back(0.5 * (times[j] + times[j - 1]));
    }
  }
  std::sort(ret.begin(), ret.end());
  return ret;
}

template <typename T, typename E>
void SCM<T, E>::refine_schedule() {
  collect_refinement_errors = true;
  const double eps = control.schedule_eps;

  for (size_t step = 0; step < control.schedule_nsteps; ++step) {
    run(); // resets, then runs with collect_refinement_errors set

    std::vector<std::vector<double>> node_error = refinement_error_by_node();

    // Flag nodes whose refinement error exceeds the threshold.
    std::vector<std::vector<bool>> split(node_error.size());
    bool any = false;
    for (size_t i = 0; i < node_error.size(); ++i) {
      split[i].assign(node_error[i].size(), false);
      for (size_t j = 0; j < node_error[i].size(); ++j) {
        if (node_error[i][j] > eps) {
          split[i][j] = true;
          any = true;
        }
      }
    }
    if (!any) {
      break; // converged: no interval needs refining
    }

    // Bisect flagged intervals and install the denser schedule.
    std::vector<std::vector<double>> times = node_schedule.get_times();
    for (size_t i = 0; i < times.size(); ++i) {
      times[i] = bisect_flagged_intervals(times[i], split[i]);
    }
    node_schedule.set_times(times);
  }

  // Leave Parameters self-describing: record the refined schedule and the
  // ode times from the final run (mirrors build_schedule.R).
  parameters.node_schedule_times = node_schedule.get_times();
  parameters.ode_times = r_ode_times();
}

// NOTE: solver.reset() sets the solver's internal time to zero. There is
// currently no other way to set that time; it might be cleaner to add an
// ode::Solver::set_time and call set_time(0) explicitly here.
template <typename T, typename E> void SCM<T, E>::reset() {
  patch.reset();
  node_schedule.reset();
  solver.reset(patch);
  history.clear();
}

template <typename T, typename E> bool SCM<T, E>::complete() const {
  return node_schedule.remaining() == 0;
}

template <typename T, typename E> double SCM<T, E>::time() const {
  return patch.time();
}

// ---- R interface ---------------------------------------------------------
//
// The fitness/offspring and per-node error computations live on the patch
// (patch.h); the SCM methods below are thin facades that preserve the R API.
//
// Several of these are diagnostic/inspection hooks rather than part of the
// production run path: outside the C++ refinement loop they are only called
// from the test suite and the node_spacing vignette (noted per method below).

// Per-species fitness: the net reproduction ratio (expected offspring per seed)
// for one species. A genuine biological quantity, not just a diagnostic.
template <typename T, typename E>
double SCM<T, E>::r_net_reproduction_ratio_for_species(
    util::index species_index) const {
  const size_t idx = species_index.check_bounds(patch.size());
  auto scalars = std::vector<double>(patch.at_species(idx).size(), 1.0);
  return patch.net_reproduction_ratio_for_species(idx, scalars);
}

// Diagnostic: per-node discretisation error in the reproduction integral. One
// of the two components of refinement_error_by_node; exposed for inspection
// and validation (tests / vignette).
template <typename T, typename E>
std::vector<std::vector<double>>
SCM<T, E>::r_net_reproduction_ratio_errors() const {
  return patch.net_reproduction_ratio_errors();
}

// The combined per-node refinement error. Used internally by refine_schedule();
// also exposed to R so tests / the vignette can inspect the signal that drives
// schedule refinement.
template <typename T, typename E>
std::vector<std::vector<double>> SCM<T, E>::refinement_error_by_node() const {
  return patch.refinement_error_by_node();
}

// Diagnostic probe: the per-node competition (light) error for one species --
// the per-step sample that collect_competition_errors() accumulates. Exposed
// mainly so tests / the vignette can reconstruct the error signal by hand.
template <typename T, typename E>
std::vector<double>
SCM<T, E>::r_compute_competition_effect_error_by_node_for_species_i(util::index species_index) const {
  // TODO(#478): I think we need to scale this by total area; that should be
  // computed for everything so will get passed in as an argument.
  // const double tot_competition_effect  = patch.compute_competition(0.0);
  const size_t idx = species_index.check_bounds(patch.size());
  return patch.r_compute_competition_effect_error_by_node_for_species_i(idx);
}

template <typename T, typename E>
void SCM<T, E>::r_set_node_schedule(NodeSchedule x) {
  if (patch.node_ode_size() > 0) {
    util::stop("Cannot set schedule without resetting first");
  }
  util::check_length(x.get_n_species(), patch.size());
  node_schedule = x;

  // Update here so that extracting Parameters reflects the new schedule,
  // keeping Parameters self-sufficient.
  parameters.node_schedule_times = node_schedule.get_times();
}

template <typename T, typename E>
void SCM<T, E>::r_set_node_schedule_times(
    std::vector<std::vector<double>> x) {
  if (patch.node_ode_size() > 0) {
    util::stop("Cannot set schedule without resetting first");
  }
  node_schedule.set_times(x);
  parameters.node_schedule_times = x;
}

template <typename T, typename E>
std::vector<double> SCM<T, E>::r_ode_times() const {
  return solver.get_times();
}

} // namespace plant

#endif
