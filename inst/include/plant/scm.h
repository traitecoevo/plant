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

template <typename T, typename E> class SCM {
public:
  typedef T                strategy_type;
  typedef E                environment_type;
  typedef Individual<T, E> individual_type;
  typedef Node<T, E>       node_type;
  typedef Species<T, E>    species_type;
  typedef Patch<T, E>      patch_type;
  typedef Parameters<T, E> parameters_type;

  SCM(parameters_type p, environment_type e, plant::Control c);

  void run();
  void run_mutant(parameters_type p);
  std::vector<size_t> run_next();

  // Adaptively refine the node-introduction schedule entirely in C++:
  // repeatedly run, flag nodes whose combined error exceeds schedule_eps,
  // and bisect the interval below each flagged node (upwind scheme), up to
  // schedule_nsteps times. Replaces the R build_schedule loop.
  void refine_schedule();

  double time() const;
  void reset();
  bool complete() const;

  // * Output total offspring calculation (not per capita)
  // These delegate to the patch, which owns the fitness/offspring computations.
  std::vector<double> net_reproduction_ratios() const { return patch.net_reproduction_ratios(); }
  std::vector<double> offspring_production() const { return patch.offspring_production(); }

  // * R interface
  std::vector<util::index> r_run_next();
  parameters_type r_parameters() const { return parameters; }
  const patch_type &r_patch() const { return patch; }
  const std::vector <patch_type> &r_history() const { return history; }

  double r_net_reproduction_ratio_for_species(util::index species_index) const;
  std::vector<std::vector<double>> r_net_reproduction_ratio_errors() const;

  // Per-node refinement error: element-wise max of the competition error
  // (sampled during the run) and the reproduction error (computed at the end).
  // This is the signal that drives schedule refinement.
  std::vector<std::vector<double>> combined_node_errors() const;
  bool r_get_collect_errors() const { return collect_errors; }
  void r_set_collect_errors(bool x) { collect_errors = x; }
  std::vector<double>
  r_compute_competition_effect_error_by_node_for_species_i(util::index species_index) const;
  std::vector<double> r_ode_times() const;
  
  bool r_use_ode_times() const;
  void r_set_use_ode_times(bool x);

  bool r_get_collect() const;
  void r_set_collect(bool x);

  NodeSchedule r_node_schedule() const { return node_schedule; }
  void r_set_node_schedule(NodeSchedule x);
  void r_set_node_schedule_times(std::vector<std::vector<double>> x);
  
  bool collect;
  bool collect_errors;
  std::vector<patch_type> history;

  Rcpp::List r_get_state() const { return patch.r_get_state(); };

private:
  // Upwind bisection: insert the midpoint of the interval below each flagged
  // node. Mirrors split_times() in build_schedule.R.
  static std::vector<double> split_times(const std::vector<double>& times,
                                         const std::vector<bool>& split);

  parameters_type parameters;
  Control control;
  patch_type patch;
  NodeSchedule node_schedule;
  ode::Solver<patch_type> solver;
};

template <typename T, typename E>
SCM<T, E>::SCM(parameters_type p, environment_type e, Control c)
    : parameters(p), control(c), patch(parameters, e, c),
      node_schedule(make_node_schedule(parameters)),
      solver(patch, make_ode_control(c)) {

  parameters.validate();

  collect = false;
  collect_errors = false;

  if (!util::identical(parameters.patch_area, 1.0)) {
    util::warning("We recommened keeping patch_area = 1 for the SCM, as need to check units for all other sizes");
  }
}

template <typename T, typename E> void SCM<T, E>::run() {
  reset();
  if (collect)
  {
    history.push_back(patch);
  }

  while (!complete()) {
    std::vector<size_t> added = run_next();
    if (collect_errors) {
      patch.collect_competition_errors(added);
    }
    // store
    if(collect)
    {
      history.push_back(patch);
    }
  }
}

template <typename T, typename E> std::vector<size_t> SCM<T, E>::run_next() {
  std::vector<size_t> ret;
  const double t0 = time();

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
  
  // some schedules have fixed integration points
  const bool use_ode_times = node_schedule.using_ode_times();
  
  if (use_ode_times) {
    solver.advance_fixed(patch, e.times);
  } else {
    solver.advance_adaptive(patch, e.time_end());
  }

  return ret;
}

// Upwind bisection of flagged intervals. For each flagged node j (j >= 1; the
// first and last nodes are never flagged), insert the midpoint of the interval
// (t[j-1], t[j]). Equivalent to sort(c(times, times[i] - dt[i-1]/2)) in R.
template <typename T, typename E>
std::vector<double> SCM<T, E>::split_times(const std::vector<double>& times,
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
  collect_errors = true;
  const double eps = control.schedule_eps;

  for (size_t step = 0; step < control.schedule_nsteps; ++step) {
    run(); // resets, then runs with collect_errors set

    std::vector<std::vector<double>> total = combined_node_errors();

    // Flag nodes whose combined error exceeds the threshold.
    std::vector<std::vector<bool>> split(total.size());
    bool any = false;
    for (size_t i = 0; i < total.size(); ++i) {
      split[i].assign(total[i].size(), false);
      for (size_t j = 0; j < total[i].size(); ++j) {
        if (total[i][j] > eps) {
          split[i][j] = true;
          any = true;
        }
      }
    }
    if (!any) {
      break;
    }

    // Bisect flagged intervals and install the denser schedule.
    std::vector<std::vector<double>> times = node_schedule.get_times();
    for (size_t i = 0; i < times.size(); ++i) {
      times[i] = split_times(times[i], split[i]);
    }
    node_schedule.set_times(times);
  }

  // Leave Parameters self-describing: record the refined schedule and the
  // ode times from the final run (mirrors build_schedule.R).
  parameters.node_schedule_times = node_schedule.get_times();
  parameters.ode_times = r_ode_times();
}

template <typename T, typename E>
void SCM<T, E>::run_mutant(parameters_type p) {
  
  // switch to cached environment
  patch.set_mutant();

  // destructive operation; overwrites resident params.
  parameters = p;

  // add strategies
  patch.overwrite_strategies(parameters.strategies);

  // resize schedule
  node_schedule = make_node_schedule(parameters);
  
  // then set ode_times to patch history
  node_schedule.r_set_ode_times(patch.step_history);
  node_schedule.r_set_use_ode_times(true);
  node_schedule.reset();

  // re-initialise solver
  reset();

  run();
}

template <typename T, typename E> double SCM<T, E>::time() const {
  return patch.time();
}

// NOTE: solver.reset() will set time within the solver to zero.
// However, there is no other current way of setting the time within
// the solver.  It might be better to add a set_time method within
// ode::Solver, and then here do explicitly ode_solver.set_time(0)?
template <typename T, typename E> void SCM<T, E>::reset() {
  patch.reset();
  node_schedule.reset();
  solver.reset(patch);
  history.clear();
}

template <typename T, typename E> bool SCM<T, E>::complete() const {
  return node_schedule.remaining() == 0;
}

template <typename T, typename E>
std::vector<util::index> SCM<T, E>::r_run_next() {
  return util::index_vector(run_next());
}

template <typename T, typename E>
std::vector<double>
SCM<T, E>::r_compute_competition_effect_error_by_node_for_species_i(util::index species_index) const {
  // TODO: I think we need to scale this by total area; that should be
  // computed for everything so will get passed in as an argument.
  // const double tot_competition_effect  = patch.compute_competition(0.0);
  const size_t idx = species_index.check_bounds(patch.size());
  return patch.r_compute_competition_effect_error_by_node_for_species_i(idx);
}

template <typename T, typename E>
std::vector<double> SCM<T, E>::r_ode_times() const {
  return solver.get_times();
}

template <typename T, typename E> bool SCM<T, E>::r_use_ode_times() const {
  return node_schedule.using_ode_times();
}

template <typename T, typename E> void SCM<T, E>::r_set_use_ode_times(bool x) {
  node_schedule.r_set_use_ode_times(x);
}


template <typename T, typename E> bool SCM<T, E>::r_get_collect() const {
  return collect;
}

template <typename T, typename E> void SCM<T, E>::r_set_collect(bool x) {
    collect = x;
}



template <typename T, typename E>
void SCM<T, E>::r_set_node_schedule(NodeSchedule x) {
  if (patch.node_ode_size() > 0) {
    util::stop("Cannot set schedule without resetting first");
  }
  util::check_length(x.get_n_species(), patch.size());
  node_schedule = x;

  // Update these here so that extracting Parameters would give the
  // new schedule, this making Parameters sufficient.
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


// The fitness/offspring and per-node error computations live on the patch
// (patch.h); the SCM methods below are thin facades that preserve the R API.

template <typename T, typename E>
double SCM<T, E>::r_net_reproduction_ratio_for_species(
    util::index species_index) const {
  const size_t idx = species_index.check_bounds(patch.size());
  auto scalars = std::vector<double>(patch.at_species(idx).size(), 1.0);
  return patch.net_reproduction_ratio_for_species(idx, scalars);
}

template <typename T, typename E>
std::vector<std::vector<double>> SCM<T, E>::combined_node_errors() const {
  return patch.combined_node_errors();
}

template <typename T, typename E>
std::vector<std::vector<double>>
SCM<T, E>::r_net_reproduction_ratio_errors() const {
  return patch.net_reproduction_ratio_errors();
}

} // namespace plant

#endif
