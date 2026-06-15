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

  double time() const;
  void reset();
  bool complete() const;

  // * Output total offspring calculation (not per capita)
  std::vector<double> net_reproduction_ratio_by_node_weighted(size_t species_index) const;
  double net_reproduction_ratio_for_species(size_t species_index, std::vector<double> const& scalars) const;
  std::vector<double> net_reproduction_ratios() const;
  std::vector<double> offspring_production() const;

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
  double total_offspring_production() const;
  // Update the running per-node competition error for the species introduced
  // this step (mirrors the per-step sampling the R refinement loop did).
  void collect_competition_errors(const std::vector<size_t>& added);

  parameters_type parameters;
  patch_type patch;
  NodeSchedule node_schedule;
  ode::Solver<patch_type> solver;

  // Per-species running max of the competition error per node, accumulated
  // across the run when collect_errors is set. Entries start at -Inf and
  // ignore NA contributions, matching apply(., 2, max, na.rm=TRUE) in R.
  std::vector<std::vector<double>> competition_error_by_node;
};

template <typename T, typename E>
SCM<T, E>::SCM(parameters_type p, environment_type e, Control c)
    : parameters(p), patch(parameters, e, c),
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
      collect_competition_errors(added);
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
  competition_error_by_node.assign(patch.size(), {});
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


// Offspring production, equal to overall fitness scaled by the birth rate
template <typename T, typename E>
std::vector<double> SCM<T, E>::offspring_production() const {
	auto ret = std::vector<double>(patch.size());
  for (size_t i = 0; i < patch.size(); ++i) {
		// scale by birth rate function over time
		auto const times = patch.at_species(i).node_times();
		auto scalars = std::vector<double>(times.size());
		for (size_t j = 0; j < times.size(); ++j) {
			scalars[j] = patch.at_species(i).extrinsic_drivers().evaluate("birth_rate", times[j]);
		}
		ret[i] = net_reproduction_ratio_for_species(i, scalars);
  }
  return ret;
}

// Overall fitness
template <typename T, typename E>
std::vector<double> SCM<T, E>::net_reproduction_ratios() const {
	auto ret = std::vector<double>(patch.size());
  for (size_t i = 0; i < patch.size(); ++i) {
		// no scaling, ie set scalars to 1.0
		auto scalars = std::vector<double>(patch.at_species(i).size(), 1.0);
		ret[i] = net_reproduction_ratio_for_species(i, scalars);
  }
  return ret;
}

// Integrate over lifetime fitness of individual nodes
template <typename T, typename E>
double
SCM<T, E>::net_reproduction_ratio_for_species(size_t species_index, std::vector<double> const& scalars) const {
	auto net_prod = net_reproduction_ratio_by_node_weighted(species_index);
	auto const times = patch.at_species(species_index).node_times();
	auto net_prod_scaled = std::vector<double>(times.size());
	for (size_t i = 0; i < times.size(); ++i) {
			net_prod_scaled[i] = net_prod[i] * scalars[i];
	}
  return util::trapezium(
      times,
      net_prod_scaled
	);
}

// R interface method
template <typename T, typename E>
double SCM<T, E>::r_net_reproduction_ratio_for_species(
    util::index species_index) const {
	const size_t idx = species_index.check_bounds(patch.size());
	auto scalars = std::vector<double>(patch.at_species(idx).size(), 1.0);
  return net_reproduction_ratio_for_species(idx, scalars);
}

// Node fitness within a meta-population of patches.
// The patch-age density weighting and S_D are now recorded on each node at
// introduction, so this is just a passthrough to the species.
template <typename T, typename E>
std::vector<double> SCM<T, E>::net_reproduction_ratio_by_node_weighted(
    size_t species_index) const {
  return patch.at_species(species_index).net_reproduction_ratio_by_node_weighted();
}

// Sum up all offspring produced
template <typename T, typename E>
double SCM<T, E>::total_offspring_production() const {
  double total = 0.0;
  std::vector<double> offspring = offspring_production();
  for (size_t i = 0; i < patch.size(); ++i) {
    total += offspring[i];
  }
  return total;
}

// Sample the competition error for each species introduced this step and fold
// it into the running per-node max (ignoring NA, matching na.rm=TRUE in R).
template <typename T, typename E>
void SCM<T, E>::collect_competition_errors(const std::vector<size_t>& added) {
  for (size_t idx : added) {
    std::vector<double> v =
        patch.r_compute_competition_effect_error_by_node_for_species_i(idx);
    std::vector<double>& acc = competition_error_by_node[idx];
    if (acc.size() < v.size()) {
      acc.resize(v.size(), -std::numeric_limits<double>::infinity());
    }
    for (size_t j = 0; j < v.size(); ++j) {
      if (!ISNAN(v[j])) {
        acc[j] = std::max(acc[j], v[j]);
      }
    }
  }
}

// Combine the competition error (sampled during the run) with the reproduction
// error (computed now) into a single per-node error vector per species. An
// all-NA node yields -Inf, matching apply(rbind(...), 2, max, na.rm=TRUE) in R.
template <typename T, typename E>
std::vector<std::vector<double>> SCM<T, E>::combined_node_errors() const {
  std::vector<std::vector<double>> repro = r_net_reproduction_ratio_errors();
  std::vector<std::vector<double>> ret(patch.size());
  for (size_t i = 0; i < patch.size(); ++i) {
    const std::vector<double>& comp = competition_error_by_node[i];
    const std::vector<double>& rep = repro[i];
    const size_t n = patch.at_species(i).size();
    std::vector<double> tot(n, -std::numeric_limits<double>::infinity());
    for (size_t j = 0; j < n; ++j) {
      if (j < comp.size() && !ISNAN(comp[j])) {
        tot[j] = std::max(tot[j], comp[j]);
      }
      if (j < rep.size() && !ISNAN(rep[j])) {
        tot[j] = std::max(tot[j], rep[j]);
      }
    }
    ret[i] = tot;
  }
  return ret;
}

// Check integration errors
template <typename T, typename E>
std::vector<std::vector<double>>
SCM<T, E>::r_net_reproduction_ratio_errors() const {
  std::vector<std::vector<double>> ret;
  double total_offspring = total_offspring_production();
  for (size_t i = 0; i < patch.size(); ++i) {
    ret.push_back(util::local_error_integration(
        patch.at_species(i).node_times(), net_reproduction_ratio_by_node_weighted(i),
        total_offspring));
  }
  return ret;
}

} // namespace plant

#endif
