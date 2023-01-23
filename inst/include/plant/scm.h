// -*-c++-*-
#ifndef PLANT_PLANT_SCM_H_
#define PLANT_PLANT_SCM_H_

#include <plant/node_schedule.h>
#include <plant/ode_solver.h>
#include <plant/patch.h>
#include <plant/scm_utils.h>

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
  void run_mutant();
  std::vector<size_t> run_next();\

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
  // TODO: These are liable to change to return all species at once by
  // default.  The pluralisation difference between
  // SCM::r_competition_effect_error and Species::r_competition_effects_error
  // will get dealt with then.
  double r_net_reproduction_ratio_for_species(util::index species_index) const;
  std::vector<std::vector<double>> r_net_reproduction_ratio_errors() const;
  std::vector<double>
  r_competition_effect_error(util::index species_index) const;
  std::vector<double> r_ode_times() const;
  bool r_use_ode_times() const;
  void r_set_use_ode_times(bool x);

  NodeSchedule r_node_schedule() const { return node_schedule; }
  void r_set_node_schedule(NodeSchedule x);
  void r_set_node_schedule_times(std::vector<std::vector<double>> x);

  std::vector<double> patch_step_history;
  std::vector<std::vector<environment_type>> environment_history;
  
private:
  double total_offspring_production() const;

  parameters_type parameters;
  patch_type patch;
  NodeSchedule node_schedule;
  ode::Solver<patch_type> solver;
};

template <typename T, typename E>
SCM<T, E>::SCM(parameters_type p, environment_type e, Control c)
    : parameters(p), patch(parameters, e, c),
      node_schedule(make_node_schedule(parameters)),
      solver(patch, make_ode_control(c)) {

  parameters.validate();
  if (!util::identical(parameters.patch_area, 1.0)) {
    util::stop("Patch area must be exactly 1 for the SCM");
  }
}

template <typename T, typename E> void SCM<T, E>::run() {
  reset();
  while (!complete()) {
    run_next();
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

  const bool use_ode_times = node_schedule.using_ode_times();
  solver.set_state_from_system(patch);
  if (use_ode_times) {
    solver.advance_fixed(patch, e.times);
  } else {
    solver.advance(patch, e.time_end());
  }

  // after each full RK45 step, save the patch cache
  if(patch.save_RK45_cache) {
    patch_step_history.push_back(t0);
    environment_history.push_back(patch.environment_cache);
  }

  return ret;
}

template <typename T, typename E> void SCM<T, E>::run_mutant() {
  patch.save_RK45_cache = false;
  patch.use_cached_environment = true;

  // need some way to add mutant strategy

  for(size_t i; i < patch_step_history.size(); ++i) {
    patch.environment_cache.clear();
    patch.environment_cache = environment_history[i];

    run_next(); 
  }
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
SCM<T, E>::r_competition_effect_error(util::index species_index) const {
  // TODO: I think we need to scale this by total area; that should be
  // computed for everything so will get passed in as an argument.
  // const double tot_competition_effect  = patch.compute_competition(0.0);
  const size_t idx = species_index.check_bounds(patch.size());
  return patch.r_competition_effect_error(idx);
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
		auto const& times = node_schedule.times(i);
		auto scalars = std::vector<double>(times.size());
		for (size_t j = 0; j < times.size(); ++j) {
			scalars[j] = patch.at(i).extrinsic_drivers().evaluate("birth_rate", times[j]);
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
		auto const& times = node_schedule.times(i);
		auto scalars = std::vector<double>(times.size(), 1.0);
		ret[i] = net_reproduction_ratio_for_species(i, scalars);
  }
  return ret;
}

// Integrate over lifetime fitness of individual nodes
template <typename T, typename E>
double
SCM<T, E>::net_reproduction_ratio_for_species(size_t species_index, std::vector<double> const& scalars) const {
	auto net_prod = net_reproduction_ratio_by_node_weighted(species_index);
	auto const& times = node_schedule.times(species_index);
	auto net_prod_scaled = std::vector<double>(times.size());
	// should be showing compiler warning for int (auto) comparison, but isn't anymore...
	for (auto i = 0; i < times.size(); ++i) {
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
	auto const& times = node_schedule.times(species_index.check_bounds(patch.size()));
	auto scalars = std::vector<double>(times.size(), 1.0);
  return net_reproduction_ratio_for_species(
      species_index.x, scalars);
}

// Node fitness within a meta-population of patches
template <typename T, typename E>
std::vector<double> SCM<T, E>::net_reproduction_ratio_by_node_weighted(
    size_t species_index) const {
  // node introduction times
  const std::vector<double> times = node_schedule.times(species_index);

  // retrieve lifetime fitness for each node
  std::vector<double> net_reproduction_ratio_by_node_weighted =
      patch.at(species_index).net_reproduction_ratio_by_node();

  // weight by probabilty of reproduction
  for (size_t i = 0; i < net_reproduction_ratio_by_node_weighted.size();
       ++i) {
    net_reproduction_ratio_by_node_weighted[i] *=
        patch.survival_weighting->density(
            times[i]) * // probability of landing in patch of a given age
        parameters.strategies[species_index]
            .S_D; // probability of survival during dispersal (assumed constant)
  }

  return net_reproduction_ratio_by_node_weighted;
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

// Check integration errors
template <typename T, typename E>
std::vector<std::vector<double>>
SCM<T, E>::r_net_reproduction_ratio_errors() const {
  std::vector<std::vector<double>> ret;
  double total_offspring = total_offspring_production();
  for (size_t i = 0; i < patch.size(); ++i) {
    ret.push_back(util::local_error_integration(
        node_schedule.times(i), net_reproduction_ratio_by_node_weighted(i),
        total_offspring));
  }
  return ret;
}

} // namespace plant

#endif
