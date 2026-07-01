// -*-c++-*-
#ifndef PLANT_PLANT_PATCH_H_
#define PLANT_PLANT_PATCH_H_

#include <plant/parameters.h>
#include <plant/species.h>
#include <plant/util.h>
#include <odelia/ode_interface.hpp>

#include <plant/disturbance_regime.h>

#include <algorithm>
#include <limits>

using namespace Rcpp;

namespace plant {

// Templated on the scalar type S (#472 scope B / #537, Milestone C) to mirror
// Species<T,E,S> / Node<T,E,S>: S defaults to double so every existing
// Patch<T,E> is Patch<T,E,double> and bit-identical. The species storage is
// Species<T,E,S>, so a Patch<...,ad> holds ad-typed individual state. NOTE the
// environment member stays type E (double): a Patch<...,ad> is MIXED -- ad
// species over a frozen-double resident environment, matching the two-pass
// replay design (freeze the resident light schedule, replay the node ODE with
// the active scalar). An ad-valued resident light spline (self-shading through
// traits) is a later piece -- the odelia AD interpolator (odelia PR #32).
template <typename T, typename E, typename S = double>
class Patch {
public:
  using value_type = double;

  typedef T                 strategy_type;
  typedef E                 environment_type;
  typedef Individual<T,E,S> individual_type;
  typedef Node<T,E,S>       node_type;
  typedef Species<T,E,S>    species_type;
  typedef Parameters<T,E>   parameters_type;

  Patch(parameters_type p, environment_type e, plant::Control c);
  void reset();
  size_t size() const {return species.size();}

  //Try using pointer in place of object itself
  double time() const {return environment.time;}
  double get_area() const { return area;}
  double height_max() const;

  double compute_competition(double height) const;

  // * Lifetime fitness / offspring production
  // These are patch-level quantities: each integrates the per-node weighted
  // net reproduction over a species' node-introduction times.
  // Integrate lifetime fitness of a species' nodes, scaled per node.
  double net_reproduction_ratio_for_species(size_t species_index,
                                            std::vector<double> const& scalars) const;
  // Offspring production: fitness scaled by the birth rate over time.
  std::vector<double> offspring_production() const;
  // Overall fitness (unscaled, scalars == 1).
  std::vector<double> net_reproduction_ratios() const;
  // Sum of offspring produced across all species.
  double total_offspring_production() const;
  // Per-node reproduction integration error for each species.
  std::vector<std::vector<double>> net_reproduction_ratio_errors() const;

  // * Schedule-refinement error collection
  // Sample the competition error for each species introduced this step and
  // fold it into the running per-node max (ignoring NA, matching na.rm=TRUE in
  // R). Accumulated across the run; cleared by reset().
  void collect_competition_errors(const std::vector<size_t>& added);
  // Combine the competition error (sampled during the run) with the
  // reproduction error (computed now) into a single per-node error vector per
  // species. An all-NA node yields -Inf. Drives schedule refinement.
  std::vector<std::vector<double>> refinement_error_by_node() const;

  void introduce_new_node(size_t species_index);
  void introduce_new_nodes(const std::vector<size_t>& species_index);

  // Open to better ways to test whether nodes have been introduced
  int node_ode_size() const {
    int node_ode_size = ode_size() - environment.ode_size();
    return(node_ode_size);
  }

  const species_type& at_species(size_t species_index) const {
    return species[species_index];
  }

  // Patch disturbance
  Disturbance_Regime* survival_weighting;

  // * ODE interface
  // Caluclate size of ode system (number of equations). Is constantly changing as 
  // new nodes are introduced into size-density distreibution 
  size_t ode_size() const;
  // How many auxiallary variables are we tracking. These are being collected but
  // are not part of core ode system
  size_t aux_size() const;
  double ode_time() const;

  // Retrieve ode state from patch and save into the ode solver
  odelia::ode::iterator ode_state(odelia::ode::iterator it) const;
  // Retrieve ode rates from patch and save into the ode solver
  odelia::ode::iterator ode_rates(odelia::ode::iterator it) const;
  // Retrieve auxillary variables and save into the ode solver
  odelia::ode::iterator ode_aux(odelia::ode::iterator it) const;

  // Returns state in structure format as opposed to single 
  // vector as given by ode_state
  Rcpp::List r_get_state() const;

  // Set state of patch, based on estimate of future state estimated by the solver
  // There are two implementations.
  //   - first function is for resident runs.
  //   - second is for mutant runs.
  // The second does not calculate environment when states are updated, as mutants only experience the environment
  // The decision which to use is determined by `use_cached_environment` below
  odelia::ode::const_iterator set_ode_state(odelia::ode::const_iterator it, double time);
  odelia::ode::const_iterator set_ode_state(odelia::ode::const_iterator it, int index);

  // * R interface
  // Data accessors:
  double r_density(double time) const {return survival_weighting->r_density(time);}
  double r_pr_survival(double time) const {return survival_weighting->pr_survival(time);}
  double r_disturbance_mean_interval() const {return survival_weighting->r_mean_interval();}
  double r_survival_weighting_cdf(double time) const {return survival_weighting->cdf(time);}
  double r_survival_weighting_icdf(double prob) const {return survival_weighting->icdf(prob);}

  parameters_type r_parameters() const {return parameters;}
  environment_type r_environment() const {return environment;}
  std::vector<species_type> r_species() const {return species;}
  std::vector<double> r_compute_competition_effect_error_by_node_for_species_i(size_t species_index) const;
  void r_set_time(double time);
  void r_set_state(double time,
                   const std::vector<double>& state,
                   const std::vector<size_t>& n,
                   const std::vector<double>& light_availability);
  void r_introduce_new_node(util::index species_index) {
    introduce_new_node(species_index.check_bounds(size()));
  }
  species_type r_at(util::index species_index) const {
    at(species_index.check_bounds(size()));
  }
  // These are only here because they wrap private functions.
  void r_compute_environment() {compute_environment(false);}
  void r_compute_rates() {
    environment_ptr = &environment;
    compute_rates();
  }

  // env. cache for assembly
  std::vector<double> step_history{0.0};  // always start at zero
  std::vector<std::vector<environment_type>> environment_history;
  std::vector<environment_type> environment_cache;

  // Per-ODE-step resident STAND state for species 0, captured alongside
  // environment_history during a save_RK45_cache run: each entry is the node
  // heights / per-node competition effects (node.compute_competition(0) =
  // k_I * area_leaf, so resident competition(z) = trapezium_i(ce_i * Q(z/h_i))).
  // Exposed so the active-knot self-shading AD driver (#472 scope B) can
  // reconstruct the resident light DIFFERENTIABLY per step (light responds to an
  // allometric trait via area_leaf) instead of holding the cached env frozen.
  // Single-species for now (the FF16 self-shading demo); multi-species is additive.
  std::vector<std::vector<double>> stand_height_history;
  std::vector<std::vector<double>> stand_competition_history;

  // Per-RK-STAGE resident stand state for species 0, captured alongside the 6
  // per-stage environments in environment_history (one inner vector per Cash-Karp
  // stage, mirroring environment_cache exactly). This is the faithful-to-the-SCM
  // harvest the RESIDENT total-gradient build (#472 scope B, R0) needs: at each
  // RK stage the active light is reconstructed from the stand that PRODUCED that
  // stage's frozen env (heights h_i, per-node competition effect ce_i), so the
  // resident feedback derivative tracks the SCM stage-by-stage rather than holding
  // a step-start census across stages (the cheaper per-step C-28 reconstruction).
  // [step][stage 0..5][cohort]. Single-species for now (the FF16 light demo).
  std::vector<std::vector<std::vector<double>>> stand_height_stage_history;
  std::vector<std::vector<std::vector<double>>> stand_competition_stage_history;
  // Per-step working caches (6 stages), filled by cache_RK45_step and flushed to
  // the *_stage_history by cache_ode_step, exactly as environment_cache is.
  std::vector<std::vector<double>> stand_height_stage_cache;
  std::vector<std::vector<double>> stand_competition_stage_cache;

  // Boundary (new_node) per-RK-stage state for species 0: its height + competition
  // effect. r_compute_competition_effect_by_nodes / r_heights iterate `nodes` only
  // and OMIT the boundary new_node, whose tail term Species::compute_competition
  // DOES add (it is the bottom of the last trapezium segment). Without it the
  // coupled-replay env reconstruction (#472 scope B, R0) is wrong at ground level
  // (z < the seedling top), where the boundary node's projected leaf area lands.
  // [step][stage 0..5], a single scalar per stage (species 0). Empty for an empty
  // species (a NaN-free placeholder is pushed so indices stay aligned).
  std::vector<std::vector<double>> stand_newnode_height_stage_history;
  std::vector<std::vector<double>> stand_newnode_competition_stage_history;
  std::vector<double> stand_newnode_height_stage_cache;
  std::vector<double> stand_newnode_competition_stage_cache;
  // Boundary-node state SNAPSHOT taken at the moment compute_environment builds the
  // light spline (species 0). compute_rates() then calls new_node.compute_initial_
  // conditions and MUTATES the boundary node, so querying it at cache time (after
  // derivs = compute_environment + compute_rates) reads the wrong state. These hold
  // the new_node baked into the cached env. (#472 scope B, R0 reconstruction.)
  double newnode_height_env_snapshot = 0.0;
  double newnode_competition_env_snapshot = 0.0;

  // ALL-SPECIES per-RK-stage harvest (#472 scope B, R2 -- the cross-species coupled
  // Jacobian). The species-0 fields above drive the single-species coupled replay;
  // these add the species dimension so the JOINT canopy can be reconstructed from
  // every species' re-evolved cohorts. [step][stage 0..5][species][cohort] for the
  // stand, [step][stage][species] for the boundary node. Additive: the single-species
  // path is untouched. Filled only on save_RK45_cache runs.
  std::vector<std::vector<std::vector<std::vector<double>>>> stand_height_stage_history_all;
  std::vector<std::vector<std::vector<std::vector<double>>>> stand_competition_stage_history_all;
  std::vector<std::vector<std::vector<double>>> stand_newnode_height_stage_history_all;
  std::vector<std::vector<std::vector<double>>> stand_newnode_competition_stage_history_all;
  std::vector<std::vector<std::vector<double>>> stand_height_stage_cache_all;
  std::vector<std::vector<std::vector<double>>> stand_competition_stage_cache_all;
  std::vector<std::vector<double>> stand_newnode_height_stage_cache_all;
  std::vector<std::vector<double>> stand_newnode_competition_stage_cache_all;
  std::vector<double> newnode_height_env_snapshot_all;
  std::vector<double> newnode_competition_env_snapshot_all;

  void cache_ode_step();
  void cache_RK45_step(int step);
  void load_ode_step();
  
  // used cache_ode_step for mutant runs
  bool save_RK45_cache;

  // used in load_ode_step for mutant runs
  bool use_cached_environment = false;

  bool is_mutant_run = false;

  void set_mutant();
  void add_strategies(std::vector<strategy_type> strategies);
  void overwrite_strategies(std::vector<strategy_type> strategies);

private:
  int idx = 0; // used to access environment cache for mutant runs
  void compute_environment(bool rescale);
  void compute_rates();

  // Seed the patch from parameters.initial_state (nodes + birth bookkeeping)
  // when present; called from reset(). Sets environment.time = initial_time.
  void set_initial_state();
  // Guard against initial conditions whose per-node log-density rates are so
  // large they would drive densities to non-finite values within a few steps.
  void check_initial_density_rates() const;

  parameters_type parameters;

  double area;
  environment_type environment;
  std::vector<species_type> species;

  //TODO(#476): Move into environment?
  std::vector<double> resource_depletion;

  environment_type* environment_ptr;

  Control control;

  // Per-species running max of the competition error per node, accumulated
  // across the run via collect_competition_errors(). Entries start at -Inf and
  // ignore NA contributions, matching apply(., 2, max, na.rm=TRUE) in R.
  std::vector<std::vector<double>> competition_error_by_node;
};

template <typename T, typename E, typename S>
Patch<T,E,S>::Patch(parameters_type p, environment_type e, Control c)
  : parameters(p),
    area(p.patch_area),
    environment(e),
    control(c),
    environment_cache(6) {  // length of odelia::ode::Step
  
  parameters.validate();

  save_RK45_cache = control.save_RK45_cache;
  survival_weighting = p.disturbance;

  // Configure the light profile's shading model before the first
  // compute_environment() in reset(). No-op for environments without a light
  // profile (only FF16 implements alternative shading models).
  environment.set_shading_model(control.shading_model,
                                control.ppa_layer_optical_depth,
                                control.ppa_layer_smoothing);

  add_strategies(parameters.strategies);

  reset();
}

template <typename T, typename E, typename S>
void Patch<T,E,S>::overwrite_strategies(std::vector<strategy_type> strategies) {
  species.clear();
  add_strategies(strategies);
}

template <typename T, typename E, typename S>
void Patch<T,E,S>::add_strategies(std::vector<strategy_type> strategies) {
  for (auto i = 0; i < strategies.size(); ++i) {
		auto s = strategies[i];
    s.control = control; // Overwrite to take the patch control object 
    auto spec = Species<T,E,S>(s);
    species.push_back(spec);
  }
}

template <typename T, typename E, typename S>
void Patch<T,E,S>::set_mutant() {
    if (environment_history.empty()) {
       util::stop("Run a resident first to generate a competitve landscape");
    }

    is_mutant_run = true;
    save_RK45_cache = false;
    use_cached_environment = true;
  idx = 0;
}

template <typename T, typename E, typename S>
void Patch<T,E,S>::reset() {
   for (auto& s : species) {
    s.clear();
    // allocate variables for tracking resource consumption
    s.resize_consumption_rates(environment.ode_size());
  }

  // resize to species count
  resource_depletion.reserve(environment.ode_size());

  // compute ephemeral effects like light_availability
  environment.clear();

  if (!parameters.initial_state.empty()) {
    // Seed the patch from an exported state / initial size distribution.
    // set_initial_state() does its own compute_environment(false)/compute_rates
    // (with the real node heights, no rescale), so skip the empty-patch path.
    set_initial_state();
    check_initial_density_rates();
  } else {
    compute_environment(false);

    // compute effects of resource consumption
    environment_ptr = &environment;
    compute_rates();
  }

  // clear accumulated per-node competition error
  competition_error_by_node.assign(species.size(), {});
}

// Seed the patch from parameters.initial_state. Introduces the requested number
// of nodes per species, loads the flat ODE state (nodes + environment) at the
// resume time, restores per-node birth bookkeeping (not part of the ODE state
// but feeds the rates and lifetime-fitness integrals), then computes the
// environment and rates. We load the state directly (rather than via the
// double-arg set_ode_state) so the first environment build is a full
// compute_environment(false): a rescale of the not-yet-built light spline would
// read uninitialised grid state.
template <typename T, typename E, typename S>
void Patch<T,E,S>::set_initial_state() {
  const size_t n_species = species.size();
  util::check_length(parameters.n_initial_cohorts.size(), n_species);

  size_t total_nodes = 0;
  for (size_t i = 0; i < n_species; ++i) {
    for (size_t j = 0; j < parameters.n_initial_cohorts[i]; ++j) {
      species[i].introduce_new_node();
    }
    total_nodes += parameters.n_initial_cohorts[i];
  }

  // Load the flat ODE state (all nodes, then environment).
  util::check_length(parameters.initial_state.size(), ode_size());
  odelia::ode::const_iterator it = parameters.initial_state.begin();
  it = odelia::ode::set_ode_state(species.begin(), species.end(), it);
  it = environment.set_ode_state(it);
  environment.time = parameters.initial_time;

  // Restore birth bookkeeping per node, sliced per species from the flat
  // parameter vectors. Skipped (left at defaults) when not supplied, e.g. a
  // from-scratch distribution seeded at patch age 0.
  if (!parameters.initial_node_times.empty()) {
    util::check_length(parameters.initial_node_times.size(), total_nodes);
    util::check_length(parameters.initial_patch_density.size(), total_nodes);
    util::check_length(parameters.initial_pr_patch_survival.size(), total_nodes);
    auto t_it = parameters.initial_node_times.begin();
    auto d_it = parameters.initial_patch_density.begin();
    auto s_it = parameters.initial_pr_patch_survival.begin();
    for (size_t i = 0; i < n_species; ++i) {
      const size_t n = parameters.n_initial_cohorts[i];
      species[i].set_birth_state(std::vector<double>(t_it, t_it + n),
                                 std::vector<double>(d_it, d_it + n),
                                 std::vector<double>(s_it, s_it + n));
      t_it += n;
      d_it += n;
      s_it += n;
    }
  }

  // Build the environment from the real node heights (full recompute, no
  // rescale) and compute rates for the seeded population.
  compute_environment(false);
  environment_ptr = &environment;
  compute_rates();
}

template <typename T, typename E, typename S>
void Patch<T,E,S>::check_initial_density_rates() const {
  for (const auto& s : species) {
    std::vector<double> rates = s.r_log_density_rates();
    if (std::any_of(rates.begin(), rates.end(),
                    [](double v) { return v < -100; })) {
      util::stop("Rates of initial node densities exceed ~1e43 and will likely "
                 "produce non-finite densities; provide more plausible initial "
                 "conditions (smaller sizes and/or lower densities).");
    }
  }
}

template <typename T, typename E, typename S>
double Patch<T,E,S>::height_max() const {
  double ret = 0.0;
  for (size_t i = 0; i < species.size(); ++i) {
    if (!is_mutant_run) {
      ret = std::max(ret, species[i].height_max());
    }
  }
  return ret;
}

template <typename T, typename E, typename S>
double Patch<T,E,S>::compute_competition(double height) const {
  double tot = 0.0;
  for (size_t i = 0; i < species.size(); ++i) {
    if (!is_mutant_run) {
      tot += species[i].compute_competition(height) / area;
    }
  }
  return tot;
}

template <typename T, typename E, typename S>
std::vector<double> Patch<T,E,S>::r_compute_competition_effect_error_by_node_for_species_i(size_t species_index) const {
  const double tot_competition_effect = compute_competition(0.0);
  return species[species_index].r_compute_competition_effect_by_nodes_error(tot_competition_effect);
}

// Integrate over lifetime fitness of individual nodes, scaled per node.
template <typename T, typename E, typename S>
double Patch<T,E,S>::net_reproduction_ratio_for_species(
    size_t species_index, std::vector<double> const& scalars) const {
  auto net_prod = species[species_index].net_reproduction_ratio_by_node_weighted();
  auto const times = species[species_index].node_times();
  auto net_prod_scaled = std::vector<double>(times.size());
  for (size_t i = 0; i < times.size(); ++i) {
    net_prod_scaled[i] = net_prod[i] * scalars[i];
  }
  return util::trapezium(times, net_prod_scaled);
}

// Offspring production, equal to overall fitness scaled by the birth rate.
template <typename T, typename E, typename S>
std::vector<double> Patch<T,E,S>::offspring_production() const {
  auto ret = std::vector<double>(species.size());
  for (size_t i = 0; i < species.size(); ++i) {
    // scale by birth rate function over time
    auto const times = species[i].node_times();
    auto scalars = std::vector<double>(times.size());
    for (size_t j = 0; j < times.size(); ++j) {
      scalars[j] = species[i].extrinsic_drivers().evaluate("birth_rate", times[j]);
    }
    ret[i] = net_reproduction_ratio_for_species(i, scalars);
  }
  return ret;
}

// Overall fitness (no scaling, ie scalars set to 1.0).
template <typename T, typename E, typename S>
std::vector<double> Patch<T,E,S>::net_reproduction_ratios() const {
  auto ret = std::vector<double>(species.size());
  for (size_t i = 0; i < species.size(); ++i) {
    auto scalars = std::vector<double>(species[i].size(), 1.0);
    ret[i] = net_reproduction_ratio_for_species(i, scalars);
  }
  return ret;
}

// Sum up all offspring produced.
template <typename T, typename E, typename S>
double Patch<T,E,S>::total_offspring_production() const {
  double total = 0.0;
  std::vector<double> offspring = offspring_production();
  for (size_t i = 0; i < species.size(); ++i) {
    total += offspring[i];
  }
  return total;
}

// Check integration errors for each species' reproduction integral.
template <typename T, typename E, typename S>
std::vector<std::vector<double>> Patch<T,E,S>::net_reproduction_ratio_errors() const {
  std::vector<std::vector<double>> ret;
  double total_offspring = total_offspring_production();
  for (size_t i = 0; i < species.size(); ++i) {
    ret.push_back(util::local_error_integration(
        species[i].node_times(),
        species[i].net_reproduction_ratio_by_node_weighted(),
        total_offspring));
  }
  return ret;
}

// Sample the competition error for each species introduced this step and fold
// it into the running per-node max (ignoring NA, matching na.rm=TRUE in R).
template <typename T, typename E, typename S>
void Patch<T,E,S>::collect_competition_errors(const std::vector<size_t>& added) {
  for (size_t idx : added) {
    std::vector<double> v =
        r_compute_competition_effect_error_by_node_for_species_i(idx);
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
template <typename T, typename E, typename S>
std::vector<std::vector<double>> Patch<T,E,S>::refinement_error_by_node() const {
  std::vector<std::vector<double>> repro = net_reproduction_ratio_errors();
  std::vector<std::vector<double>> ret(species.size());
  for (size_t i = 0; i < species.size(); ++i) {
    const std::vector<double>& comp = competition_error_by_node[i];
    const std::vector<double>& rep = repro[i];
    const size_t n = species[i].size();
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

// Pre-compute environment, as shaped by residents
// Creates splines of resource availability
template <typename T, typename E, typename S>
void Patch<T,E,S>::compute_environment(bool rescale) {
  
  // Define an anonymous function to use in creation of environment
  auto f = [&](double x) -> double { return compute_competition(x); };

  if (size() > 0 & !is_mutant_run) {
    environment.compute_environment(f, height_max(), rescale);
    // Snapshot the boundary node baked into THIS spline, before compute_rates()
    // mutates it (species 0; the per-RK-stage harvest reads these in cache).
    if (!species.empty()) {
      const auto& nn = species[0].r_new_node();
      newnode_height_env_snapshot = nn.height();
      newnode_competition_env_snapshot = nn.compute_competition(0.0);
    }
    // Per-species boundary-node snapshot for the all-species coupled harvest (R2).
    newnode_height_env_snapshot_all.assign(species.size(), 0.0);
    newnode_competition_env_snapshot_all.assign(species.size(), 0.0);
    for (size_t k = 0; k < species.size(); ++k) {
      const auto& nnk = species[k].r_new_node();
      newnode_height_env_snapshot_all[k] = nnk.height();
      newnode_competition_env_snapshot_all[k] = nnk.compute_competition(0.0);
    }
  }
}


template <typename T, typename E, typename S>
void Patch<T,E,S>::compute_rates() {

  // Computes rates of change for the patch, including all the component species
  // While the patch has an `environment`, the rates here are calculated from
  // the env_ptr, which is a pointer to an environment object
  //  -- for the resident the pointer points to the internal environment object
  //  -- for a mutant, the pointer points to a cached environment object
  double time_ = environment_ptr->time;

  double pr_patch_survival = survival_weighting->pr_survival(time_);
  for (size_t i = 0; i < size(); ++i) {
    double birth_rate = species[i].extrinsic_drivers().evaluate("birth_rate", time_);

    // Pass the environment that pointer is tracking into compute rates.
    species[i].compute_rates(*environment_ptr, pr_patch_survival, birth_rate);
  }

  resource_depletion.reserve(environment_ptr->ode_size());
  for(size_t i = 0; i < environment_ptr->ode_size(); i++) {
    double resource_consumed = std::accumulate(species.begin(), species.end(), 0.0, [i](double r, const species_type& s) {
      return r + s.consumption_rate(i); // accumulates r from zero
    });

    resource_depletion.push_back(resource_consumed/area);
  }
  

  environment_ptr->compute_rates(resource_depletion);

  //todo do we need to clear this every step?
  resource_depletion.clear();

}

// TODO(#478): We should only be recomputing the light environment for the
// points that are below the height of the seedling -- not the entire
// light environment; probably worth just doing a rescale there?
template <typename T, typename E, typename S>
void Patch<T,E,S>::introduce_new_node(size_t species_index) {
  
  species[species_index].introduce_new_node();

  compute_environment(false);
}

template <typename T, typename E, typename S>
void Patch<T,E,S>::introduce_new_nodes(const std::vector<size_t>& species_index) {
  // Record introduction time and patch-age density on each node as it is
  // introduced, so lifetime-fitness calcs need not look these up later.
  const double t = time();
  const double patch_density = survival_weighting->density(t);
  for (size_t i : species_index) {
    species[i].introduce_new_node(t, patch_density);
  }

  compute_environment(false);
}

template <typename T, typename E, typename S>
void Patch<T,E,S>::r_set_time(double time) {
  environment.time = time;
}

// Arguments here are:
//   time: time
//   state: vector of ode state; we'll pass an iterator with that in
//   n: number of *individuals* of each species
template <typename T, typename E, typename S>
void Patch<T,E,S>::r_set_state(double time,
                           const std::vector<double>& state,
                           const std::vector<size_t>& n,
                           const std::vector<double>& light_availability) {
  const size_t n_species = species.size();
  util::check_length(n.size(), n_species);
  reset();
  for (size_t i = 0; i < n_species; ++i) {
    for (size_t j = 0; j < n[i]; ++j) {
      species[i].introduce_new_node();
    }
  }
  util::check_length(state.size(), ode_size());
  set_ode_state(state.begin(), time);
  environment.r_init_interpolators(light_availability);
}

// ODE interface
template <typename T, typename E, typename S>
size_t Patch<T,E,S>::ode_size() const {
  return odelia::ode::ode_size(species.begin(), species.end()) + environment.ode_size();
}

template <typename T, typename E, typename S>
size_t Patch<T,E,S>::aux_size() const {
  // TODO(#478): Is this useful for environment vectors?
  // no use for auxiliary environment variables (yet)
  return odelia::ode::aux_size(species.begin(), species.end());// + environment.ode_size();
}

template <typename T, typename E, typename S>
double Patch<T,E,S>::ode_time() const {
  return time();
}

// First set_ode_state function is for resident runs. Second is for mutant runs
template <typename T, typename E, typename S>
odelia::ode::const_iterator Patch<T,E,S>::set_ode_state(odelia::ode::const_iterator it,
                                              double time) {
  
  // Set ode states
  it = odelia::ode::set_ode_state(species.begin(), species.end(), it);
  it = environment.set_ode_state(it);

  // update time
  environment.time = time;

  // Pre-compute environment, as shaped by residents
  compute_environment(true);
  environment_ptr = &environment;

  // Compute rates of change
  compute_rates();
  return it;
}

// used for mutant runs
// -- differs from above in that an index is passed in as argument
// -- environments are loaded from ODE history, instead of being calculated 
template <typename T, typename E, typename S>
odelia::ode::const_iterator Patch<T,E,S>::set_ode_state(odelia::ode::const_iterator it,
                                              int index) {

  it = odelia::ode::set_ode_state(species.begin(), species.end(), it);

  // using a pointer here to avoid copying environment object
  // just point the pointer, used inside compute rates to get env, to relevant env object
  environment_ptr = &(environment_history[idx][index]);
  environment.time = environment_ptr->time;

  // increment the iterator by an appropriate amount, but don't actually do anything in the env
  for (size_t i = 0; i < environment_ptr->ode_size(); i++) {*it++;}
 
  compute_rates();
  return it;
}

// called from ode_solver->cache
// saves cached set of environments(6) from each ODE step to the step history
template <typename T, typename E, typename S>
void Patch<T,E,S>::cache_ode_step() {
  if(save_RK45_cache) {
    step_history.push_back(time());
    environment_history.push_back(environment_cache);
    // Capture the species-0 stand (heights + per-node competition effect) at this
    // step boundary, for the active-knot resident-light reconstruction.
    if (!species.empty()) {
      stand_height_history.push_back(species[0].r_heights());
      stand_competition_history.push_back(species[0].r_compute_competition_effect_by_nodes());
    }
    // Flush the per-RK-stage stand caches accumulated over this step's 6 stages.
    stand_height_stage_history.push_back(stand_height_stage_cache);
    stand_competition_stage_history.push_back(stand_competition_stage_cache);
    stand_newnode_height_stage_history.push_back(stand_newnode_height_stage_cache);
    stand_newnode_competition_stage_history.push_back(stand_newnode_competition_stage_cache);
    // All-species harvest (R2).
    stand_height_stage_history_all.push_back(stand_height_stage_cache_all);
    stand_competition_stage_history_all.push_back(stand_competition_stage_cache_all);
    stand_newnode_height_stage_history_all.push_back(stand_newnode_height_stage_cache_all);
    stand_newnode_competition_stage_history_all.push_back(stand_newnode_competition_stage_cache_all);
  }
}

// called from ode_step->cache
// saves environment at each RK45 step to the environment cache
template <typename T, typename E, typename S>
void Patch<T,E,S>::cache_RK45_step(int step) {
  if(save_RK45_cache) {
    if(step == 0) {
      environment_cache.clear();
      stand_height_stage_cache.clear();
      stand_competition_stage_cache.clear();
      stand_newnode_height_stage_cache.clear();
      stand_newnode_competition_stage_cache.clear();
      stand_height_stage_cache_all.clear();
      stand_competition_stage_cache_all.clear();
      stand_newnode_height_stage_cache_all.clear();
      stand_newnode_competition_stage_cache_all.clear();
    }
    environment_cache.push_back(environment);
    // All-species stand for THIS stage: [species][cohort] heights + per-node effect,
    // and the per-species boundary node. Mirrors the species-0 capture below.
    {
      std::vector<std::vector<double>> h_all(species.size()), c_all(species.size());
      std::vector<double> nnh_all(species.size(), 0.0), nnc_all(species.size(), 0.0);
      for (size_t k = 0; k < species.size(); ++k) {
        h_all[k] = species[k].r_heights();
        c_all[k] = species[k].r_compute_competition_effect_by_nodes();
        if (k < newnode_height_env_snapshot_all.size()) {
          nnh_all[k] = newnode_height_env_snapshot_all[k];
          nnc_all[k] = newnode_competition_env_snapshot_all[k];
        }
      }
      stand_height_stage_cache_all.push_back(h_all);
      stand_competition_stage_cache_all.push_back(c_all);
      stand_newnode_height_stage_cache_all.push_back(nnh_all);
      stand_newnode_competition_stage_cache_all.push_back(nnc_all);
    }
    // Capture the species-0 stand that produced THIS stage's environment (the ODE
    // state was just set to the stage trial point and compute_environment ran in
    // derivs, immediately before this cache call), aligned 1:1 with environment_cache.
    if (!species.empty()) {
      stand_height_stage_cache.push_back(species[0].r_heights());
      stand_competition_stage_cache.push_back(species[0].r_compute_competition_effect_by_nodes());
      // Boundary node tail term (compute_competition adds it beyond the `nodes` loop),
      // snapshotted in compute_environment before compute_rates mutated it.
      stand_newnode_height_stage_cache.push_back(newnode_height_env_snapshot);
      stand_newnode_competition_stage_cache.push_back(newnode_competition_env_snapshot);
    } else {
      stand_height_stage_cache.emplace_back();
      stand_competition_stage_cache.emplace_back();
      stand_newnode_height_stage_cache.push_back(0.0);
      stand_newnode_competition_stage_cache.push_back(0.0);
    }
  }
}

// called from ode_solver->load, only gets called for mutant runs
template <typename T, typename E, typename S>
void Patch<T,E,S>::load_ode_step() {
  if (use_cached_environment)
  {
    // Minor optimization to check the current and next index before doing a search, as the most common case is that the ODE solver is stepping through the cached environments in order. If the call sequence was not strictly sequential, we fallback to a search through the step history to find the correct environment.

    const double t = time();
    const size_t n = step_history.size();

    // Fast path: step_to() advances through ode_times in order, so this is
    // usually either the current cached step index or the next one.
    if (static_cast<size_t>(idx) < n && util::identical(step_history[idx], t)) {
      return;
    }
    if (static_cast<size_t>(idx + 1) < n &&
        util::identical(step_history[idx + 1], t)) {
      ++idx;
      return;
    }

    // Fallback to search if the call sequence was not strictly sequential.
    auto step = std::find(step_history.begin(), step_history.end(), t);
    if (step == step_history.end()) {
      util::stop("ODE time not found in step history");
    }
    idx = static_cast<int>(std::distance(step_history.begin(), step));
  }
}

template <typename T, typename E, typename S>
odelia::ode::iterator Patch<T,E,S>::ode_state(odelia::ode::iterator it) const {
  it = odelia::ode::ode_state(species.begin(), species.end(), it);
  it = environment.ode_state(it);
  return it;
}

template <typename T, typename E, typename S>
Rcpp::List Patch<T,E,S>::r_get_state() const
{

  // Aseemble commkunity state, icnluding auxiallry variables
  Rcpp::List community_state;
  for (size_t i = 0; i < species.size(); ++i)
  {
    community_state.push_back(species[i].r_get_state());
  }

  return Rcpp::List::create(_["time"] = time(),
                            _["patch_density"] = r_density(time()),
                            _["species"] = community_state,
                            _["env"] = environment.r_get_state());
}

template <typename T, typename E, typename S>
odelia::ode::iterator Patch<T,E,S>::ode_rates(odelia::ode::iterator it) const {
  it = odelia::ode::ode_rates(species.begin(), species.end(), it);
  it = environment.ode_rates(it);
  return it;
}

template <typename T, typename E, typename S>
odelia::ode::iterator Patch<T,E,S>::ode_aux(odelia::ode::iterator it) const {
  it = odelia::ode::ode_aux(species.begin(), species.end(), it);
  return it;
}

}

#endif
