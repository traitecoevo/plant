// -*-c++-*-
#ifndef PLANT_PLANT_PATCH_H_
#define PLANT_PLANT_PATCH_H_

#include <plant/parameters.h>
#include <plant/species.h>
#include <plant/util.h>
#include <plant/adaptive_interpolator.h> // interpolator::refinement_failure
#include <odelia/ode_interface.hpp>
#include <odelia/ode_util.hpp> // odelia::util::stop_domain

#include <plant/disturbance_regime.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility> // std::pair

using namespace Rcpp;

namespace plant {

template <typename T, typename E>
class Patch {
public:
  using value_type = double;

  typedef T                 strategy_type;
  typedef E                 environment_type;
  typedef Individual<T,E>   individual_type;
  typedef Node<T,E>         node_type;
  typedef Species<T,E>      species_type;
  typedef Parameters<T,E>   parameters_type;

  Patch(parameters_type p, environment_type e, plant::Control c);
  void reset();
  size_t size() const {return species.size();}

  //Try using pointer in place of object itself
  double time() const {return environment.time;}
  double get_area() const { return area;}
  double height_max() const;

  double compute_competition(double height) const;

  // Describe the size distribution around `height`, for error messages. The
  // light spline is built from the cohort heights, so when its refinement fails
  // the useful thing to report is what the cohorts are doing there.
  std::string describe_nodes_near(double height) const;

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

  // Apply one non-introduction scheduled event (issue #628). Called between
  // solver legs, so it may change state and ODE size but must not move the
  // clock -- the caller re-reads both with set_state_from_system(). An action
  // is free to reach the answer however it likes, including by integrating its
  // own fast sub-model over the event's nominal duration with the patch's
  // demography frozen; to this solver that is still one instantaneous jump.
  EventRecord apply_event(const NodeScheduleEvent& event);

  // Scale every selected node's density by phi, and report how many nodes were
  // touched. Shared by harvest and climate extremes, which differ only in how phi
  // is chosen per node -- so there is one place where a discrete removal meets
  // the state, and one place to get it right.
  //
  // Both accountings have to move together. On the birth-date coordinate the
  // model maintains log_density(t) = log(birth_rate) - mortality(t), because
  // log_density integrates -mortality_rate from log(birth_rate * pr_estab)
  // while mortality integrates +mortality_rate from -log(pr_estab). Moving
  // only log_density leaves fecundity undiscounted (it is weighted by
  // exp(-mortality) independently); moving only mortality leaves the standing
  // density untouched. Either way the two views of the same cohort disagree.
  // Returns {nodes_affected, density_removed} -- the count is bookkeeping, the
  // density is the quantity actually taken out.
  template <typename Select>
  std::pair<size_t, double> scale_node_densities(size_t species_index,
                                                 Select select);

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
  // Compute rates of change at the state currently loaded and save them into the
  // ode solver. Computing here rather than in set_ode_state is what makes the
  // rates always those of that state, however the caller arrived at it.
  odelia::ode::iterator ode_rates(odelia::ode::iterator it);
  // Retrieve auxillary variables and save into the ode solver
  odelia::ode::iterator ode_aux(odelia::ode::iterator it) const;

  // The strategy's own bound on its states, read by the solver: a step landing
  // on a state it refuses is rejected and retried smaller, rather than
  // committed and clamped by whoever reads it next. Walks the state vector in
  // the order ode_state writes it.
  bool ode_state_valid(const std::vector<double>& y) const;

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
  // This is only here because it wraps a private function.
  void r_compute_environment() {compute_environment(false);}

  // env. cache for assembly
  std::vector<double> step_history{0.0};  // always start at zero
  std::vector<std::vector<environment_type>> environment_history;
  std::vector<environment_type> environment_cache;

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

  // The environment the rates are computed against: the patch's own on a
  // resident run, and on a mutant step the one recorded for that step, which the
  // mutant experiences rather than shapes. Derived from what the patch owns, so
  // it still refers into the patch after the patch is copied.
  environment_type& rate_environment() {
    return use_cached_environment
      ? environment_history.at(idx).at(cached_environment_index)
      : environment;
  }

  // Seed the patch from parameters.initial_state (nodes + birth bookkeeping)
  // when present; called from reset(). Sets environment.time = initial_time.
  void set_initial_state();
  // Guard against initial conditions whose per-node log-density rates are so
  // large they would drive densities to non-finite values within a few steps.
  void check_initial_density_rates() const;
  // Guard against the SCM equations running away mid-integration: a cohort
  // density (exp(log_density)) or an environment state (e.g. TF24 soil water)
  // going non-finite. Called each derivs evaluation, before the non-finite
  // value can propagate into competition, resource uptake, or physiology and
  // surface as an opaque downstream error.
  void check_finite_ode_state() const;
  // Guard the birth-date coordinate's quadrature grid: it is the per-node
  // introduction times, so nodes sharing one give zero-width intervals that drop
  // silently out of the integral. Only reachable for a patch whose nodes were
  // created outside the schedule (seeded or imported without per-node times).
  void check_birth_dates_distinct() const;

  parameters_type parameters;

  double area;
  environment_type environment;
  std::vector<species_type> species;

  //TODO(#476): Move into environment?
  std::vector<double> resource_depletion;

  // Which recorded environment a mutant step is evaluated against; see
  // rate_environment(). Unused on a resident run.
  int cached_environment_index = 0;

  Control control;

  // Per-species running max of the competition error per node, accumulated
  // across the run via collect_competition_errors(). Entries start at -Inf and
  // ignore NA contributions, matching apply(., 2, max, na.rm=TRUE) in R.
  std::vector<std::vector<double>> competition_error_by_node;
};

template <typename T, typename E>
Patch<T,E>::Patch(parameters_type p, environment_type e, Control c)
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

template <typename T, typename E>
void Patch<T,E>::overwrite_strategies(std::vector<strategy_type> strategies) {
  species.clear();
  add_strategies(strategies);
}

template <typename T, typename E>
void Patch<T,E>::add_strategies(std::vector<strategy_type> strategies) {
  for (auto i = 0; i < strategies.size(); ++i) {
		auto s = strategies[i];
    s.control = control; // Overwrite to take the patch control object 
    auto spec = Species<T,E>(s);
    species.push_back(spec);
  }
}

template <typename T, typename E>
void Patch<T,E>::set_mutant() {
    if (environment_history.empty()) {
       util::stop("Run a resident first to generate a competitve landscape");
    }

    is_mutant_run = true;
    save_RK45_cache = false;
    use_cached_environment = true;
  idx = 0;
}

template <typename T, typename E>
void Patch<T,E>::reset() {
   for (auto& s : species) {
    s.clear();
    // allocate variables for tracking resource consumption
    s.resize_consumption_rates(environment.n_resources());
  }

  // resize to species count
  resource_depletion.reserve(environment.n_resources());

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
template <typename T, typename E>
void Patch<T,E>::set_initial_state() {
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

  check_birth_dates_distinct();

  // Build the environment from the real node heights (full recompute, no
  // rescale) and compute rates for the seeded population.
  compute_environment(false);
  compute_rates();
}

template <typename T, typename E>
void Patch<T,E>::check_initial_density_rates() const {
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

template <typename T, typename E>
void Patch<T,E>::check_birth_dates_distinct() const {
  for (size_t i = 0; i < species.size(); ++i) {
    // Only the birth-date coordinate integrates over these times. On the height
    // path they feed the lifetime-fitness integral after the run, where a
    // repeated time has always been tolerated.
    if (!species[i].density_in_birth_date()) {
      continue;
    }
    if (!species[i].birth_dates_are_distinct()) {
      util::stop("Species " + util::to_string(i + 1) + " has nodes sharing an "
                 "introduction time, which the birth-date size-density "
                 "coordinate integrates over: the repeated nodes span zero width "
                 "and drop out of the competition and resource integrals. Supply "
                 "per-node introduction times (parameters$initial_node_times) "
                 "with the initial state, or run with "
                 "control$node_density_in_birth_date = FALSE.");
    }
  }
}

template <typename T, typename E>
void Patch<T,E>::check_finite_ode_state() const {
  // The two failure modes of the same runaway (issue #550), caught here so they
  // fail with an actionable message instead of an opaque downstream one:
  //
  //  (1) A cohort density overflowing to +Inf, which then poisons the
  //      competition integral ("Detected non-finite contribution").
  //  (2) The coupled environment state (TF24 soil water) going non-finite,
  //      because the density-weighted resource uptake diverges as density
  //      grows; the RK stage arithmetic can produce a non-finite soil state a
  //      step before density itself overflows, surfacing as a non-finite soil
  //      potential ("non-finite psi_soil").
  //
  // Density ceiling (~1e43): the same "already unphysical, will drive
  // non-finite values" bar used by check_initial_density_rates. Catching the
  // runaway at this magnitude -- before density reaches a literal +Inf --
  // widens coverage of mode (1) and heads off some instances of mode (2).
  const double log_density_ceiling = 50.0;
  for (size_t i = 0; i < species.size(); ++i) {
    for (auto it = species[i].node_begin(); it != species[i].node_end(); ++it) {
      if (!util::is_finite(it->get_density()) ||
          it->get_log_density() > log_density_ceiling) {
        // The size-density characteristic equation integrates
        //   d(log density)/dt = -d(growth)/d(height) - mortality
        // (see Node::compute_rates). When growth rate falls steeply with size
        // -- e.g. a cohort hitting an extreme drought trough -- the gradient
        // term can spike sharply positive, so log_density integrates upward
        // until density = exp(log_density) overflows to +Inf. The individual's
        // physiology (height, leaf area, per-individual competition) stays
        // finite; it is the cohort *density* that blows up. Smaller ODE steps
        // do not help (the divergence is in the equations, not the stepper), so
        // fail here with an actionable message rather than letting the +Inf
        // propagate into the competition integral or density-weighted resource
        // uptake, where it surfaces as an opaque downstream error.
        util::stop("Non-finite cohort density in the SCM size-density "
                   "(characteristic) equations: species " +
                   util::to_string(i + 1) + " has a node with density=" +
                   util::to_string(it->get_density()) + " (log_density=" +
                   util::to_string(it->get_log_density()) + ", height=" +
                   util::to_string(it->height()) + ") at time=" +
                   util::to_string(environment.time) +
                   ". The density derivative -d(growth)/d(height) - mortality "
                   "can grow without bound when growth rate falls steeply with "
                   "size under rapidly changing or extreme environmental "
                   "forcing, driving density to overflow. Try a shorter "
                   "max_patch_lifetime or less extreme environmental drivers.");
      }
    }
  }

  // Mode (2): a non-finite environment state. `vars.states` is generic across
  // environments (size 0 for light-only environments like FF16, so this loop is
  // a no-op there and cannot false-positive). For TF24 these are the per-depth
  // soil-water states.
  const Internals& env_vars = environment.vars;
  for (size_t i = 0; i < env_vars.state_size; ++i) {
    if (!util::is_finite(env_vars.states[i])) {
      // A rejection, not a failure. #608 measured every observed TF24 soil
      // excursion to be explicit-integrator overshoot -- the step a leg
      // inherits across a discrete change, not a defect in the water balance --
      // and the same case integrates cleanly from a smaller step. So hand this
      // to the stepper to shrink and retry (odelia #55/#56). If the minimum
      // step still lands here, odelia stops and reports this message, so
      // nothing is lost when it really is divergent.
      //
      // Contrast mode (1) above, which stays fatal: that divergence is in the
      // equations rather than the stepper, so shrinking cannot recover it and
      // trying would only burn steps before failing anyway.
      odelia::util::stop_domain("Non-finite environment state (index " + util::to_string(i) +
                 " = " + util::to_string(env_vars.states[i]) + ") at time=" +
                 util::to_string(environment.time) +
                 ". For TF24 this is a soil-water state driven non-finite by the "
                 "density-weighted resource uptake as a cohort density runs away "
                 "in the SCM size-density equations (the same failure that "
                 "otherwise overflows density to +Inf; see #550). Try a shorter "
                 "max_patch_lifetime or less extreme environmental drivers.");
    }
  }
}

template <typename T, typename E>
double Patch<T,E>::height_max() const {
  double ret = 0.0;
  for (size_t i = 0; i < species.size(); ++i) {
    if (!is_mutant_run) {
      ret = std::max(ret, species[i].height_max());
    }
  }
  return ret;
}

template <typename T, typename E>
double Patch<T,E>::compute_competition(double height) const {
  double tot = 0.0;
  for (size_t i = 0; i < species.size(); ++i) {
    if (!is_mutant_run) {
      tot += species[i].compute_competition(height) / area;
    }
  }
  return tot;
}

template <typename T, typename E>
std::vector<double> Patch<T,E>::r_compute_competition_effect_error_by_node_for_species_i(size_t species_index) const {
  const double tot_competition_effect = compute_competition(0.0);
  return species[species_index].r_compute_competition_effect_by_nodes_error(tot_competition_effect);
}

// Integrate over lifetime fitness of individual nodes, scaled per node.
template <typename T, typename E>
double Patch<T,E>::net_reproduction_ratio_for_species(
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
template <typename T, typename E>
std::vector<double> Patch<T,E>::offspring_production() const {
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
template <typename T, typename E>
std::vector<double> Patch<T,E>::net_reproduction_ratios() const {
  auto ret = std::vector<double>(species.size());
  for (size_t i = 0; i < species.size(); ++i) {
    auto scalars = std::vector<double>(species[i].size(), 1.0);
    ret[i] = net_reproduction_ratio_for_species(i, scalars);
  }
  return ret;
}

// Sum up all offspring produced.
template <typename T, typename E>
double Patch<T,E>::total_offspring_production() const {
  double total = 0.0;
  std::vector<double> offspring = offspring_production();
  for (size_t i = 0; i < species.size(); ++i) {
    total += offspring[i];
  }
  return total;
}

// Check integration errors for each species' reproduction integral.
template <typename T, typename E>
std::vector<std::vector<double>> Patch<T,E>::net_reproduction_ratio_errors() const {
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
template <typename T, typename E>
void Patch<T,E>::collect_competition_errors(const std::vector<size_t>& added) {
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
template <typename T, typename E>
std::vector<std::vector<double>> Patch<T,E>::refinement_error_by_node() const {
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
template <typename T, typename E>
void Patch<T,E>::compute_environment(bool rescale) {

  // The boundary node is one end of the birth-date quadrature and its birth date
  // is the current time, so refresh it before the profile is built. Its own
  // stamp is set in compute_rates(), which the stepper calls *after* the
  // set_ode_state() that brings us here, so reading that stamp would use the
  // previous derivs call's time and shorten the boundary interval by a
  // Runge-Kutta stage. Measured effect on FF16 offspring production is below
  // 1e-6 -- the boundary node carries almost no leaf area, so the segment it
  // ends contributes little -- but the interval is then a function of the step
  // size, which the spatial quadrature has no business depending on, and the
  // whole integral *is* that one segment while a species has a single node.
  // No-op for the height coordinate, where this abscissa is the constant
  // initial height.
  for (auto& s : species) {
    s.set_new_node_birth_date(environment.time);
  }

  // Define an anonymous function to use in creation of environment
  auto f = [&](double x) -> double { return compute_competition(x); };

  if (size() > 0 & !is_mutant_run) {
    try {
      environment.compute_environment(f, height_max(), rescale);
    } catch (const interpolator::refinement_failure& e) {
      // The refiner can only say that the profile has a feature it cannot
      // resolve, at some height. We know what is at that height, and it is
      // usually the actual problem: cohorts stacked at one size put a step in
      // the competition profile (#571). Say so rather than making the reader
      // reconstruct the patch state by hand.
      util::stop(std::string(e.what()) + " " +
                 describe_nodes_near(e.report.x_lo));
    }
  }
}

// Report what the size distribution is doing around `height`: how many cohorts
// sit within a narrow window of it, and -- the usual culprit -- whether the node
// list is still ordered by decreasing height. Species::compute_competition() and
// Species::height_max() both take that ordering as given (see the invariant noted
// on both), so once it is violated the competition profile they build has
// fictitious steps in it, and the refiner fails on one of them (#571).
//
// A window rather than a single point because the unresolved interval is ~1e-5 m
// wide, while what matters is whether cohorts share effectively the same size.
template <typename T, typename E>
std::string Patch<T,E>::describe_nodes_near(double height) const {
  const double window = 1e-3; // m

  std::string ret = "Patch state at that height (time " +
                    util::format_double(environment.time) + "):";

  for (size_t i = 0; i < species.size(); ++i) {
    size_t n_near = 0, n_inversions = 0, n_inversions_live = 0, n_zero_density = 0;
    double h_near_min = std::numeric_limits<double>::infinity(),
           h_near_max = -std::numeric_limits<double>::infinity(),
           h_max = -std::numeric_limits<double>::infinity(),
           h_front = NA_REAL, h_prev = NA_REAL;
    bool prev_live = false;

    for (auto it = species[i].node_begin(); it != species[i].node_end(); ++it) {
      const double h = it->height();
      const bool live = it->get_density() > 0.0;
      if (!live) {
        n_zero_density++;
      }
      if (!util::is_finite(h_front)) {
        h_front = h;
      } else if (h > h_prev) {
        n_inversions++;
        // Whether the *live* cohorts are still ordered is the question that
        // matters: the method of characteristics guarantees it for them (growth
        // trajectories sharing an environment cannot cross), so inversions among
        // zero-density nodes are bookkeeping debris in the quadrature grid,
        // whereas an inversion between two live cohorts would mean the
        // characteristics themselves had crossed.
        if (live && prev_live) {
          n_inversions_live++;
        }
      }
      h_prev = h;
      prev_live = live;
      h_max = std::max(h_max, h);
      if (fabs(h - height) <= window) {
        n_near++;
        h_near_min = std::min(h_near_min, h);
        h_near_max = std::max(h_near_max, h);
      }
    }

    ret += " species " + util::to_string(i + 1) + ": " +
           util::to_string(n_near) + " of " +
           util::to_string(species[i].size()) + " cohorts within " +
           util::format_double(window) + " m";
    if (n_near > 0) {
      ret += ", spanning [" + util::format_double(h_near_min) + ", " +
             util::format_double(h_near_max) + "]";
      // Several cohorts at one size means growth has stalled and the schedule is
      // introducing new cohorts into a patch with nothing to separate them.
      if (n_near > 1 && (h_near_max - h_near_min) < window / 100) {
        ret += " -- effectively a single size, so the size distribution has"
               " collapsed here";
      }
    }
    if (n_inversions > 0) {
      ret += "; node heights are NOT decreasing (" +
             util::to_string(n_inversions) + " inversions, " +
             util::to_string(n_inversions_live) +
             " of them between two cohorts of non-zero density; " +
             util::to_string(n_zero_density) + " of " +
             util::to_string(species[i].size()) +
             " nodes have zero density), which breaks the ordering that"
             " Species::compute_competition() and height_max() assume: the"
             " tallest cohort is " + util::format_double(h_max) +
             " but height_max() reports " + util::format_double(h_front) +
             " (the front node), so the competition profile is wrong and its"
             " steps are an artefact rather than a feature of the model";
      if (n_inversions_live == 0) {
        ret += " (the live cohorts are still correctly ordered, so this is the"
               " quadrature grid being scrambled by zero-density nodes rather"
               " than growth trajectories crossing)";
      }
    }
    ret += ";";
  }

  return ret;
}


template <typename T, typename E>
void Patch<T,E>::compute_rates() {

  // Computes rates of change for the patch, including all the component species,
  // against the environment the patch experiences (see rate_environment()).
  environment_type& env = rate_environment();
  double time_ = env.time;

  double pr_patch_survival = survival_weighting->pr_survival(time_);
  for (size_t i = 0; i < size(); ++i) {
    double birth_rate = species[i].extrinsic_drivers().evaluate("birth_rate", time_);

    species[i].compute_rates(env, pr_patch_survival, birth_rate);
  }

  resource_depletion.reserve(env.n_resources());
  for(size_t i = 0; i < env.n_resources(); i++) {
    double resource_consumed = std::accumulate(species.begin(), species.end(), 0.0, [i](double r, const species_type& s) {
      return r + s.consumption_rate(i); // accumulates r from zero
    });

    resource_depletion.push_back(resource_consumed/area);
  }
  

  env.compute_rates(resource_depletion);

  //todo do we need to clear this every step?
  resource_depletion.clear();

}

// TODO(#478): We should only be recomputing the light environment for the
// points that are below the height of the seedling -- not the entire
// light environment; probably worth just doing a rescale there?
template <typename T, typename E>
void Patch<T,E>::introduce_new_node(size_t species_index) {
  
  species[species_index].introduce_new_node();

  compute_environment(false);
}

template <typename T, typename E>
void Patch<T,E>::introduce_new_nodes(const std::vector<size_t>& species_index) {
  // Record introduction time and patch-age density on each node as it is
  // introduced, so lifetime-fitness calcs need not look these up later.
  const double t = time();
  const double patch_density = survival_weighting->density(t);
  for (size_t i : species_index) {
    species[i].introduce_new_node(t, patch_density);
  }

  compute_environment(false);
}

template <typename T, typename E>
template <typename Select>
std::pair<size_t, double>
Patch<T,E>::scale_node_densities(size_t species_index, Select select) {
  species_type& sp = species[species_index];
  size_t n_affected = 0;
  double density_removed = 0.0;
  for (auto n = sp.node_begin(); n != sp.node_end(); ++n) {
    const double phi = select(n->height());
    if (phi >= 1.0) {
      continue;
    }
    if (!(phi > 0.0)) {
      util::stop("An event must leave a positive fraction of each cohort: "
                 "removing all of one would take its density to -Inf, which "
                 "the density transport cannot carry back.");
    }
    const double log_phi = std::log(phi);
    const double before = n->get_density();
    n->set_log_density(n->get_log_density() + log_phi);
    if (util::is_finite(before)) {
      density_removed += before - n->get_density();
    }
    n->individual.set_state("mortality",
                            n->individual.state(MORTALITY_INDEX) - log_phi);
    ++n_affected;
  }
  return {n_affected, density_removed};
}

template <typename T, typename E>
EventRecord Patch<T,E>::apply_event(const NodeScheduleEvent& event) {
  EventRecord rec;
  rec.time = time();
  rec.type = event.type;
  rec.target = event.target;
  rec.target_index = event.target_index;
  rec.requested = event.params;

  // Which species the event reaches: one, or all of them.
  std::vector<size_t> targets;
  if (event.target == EventTarget::Species) {
    targets.push_back(event.target_index);
  } else {
    for (size_t i = 0; i < species.size(); ++i) {
      targets.push_back(i);
    }
  }

  switch (event.type) {
  case EventType::NodeIntroduction:
    // Introductions are batched by the caller so that a run of them costs one
    // environment recompute rather than one each; they never arrive here.
    util::stop("Node introductions are applied via introduce_new_nodes()");
    break;

  case EventType::ResourcePulse: {
    // The environment only: no individual changes, so the competition profile
    // is untouched and the nodes keep their state. The rates are stale
    // afterwards, but the solver recomputes them when it re-reads the system
    // (Patch::ode_rates computes), so there is nothing to do here.
    if (event.target_index >= environment.n_resources()) {
      util::stop("Resource " + util::to_string(event.target_index + 1) +
                 " does not exist: this environment has " +
                 util::to_string(environment.n_resources()) + " resources");
    }
    rec.applied = environment.add_resource_pulse(event.target_index,
                                                 event.params.at(0));
    break;
  }

  case EventType::Harvest: {
    // Remove a fraction of the individuals in a size band. One rule covers the
    // cases that matter: a whole-population removal (the default band), taking
    // everything above a size, and the size-selective removal the restoration
    // work asks for (#627).
    //
    // Individuals are removed rather than shrunk. Changing sizes instead would
    // have to keep them in decreasing order, which Species' competition
    // integral relies on, and the unordered fallback exists only on the size
    // coordinate. Removing part of an individual rather than the whole of it
    // is a real and different thing, and needs that question answered first.
    const double fraction = event.params.at(0);
    const double size_min = event.params.at(1);
    const double size_max = event.params.at(2);
    if (!(fraction >= 0.0) || fraction >= 1.0) {
      util::stop("Harvest fraction must be in [0, 1)");
    }
    const double phi = 1.0 - fraction;
    size_t n_affected = 0;
    double removed = 0.0;
    for (size_t i : targets) {
      const auto r = scale_node_densities(i, [&](double size) {
        return (size >= size_min && size <= size_max) ? phi : 1.0;
      });
      n_affected += r.first;
      removed += r.second;
    }
    rec.applied = {fraction, static_cast<double>(n_affected), removed};
    break;
  }

  case EventType::ClimateExtreme: {
    // An episode of extreme conditions -- heat, cold, salinity, whatever the
    // model's `intensity` means -- that kills in proportion to the dose
    // accumulated above a threshold.
    //
    // This is a worked example of the sub-integration pattern rather than
    // defensible physiology. The event has a nominal duration in the world;
    // the solver's clock does not move across it, so the action steps its own
    // damage variable over that duration under a simple daily cycle and hands
    // back a single survival fraction. The hook is the part worth keeping: a
    // model with a real damage state (#566 for TF24's leaf) plugs in here, and
    // until one exists a dose-dependent mortality is the honest crude stand-in.
    const double intensity = event.params.at(0);
    const double duration = event.params.at(1);
    const double threshold = event.params.at(2);
    const double sensitivity = event.params.at(3);
    if (!util::is_finite(duration) || duration < 0.0) {
      util::stop("Climate extreme duration must be finite and non-negative");
    }

    // Half-hourly sub-steps through a sinusoidal daily cycle peaking at
    // `intensity`, with amplitude set by its excess over `threshold`. Dose
    // accrues only above the threshold, so a mild episode accrues none at all.
    const double dt_days = 0.5 / 24.0;
    const size_t n_steps =
      static_cast<size_t>(std::max(1.0, std::ceil(duration * 365.0 / dt_days)));
    const double amplitude = std::max(0.0, intensity - threshold);
    double damage = 0.0;
    for (size_t k = 0; k < n_steps; ++k) {
      const double day_fraction =
        std::fmod(static_cast<double>(k) * dt_days, 1.0);
      // Peak mid-cycle; the trough sits one amplitude below the threshold.
      const double level =
        threshold + amplitude * (2.0 * std::sin(M_PI * day_fraction) - 1.0);
      const double excess = std::max(0.0, level - threshold);
      damage += sensitivity * excess * dt_days / 365.0;
    }

    const double phi = std::exp(-damage);
    size_t n_affected = 0;
    double removed = 0.0;
    for (size_t i : targets) {
      const auto r = scale_node_densities(i, [&](double) { return phi; });
      n_affected += r.first;
      removed += r.second;
    }
    rec.applied = {1.0 - phi, static_cast<double>(n_affected), removed};
    break;
  }
  }

  // Anything that changed the vegetation changes the light profile too, and
  // every cohort's rates are computed against it.
  if (event.type != EventType::ResourcePulse) {
    compute_environment(false);
  }
  return rec;
}

template <typename T, typename E>
void Patch<T,E>::r_set_time(double time) {
  environment.time = time;
}

// Arguments here are:
//   time: time
//   state: vector of ode state; we'll pass an iterator with that in
//   n: number of *individuals* of each species
template <typename T, typename E>
void Patch<T,E>::r_set_state(double time,
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
  // Every node here is a copy of the boundary node, so they all carry the same
  // birth date and the ODE state does not restore it (it is bookkeeping, not a
  // state variable). Fine for the height coordinate; fatal for the birth-date
  // one, which uses these times as its quadrature grid.
  check_birth_dates_distinct();
  set_ode_state(state.begin(), time);
  environment.r_init_interpolators(light_availability);
}

// ODE interface
template <typename T, typename E>
size_t Patch<T,E>::ode_size() const {
  return odelia::ode::ode_size(species.begin(), species.end()) + environment.ode_size();
}

template <typename T, typename E>
size_t Patch<T,E>::aux_size() const {
  // TODO(#478): Is this useful for environment vectors?
  // no use for auxiliary environment variables (yet)
  return odelia::ode::aux_size(species.begin(), species.end());// + environment.ode_size();
}

template <typename T, typename E>
double Patch<T,E>::ode_time() const {
  return time();
}

// First set_ode_state function is for resident runs. Second is for mutant runs
template <typename T, typename E>
odelia::ode::const_iterator Patch<T,E>::set_ode_state(odelia::ode::const_iterator it,
                                              double time) {
  
  // Set ode states
  it = odelia::ode::set_ode_state(species.begin(), species.end(), it);
  it = environment.set_ode_state(it);

  // update time
  environment.time = time;

  // Catch a runaway size-density equation (non-finite cohort density or
  // environment state) before it feeds into competition, resource uptake, or
  // physiology below and surfaces as an opaque downstream error (issue #550).
  check_finite_ode_state();

  // Pre-compute environment, as shaped by residents
  compute_environment(true);

  return it;
}

// used for mutant runs
// -- differs from above in that an index is passed in as argument
// -- environments are loaded from ODE history, instead of being calculated 
template <typename T, typename E>
odelia::ode::const_iterator Patch<T,E>::set_ode_state(odelia::ode::const_iterator it,
                                              int index) {

  it = odelia::ode::set_ode_state(species.begin(), species.end(), it);

  // Record which of this step's cached environments the rates are evaluated
  // against; rate_environment() reads it back without copying the object.
  cached_environment_index = index;
  const environment_type& cached = rate_environment();
  environment.time = cached.time;

  // increment the iterator by an appropriate amount, but don't actually do anything in the env
  for (size_t i = 0; i < cached.ode_size(); i++) {*it++;}

  return it;
}

// called from ode_solver->cache
// saves cached set of environments(6) from each ODE step to the step history
template <typename T, typename E>
void Patch<T,E>::cache_ode_step() {
  if(save_RK45_cache) { 
    step_history.push_back(time());
    environment_history.push_back(environment_cache);
  }
}

// called from ode_step->cache
// saves environment at each RK45 step to the environment cache
template <typename T, typename E>
void Patch<T,E>::cache_RK45_step(int step) {
  if(save_RK45_cache) {  
    if(step == 0) {
      environment_cache.clear();
    }
    environment_cache.push_back(environment);
  }
}

// called from ode_solver->load, only gets called for mutant runs
template <typename T, typename E>
void Patch<T,E>::load_ode_step() {
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

template <typename T, typename E>
odelia::ode::iterator Patch<T,E>::ode_state(odelia::ode::iterator it) const {
  it = odelia::ode::ode_state(species.begin(), species.end(), it);
  it = environment.ode_state(it);
  return it;
}

template <typename T, typename E>
Rcpp::List Patch<T, E>::r_get_state() const
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

template <typename T, typename E>
odelia::ode::iterator Patch<T,E>::ode_rates(odelia::ode::iterator it) {
  compute_rates();
  it = odelia::ode::ode_rates(species.begin(), species.end(), it);
  it = environment.ode_rates(it);
  return it;
}

template <typename T, typename E>
odelia::ode::iterator Patch<T,E>::ode_aux(odelia::ode::iterator it) const {
  it = odelia::ode::ode_aux(species.begin(), species.end(), it);
  return it;
}

template <typename T, typename E>
bool Patch<T,E>::ode_state_valid(const std::vector<double>& y) const {
  // The environment block first, and separately: it is the trailing part of the
  // ODE vector, it is where integrator overshoot shows up (TF24's soil water,
  // driven past its bounds by a step inherited across a discrete change), and
  // the completed step's state never reaches check_finite_ode_state(), which
  // runs per stage. Finiteness only -- a node's log_density is legitimately
  // -Inf, so this must not be extended over the species block.
  const size_t n_env = environment.ode_size();
  if (n_env > 0 && y.size() >= n_env) {
    for (size_t i = y.size() - n_env; i < y.size(); ++i) {
      if (!util::is_finite(y[i])) {
        return false;
      }
    }
  }

  // Resolved on the concrete strategy, so a model that appends or inserts a
  // state gets its own positions rather than its base's. Named rather than
  // indexed because a strategy that declared the index would be declaring
  // something its own state_names() already says.
  static const std::vector<size_t> bounded = [] {
    const std::vector<std::string> names = T::state_names();
    std::vector<size_t> ret;
    for (const std::string& name : T::non_negative_states()) {
      const auto found = std::find(names.begin(), names.end(), name);
      if (found == names.end()) {
        util::stop("Strategy declares '" + name + "' as a non-negative state "
                   "and does not carry it");
      }
      ret.push_back(static_cast<size_t>(found - names.begin()));
    }
    return ret;
  }();
  if (bounded.empty()) {
    return true;
  }
  const size_t stride = node_type::ode_size();
  size_t at = 0;
  for (const auto& sp : species) {
    for (size_t k = 0; k < sp.size(); ++k, at += stride) {
      for (const size_t i : bounded) {
        if (y[at + i] < 0.0) {
          return false;
        }
      }
    }
  }
  return true;
}

}

#endif
