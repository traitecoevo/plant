// -*-c++-*-
#ifndef PLANT_PLANT_PATCH_H_
#define PLANT_PLANT_PATCH_H_

#include <memory> // std::shared_ptr

#include <plant/parameters.h>
#include <plant/species.h>
#include <plant/ode_interface.h>

#include <plant/disturbance_regime.h>

using namespace Rcpp;

namespace plant {

template <typename T, typename E>
class Patch {
public:
  typedef T                 strategy_type;
  typedef E                 environment_type;
  typedef Individual<T,E>   individual_type;
  typedef Node<T,E>         node_type;
  typedef Species<T,E>      species_type;
  typedef Parameters<T,E>   parameters_type;

  // typedef std::shared_ptr<environment_type> environment_type_ptr;

  Patch(parameters_type p, environment_type e, plant::Control c);
  void reset();
  size_t size() const {return species.size();}
  //Try using pointer in place of object itself
  double time() const {return environment.time;}

  double height_max() const;

  // [eqn 11] Canopy openness at `height`
  double compute_competition(double height) const;

  void introduce_new_node(size_t species_index);
  void introduce_new_nodes(const std::vector<size_t>& species_index);

  // Open to better ways to test whether nodes have been introduced
  int node_ode_size() const {
    int node_ode_size = ode_size() - environment.ode_size();
    return(node_ode_size);
  }

  const species_type& at(size_t species_index) const {
    return species[species_index];
  }

  // Patch disturbance
  Disturbance_Regime* survival_weighting;

  // * ODE interface
  size_t ode_size() const;
  size_t aux_size() const;
  double ode_time() const;
  // set_ode_state is the main interface for stepping the patch as an ode system
  // There are two implementations.
  // First set_ode_state function is for resident runs.
  // Second is for mutant runs
  // the decision which to use is based on whether we are using the stored cache of
  // environment (determined by `use_cached_environment` in ode_interface.h)
  ode::const_iterator set_ode_state(ode::const_iterator it, double time);
  ode::const_iterator set_ode_state(ode::const_iterator it, int index);
  ode::iterator       ode_state(ode::iterator it) const;
  ode::iterator       ode_rates(ode::iterator it) const;
  ode::iterator       ode_aux(ode::iterator it) const;

  // * R interface
  // Data accessors:

  std::vector<double> r_density(std::vector<double> time) const {return survival_weighting->r_density(time);}
  double r_pr_survival(double time) const {return survival_weighting->pr_survival(time);}
  double r_disturbance_mean_interval() const {return survival_weighting->r_mean_interval();}
  double r_survival_weighting_cdf(double time) const {return survival_weighting->cdf(time);}
  double r_survival_weighting_icdf(double prob) const {return survival_weighting->icdf(prob);}

  parameters_type r_parameters() const {return parameters;}
  environment_type r_environment() const {return environment;}
  std::vector<species_type> r_species() const {return species;}
  std::vector<double> r_competition_effect_error(size_t species_index) const;
  void r_set_time(double time);
  void r_set_state(double time,
                   const std::vector<double>& state,
                   const std::vector<size_t>& n,
                   const std::vector<double>& canopy);
  void r_introduce_new_node(util::index species_index) {
    introduce_new_node(species_index.check_bounds(size()));
  }
  species_type r_at(util::index species_index) const {
    at(species_index.check_bounds(size()));
  }
  // These are only here because they wrap private functions.
  void r_compute_environment() {compute_environment();}
  void r_compute_rates() {compute_rates();}

  // Prototype env. cache for assembly
  std::vector<double> step_history{0.0};  // always start at zero
  std::vector<std::vector<environment_type>> environment_history;

  void cache_ode_step();
  void cache_RK45_step(int step);
  void load_ode_step();
  std::vector<environment_type> environment_cache;

  // use cache for mutant runs
  // todo: document more -- what variable does save_RK45_cache relate to in ode_interface
  bool save_RK45_cache;
  bool use_cached_environment = false;

  bool is_mutant_run = false;

  void set_mutant();
  void add_strategies(std::vector<strategy_type> strategies);
  void overwrite_strategies(std::vector<strategy_type> strategies);

private:
  int idx;
  void compute_environment();
  void rescale_environment();
  void compute_rates();

  parameters_type parameters;
  environment_type environment;
  std::vector<species_type> species;
  std::vector<double> resource_depletion;

  environment_type* env_ptr1;
  std::shared_ptr<environment_type> env_ptr2;

  Control control;
};

template <typename T, typename E>
Patch<T,E>::Patch(parameters_type p, environment_type e, Control c)
  : parameters(p),
    environment(e),
// Ideally we'd set pointer here, code below runs but the connection between pointer and env get's lost somewhere
//   env_ptr1(&environment),
    control(c),
    environment_cache(6) {  // length of ode::Step
  parameters.validate();

  // todo: document this.
  // save_RK45_cache is not in ode interface
  save_RK45_cache = control.save_RK45_cache;
  survival_weighting = p.disturbance;

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
    is_mutant_run = true;
    save_RK45_cache = false;
    use_cached_environment = true;

    // check cache is valid
    // if cache = invalid
    //  util::stop("Run a resident first to generate a competitve landscape")
}

template <typename T, typename E>
void Patch<T,E>::reset() {
 
    for (auto& s : species)
    {
    s.clear();

    // allocate variables for tracking resource consumption
    s.resize_consumption_rates(environment.ode_size());
  }

  // resize to species count
  resource_depletion.reserve(environment.ode_size());

  // compute ephemeral effects like canopy
  environment.clear();
  compute_environment();

  // compute effects of resource consumption
  compute_rates();
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
      tot += species[i].compute_competition(height);
    }
  }
  return tot;
}

template <typename T, typename E>
std::vector<double> Patch<T,E>::r_competition_effect_error(size_t species_index) const {
  const double tot_competition_effect = compute_competition(0.0);
  return species[species_index].r_competition_effects_error(tot_competition_effect);
}

template <typename T, typename E>
void Patch<T,E>::compute_environment() {
  if (size() > 0)
  {
    if (!is_mutant_run)
    {
      auto f = [&] (double x) -> double
      {return compute_competition(x);};
      environment.compute_environment(f, height_max());
    }
  }
}

template <typename T, typename E>
void Patch<T,E>::rescale_environment() {
  if (size() > 0) {
    if (!is_mutant_run ) {
      auto f = [&] (double x) -> double
      {return compute_competition(x);};
    environment.rescale_environment(f, height_max());
    }
  }
}

template <typename T, typename E>
void Patch<T,E>::compute_rates() {

  // Setting the pointer here, as a proof of concept that pointers work
  // env_ptr1 and environment should access the same data
  // Below we use the pointer when computing rates, this works
  // next need to figure out how to redirect to correct env object
  env_ptr1 = &environment;

  double pr_patch_survival = survival_weighting->pr_survival(time());

//  std::cout << "E " << environment.time << " " << env_ptr1->time << std::endl;

  for (size_t i = 0; i < size(); ++i) {
    double pr_patch_survival = survival_weighting->pr_survival(time());
		double birth_rate = species[i].extrinsic_drivers().evaluate("birth_rate", time());

   // species[i].compute_rates(environment, pr_patch_survival, birth_rate);
   //Pass pointer into compute rates. This works
    species[i].compute_rates(*env_ptr1, pr_patch_survival, birth_rate);
  }

  resource_depletion.reserve(env_ptr1->ode_size());
  for(size_t i = 0; i < env_ptr1->ode_size(); i++) {
    double resource_consumed = std::accumulate(species.begin(), species.end(), 0.0, [i](double r, const species_type& s) {
      return r + s.consumption_rate(i); // accumulates r from zero
    });

    resource_depletion.push_back(resource_consumed);
  }

  env_ptr1->compute_rates(resource_depletion);
  //todo do we need to clear this everyt step?
  resource_depletion.clear();
}

// TODO: We should only be recomputing the light environment for the
// points that are below the height of the seedling -- not the entire
// light environment; probably worth just doing a rescale there?
template <typename T, typename E>
void Patch<T,E>::introduce_new_node(size_t species_index) {
  species[species_index].introduce_new_node();
  if (!is_mutant_run) {
    compute_environment();
  }
}

template <typename T, typename E>
void Patch<T,E>::introduce_new_nodes(const std::vector<size_t>& species_index) {
  // std::cout << "introducing nodes for species: ";
  for (size_t i : species_index) {
    // std::cout << i << " ";
    species[i].introduce_new_node();
  }
  // std::cout << std::endl;
  if (!is_mutant_run) {
    compute_environment();
  }
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
                           const std::vector<double>& canopy) {
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
  environment.r_init_interpolators(canopy);
}

// ODE interface
template <typename T, typename E>
size_t Patch<T,E>::ode_size() const {
  return ode::ode_size(species.begin(), species.end()) + environment.ode_size();
}

template <typename T, typename E>
size_t Patch<T,E>::aux_size() const {
  // no use for auxiliary environment variables (yet)
  return ode::aux_size(species.begin(), species.end());// + environment.ode_size();
}

template <typename T, typename E>
double Patch<T,E>::ode_time() const {
  return time();
}

// First set_ode_state function is for resident runs. Second is for mutant runs
template <typename T, typename E>
ode::const_iterator Patch<T,E>::set_ode_state(ode::const_iterator it,
                                              double time) {
  
  it = ode::set_ode_state(species.begin(), species.end(), it);
  it = environment.set_ode_state(it);

  environment.time = time;
  // std::cout << "resident: etime " << environment.time << std::endl;
  
  if (environment.canopy_rescale_usually) {
    rescale_environment();
  } else {
    compute_environment();
  }
  compute_rates();
  return it;
}

// pre-cached environments, used for mutant runs
// differs from above in that an index is passed in as argument
template <typename T, typename E>
ode::const_iterator Patch<T,E>::set_ode_state(ode::const_iterator it,
                                              int index) {

  it = ode::set_ode_state(species.begin(), species.end(), it);

  // todo: use pointer here
  environment = environment_history[idx][index];

  // todo: Possibly should skip next step. We don't want to write to env
  // Instead could increment the iterator by an appropriate amount
  // it = environment.set_ode_state(it);
  for (size_t i = 0; i < environment.ode_size(); i++) {*it++;}
 
  // std::cout << "mutant: etime " << environment.time << "; index " << index << std::endl;

  compute_rates();
  return it;
}

template <typename T, typename E>
void Patch<T,E>::cache_ode_step() {
  if(save_RK45_cache) { 
    step_history.push_back(time());
    environment_history.push_back(environment_cache);
    // environment_cache.clear()
  }
}

// This function used in ode_solver.h
// only gets called for mutant runs
template <typename T, typename E>
void Patch<T,E>::load_ode_step() {
  if(use_cached_environment) { 
    std::vector<double>::iterator step;
    
    // find where we are in the cache
    step = std::find(step_history.begin(), step_history.end(), time());

    if(*step != time()) {
      util::stop("ODE time not found in step history");
    }

    // index into the right environment
    // todo: pointer here
    // loads a set of 6 environments for the current step
    idx = std::distance(step_history.begin(), step);
   // environment_cache.clear();
   // environment_cache = environment_history[idx];
  }
}

template <typename T, typename E>
void Patch<T,E>::cache_RK45_step(int step) {
  if(save_RK45_cache) {  
    if(step == 0) {
      environment_cache.clear();
    }
    environment_cache.push_back(environment);
  }
}

template <typename T, typename E>
ode::iterator Patch<T,E>::ode_state(ode::iterator it) const {
  it = ode::ode_state(species.begin(), species.end(), it);
  it = environment.ode_state(it);
  return it;
}

template <typename T, typename E>
ode::iterator Patch<T,E>::ode_rates(ode::iterator it) const {
  it = ode::ode_rates(species.begin(), species.end(), it);
  it = environment.ode_rates(it);
  return it;
}

template <typename T, typename E>
ode::iterator Patch<T,E>::ode_aux(ode::iterator it) const {
  it = ode::ode_aux(species.begin(), species.end(), it);
  //it = environment.ode_rates(it);
  return it;
}

}

#endif
