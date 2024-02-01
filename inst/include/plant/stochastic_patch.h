// -*-c++-*-
#ifndef PLANT_PLANT_STOCHASTIC_PATCH_H_
#define PLANT_PLANT_STOCHASTIC_PATCH_H_

namespace plant {

// NOTE: compute_environment() here might fail (especially for
// rare seed arrivals) because the adaptive refinement can't deal with
// the sharp corners that are implied.  The simplest thing to do is to
// tone down the tolerance (fast_control() seems good enough) but that
// might not be enough.  It might be best to do something more clever
// than just fail on refining, but doing that will require
// confirmation that the issue is simply in a couple of places rather
// than throughout.  Running the spline piecewise would be the best
// bet there.
template <typename T, typename E>
class StochasticPatch {
public:
  typedef T                      strategy_type;
  typedef E                      environment_type;
  typedef Individual<T,E>        individual_type;
  typedef StochasticSpecies<T,E> species_type;
  typedef Parameters<T,E>        parameters_type;
  StochasticPatch(parameters_type p, environment_type e, Control c);
  void reset();

  size_t size() const {return species.size();}
  double time() const {return environment.time;}

  double height_max() const;

  // [eqn 11] Canopy openness at `height`
  double compute_competition(double height) const;

  bool introduce_new_node(size_t species_index);
  void introduce_new_node_and_update(size_t species_index);

  std::vector<size_t> deaths();

  const species_type& at(size_t species_index) const {
    return species[species_index];
  }

  // * ODE interface
  size_t ode_size() const;
  double ode_time() const;
  ode::const_iterator set_ode_state(ode::const_iterator it, double time);
  ode::iterator       ode_state(ode::iterator it) const;
  ode::iterator       ode_rates(ode::iterator it) const;

  // * R interface
  // Data accessors:
  parameters_type r_parameters() const {return parameters;}
  E r_environment() const {return environment;}
  std::vector<species_type> r_species() const {return species;}
  void r_set_state(double time,
                   const std::vector<double>& state,
                   const std::vector<size_t>& n);
  // TODO: No support here for setting *vectors* of species.  Might
  // want to supoprt that?
  bool r_introduce_new_node(util::index species_index) {
    return introduce_new_node(species_index.check_bounds(size()));
  }
  void r_introduce_new_node_and_update(util::index species_index) {
    introduce_new_node_and_update(species_index.check_bounds(size()));
  }

  species_type r_at(util::index species_index) const {
    at(species_index.check_bounds(size()));
  }
  // These are only here because they wrap private functions.
  void r_compute_environment() {compute_environment(false);}
  void r_compute_rates() {compute_rates();}
private:
  void compute_environment(bool rescale);
  void compute_rates();

  parameters_type parameters;

  E environment;
  std::vector<species_type> species;
  Control control;
};

template <typename T, typename E>
StochasticPatch<T,E>::StochasticPatch(parameters_type p, environment_type e, Control c)
  : parameters(p),

    environment(e),
    control(c) {
  parameters.validate();
  for (auto s : parameters.strategies) {
    s.control = control;
    species.push_back(species_type(s));
  }
  reset();
}

template <typename T, typename E>
void StochasticPatch<T,E>::reset() {
  for (auto& s : species) {
    s.clear();
  }
  environment.clear();
  compute_environment(false);
  compute_rates();
}

template <typename T, typename E>
double StochasticPatch<T,E>::height_max() const {
  double ret = 0.0;
  for (size_t i = 0; i < species.size(); ++i) {
      ret = std::max(ret, species[i].height_max());
  }
  return ret;
}

template <typename T, typename E>
double StochasticPatch<T,E>::compute_competition(double height) const {
  double tot = 0.0;
  for (size_t i = 0; i < species.size(); ++i) {
      tot += species[i].compute_competition(height);
  }
  return tot;
}

template <typename T, typename E>
void StochasticPatch<T,E>::compute_environment(bool rescale) {
  if (height_max() > 0.0) {
    auto f = [&] (double x) -> double {return compute_competition(x);};
    environment.compute_environment(f, height_max(), rescale);
  } else {
    environment.clear_environment();
  }
}


template <typename T, typename E>
void StochasticPatch<T,E>::compute_rates() {
  for (size_t i = 0; i < size(); ++i) {
    species[i].compute_rates(environment);
  }
}

// In theory, this could be done more efficiently by, in the introdudce_new_node
// case, using the values stored in the species offspring.  But we don't
// really get that here.  It might be better to move introdudce_new_node /
// introduce_new_node_and_update within Species, given this.
template <typename T, typename E>
void StochasticPatch<T,E>::introduce_new_node_and_update(size_t species_index) {
  // Add a offspring, setting ODE variables based on the *current* light environment
  species[species_index].introduce_new_node(environment);
  // Then we update the light environment.
  compute_environment(false);
}

template <typename T, typename E>
bool StochasticPatch<T,E>::introduce_new_node(size_t species_index) {
  const double pr_germinate =
    species[species_index].establishment_probability(environment);
  const bool added = unif_rand() < pr_germinate;
  if (added) {
    introduce_new_node_and_update(species_index);
  }
  return added;
}

template <typename T, typename E>
std::vector<size_t> StochasticPatch<T,E>::deaths() {
  std::vector<size_t> ret;
  ret.reserve(size());
  bool recompute = false;
  for (auto& s : species) {
    const size_t n_deaths = s.deaths();
    ret.push_back(n_deaths);
    recompute = recompute || n_deaths > 0;
  }
  if (recompute) {
    compute_environment(false);
    compute_rates();
  }
  return ret;
}

// Arguments here are:
//   time: time
//   state: vector of ode state; we'll pass an iterator with that in
//   n: number of *individuals* of each species
template <typename T, typename E>
void StochasticPatch<T,E>::r_set_state(double time,
                           const std::vector<double>& state,
                           const std::vector<size_t>& n) {
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
}

// ODE interface
template <typename T, typename E>
size_t StochasticPatch<T,E>::ode_size() const {
  return ode::ode_size(species.begin(), species.end());
}

template <typename T, typename E>
double StochasticPatch<T,E>::ode_time() const {
  return time();
}

template <typename T, typename E>
ode::const_iterator StochasticPatch<T,E>::set_ode_state(ode::const_iterator it,
                                                      double time) {
  
  // set ode sates
  it = ode::set_ode_state(species.begin(), species.end(), it);
  environment.time = time;

  // pre-compute resources avaialability and competion, as defined by residents
  compute_environment(true);

  // compute rates of changes
  compute_rates();
  return it;
}

template <typename T, typename E>
ode::iterator StochasticPatch<T,E>::ode_state(ode::iterator it) const {
  return ode::ode_state(species.begin(), species.end(), it);
}

template <typename T, typename E>
ode::iterator StochasticPatch<T,E>::ode_rates(ode::iterator it) const {
  return ode::ode_rates(species.begin(), species.end(), it);
}

}

#endif
