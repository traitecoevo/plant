// -*-c++-*-
#ifndef PLANT_PLANT_STOCHASTIC_PATCH_H_
#define PLANT_PLANT_STOCHASTIC_PATCH_H_

#include <plant/parameters.h>
#include <plant/stochastic_species.h>
#include <plant/util.h>
#include <plant/with_slope.h>
#include <numeric> // std::accumulate, in compute_rates

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
  using value_type = typename T::value_type;

  typedef T                      strategy_type;
  typedef E                      environment_type;
  typedef Individual<T,E>        individual_type;
  typedef StochasticSpecies<T,E> species_type;
  typedef Parameters<T,E>        parameters_type;
  StochasticPatch(parameters_type p, environment_type e, Control c);
  void reset();

  size_t size() const {return species.size();}
  double time() const {return environment.time;}
  double get_area() const { return area;}

  value_type height_max() const;

  // [eqn 11] Canopy openness at `height`
  value_type compute_competition(double height) const;
  // That profile and its vertical derivative, from one pass over the species.
  with_slope<value_type>
  compute_competition_and_slope(double height) const;

  bool introduce_new_node(size_t species_index);
  void introduce_new_node_and_update(size_t species_index);

  std::vector<size_t> deaths();

  const species_type& at_species(size_t species_index) const {
    return species[species_index];
  }

  // * ODE interface
  size_t ode_size() const;
  // The individuals' share of the ODE system, i.e. zero on an empty patch
  // whatever state the environment carries.
  size_t node_ode_size() const { return ode_size() - environment.ode_size(); }
  double ode_time() const;
  double area;

  template <typename It> It set_ode_state(It it, double time);
  template <typename It> It ode_state(It it) const;
  template <typename It> It ode_rates(It it) const;

  // * R interface
  // Data accessors:
  parameters_type r_parameters() const {return parameters;}
  E r_environment() const {return environment;}
  std::vector<species_type> r_species() const {return species;}
  void r_set_state(double time,
                   const std::vector<double>& state,
                   const std::vector<size_t>& n);
  Rcpp::List r_get_state() const;
  // TODO(#479): No support here for setting *vectors* of species.  Might
  // want to supoprt that?
  bool r_introduce_new_node(util::index species_index) {
    return introduce_new_node(species_index.check_bounds(size()));
  }
  void r_introduce_new_node_and_update(util::index species_index) {
    introduce_new_node_and_update(species_index.check_bounds(size()));
  }

  species_type r_at(util::index species_index) const {
    return species[species_index.check_bounds(size())];
  }
  // These are only here because they wrap private functions.
  void r_compute_environment() {compute_environment();}
  void r_compute_rates() {compute_rates();}
private:
  void compute_environment();
  void compute_rates();

  parameters_type parameters;

  E environment;
  std::vector<species_type> species;
  Control control;

  // Scratch for compute_rates(), which the solver calls on every derivatives
  // evaluation. A member reserved once in reset() and cleared (not freed) after
  // use, so the RHS does not allocate; mirrors Patch.
  std::vector<double> resource_depletion;
};

template <typename T, typename E>
StochasticPatch<T,E>::StochasticPatch(parameters_type p, environment_type e, Control c)
  : parameters(p),
    area(p.patch_area),
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
  resource_depletion.reserve(environment.n_resources());
  environment.clear();
  compute_environment();
  compute_rates();
}

template <typename T, typename E>
typename StochasticPatch<T,E>::value_type
StochasticPatch<T,E>::height_max() const {
  value_type ret = 0.0;
  for (size_t i = 0; i < species.size(); ++i) {
      const value_type h = species[i].height_max();
      if (h > ret) {
        ret = h;
      }
  }
  return ret;
}

template <typename T, typename E>
typename StochasticPatch<T,E>::value_type
StochasticPatch<T,E>::compute_competition(double height) const {
  value_type tot = 0.0;
  for (size_t i = 0; i < species.size(); ++i) {
    tot += species[i].compute_competition(height) / area;
  }
  return tot;
}

template <typename T, typename E>
with_slope<typename StochasticPatch<T,E>::value_type>
StochasticPatch<T,E>::compute_competition_and_slope(double height) const {
  with_slope<value_type> sum{0.0, 0.0};
  for (size_t i = 0; i < species.size(); ++i) {
    const with_slope<value_type> fs =
      species[i].compute_competition_and_slope(height);
    sum.value += fs.value / area;
    sum.slope += fs.slope / area;
  }
  return sum;
}

template <typename T, typename E>
void StochasticPatch<T,E>::compute_environment() {
  if (height_max() > 0.0) {
    // Written as std::vector<double> this still compiles, taking the value of an
    // active profile, and the field's knot values and slopes would then be
    // constants with nothing raised to say so.
    //
    // Individuals rather than nodes here, so there is no trapezium over a grid and
    // nothing to accumulate across knots: the per-knot sum IS the reduction, and
    // this walks the knots itself.
    auto f = [&] (const std::vector<double>& x, std::vector<value_type>& y,
                  std::vector<value_type>& m) -> void {
      for (size_t k = 0; k < x.size(); ++k) {
        const with_slope<value_type> fs = compute_competition_and_slope(x[k]);
        y[k] = fs.value;
        m[k] = fs.slope;
      }
    };
    environment.compute_environment(f, height_max());
  } else {
    environment.clear_environment();
  }
}


template <typename T, typename E>
void StochasticPatch<T,E>::compute_rates() {
  for (size_t i = 0; i < size(); ++i) {
    species[i].compute_rates(environment);
  }

  // Resource uptake per unit area, which closes the environment's own balance
  // (for TF24, soil water). Sized by the number of resources the environment
  // consumes, which is not its ODE width: TF24 has five soil layers plus four
  // cumulative-flux states that nothing draws on.
  for (size_t i = 0; i < environment.n_resources(); i++) {
    double resource_consumed = std::accumulate(
        species.begin(), species.end(), 0.0,
        [i](double r, const species_type& s) {return r + s.consumption_rate(i);});
    resource_depletion.push_back(resource_consumed / area);
  }

  environment.compute_rates(resource_depletion);
  resource_depletion.clear();
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
  compute_environment();
}

template <typename T, typename E>
bool StochasticPatch<T,E>::introduce_new_node(size_t species_index) {
  // The draw decides whether an individual is added, so the probability is read
  // at its value here.
  const double pr_germinate = odelia::util::to_passive(
    species[species_index].establishment_probability(environment));
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
    compute_environment();
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

template <typename T, typename E>
Rcpp::List StochasticPatch<T, E>::r_get_state() const
{

  // Aseemble commkunity state, icnluding auxiallry variables
  Rcpp::List community_state;
  for (size_t i = 0; i < species.size(); ++i)
  {
    community_state.push_back(species[i].r_get_state());
  }

  return Rcpp::List::create(_["time"] = time(),
                            _["species"] = community_state,
                            _["env"] = environment.r_get_state());
}

// ODE interface
template <typename T, typename E>
size_t StochasticPatch<T,E>::ode_size() const {
  return odelia::ode::ode_size(species.begin(), species.end()) + environment.ode_size();
}

template <typename T, typename E>
double StochasticPatch<T,E>::ode_time() const {
  return time();
}

template <typename T, typename E>
template <typename It>
It StochasticPatch<T,E>::set_ode_state(It it, double time) {
  
  // set ode sates
  it = odelia::ode::set_ode_state(species.begin(), species.end(), it);
  it = environment.set_ode_state(it);
  environment.time = time;

  // pre-compute resources avaialability and competion, as defined by residents
  compute_environment();

  // compute rates of changes
  compute_rates();
  return it;
}

template <typename T, typename E>
template <typename It>
It StochasticPatch<T,E>::ode_state(It it) const {
  it = odelia::ode::ode_state(species.begin(), species.end(), it);
  return environment.ode_state(it);
}

template <typename T, typename E>
template <typename It>
It StochasticPatch<T,E>::ode_rates(It it) const {
  it = odelia::ode::ode_rates(species.begin(), species.end(), it);
  return environment.ode_rates(it);
}

}

#endif
