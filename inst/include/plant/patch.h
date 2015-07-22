// -*-c++-*-
#ifndef PLANT_PLANT_PATCH_H_
#define PLANT_PLANT_PATCH_H_

#include <plant/parameters.h>
#include <plant/species.h>
#include <plant/ode_interface.h>

namespace plant {

template <typename T>
class Patch {
public:
  typedef T             strategy_type;
  typedef Plant<T>      plant_type;
  typedef Cohort<T>     cohort_type;
  typedef Species<T>    species_type;
  typedef Parameters<T> parameters_type;

  Patch(parameters_type p);

  void reset();
  size_t size() const {return species.size();}
  double time() const {return environment.time;}

  double height_max() const;

  // [eqn 11] Canopy openness at `height`
  double area_leaf_above(double height) const;
  double canopy_openness(double height) const;

  void add_seed(size_t species_index);
  void add_seeds(const std::vector<size_t>& species_index);

  const species_type& at(size_t species_index) const {
    return species[species_index];
  }
  const Disturbance& disturbance_regime() const {
    return environment.disturbance_regime;
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
  Environment r_environment() const {return environment;}
  std::vector<species_type> r_species() const {return species;}
  std::vector<double> r_area_leaf_error(size_t species_index) const;
  void r_set_state(double time,
                   const std::vector<double>& state,
                   const std::vector<size_t>& n);
  void r_add_seed(util::index species_index) {
    add_seed(species_index.check_bounds(size()));
  }
  species_type r_at(util::index species_index) const {
    at(species_index.check_bounds(size()));
  }
  // These are only here because they wrap private functions.
  void r_compute_light_environment() {compute_light_environment();}
  void r_compute_vars_phys() {compute_vars_phys();}

private:
  void compute_light_environment();
  void rescale_light_environment();
  void compute_vars_phys();

  parameters_type parameters;
  std::vector<bool> is_resident;
  Environment environment;
  std::vector<species_type> species;
};

template <typename T>
Patch<T>::Patch(parameters_type p)
  : parameters(p),
    is_resident(p.is_resident),
    environment(make_environment(parameters)) {
  parameters.validate();
  for (auto s : parameters.strategies) {
    species.push_back(Species<T>(s));
  }
  reset();
}

template <typename T>
void Patch<T>::reset() {
  for (auto& s : species) {
    s.clear();
  }
  environment.clear();
  compute_light_environment();
  compute_vars_phys();
}

template <typename T>
double Patch<T>::height_max() const {
  double ret = 0.0;
  for (size_t i = 0; i < species.size(); ++i) {
    if (is_resident[i]) {
      ret = std::max(ret, species[i].height_max());
    }
  }
  return ret;
}

template <typename T>
double Patch<T>::area_leaf_above(double height) const {
  double tot = 0.0;
  for (size_t i = 0; i < species.size(); ++i) {
    if (is_resident[i]) {
      tot += species[i].area_leaf_above(height);
    }
  }
  return tot;
}

template <typename T>
double Patch<T>::canopy_openness(double height) const {
  // NOTE: patch_area does not appear in the EBT model formulation;
  // really we should require that it is 1.0, or drop it entirely.
  return exp(-parameters.c_ext * area_leaf_above(height) /
             parameters.patch_area);
}

template <typename T>
std::vector<double> Patch<T>::r_area_leaf_error(size_t species_index) const {
  const double tot_area_leaf = area_leaf_above(0.0);
  return species[species_index].r_area_leafs_error(tot_area_leaf);
}

template <typename T>
void Patch<T>::compute_light_environment() {
  if (parameters.n_residents() > 0) {
    auto f = [&] (double x) -> double {return canopy_openness(x);};
    environment.compute_light_environment(f, height_max());
  }
}

template <typename T>
void Patch<T>::rescale_light_environment() {
  if (parameters.n_residents() > 0) {
    auto f = [&] (double x) -> double {return canopy_openness(x);};
    environment.rescale_light_environment(f, height_max());
  }
}

template <typename T>
void Patch<T>::compute_vars_phys() {
  for (size_t i = 0; i < size(); ++i) {
    environment.set_seed_rain_index(i);
    species[i].compute_vars_phys(environment);
  }
}

// TODO: We should only be recomputing the light environment for the
// points that are below the height of the seedling -- not the entire
// light environment; probably worth just doing a rescale there?
template <typename T>
void Patch<T>::add_seed(size_t species_index) {
  species[species_index].add_seed();
  if (parameters.is_resident[species_index]) {
    compute_light_environment();
  }
}

template <typename T>
void Patch<T>::add_seeds(const std::vector<size_t>& species_index) {
  bool recompute = false;
  for (size_t i : species_index) {
    species[i].add_seed();
    recompute = recompute || parameters.is_resident[i];
  }
  if (recompute) {
    compute_light_environment();
  }
}

// Arguments here are:
//   time: time
//   state: vector of ode state; we'll pass an iterator with that in
//   n: number of *individuals* of each species
template <typename T>
void Patch<T>::r_set_state(double time,
                           const std::vector<double>& state,
                           const std::vector<size_t>& n) {
  const size_t n_species = species.size();
  util::check_length(n.size(), n_species);
  reset();
  for (size_t i = 0; i < n_species; ++i) {
    for (size_t j = 0; j < n[i]; ++j) {
      species[i].add_seed();
    }
  }
  util::check_length(state.size(), ode_size());
  set_ode_state(state.begin(), time);
}

// ODE interface
template <typename T>
size_t Patch<T>::ode_size() const {
  return ode::ode_size(species.begin(), species.end());
}

template <typename T>
double Patch<T>::ode_time() const {
  return time();
}

template <typename T>
ode::const_iterator Patch<T>::set_ode_state(ode::const_iterator it,
                                            double time) {
  it = ode::set_ode_state(species.begin(), species.end(), it);
  environment.time = time;
  if (parameters.control.environment_light_rescale_usually) {
    rescale_light_environment();
  } else {
    compute_light_environment();
  }
  compute_vars_phys();
  return it;
}

template <typename T>
ode::iterator Patch<T>::ode_state(ode::iterator it) const {
  return ode::ode_state(species.begin(), species.end(), it);
}

template <typename T>
ode::iterator Patch<T>::ode_rates(ode::iterator it) const {
  return ode::ode_rates(species.begin(), species.end(), it);
}

}

#endif
