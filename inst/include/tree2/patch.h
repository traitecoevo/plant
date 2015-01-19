// -*-c++-*-
#ifndef TREE_PATCH_H_
#define TREE_PATCH_H_

#include <tree2/parameters.h>
#include <tree2/species.h>

// EBT uses:
//   .size
//   .add_seedling(int)
//   .get_disturbance_regime() (!) // used in seed_rain_cohort
//   .at // also in seed_rain_cohort -- probably needs work there...
//   OK .leaf_area_above(double)
//   .leaf_area_error(int)
//   .reset // hmm...
//   .{various ode functions}
//   // access cohort_schedule (will do via Parameters now?)
//   .seed_rains (all of them)
//   .set_times (schedule crap)

namespace tree2 {

template <typename T>
class Patch {
public:
  typedef T cohort_type;
  Patch(Parameters p);

  size_t size() const {return species.size();}
  void reset();

  // Not sure how I tended to add seeds/seedlings to these?  I might
  // not have really, outside of the main model.  Wait and see...

  double height_max() const;

  // [eqn 11] Canopy openness at `height`
  double leaf_area_above(double height) const;
  double canopy_openness(double height) const;

  void add_seed(size_t species_index);
  void add_seeds(const std::vector<size_t>& species_index);

  // * ODE interface
  size_t ode_size() const;
  ode::const_iterator set_ode_values(ode::const_iterator it, double time);
  ode::iterator       ode_values(ode::iterator it) const;
  ode::iterator       ode_rates(ode::iterator it) const;

  // * R interface
  // Data accessors:
  Parameters r_parameters() const {return parameters;}
  Environment r_environment() const {return environment;}
  std::vector<Species<T> > r_species() const {return species;}

  void r_add_seed(util::index species_index) {
    add_seed(species_index.check_bounds(size()));
  }
  void r_compute_light_environment() {compute_light_environment();}
  void r_compute_vars_phys() {compute_vars_phys();}

private:
  void compute_light_environment();
  void rescale_light_environment();
  void compute_vars_phys();

  Parameters parameters;
  std::vector<bool> is_resident;
  Environment environment;
  std::vector<Species<T> > species;
};

template <typename T>
Patch<T>::Patch(Parameters p)
  : parameters(p),
    is_resident(p.is_resident),
    environment(parameters) {
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

// TODO: use this https://github.com/klmr/cpp11-range for ranged loops:
//   for (size_t i : range(0, species.size())) {
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
double Patch<T>::leaf_area_above(double height) const {
  double tot = 0.0;
  for (size_t i = 0; i < species.size(); ++i) {
    if (is_resident[i]) {
      tot += species[i].leaf_area_above(height);
    }
  }
  return tot;
}

template <typename T>
double Patch<T>::canopy_openness(double height) const {
  // NOTE: patch_area does not appear in the EBT model formulation;
  // really we should require that it is 1.0, or drop it entirely.
  return exp(-parameters.c_ext * leaf_area_above(height) /
             parameters.patch_area);
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
  for (size_t species_index = 0; species_index < size(); ++species_index) {
    environment.set_seed_rain_index(species_index);
    species[species_index].compute_vars_phys(environment);
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

// ODE interface
template <typename T>
size_t Patch<T>::ode_size() const {
  return ode::ode_size(species.begin(), species.end());
}

template <typename T>
ode::const_iterator Patch<T>::set_ode_values(ode::const_iterator it,
                                             double time) {
  it = ode::set_ode_values(species.begin(), species.end(), it);
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
ode::iterator Patch<T>::ode_values(ode::iterator it) const {
  return ode::ode_values(species.begin(), species.end(), it);
}

template <typename T>
ode::iterator Patch<T>::ode_rates(ode::iterator it) const {
  return ode::ode_rates(species.begin(), species.end(), it);
}

}

#endif
