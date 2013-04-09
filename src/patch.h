// -*-c++-*-
#ifndef TREE_PATCH_
#define TREE_PATCH_

#include <vector>

#include "ode_target.h"
#include "adaptive_spline.h"
#include "parameters.h"
#include "species.h"

namespace model {

class Patch : public ode::OdeTarget {
public:
  Patch(Parameters p);
  Patch(Parameters *p);

  Patch(const Patch &other);
  Patch& operator=(Patch rhs);
  ~Patch();

  // Advance the system through one complete time step.
  void step();

  // Advance the system through one time step deterministically
  // (plant growth, physiological accounting)
  void step_deterministic();
  // Advance the system through the stochastic life cycle stages
  // (producing seeds and dying).
  void step_stochastic();

  // * ODE interface.
  void derivs(double time,
	      std::vector<double>::const_iterator y,
	      std::vector<double>::iterator dydt);
  size_t ode_size() const;
  ode::iter_const ode_values_set(ode::iter_const it, bool &changed);
  ode::iter       ode_values(ode::iter it) const;
  ode::iter       ode_rates(ode::iter it)  const;

  // * R interface.

  // Actually public functions for interrogating & modifying
  Rcpp::List r_get_plants(int idx) const; // TODO: return all?
  spline::Spline r_light_environment() const;
  void r_add_seed(int idx); // TODO: May change?
  void r_step();
  void r_step_stochastic(); // step_stochastic, plus RNG control

  // Wrapper functions for testing
  size_t r_size() const;
  double r_height_max() const;
  double r_canopy_openness(double height);
  void r_compute_light_environment();
  void r_compute_vars_phys();
  double r_age() const;

  // TODO: This is likely to change as more is written.
  std::vector<double> r_get_mass_leaf(int idx) const;
  void r_set_mass_leaf(std::vector<double> x, int idx);
  
  // Also this
  void r_clear();
  
private:
  void swap(Patch &a, Patch &b);

  void initialise();

  // Number of species
  size_t size() const;

  // Maximum height for any species in the Patch
  double height_max() const;

  // [eqn 11] Canopy openness at `height`
  double canopy_openness(double height);

  void compute_light_environment();
  void compute_vars_phys();

  bool standalone;
  Parameters *parameters;
  double age;

  spline::AdaptiveSpline light_environment;

  std::vector< Species > species;
  ode::Solver<Patch> ode_solver;
};

}

#endif
