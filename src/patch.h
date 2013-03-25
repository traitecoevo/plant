// -*-c++-*-
#ifndef TREE_PATCH_
#define TREE_PATCH_

#include <vector>

#include "parameters.h"
#include "species.h"
#include "adaptive_spline.h"

namespace model {

class Patch {
public:
  Patch(Parameters p);
  Patch(Parameters *p);

  Patch(const Patch &other);
  Patch& operator=(const Patch &rhs);
  ~Patch();

  size_t size() const;
  double height_max() const;

  // [eqn 11] Canopy openness at `height`
  // NOTE: I'd rather this was const, but that interferes with the
  // functor code for now (TODO?)
  double canopy_openness(double height);

  void compute_light_environment();
  spline::Spline get_light_environment() const;

  void compute_vars_phys();

  void derivs(double time,
	      std::vector<double>::const_iterator y,
	      std::vector<double>::iterator dydt);
  size_t ode_size() const;
  bool set_values(ode::iter_const it);
  void get_values(ode::iter it) const;
  void get_rates(ode::iter it)  const;

  // Below here is a jumble of R interface stuff, mostly.
  Rcpp::List get_plants(int idx) const;

  void add_seed(int idx);
  void r_add_seed(int idx);

  // This is likely to change as more is written.
  std::vector<double> r_get_mass_leaf(int idx) const;
  void r_set_mass_leaf(std::vector<double> x, int idx);

  std::vector<double> r_derivs(std::vector<double> y);
  void r_ode_values_set(std::vector<double> y);
  std::vector<double> r_ode_values() const;
  std::vector<double> r_ode_rates() const;

private:
  void set_strategies();

  bool standalone;
  Parameters *parameters;

  spline::AdaptiveSpline light_environment;

  std::vector< Species > species;
};

}

#endif
