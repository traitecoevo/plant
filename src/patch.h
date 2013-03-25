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

  Rcpp::List get_plants(int idx) const;

  void add_seed(int idx);
  void r_add_seed(int idx);

  // This is likely to change as more is written.
  std::vector<double> r_get_mass_leaf(int idx) const;
  void r_set_mass_leaf(std::vector<double> x, int idx);

private:
  void set_strategies();

  bool standalone;
  Parameters *parameters;

  spline::AdaptiveSpline light_environment;

  std::vector< Species > species;
};

}

#endif
