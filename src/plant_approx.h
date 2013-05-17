// -*-c++-*-
#ifndef TREE_PLANT_APPROX_H_
#define TREE_PLANT_APPROX_H_

#include "plant_spline.h"

namespace model {

// This is a bit ugly as it requires a PlantSpline and a strategy --
// however, at the moment this is the least bad thing I can think of
// quickly.  It's possible that the C-indended constructor could skip
// this though.  And then the R-intended constructor could go through
// and check that the Strategy is actually OK.
class PlantApprox : public Plant {
public:
  PlantApprox(Strategy  s, PlantSpline  ps);
  PlantApprox(Strategy *s, PlantSpline *ps);

  void compute_vars_phys(const Environment& env);

  // * ODE interface
  ode::iter ode_rates(ode::iter it) const;

  // * R interface
  void r_compute_vars_phys_spline(const Environment& env);

private:
  bool large_plant_do_exact() const;

  PlantSpline::ptr plant_spline;
};

}

#endif
