// -*-c++-*-
#ifndef TREE_PLANT_APPROX_
#define TREE_PLANT_APPROX_

#include "plant_spline.h"

namespace model {

class PlantApprox : public Plant {
public:
  PlantApprox(PlantSpline ps);
  PlantApprox(PlantSpline *ps);

  // get_height -- still OK
  // get_mass_leaf -- still OK
  // set_mass_leaf -- still OK
  // leaf_area_above -- still OK [cf CohortDiscrete]
  void compute_vars_phys(spline::Spline *env);

  // * ODE interface
  ode::iter ode_rates(ode::iter it) const;

  // * R interface:
  void r_compute_vars_phys(spline::Spline env);

private:
  bool large_plant_do_exact() const;

  PlantSpline::ptr plant_spline;
  Strategy *strategy;
};

}

#endif
