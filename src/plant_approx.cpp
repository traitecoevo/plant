#include "plant_approx.h"

namespace model {

PlantApprox::PlantApprox(PlantSpline ps)
  : Plant(ps.get_strategy()),
    plant_spline(ps), 
    strategy(plant_spline->get_strategy()) {
}

PlantApprox::PlantApprox(PlantSpline *ps)
  : Plant(ps->get_strategy()),
    plant_spline(ps), 
    strategy(plant_spline->get_strategy()) {
}

void PlantApprox::compute_vars_phys(spline::Spline *env) {
  if ( large_plant_do_exact() )
    Plant::compute_vars_phys(env);
}

ode::iter PlantApprox::ode_rates(ode::iter it) const {
  if ( large_plant_do_exact() )
    return ode_rates(it);
  else
    return plant_spline->ode_rates(get_mass_leaf(), it);
}

void PlantApprox::r_compute_vars_phys(spline::Spline env) {
  plant_spline->compute_vars_phys(&env);
  compute_vars_phys(&env);
}

bool PlantApprox::large_plant_do_exact() const {
  return get_mass_leaf() > plant_spline->mass_leaf_max();
}

}
