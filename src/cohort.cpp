#include "cohort.h"

namespace model {

// On initialisation, both the "mean" and "top" plant will have the
// same initial mass, which will be the seed mass.
Cohort::Cohort(Strategy s)
  : strategy(s),
    plant_mean(strategy.get()),
    plant_top(strategy.get()) {
}

Cohort::Cohort(Strategy *s)
  : strategy(s),
    plant_mean(strategy.get()),
    plant_top(strategy.get()) {
}

// * ODE interface
size_t Cohort::ode_size() const { 
  return ode_dimension; 
}

ode::iter_const Cohort::ode_values_set(ode::iter_const it, bool &changed) {
  *it++; // TODO: Cohort-level variables to add.
  it = plant_mean.ode_values_set(it, changed);
  it = plant_top.ode_values_set(it, changed);
  return it;
}

ode::iter Cohort::ode_values(ode::iter it) const {
  *it++; // TODO: Cohort-level variables to add.
  it = plant_mean.ode_values(it);
  it = plant_top.ode_values(it);
  return it;
}

ode::iter Cohort::ode_rates(ode::iter it) const {
  *it++; // TODO: Cohort-level variables to add.
  it = plant_mean.ode_rates(it);
  it = plant_top.ode_rates(it);
  return it;
}

}
