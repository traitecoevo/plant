#include "cohort.h"

namespace model {

// On initialisation, both the "mean" and "top" plant will have the
// same initial mass, which will be the seed mass.

// What I had originally tried was 
//   Cohort::Cohort(Strategy s)
//     : standalone(true),
//       strategy(new Strategy(s)),
//       plant_mean(strategy),
//       plant_top(strategy) {
//   }
// but that doesn't work because by the time strategy is passed to
// initialse `plant_mean` and `plant_top` the needed variables are not
// ready.  A better might be to do the prepare_strategy elsewhere
// entirely.
// 
// Another option would be to to just pass `s` to `plant_mean` and
// `plant_top`, which will initialse them there.  In that case though,
// we then don't have the same pointer in the two different plants
// (not a huge deal, but seems not that nice).
// 
// A different option would be to *require* that the strategy is
// initialised before making a plant out of it.  This only needs
// checking on initialisation of Plant, and can be done just before
// setting the leaf mass.
//
// For now, I'm just ignoring this, even though it means that the
// cohorts are set to a stupid size in standalone now.
Cohort::Cohort(Strategy s)
  : standalone(true),
    strategy(new Strategy(s)),
    plant_mean(strategy),
    plant_top(strategy) {
}

Cohort::Cohort(Strategy *s)
  : standalone(false),
    strategy(s),
    plant_mean(strategy),
    plant_top(strategy) {
}

Cohort::Cohort(const Cohort &other)
  : standalone(other.standalone),
    strategy(standalone ? new Strategy(*other.strategy) : other.strategy),
    plant_mean(strategy),
    plant_top(strategy) {
}

Cohort& Cohort::operator=(Cohort rhs) {
  swap(*this, rhs);
  return *this;
}

void swap(Cohort &a, Cohort &b) {
  using std::swap;
  swap(a.standalone, b.standalone);
  swap(a.strategy,   b.strategy);
  swap(a.plant_mean, b.plant_mean);
  swap(a.plant_top,  b.plant_top);
}

Cohort::~Cohort() {
  if ( standalone )
    delete strategy;
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
