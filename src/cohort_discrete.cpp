#include "cohort_discrete.h"

namespace model {

CohortDiscrete::CohortDiscrete(Strategy s, int n_individuals)
  : standalone(true),
    strategy(new Strategy(s)),
    plant(strategy),
    n_individuals(n_individuals) {
}

CohortDiscrete::CohortDiscrete(Strategy *s, int n_individuals)
  : standalone(false),
    strategy(s),
    plant(strategy),
    n_individuals(n_individuals) {
}

CohortDiscrete::CohortDiscrete(const CohortDiscrete &other)
  : standalone(other.standalone),
    strategy(standalone ? new Strategy(*other.strategy) : other.strategy),
    plant(other.plant),
    n_individuals(other.n_individuals) {
}

CohortDiscrete& CohortDiscrete::operator=(CohortDiscrete rhs) {
  swap(*this, rhs);
  return *this;
}

void swap(CohortDiscrete &a, CohortDiscrete &b) {
  using std::swap;
  swap(a.standalone,    b.standalone);
  swap(a.strategy,      b.strategy);
  swap(a.plant,         b.plant);
  swap(a.n_individuals, b.n_individuals);
}

CohortDiscrete::~CohortDiscrete() {
  if ( standalone )
    delete strategy;
}

// * Wrappers to Plant
double CohortDiscrete::get_height() const {
  return plant.get_height();
}

void CohortDiscrete::set_mass_leaf(double mass_leaf_) {
  plant.set_mass_leaf(mass_leaf_);
}

double CohortDiscrete::leaf_area_above(double z) const {
  return n_individuals * plant.leaf_area_above(z);
}

void CohortDiscrete::compute_vars_phys(spline::Spline *env) {
  plant.compute_vars_phys(env);
}

int CohortDiscrete::offspring() {
  return n_individuals * plant.offspring();
}

bool CohortDiscrete::died() {
  n_individuals -= Rf_rbinom(n_individuals, plant.mortality());
  plant.mortality_reset();
  if ( n_individuals < 0 )
    Rf_error("Somehow we have negative individuals!");
  return n_individuals == 0;
}

double CohortDiscrete::germination_probability(spline::Spline *env) {
  return plant.germination_probability(env);
}

// * ODE interface
size_t CohortDiscrete::ode_size() const { 
  return plant.ode_size();
}

ode::iter_const CohortDiscrete::ode_values_set(ode::iter_const it, 
					       bool &changed) {
  return plant.ode_values_set(it, changed);
}

ode::iter CohortDiscrete::ode_values(ode::iter it) const {
  return plant.ode_values(it);
}

ode::iter CohortDiscrete::ode_rates(ode::iter it) const {
  return plant.ode_rates(it);
}

}
