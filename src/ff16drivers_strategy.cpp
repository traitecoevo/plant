#include <plant/models/ff16drivers_strategy.h>

namespace plant {

FF16drivers_Strategy::FF16drivers_Strategy() {
  collect_all_auxiliary = false;
  // build the string state/aux name to index map
  refresh_indices();
  name = "FF16";
}

// one-shot update of the scm variables
// i.e. setting rates of ode vars from the state and updating aux vars
void FF16drivers_Strategy::compute_rates(const FF16_Environment& environment,
                              bool reuse_intervals,
                              Internals& vars) {

  double height = vars.state(HEIGHT_INDEX);
  double area_leaf_ = vars.aux(aux_index.at("competition_effect"));

  const double net_mass_production_dt_ =
    net_mass_production_dt(environment, height, area_leaf_, environment.time, reuse_intervals);

  // store the aux sate
  vars.set_aux(aux_index.at("net_mass_production_dt"), net_mass_production_dt_);

  if (net_mass_production_dt_ > 0) {

    const double fraction_allocation_reproduction_ = fraction_allocation_reproduction(height);
    const double darea_leaf_dmass_live_ = darea_leaf_dmass_live(area_leaf_);
    const double fraction_allocation_growth_ = fraction_allocation_growth(height);
    const double area_leaf_dt = net_mass_production_dt_ * fraction_allocation_growth_ * darea_leaf_dmass_live_;

    vars.set_rate(HEIGHT_INDEX, dheight_darea_leaf(area_leaf_) * area_leaf_dt);
    vars.set_rate(FECUNDITY_INDEX,
      fecundity_dt(net_mass_production_dt_, fraction_allocation_reproduction_));

    vars.set_rate(state_index.at("area_heartwood"), area_heartwood_dt(area_leaf_));
    const double area_sapwood_ = area_sapwood(area_leaf_);
    const double mass_sapwood_ = mass_sapwood(area_sapwood_, height);
    vars.set_rate(state_index.at("mass_heartwood"), mass_heartwood_dt(mass_sapwood_));

    if (collect_all_auxiliary) {
      vars.set_aux(aux_index.at("area_sapwood"), area_sapwood_);
    }
  } else {
    vars.set_rate(HEIGHT_INDEX, 0.0);
    vars.set_rate(FECUNDITY_INDEX, 0.0);
    vars.set_rate(state_index.at("area_heartwood"), 0.0);
    vars.set_rate(state_index.at("mass_heartwood"), 0.0);
  }
  // [eqn 21] - Instantaneous mortality rate
  vars.set_rate(MORTALITY_INDEX,
      mortality_dt(net_mass_production_dt_ / area_leaf_, vars.state(MORTALITY_INDEX), environment.time));
}

// One shot calculation of net_mass_production_dt
// Used by establishment_probability() and compute_rates().
double FF16drivers_Strategy::net_mass_production_dt(const FF16_Environment& environment,
                                              double height, double area_leaf_, double time,
                                              bool reuse_intervals) {

  const double mass_leaf_    = mass_leaf(area_leaf_);
  const double area_sapwood_ = area_sapwood(area_leaf_);
  const double mass_sapwood_ = mass_sapwood(area_sapwood_, height);
  const double area_bark_    = area_bark(area_leaf_);
  const double mass_bark_    = mass_bark(area_bark_, height);
  const double mass_root_    = mass_root(area_leaf_);

  const double assimilation_ = assimilator.assimilate(environment, 
                                                      height,
                                                      area_leaf_, 
                                                      reuse_intervals) *
                                extrinsic_drivers.evaluate("growth_multiplier", time);

  const double respiration_ =
    respiration(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
  const double turnover_ =
    turnover(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
  return net_mass_production_dt_A(assimilation_, respiration_, turnover_);
}


double FF16drivers_Strategy::mortality_dt(double productivity_area,
                                   double cumulative_mortality,
                                   double time) const {

  // NOTE: When plants are extremely inviable, the rate of change in
  // mortality can be Inf, because net production is negative, leaf
  // area is small and so we get exp(big number).  However, most of
  // the time that happens we should get infinite mortality variable
  // levels and the rate of change won't matter.  It is possible that
  // we will need to trim this to some large finite value, but for
  // now, just checking that the actual mortality rate is finite.
  if (R_FINITE(cumulative_mortality)) {
    return
      mortality_growth_independent_dt(time) +
      mortality_growth_dependent_dt(productivity_area);
 } else {
    // If mortality probability is 1 (latency = Inf) then the rate
    // calculations break.  Setting them to zero gives the correct
    // behaviour.
    return 0.0;
  }
}

double FF16drivers_Strategy::mortality_growth_independent_dt(double time) const {
  return d_I * extrinsic_drivers.evaluate("growth_independent_mortality_multiplier", time);
}

// [eqn 20] Survival of seedlings during establishment
double FF16drivers_Strategy::establishment_probability(const FF16_Environment& environment) {
  
  double decay_over_time = exp(-recruitment_decay * environment.time);
  
  const double net_mass_production_dt_ =
    net_mass_production_dt(environment, height_0, area_leaf_0, environment.time);
  if (net_mass_production_dt_ > 0) {
    const double tmp = a_d0 * area_leaf_0 / net_mass_production_dt_;
    return 1.0 / (tmp * tmp + 1.0) * decay_over_time;
  } else {
    return 0.0;
  }
}

void FF16drivers_Strategy::prepare_strategy() {
  // Set up the integrator
  assimilator.initialize(a_p1, a_p2, eta,
                         control.assimilator_adaptive_integration,
                         control.assimilator_integration_rule,
                         control.assimilator_integration_iterations,
                         control.assimilator_integration_tol);

  // NOTE: this pre-computes something to save a very small amount of time
  eta_c = 1 - 2/(1 + eta) + 1/(1 + 2*eta);
  // NOTE: Also pre-computing, though less trivial
  height_0 = height_seed();
  area_leaf_0 = area_leaf(height_0);

  if (is_variable_birth_rate) {
    extrinsic_drivers.set_variable("birth_rate", birth_rate_x, birth_rate_y);
  } else {
    extrinsic_drivers.set_constant("birth_rate", birth_rate_y[0]);
  }

  if (is_variable_growth_rate) {
    extrinsic_drivers.set_variable("growth_multiplier", growth_rate_x, growth_rate_y);
  } else {
    extrinsic_drivers.set_constant("growth_multiplier", growth_rate_y[0]);
  }

  if (is_variable_mortality_rate) {
    extrinsic_drivers.set_variable("growth_independent_mortality_multiplier", mortality_rate_x, mortality_rate_y);
  } else {
    extrinsic_drivers.set_constant("growth_independent_mortality_multiplier", mortality_rate_y[0]);
  }
}

FF16drivers_Strategy::ptr make_strategy_ptr(FF16drivers_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<FF16drivers_Strategy>(s);
}
}
