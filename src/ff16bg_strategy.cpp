// Built from  src/ff16_strategy.cpp on Fri Oct 30 11:43:30 2020 using the
// scaffolder, from the strategy:  FF16
#include <plant/models/ff16bg_strategy.h>

namespace plant {

// TODO: Document consistent argument order: l, b, s, h, r
// TODO: Document ordering of different types of variables (size
// before physiology, before compound things?)
// TODO: Consider moving to activating as an initialisation list?
FF16bg_Strategy::FF16bg_Strategy() {
  collect_all_auxiliary = false;
  // build the string state/aux name to index map
  refresh_indices();
  name = "FF16bg";
}

double FF16bg_Strategy::below_ground_influence(const FF16_Environment &environment) {
  double e_0 = environment.get_environment_at_height(height_0 + 1e-6);
  return (pow(e_0, k_2 / k_I));
}

// Add rate limit based on total stand leaf area
double FF16bg_Strategy::net_mass_production_dt(const FF16_Environment &environment,
                                        double height, double area_leaf_,
                                        bool reuse_intervals) {
  const double mass_leaf_ = mass_leaf(area_leaf_);
  const double area_sapwood_ = area_sapwood(area_leaf_);
  const double mass_sapwood_ = mass_sapwood(area_sapwood_, height);
  const double area_bark_ = area_bark(area_leaf_);
  const double mass_bark_ = mass_bark(area_bark_, height);
  const double mass_root_ = mass_root(area_leaf_);

  double assimilation_ =
      assimilator.assimilate(environment, height, area_leaf_, reuse_intervals);

  // Reduce net assimilation of all individuals
  assimilation_ *= below_ground_influence(environment);

  const double respiration_ =
      respiration(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
  const double turnover_ =
      turnover(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
  return net_mass_production_dt_A(assimilation_, respiration_, turnover_);
}

// one - shot update of the scm variables
// i.e. setting rates of ode vars from the state and updating aux vars
void FF16bg_Strategy::compute_rates(const FF16_Environment &environment,
                                   bool reuse_intervals, Internals &vars) {

  double height = vars.state(HEIGHT_INDEX);
  double area_leaf_ = vars.aux(aux_index.at("competition_effect"));

  const double net_mass_production_dt_ =
      net_mass_production_dt(environment, height, area_leaf_, reuse_intervals);

  // store the aux sate
  vars.set_aux(aux_index.at("net_mass_production_dt"), net_mass_production_dt_);

  if (net_mass_production_dt_ > 0) {

    const double fraction_allocation_reproduction_ =
        fraction_allocation_reproduction(height);
    const double darea_leaf_dmass_live_ = darea_leaf_dmass_live(area_leaf_);
    const double fraction_allocation_growth_ =
        fraction_allocation_growth(height);
    const double area_leaf_dt = net_mass_production_dt_ *
                                fraction_allocation_growth_ *
                                darea_leaf_dmass_live_;

    vars.set_rate(HEIGHT_INDEX, dheight_darea_leaf(area_leaf_) * area_leaf_dt);
    vars.set_rate(FECUNDITY_INDEX,
                  fecundity_dt(net_mass_production_dt_,
                               fraction_allocation_reproduction_));

    vars.set_rate(state_index.at("area_heartwood"),
                  area_heartwood_dt(area_leaf_));
    const double area_sapwood_ = area_sapwood(area_leaf_);
    const double mass_sapwood_ = mass_sapwood(area_sapwood_, height);
    vars.set_rate(state_index.at("mass_heartwood"),
                  mass_heartwood_dt(mass_sapwood_));

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
                mortality_dt(net_mass_production_dt_ / area_leaf_,
                             vars.state(MORTALITY_INDEX), environment.time));
}

double FF16bg_Strategy::mortality_dt(double productivity_area,
                                    double cumulative_mortality,
                                    double time) const {
  if (R_FINITE(cumulative_mortality)) {
    return mortality_growth_independent_dt() +
           mortality_growth_dependent_dt(productivity_area) +
           extrinsic_drivers.evaluate("mortality_rate", time);
  } else {
    // If mortality probability is 1 (latency = Inf) then the rate
    // calculations break.  Setting them to zero gives the correct
    // behaviour.
    return 0.0;
  }
}

void FF16bg_Strategy::prepare_strategy() {
  // Set up the integrator
  assimilator.initialize(a_p1, a_p2, eta,
                         control.assimilator_adaptive_integration,
                         control.assimilator_integration_rule,
                         control.assimilator_integration_iterations,
                         control.assimilator_integration_tol);

  // NOTE: this pre-computes something to save a very small amount of time
  eta_c = 1 - 2 / (1 + eta) + 1 / (1 + 2 * eta);
  // NOTE: Also pre-computing, though less trivial
  height_0 = height_seed();
  area_leaf_0 = area_leaf(height_0);

  if (is_variable_birth_rate) {
    extrinsic_drivers.set_variable("birth_rate", birth_rate_x, birth_rate_y);
  } else {
    extrinsic_drivers.set_constant("birth_rate", birth_rate_y[0]);
  }

  if (is_variable_mortality_rate) {
    extrinsic_drivers.set_variable("mortality_rate", mortality_rate_x,
                                   mortality_rate_y);
  } else {
    extrinsic_drivers.set_constant("mortality_rate", mortality_rate_y[0]);
  }
}

FF16bg_Strategy::ptr make_strategy_ptr(FF16bg_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<FF16bg_Strategy>(s);
}
} // namespace plant
