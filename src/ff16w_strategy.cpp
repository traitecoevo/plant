// Inherit from FF16, like FF16r
#include <plant/models/ff16w_strategy.h>

namespace plant {
FF16w_Strategy::FF16w_Strategy() {
  collect_all_auxiliary = false;

  // initialise leaf traits
  // leaf = Leaf(vcmax, p_50, c, b, psi_crit, beta, beta_2, huber_value);

  // build the string state/aux name to index map
  refresh_indices();
  name = "FF16w";
}

double
FF16w_Strategy::compute_assimilation(double z, double height,
                                     const FF16_Environment &environment) {

  double radiation =
      environment.PPFD * environment.get_environment_at_height(z);

  const double psi_soil = environment.get_psi_soil();

  const double leaf_profit =
      leaf.optimise_profit_gss_bartlett(radiation, psi_soil, height);

  return leaf_profit * assimilator.q(z, height);
}

// this is super janky - we integrate assimilation twice and return different
// integrals - to be replaced with something less dumb.
double FF16w_Strategy::compute_E(double z, double height,
                                 const FF16_Environment &environment) {

  double radiation =
      environment.PPFD * environment.get_environment_at_height(z);

  const double psi_soil = environment.get_psi_soil();

  // E gets set in leaf.calc_assim_gross()
  const double E = leaf.optimise_E_gss_bartlett(radiation, psi_soil, height);

  return E * assimilator.q(z, height);
}



double FF16w_Strategy::evapotranspiration_dt(const FF16_Environment &environment,
                                      double height, double area_leaf,
                                      bool reuse_intervals) {

  // integrate over x from zero to `height`, with fixed canopy openness
  auto f = [&](double x) -> double {
    return compute_E(x, height, environment);
  };

  double E_integrated =
      assimilator.assimilate(f, height, area_leaf, reuse_intervals);

  return E_integrated * area_leaf;
}

// one-shot update of the scm variables
// i.e. setting rates of ode vars from the state and updating aux vars
void FF16w_Strategy::compute_rates(const FF16_Environment &environment,
                                   bool reuse_intervals, Internals &vars) {

  double height = vars.state(HEIGHT_INDEX);
  double area_leaf_ = vars.aux(aux_index.at("competition_effect"));

  const double net_mass_production_dt_ =
      net_mass_production_dt(environment, height, area_leaf_, reuse_intervals);

  // store the aux sate
  vars.set_aux(aux_index.at("net_mass_production_dt"), net_mass_production_dt_);

  // stubbing out E_p for integration
  for (size_t i = 0; i < environment.ode_size(); i++) {
    vars.set_consumption_rate(i, evapotranspiration_dt(environment, height,
                                                       area_leaf_,
                                                       reuse_intervals));
  }

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
                             vars.state(MORTALITY_INDEX)));
}

FF16w_Strategy::ptr make_strategy_ptr(FF16w_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<FF16w_Strategy>(s);
}
} // namespace plant
