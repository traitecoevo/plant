// Inherit from FF16, like FF16r
#include <plant/models/ff16w_strategy.h>

namespace plant {
FF16w_Strategy::FF16w_Strategy() {
  collect_all_auxiliary = false;

  // initialise leaf traits
  // leaf = Leaf(vcmax, p_50, c, b, psi_crit, huber_value, K_s);


  // build the string state/aux name to index map
  refresh_indices();
  name = "FF16w";
}

// not sure 'average' is the right term here..
double FF16w_Strategy::compute_average_light_environment(
    double z, double height, const FF16_Environment &environment) {
//NOTE: this function is currently being constrained at 0 because 

     return std::max(environment.get_environment_at_height(z), 0.0001) * q(z, height);
}

// assumes calc_profit_bartlett has been run for optimal psi_stem
double FF16w_Strategy::evapotranspiration_dt(double area_leaf_) {
  return leaf.E * area_leaf_;
}

// One shot calculation of net_mass_production_dt
// Used by establishment_probability() and compute_rates().
double FF16w_Strategy::net_mass_production_dt(const FF16_Environment &environment,
                                       double height, double area_leaf_, Internals &vars, 
                                       bool reuse_intervals) {
  const double mass_leaf_ = mass_leaf(area_leaf_);
  const double area_sapwood_ = area_sapwood(area_leaf_);
  const double mass_sapwood_ = mass_sapwood(area_sapwood_, height);
  const double area_bark_ = area_bark(area_leaf_);
  const double mass_bark_ = mass_bark(area_bark_, height);
  const double mass_root_ = mass_root(area_leaf_);

  // integrate over x from zero to `height`, with fixed canopy openness
  auto f = [&](double x) -> double {
    return compute_average_light_environment(x, height, environment);
 
  };

  double average_light_environment =
      integrator.is_adaptive() && reuse_intervals
          ? integrator.integrate_with_last_intervals(f, 0.0, height)
          : integrator.integrate(f, 0.0, height);

  // std::cout << "is_adaptive " << integrator.is_adaptive() << reuse_intervals << "reuse_intervals" << std::endl;

  if(average_light_environment < 0){
    for (double height_test = 0; height_test <= height; height_test += height/control.assimilator_integration_iterations){
      double compute_test = compute_average_light_environment(height_test, height, environment);

      std::cout  << "compute_test"  << compute_test  << "height_test" << height_test << "is_adaptive " << integrator.is_adaptive()  << "reuse_intervals" << reuse_intervals << std::endl;
    }

    util::stop("Error");
    }
    
    // k_I adds self shading in addition to patch competition
  const double average_radiation = k_I * average_light_environment * environment.PPFD;
    
  const double psi_soil = environment.get_psi_soil() / 1000000;
    
  // height * eta_c = height of average leaf area
  const double k_l_max = K_s * huber_value / (height * eta_c);
    
  leaf.set_physiology(average_radiation, psi_soil, k_l_max, leaf.atm_vpd, leaf.ca);
    
  // double ci_guess = vars.aux(aux_index.at("opt_ci"));
    
  // leaf.optimise_ci_Sperry_Newton_analytical(ci_guess);
    // std::cout << "started a strategy" << std::endl;
  double psi_guess = vars.aux(aux_index.at("opt_psi_stem"));

  // leaf.optimise_psi_stem_Sperry_Newton_analytical();


  leaf.optimise_psi_stem_Sperry_Newton_analytical(psi_guess);
// std::cout << "finished a strategy" << std::endl;


  vars.set_aux(aux_index.at("opt_psi_stem"), leaf.opt_psi_stem);
    
  vars.set_aux(aux_index.at("count"), leaf.count);

  vars.set_aux(aux_index.at("profit"), leaf.profit);

  vars.set_aux(aux_index.at("A_lim_"), leaf.A_lim_);


  const double assimilation_per_area = leaf.profit;
    
  const double assimilation = assimilation_per_area * area_leaf_* 60*60*12*365/1000000;
    
  const double respiration_ = 
  respiration(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
      
  const double turnover_ = 
  turnover(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
  
  return net_mass_production_dt_A(assimilation, respiration_, turnover_);

}

// one-shot update of the scm variables
// i.e. setting rates of ode vars from the state and updating aux vars
void FF16w_Strategy::compute_rates(const FF16_Environment &environment,
                                   bool reuse_intervals, Internals &vars) {

  double height = vars.state(HEIGHT_INDEX);
  double area_leaf_ = vars.aux(aux_index.at("competition_effect"));

  const double net_mass_production_dt_ =
      net_mass_production_dt(environment, height, area_leaf_, vars, reuse_intervals);


  // store the aux sate
  vars.set_aux(aux_index.at("net_mass_production_dt"), net_mass_production_dt_);

  // stubbing out E_p for integration
  for (size_t i = 0; i < environment.ode_size(); i++) {

    vars.set_consumption_rate(i, evapotranspiration_dt(area_leaf_)*60*60*12*365/1000);

    // std::cout << "area_leaf" << area_leaf_ << "water_use" << evapotranspiration_dt(area_leaf_)*60*60*12*365/1000 << std::endl;
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

void FF16w_Strategy::prepare_strategy() {
  // Set up the integrator
  initialize_integrator(control.assimilator_adaptive_integration,
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
  leaf = Leaf(vcmax, p_50, c, b, psi_crit, huber_value, K_s, epsilon_leaf);
}


FF16w_Strategy::ptr make_strategy_ptr(FF16w_Strategy s) {
  s.prepare_strategy();

  return std::make_shared<FF16w_Strategy>(s);
}
} // namespace plant
