// Inherit from FF16, like FF16r
#include <plant/models/ff16w_strategy.h>

namespace plant {
FF16w_Strategy::FF16w_Strategy() {
  collect_all_auxiliary = false;

  // build the string state/aux name to index map
  refresh_indices();
  name = "FF16w";
}



 std::vector<double> plant::FF16w_Strategy::compute_root_frac_per_layer(double soil_number_of_depths, double depth, double height){

  int n = soil_number_of_depths;

  const double taproot_depth = 0.5 * height;
  const double root_eta = 0.1; 
  
  std::vector<double> input_root_per_layer;

  for (size_t i = 0; i < n; i++) {

    double depth_1 = depth*(i+1);
    double depth_0 = depth*(i);
    double root_frac_1 = pow(1 - pow(depth_1/taproot_depth, root_eta),2);
    double root_frac_0 = pow(1 - pow(depth_0/taproot_depth, root_eta),2);

    double root_frac = root_frac_0 - root_frac_1;

    if(root_frac < 0){
      root_frac = 0;
    }

      // std::cout << "depth_1: " << depth_1 << "depth_0: " << depth_0 << "root_frac_1: " <<  root_frac_1 << "root_frac_0: " <<  root_frac_0 << "root_frac: " <<  root_frac << "taproot_depth: " << taproot_depth << std::endl;

    input_root_per_layer.push_back(root_frac);
    }
    return{input_root_per_layer};
  }

// not sure 'average' is the right term here..
double FF16w_Strategy::compute_average_light_environment(
    double z, double height, const FF16_Environment &environment) {
//NOTE: this function is currently being constrained at 0 because 

     return std::max(environment.get_environment_at_height(z), 0.0001) * q(z, height);
}

// assumes calc_profit_bartlett has been run for optimal psi_stem
double FF16w_Strategy::evapotranspiration_dt(double area_leaf_) {
  return leaf.transpiration_ * area_leaf_;
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

// calculate the average light environment based on the height of the plant
  double average_light_environment =
      integrator.is_adaptive() && reuse_intervals
          ? integrator.integrate_with_last_intervals(f, 0.0, height)
          : integrator.integrate(f, 0.0, height);
    
    // calculate average radiation by multipling average canopy openness by PPFD and accounting for self-shading k_I.
  const double average_radiation = k_I * average_light_environment * environment.PPFD;
    // calculate psi_soil (-MPa)
  // const double psi_soil = environment.get_psi_soil() / 1000000;

  //TOREMOVE

  


  std::vector<double> soil_moist = environment.get_soil_water_state();
  std::vector<double> input_root_per_layer = compute_root_frac_per_layer(environment.ode_size(), environment.delta_z, height);

  std::vector<double> input_soil_moist;

  double average_soil_moist = 0;

  for (size_t i = 0; i < environment.ode_size(); i++) {
    average_soil_moist += soil_moist[i]*input_root_per_layer[i];

    // std::cout << soil_moist[i] << "soil_moist" << std::endl; 
  }

  // std::cout << average_soil_moist << "average_soil_moist" << std::endl;
  
  
  const double psi_soil = environment.psi_from_soil_moist(average_soil_moist)/1000000;



  // find leaf specific max hydraulic conductance
  // K_s: max hydraulic conductivity (kg m^-2 s^-1 MPa^-1),
  // theta: huber value
  // eta_c: accounts for average position of leaf mass
  // height: maximum plant height
  const double leaf_specific_conductance_max = K_s * theta / (height * eta_c);

  // find sapwood volume per leaf area
  // theta: huber value
  // eta_c: accounts for average position of leaf mass

  // const double sapwood_volume_per_leaf_area = theta * (height * eta_c);

  const double sapwood_volume_per_leaf_area = (0.000157*(1-var_sapwood_volume_cost) + theta*var_sapwood_volume_cost)  * (height * eta_c);
// set strategy-level physiological parameters for the leaf-submodel.

  leaf.set_physiology(rho, a_bio, average_radiation, psi_soil, leaf_specific_conductance_max, environment.get_atm_vpd(), environment.get_ca(), sapwood_volume_per_leaf_area, environment.get_leaf_temp(), environment.get_atm_o2_kpa(), environment.get_atm_kpa());

  // optimise psi_stem, setting opt_psi_stem_, profit_, hydraulic_cost_, assim_colimited_ etc.
  // leaf.optimise_psi_stem_Bartlett_analytical();
  leaf.optimise_psi_stem_TF();

  // stomatal conductance to c02 (umol m^-2 s^-1)
  vars.set_aux(aux_index.at("transpiration_"), leaf.transpiration_);
  vars.set_aux(aux_index.at("stom_cond_CO2_"), leaf.stom_cond_CO2_);

  // optimum psi_stem (-MPa)
  vars.set_aux(aux_index.at("ci_"), leaf.ci_);

  vars.set_aux(aux_index.at("opt_psi_stem_"), leaf.opt_psi_stem_);
    
  // profit (umol m^-2 s^-1), assim_colimited_ - hydraulic_cost_
  vars.set_aux(aux_index.at("profit_"), leaf.profit_);

  // assim_colimted_(umol m^-2 s^-1), per leaf area
  vars.set_aux(aux_index.at("assim_colimited_"), leaf.assim_colimited_);
  
  // cost (umol m^-2 s^-1), hydraulic_cost_
  vars.set_aux(aux_index.at("hydraulic_cost_"), leaf.hydraulic_cost_);

  // convert assimilation per leaf area per second (umol m^-2 s^-1) to canopy-level total yearly assimilation (mol yr^-1)

  const double assimilation = leaf.profit_ * area_leaf_* 60*60*12*365/1e6;
    
  const double respiration_ = 
  respiration(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
      

  const double turnover_ = 
  turnover(mass_leaf_, mass_bark_, mass_sapwood_, mass_root_);

  vars.set_aux(aux_index.at("respiration_"), respiration_);
  vars.set_aux(aux_index.at("turnover_"), turnover_);


  return net_mass_production_dt_A(assimilation, respiration_, turnover_);

}

 double FF16w_Strategy::darea_leaf_dmass_live(double area_leaf) const {

  return 1.0 /
         (dmass_leaf_darea_leaf(area_leaf) * (1 + dmass_dN * nmass_l) +
          dmass_sapwood_darea_leaf(area_leaf) * (1 + dmass_dN * nmass_s)  +
          dmass_bark_darea_leaf(area_leaf) * (1 + dmass_dN * nmass_b) + 
          dmass_root_darea_leaf(area_leaf) * (1 + dmass_dN * nmass_r));

}

  std::vector<double> FF16w_Strategy::compute_root_frac(double soil_number_of_depths, double delta_z, double height) {

  int n = soil_number_of_depths;

  const double taproot_depth = 0.5 * height;
  const double root_eta = 0.1; 
  
  std::vector<double> input_root;

  for (size_t i = 0; i < n; i++) {

    double depth = delta_z*(i+1);
    double root_prob_dens = 2*root_eta*(1-pow(depth,root_eta)*pow(taproot_depth,-root_eta))*pow(depth,root_eta-1)*pow(taproot_depth,-root_eta);

    if(root_prob_dens < 0){
      root_prob_dens = 0;
    }

      // std::cout << root_prob_dens << "input_root_ex" << taproot_depth << "t_d" << delta_z << "d_z" << depth << "depth" <<  std::endl;

    input_root.push_back(root_prob_dens);
    }
    return{input_root};
  }

 // one-shot update of the scm variables
  // i.e. setting rates of ode vars from the state and updating aux vars
  void FF16w_Strategy::compute_rates(const FF16_Environment &environment,
                                     bool reuse_intervals, Internals &vars)
  {

      double height = vars.state(HEIGHT_INDEX);

      double area_leaf_ = vars.aux(aux_index.at("competition_effect"));

      const double net_mass_production_dt_ =
          net_mass_production_dt(environment, height, area_leaf_, vars, reuse_intervals);

      // store the aux sate
      vars.set_aux(aux_index.at("net_mass_production_dt"), net_mass_production_dt_);

      std::vector<double> input_root = compute_root_frac(environment.ode_size(), environment.delta_z, height);

      // convert evapotranspiration per leaf area (kg H20 m^-2 s^-1) to canopy-level total yearly assimilation (m yr^-1)
      // stubbing out E_p for integration
      for (size_t i = 0; i < environment.ode_size(); i++)
      {

          vars.set_consumption_rate(i, evapotranspiration_dt(area_leaf_) * 60 * 60 * 12 * 365 / 1000 * input_root[i]);
      }

      if (net_mass_production_dt_ > 0)
      {
          const double fraction_allocation_reproduction_ =
              fraction_allocation_reproduction(height);
          const double darea_leaf_dmass_live_ = darea_leaf_dmass_live(area_leaf_);
          const double fraction_allocation_growth_ =
              fraction_allocation_growth(height);
          const double area_leaf_dt = net_mass_production_dt_ *
                                      fraction_allocation_growth_ *
                                      darea_leaf_dmass_live_;

          vars.set_aux(aux_index.at("darea_leaf_dmass_live_"), darea_leaf_dmass_live_);

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

          if (collect_all_auxiliary)
          {
              vars.set_aux(aux_index.at("area_sapwood"), area_sapwood_);
          }
      }
      else
      {
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


void FF16w_Strategy::prepare_strategy()
{
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

    if (is_variable_birth_rate)
    {
        extrinsic_drivers.set_variable("birth_rate", birth_rate_x, birth_rate_y);
    } else {
    extrinsic_drivers.set_constant("birth_rate", birth_rate_y[0]);
  }
  leaf = Leaf(vcmax_25,  c,  b, psi_crit,  beta1,  beta2, jmax_25, hk_s, a, curv_fact_elec_trans,curv_fact_colim, control.newton_tol_abs, control.GSS_tol_abs,
           control.vulnerability_curve_ncontrol,
           control.ci_abs_tol,
           control.ci_niter);
}

FF16w_Strategy::ptr make_strategy_ptr(FF16w_Strategy s) {
  s.prepare_strategy();

  return std::make_shared<FF16w_Strategy>(s);
}
} // namespace plant
