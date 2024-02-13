// Built from  src/ff16_strategy.cpp on Mon Feb 12 09:52:27 2024 using the scaffolder, from the strategy:  FF16
#include <plant/models/tf24_strategy.h>

namespace plant {

// TODO: Document consistent argument order: l, b, s, h, r
// TODO: Document ordering of different types of variables (size
// before physiology, before compound things?)
// TODO: Consider moving to activating as an initialisation list?
TF24_Strategy::TF24_Strategy() {
  collect_all_auxiliary = false;
  // build the string state/aux name to index map
  refresh_indices();
  name = "TF24";
}

void TF24_Strategy::refresh_indices () {
    // Create and fill the name to state index maps
  state_index = std::map<std::string,int>();
  aux_index   = std::map<std::string,int>();
  std::vector<std::string> aux_names_vec = aux_names();
  std::vector<std::string> state_names_vec = state_names();
  for (size_t i = 0; i < state_names_vec.size(); i++) {
    state_index[state_names_vec[i]] = i;
  }
  for (size_t i = 0; i < aux_names_vec.size(); i++) {
    aux_index[aux_names_vec[i]] = i;
  }
}


std::vector<double> plant::TF24_Strategy::compute_root_frac_per_layer(double soil_number_of_depths, double depth, double height){

  int n = soil_number_of_depths;

  const double taproot_depth = 0.5 * height;
  const double root_eta = 0.5; 
  
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
    input_root_per_layer.push_back(root_frac);
    }
    return{input_root_per_layer};
  }


   std::vector<double> TF24_Strategy::compute_root_frac(double soil_number_of_depths, double delta_z, double height) {

  int n = soil_number_of_depths;

  const double taproot_depth = 0.5 * height;
  const double root_eta = 0.1; 
  
  std::vector<double> input_root;

  for (size_t i = 0; i < n; i++) {

    double depth = delta_z*(i+1);
    //Falster et al. 2011
    double root_prob_dens = 2*root_eta*(1-pow(depth,root_eta)*pow(taproot_depth,-root_eta))*pow(depth,root_eta-1)*pow(taproot_depth,-root_eta);

    if(root_prob_dens < 0){
      root_prob_dens = 0;
    }

      // std::cout << root_prob_dens << "input_root_ex" << taproot_depth << "t_d" << delta_z << "d_z" << depth << "depth" <<  std::endl;

    input_root.push_back(root_prob_dens);
    }
    return{input_root};
  }

// not sure 'average' is the right term here..
double TF24_Strategy::compute_average_light_environment(
    double z, double height, const TF24_Environment &environment) {
//NOTE: this function is currently being constrained at 0 because 

     return std::max(environment.get_environment_at_height(z), 0.0001) * q(z, height);
}

// assumes calc_profit_bartlett has been run for optimal psi_stem
double TF24_Strategy::evapotranspiration_dt(double area_leaf_) {
  return leaf.transpiration_ * area_leaf_;
}


// [eqn 2] area_leaf (inverse of [eqn 3])
double TF24_Strategy::area_leaf(double height) const {
  return pow(height / a_l1, 1.0 / a_l2);
}

// [eqn 1] mass_leaf (inverse of [eqn 2])
double TF24_Strategy::mass_leaf(double area_leaf) const {
  return area_leaf * lma;
}

// [eqn 4] area and mass of sapwood
double TF24_Strategy::area_sapwood(double area_leaf) const {
  return area_leaf * theta;
}

double TF24_Strategy::mass_sapwood(double area_sapwood, double height) const {
  return area_sapwood * height * eta_c * rho;
}

// [eqn 5] area and mass of bark
double TF24_Strategy::area_bark(double area_leaf) const {
  return a_b1 * area_leaf * theta;
}

double TF24_Strategy::mass_bark(double area_bark, double height) const {
  return area_bark * height * eta_c * rho;
}

double TF24_Strategy::area_stem(double area_bark, double area_sapwood,
                            double area_heartwood) const {
  return area_bark + area_sapwood + area_heartwood;
}

double TF24_Strategy::diameter_stem(double area_stem) const {
  return std::sqrt(4 * area_stem / M_PI);
}

// [eqn 7] Mass of (fine) roots
double TF24_Strategy::mass_root(double area_leaf) const {
  return a_r1 * area_leaf;
}

// [eqn 8] Total mass
double TF24_Strategy::mass_live(double mass_leaf, double mass_bark,
                           double mass_sapwood, double mass_root) const {
  return mass_leaf + mass_sapwood + mass_bark + mass_root;
}

double TF24_Strategy::mass_total(double mass_leaf, double mass_bark,
                            double mass_sapwood, double mass_heartwood,
                            double mass_root) const {
  return mass_leaf + mass_bark + mass_sapwood +  mass_heartwood + mass_root;
}

double TF24_Strategy::mass_above_ground(double mass_leaf, double mass_bark,
                            double mass_sapwood, double mass_root) const {
  return mass_leaf + mass_bark + mass_sapwood + mass_root;
}

// for updating auxiliary state
void TF24_Strategy::update_dependent_aux(const int index, Internals& vars) {
  if (index == HEIGHT_INDEX) {
    double height = vars.state(HEIGHT_INDEX);
    vars.set_aux(aux_index.at("competition_effect"), area_leaf(height));
  }
}


// one-shot update of the scm variables
// i.e. setting rates of ode vars from the state and updating aux vars
void TF24_Strategy::compute_rates(const TF24_Environment &environment,
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

// [eqn 13] Total maintenance respiration
// NOTE: In contrast with Falster ref model, we do not normalise by a_y*a_bio.
double TF24_Strategy::respiration(double mass_leaf, double mass_sapwood,
                             double mass_bark, double mass_root) const {
  return respiration_leaf(mass_leaf) +
         respiration_bark(mass_bark) +
         respiration_sapwood(mass_sapwood) +
         respiration_root(mass_root);
}

double TF24_Strategy::respiration_leaf(double mass) const {
  return r_l * mass;
}

double TF24_Strategy::respiration_bark(double mass) const {
  return r_b * mass;
}

double TF24_Strategy::respiration_sapwood(double mass) const {
  return r_s * mass;
}

double TF24_Strategy::respiration_root(double mass) const {
  return r_r * mass;
}

// [eqn 14] Total turnover
double TF24_Strategy::turnover(double mass_leaf, double mass_bark,
                          double mass_sapwood, double mass_root) const {
   return turnover_leaf(mass_leaf) +
          turnover_bark(mass_bark) +
          turnover_sapwood(mass_sapwood) +
          turnover_root(mass_root);
}

double TF24_Strategy::turnover_leaf(double mass) const {
  return k_l * mass;
}

double TF24_Strategy::turnover_bark(double mass) const {
  return k_b * mass;
}

double TF24_Strategy::turnover_sapwood(double mass) const {
  return k_s * mass;
}

double TF24_Strategy::turnover_root(double mass) const {
  return k_r * mass;
}

// [eqn 15] Net production
//
// NOTE: Translation of variable names from the Falster 2011.  Everything
// before the minus sign is SCM's N, our `net_mass_production_dt` is SCM's P.
double TF24_Strategy::net_mass_production_dt_A(double assimilation, double respiration,
                                double turnover) const {
  return a_bio * a_y * (assimilation - respiration) - turnover;
}

// One shot calculation of net_mass_production_dt
// Used by establishment_probability() and compute_rates().
double TF24_Strategy::net_mass_production_dt(const TF24_Environment &environment,
                                       double height, double area_leaf_, Internals& vars, 
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
      function_integrator.integrate(f, 0.0, height);
    
    // calculate average radiation by multipling average canopy openness by PPFD and accounting for self-shading k_I.
  const double average_radiation = k_I * average_light_environment * environment.PPFD;

  std::vector<double> soil_moist = environment.get_soil_water_state();
  std::vector<double> input_root_per_layer = compute_root_frac_per_layer(environment.ode_size(), environment.delta_z, height);

  std::vector<double> input_soil_moist;

  double average_soil_moist = 0;

  for (size_t i = 0; i < environment.ode_size(); i++) {
    average_soil_moist += soil_moist[i]*input_root_per_layer[i];
  }
    
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
  // vars.set_aux(aux_index.at("transpiration_"), leaf.transpiration_);
std::cout << "leaf.stom_cond_CO2_" << vars.aux_size << std::endl;

  vars.set_aux(aux_index.at("stom_cond_CO2_"), leaf.stom_cond_CO2_);
  std::cout << "progress" << std::endl;

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

// [eqn 16] Fraction of production allocated to reproduction
double TF24_Strategy::fraction_allocation_reproduction(double height) const {
  return a_f1 / (1.0 + exp(a_f2 * (1.0 - height / hmat)));
}

// Fraction of production allocated to growth
double TF24_Strategy::fraction_allocation_growth(double height) const {
  return 1.0 - fraction_allocation_reproduction(height);
}

// [eqn 17] Rate of offspring production
double TF24_Strategy::fecundity_dt(double net_mass_production_dt,
                               double fraction_allocation_reproduction) const {
  return net_mass_production_dt * fraction_allocation_reproduction /
    (omega + a_f3);
}

double TF24_Strategy::darea_leaf_dmass_live(double area_leaf) const {
  return 1.0/(  dmass_leaf_darea_leaf(area_leaf)
              + dmass_sapwood_darea_leaf(area_leaf)
              + dmass_bark_darea_leaf(area_leaf)
              + dmass_root_darea_leaf(area_leaf));
}

// TODO: Ordering below here needs working on, probably as @dfalster
// does equation documentation?
double TF24_Strategy::dheight_darea_leaf(double area_leaf) const {
  return a_l1 * a_l2 * pow(area_leaf, a_l2 - 1);
}

// Mass of leaf needed for new unit area leaf, d m_s / d a_l
double TF24_Strategy::dmass_leaf_darea_leaf(double /* area_leaf */) const {
  return lma;
}

// Mass of stem needed for new unit area leaf, d m_s / d a_l
double TF24_Strategy::dmass_sapwood_darea_leaf(double area_leaf) const {
  return rho * eta_c * a_l1 * theta * (a_l2 + 1.0) * pow(area_leaf, a_l2);
}

// Mass of bark needed for new unit area leaf, d m_b / d a_l
double TF24_Strategy::dmass_bark_darea_leaf(double area_leaf) const {
  return a_b1 * dmass_sapwood_darea_leaf(area_leaf);
}

// Mass of root needed for new unit area leaf, d m_r / d a_l
double TF24_Strategy::dmass_root_darea_leaf(double /* area_leaf */) const {
  return a_r1;
}

// Growth rate of basal diameter_stem per unit time
double TF24_Strategy::ddiameter_stem_darea_stem(double area_stem) const {
  return pow(M_PI * area_stem, -0.5);
}

// Growth rate of sapwood area at base per unit time
double TF24_Strategy::area_sapwood_dt(double area_leaf_dt) const {
  return area_leaf_dt * theta;
}

// Note, unlike others, heartwood growth does not depend on leaf area growth, but
// rather existing sapwood
double TF24_Strategy::area_heartwood_dt(double area_leaf) const {
  return k_s * area_sapwood(area_leaf);
}

// Growth rate of bark area at base per unit time
double TF24_Strategy::area_bark_dt(double area_leaf_dt) const {
  return a_b1 * area_leaf_dt * theta;
}

// Growth rate of stem basal area per unit time
double TF24_Strategy::area_stem_dt(double area_leaf,
                               double area_leaf_dt) const {
  return area_sapwood_dt(area_leaf_dt) +
    area_bark_dt(area_leaf_dt) +
    area_heartwood_dt(area_leaf);
}

// Growth rate of basal diameter_stem per unit time
double TF24_Strategy::diameter_stem_dt(double area_stem, double area_stem_dt) const {
  return ddiameter_stem_darea_stem(area_stem) * area_stem_dt;
}

// Growth rate of root mass per unit time
double TF24_Strategy::mass_root_dt(double area_leaf,
                               double area_leaf_dt) const {
  return area_leaf_dt * dmass_root_darea_leaf(area_leaf);
}

double TF24_Strategy::mass_live_dt(double fraction_allocation_reproduction,
                               double net_mass_production_dt) const {
  return (1 - fraction_allocation_reproduction) * net_mass_production_dt;
}

// TODO: Change top two to use mass_live_dt
double TF24_Strategy::mass_total_dt(double fraction_allocation_reproduction,
                                     double net_mass_production_dt,
                                     double mass_heartwood_dt) const {
  return mass_live_dt(fraction_allocation_reproduction, net_mass_production_dt) +
    mass_heartwood_dt;
}

// TODO: Do we not track root mass change?
double TF24_Strategy::mass_above_ground_dt(double area_leaf,
                                       double fraction_allocation_reproduction,
                                       double net_mass_production_dt,
                                       double mass_heartwood_dt,
                                       double area_leaf_dt) const {
  const double mass_root_dt =
    area_leaf_dt * dmass_root_darea_leaf(area_leaf);
  return mass_total_dt(fraction_allocation_reproduction, net_mass_production_dt,
                        mass_heartwood_dt) - mass_root_dt;
}

double TF24_Strategy::mass_heartwood_dt(double mass_sapwood) const {
  return turnover_sapwood(mass_sapwood);
}


double TF24_Strategy::mass_live_given_height(double height) const {
  double area_leaf_ = area_leaf(height);
  return mass_leaf(area_leaf_) +
         mass_bark(area_bark(area_leaf_), height) +
         mass_sapwood(area_sapwood(area_leaf_), height) +
         mass_root(area_leaf_);
}

double TF24_Strategy::height_given_mass_leaf(double mass_leaf) const {
  return a_l1 * pow(mass_leaf / lma, a_l2);
}

double TF24_Strategy::mortality_dt(double productivity_area,
                              double cumulative_mortality) const {

  // NOTE: When plants are extremely inviable, the rate of change in
  // mortality can be Inf, because net production is negative, leaf
  // area is small and so we get exp(big number).  However, most of
  // the time that happens we should get infinite mortality variable
  // levels and the rate of change won't matter.  It is possible that
  // we will need to trim this to some large finite value, but for
  // now, just checking that the actual mortality rate is finite.
  if (R_FINITE(cumulative_mortality)) {
    return
      mortality_growth_independent_dt() +
      mortality_growth_dependent_dt(productivity_area);
 } else {
    // If mortality probability is 1 (latency = Inf) then the rate
    // calculations break.  Setting them to zero gives the correct
    // behaviour.
    return 0.0;
  }
}

double TF24_Strategy::mortality_growth_independent_dt() const {
  return d_I;
}

double TF24_Strategy::mortality_growth_dependent_dt(double productivity_area) const {
  return a_dG1 * exp(-a_dG2 * productivity_area);
}

// [eqn 20] Survival of seedlings during establishment
double TF24_Strategy::establishment_probability(const TF24_Environment& environment) {
  
  double decay_over_time = exp(-recruitment_decay * environment.time);
  
  const double net_mass_production_dt_ =
    net_mass_production_dt(environment, height_0, area_leaf_0);
  if (net_mass_production_dt_ > 0) {
    const double tmp = a_d0 * area_leaf_0 / net_mass_production_dt_;
    return 1.0 / (tmp * tmp + 1.0) * decay_over_time;
  } else {
    return 0.0;
  }
}

double TF24_Strategy::compute_competition(double z, double height) const {
  return k_I * area_leaf(height) * Q(z, height);
}

// [eqn  9] Probability density of leaf area at height `z`
double TF24_Strategy::q(double z, double height) const {
  const double tmp = pow(z / height, eta);
  return 2 * eta * (1 - tmp) * tmp / z;
}

// [eqn 10] ... Fraction of leaf area above height 'z' for an
//              individual of height 'height'
double TF24_Strategy::Q(double z, double height) const {
  if (z > height) {
    return 0.0;
  }
  const double tmp = 1.0-pow(z / height, eta);
  return tmp * tmp;
}

// (inverse of [eqn 10]; return the height above which fraction 'x' of
// the leaf mass would be found).
double TF24_Strategy::Qp(double x, double height) const { // x in [0,1], unchecked.
  return pow(1 - sqrt(x), (1/eta)) * height;
}

// The aim is to find a plant height that gives the correct seed mass.
double TF24_Strategy::height_seed(void) const {

  // Note, these are not entirely correct bounds. Ideally we would use height
  // given *total* mass, not leaf mass, but that is difficult to calculate.
  // Using "height given leaf mass" will expand upper bound, but that's ok
  // most of time. Only issue is that could break with obscure parameter
  // values for LMA or height-leaf area scaling. Could instead use some
  // absolute maximum height for new seedling, e.g. 1m?
  const double
    h0 = height_given_mass_leaf(std::numeric_limits<double>::min()),
    h1 = height_given_mass_leaf(omega);

  const double tol = control.offspring_production_tol;
  const size_t max_iterations = control.offspring_production_iterations;

  auto target = [&] (double x) mutable -> double {
    return mass_live_given_height(x) - omega;
  };

  return util::uniroot(target, h0, h1, tol, max_iterations);
}

void TF24_Strategy::prepare_strategy()
{
  // Set up the function_integrator
  function_integrator = quadrature::QK(
      // Gauss-Kronrod quadrature integeration rule (see qkrules)
      control.function_integration_rule);

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

TF24_Strategy::ptr make_strategy_ptr(TF24_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<TF24_Strategy>(s);
}
}
