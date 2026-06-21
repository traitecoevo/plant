// Built from  src/ff16_strategy.cpp on Mon Feb 12 09:52:27 2024 using the scaffolder, from the strategy:  FF16
#include <plant/models/tf24_strategy.h>

namespace plant {

// --- Hard-coded root-distribution constants (review #9) ---------------------
// Named here for clarity; promotion to user-tunable traits (RcppR6) is a
// deliberate follow-up (see vignettes/models/code_review_leaf_tf24.qmd #9).
// rescales total fine-root mass into the per-layer carbon units expected by the
// root hydraulic network in Leaf::set_physiology.
static const double root_mass_carbon_scale = 83.26 * 0.5;
// shape exponent for the Q() root-fraction-with-depth profile.
static const double root_depth_shape_eta = 0.2;
// rooting depth cap (m), i.e. the depth of the soil column.
static const double rooting_depth_max = 1.5;

// NOTE (review #9): the per-second -> annual factor 60*60*12*365 (seconds of
// daylight per year, 12 h day x 365 d) recurs in compute_rates and
// net_mass_production_dt below. It is deliberately left inline rather than
// hoisted to a constant: collapsing the 4-step integer product into one double
// changes the floating-point rounding, and the adaptive ODE amplifies it
// (offspring_production shifts ~0.2%). Kept inline to preserve bit-identical
// results.

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

// not sure 'average' is the right term here..
double TF24_Strategy::compute_average_light_environment(
    double z, double height, const TF24_Environment &environment) {
// NOTE: the light environment is clamped to a small positive floor (1e-4)
// rather than allowed to reach 0 (original rationale was never recorded;
// preserved as-is).

     return std::max(environment.get_environment_at_height(z), 0.0001) * q(z, height);
}

// assumes optimise_psi_stem_TF has been run for optimal psi_stem
double TF24_Strategy::evapotranspiration_dt(double area_leaf_, int soil_layer) {
  return leaf.soil_consumption_[soil_layer] * area_leaf_;
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

  // Cache integer indices for the keys used in the hot compute_rates path, so
  // it no longer does a std::map<string,int> lookup per derivs evaluation.
  aux_idx_competition_effect    = aux_index.at("competition_effect");
  aux_idx_height_inverse        = aux_index.at("height_inverse");
  aux_idx_net_mass_production_dt = aux_index.at("net_mass_production_dt");
  aux_idx_root_mass             = aux_index.at("root_mass");
  aux_idx_opt_psi_stem          = aux_index.at("opt_psi_stem");
  aux_idx_opt_root_psi          = aux_index.at("opt_root_psi");
  aux_idx_transpiration         = aux_index.at("transpiration");
  aux_idx_E_up                  = aux_index.at("E_up_");
  aux_idx_profit                = aux_index.at("profit");
  aux_idx_stom_cond_CO2         = aux_index.at("stom_cond_CO2");
  // area_sapwood is only registered when collect_all_auxiliary is set.
  aux_idx_area_sapwood = aux_index.count("area_sapwood") ? aux_index.at("area_sapwood") : -1;
  state_idx_area_heartwood      = state_index.at("area_heartwood");
  state_idx_mass_heartwood      = state_index.at("mass_heartwood");
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
    vars.set_aux(aux_idx_competition_effect, area_leaf(height));
    vars.set_aux(aux_idx_height_inverse, 1.0 / height);
  }
}


// one-shot update of the scm variables
// i.e. setting rates of ode vars from the state and updating aux vars
void TF24_Strategy::compute_rates(const TF24_Environment& environment,  Internals& vars) {
  double height = vars.state(HEIGHT_INDEX);
  double area_leaf_ = vars.aux(aux_idx_competition_effect);

  const double net_mass_production_dt_ =
    net_mass_production_dt(environment, height, area_leaf_,
                           vars.aux(aux_idx_height_inverse));

  // store the aux sate
  vars.set_aux(aux_idx_net_mass_production_dt, net_mass_production_dt_);
  vars.set_aux(aux_idx_root_mass, mass_root(area_leaf_));
  vars.set_aux(aux_idx_opt_psi_stem, leaf.opt_psi_stem_);
  vars.set_aux(aux_idx_opt_root_psi, leaf.root_collar_psi_);
  vars.set_aux(aux_idx_transpiration, leaf.transpiration_);
  vars.set_aux(aux_idx_E_up, leaf.E_up_);
  vars.set_aux(aux_idx_profit, leaf.profit_);
  vars.set_aux(aux_idx_stom_cond_CO2, leaf.stom_cond_CO2_);




  // consumption rates should be emerging from net_mass_produciton_dt
  // convert evapotranspiration per leaf area per soil layer (mol H20 m^-2 s^-1) to canopy-level total 
  // yearly evapotranspiration per soil layer (m yr^-1)
  // stubbing out E_p for integration
  int soil_number_of_depths_ = environment.get_soil_number_of_depths();


  for (int i = 0; i < soil_number_of_depths_; i++) {

    // evapotranspiration (mol H20 m^-2 s^-1 layer^-1)
    // consumption rate (m yr^-1 layer ^-1)
    vars.set_consumption_rate(i, evapotranspiration_dt(area_leaf_, i)*60*60*12*365/1000*kg_per_mol_h2o);
  }

  if (net_mass_production_dt_ > 0) {

    const double fraction_allocation_reproduction_ = fraction_allocation_reproduction(height);
    const double darea_leaf_dmass_live_ = darea_leaf_dmass_live(area_leaf_);
    const double fraction_allocation_growth_ = fraction_allocation_growth(height);
    const double area_leaf_dt = net_mass_production_dt_ * fraction_allocation_growth_ * darea_leaf_dmass_live_;

    vars.set_rate(HEIGHT_INDEX, dheight_darea_leaf(area_leaf_) * area_leaf_dt);
    vars.set_rate(FECUNDITY_INDEX,
      fecundity_dt(net_mass_production_dt_, fraction_allocation_reproduction_));

    vars.set_rate(state_idx_area_heartwood, area_heartwood_dt(area_leaf_));
    const double area_sapwood_ = area_sapwood(area_leaf_);
    const double mass_sapwood_ = mass_sapwood(area_sapwood_, height);
    vars.set_rate(state_idx_mass_heartwood, mass_heartwood_dt(mass_sapwood_));

    if (collect_all_auxiliary) {
      vars.set_aux(aux_idx_area_sapwood, area_sapwood_);
    }
  } else {
    vars.set_rate(HEIGHT_INDEX, 0.0);
    vars.set_rate(FECUNDITY_INDEX, 0.0);
    vars.set_rate(state_idx_area_heartwood, 0.0);
    vars.set_rate(state_idx_mass_heartwood, 0.0);
  }
  // [eqn 21] - Instantaneous mortality rate
  vars.set_rate(MORTALITY_INDEX,
      mortality_dt(net_mass_production_dt_ / area_leaf_, vars.state(MORTALITY_INDEX)));

}

// [eqn 12] Gross annual CO2 assimilation (!!not in use for TF24 model!!)
double TF24_Strategy::assimilation(const TF24_Environment& environment,
                                    double height,
                                    double area_leaf) {


  double A = 0.0;

  // Define an anonymous function to integrate
  // For given height in crown, take photosynthesis at depth multipled by 
  //   amount of leaf at that depth
  std::function<double(double)> f = [&](double z) -> double {
    return assimilation_leaf(environment.get_environment_at_height(z)) * q(z, height);
  };

  // Integrate over crown depth using using Gauss-Kronrod quadrature.
  // The number of points used in the integration is determined by the control parameter
  // function_integration_rule. Rules defined in qk_rules.cpp
  A = function_integrator.integrate(f, 0.0, height);

  return area_leaf * A;
}

// Photosynthetic rate per leaf area
// `x` is openness, ranging from 0 to 1.
double TF24_Strategy::assimilation_leaf(double x) const {
  return a_p1 * x / (x + a_p2);
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
double TF24_Strategy::net_mass_production_dt(const TF24_Environment& environment,
                                double height, double area_leaf_,
                                double height_inverse) {
  // height_inverse (= 1/height) is supplied by the shared individual.h interface
  // (cached aux); unused here as the TF24 root-water path works in height directly.
  (void)height_inverse;
  const double mass_leaf_    = mass_leaf(area_leaf_);
  const double area_sapwood_ = area_sapwood(area_leaf_);
  const double mass_sapwood_ = mass_sapwood(area_sapwood_, height);
  const double area_bark_    = area_bark(area_leaf_);
  const double mass_bark_    = mass_bark(area_bark_, height);
  const double mass_root_    = mass_root(area_leaf_);

  int soil_number_of_depths_ = environment.get_soil_number_of_depths();
  const std::vector<double>& soil_depths_ = environment.z;



  // The radiation that drives the leaf optimisation depends on the shading
  // model and is computed below (just before the optimisation), once the
  // depth-independent inputs are ready.

  // psi_soil (-MPa), computed once per soil state and cached in environment.
  const std::vector<double>& psi_soil = environment.get_soil_water_potential_state();
  
// find leaf specific max hydraulic conductance (kg m^-2 LA s^-1 MPa ^-1)
  // K_s: max hydraulic conductivity (kg m^-2 s^-1 MPa^-1),
  // theta: huber value
  // eta_c: accounts for average position of leaf mass
  // height: maximum plant height
  const double leaf_specific_conductance_max = K_s * theta / (height * eta_c);

  // find sapwood volume per leaf area
  // theta: huber value
  // eta_c: accounts for average position of leaf mass

  const double sapwood_volume_per_leaf_area = theta * (height * eta_c);
  
  // ----------------------------------------------------------------------
  // ROOT MASS DISTRIBUTION ACROSS SOIL LAYERS
  // ----------------------------------------------------------------------
  // Total fine-root mass (mass_root_) is distributed over depth using the same
  // cumulative shape function Q() used for the leaf canopy, but parameterised
  // over soil depth instead of crown height. Q(z, rooting_depth, 0.2) gives the
  // fraction of roots *below* depth z, so the mass in layer a is
  //   root_mass_scale * (Q(z_{a-1}) - Q(z_a)).
  // rooting_depth is capped at 1.5 m (the soil column depth). The constant
  // 83.26 * 0.5 rescales mass_root_ into the per-layer carbon units expected by
  // the root hydraulic network (set_physiology). The loop breaks early once Q
  // reaches 0 (below the rooting depth) to avoid touching empty deep layers.
  //
  // Reuse the member buffer (assign refills + zeroes without reallocating when
  // the layer count is unchanged); zeroing matters because the loop below breaks
  // early below the rooting depth, leaving deep layers that must read as 0.
  // TODO (perf): the rooting depth fraction (0.2), depth cap (1.5) and scale
  // (83.26) are hard-coded and should become traits.
  mass_root_prop_.assign(soil_number_of_depths_, 0.0);



  // Use Q function with new arghument
  // std::fill(mass_root_prop_.begin(), mass_root_prop_.end(), 0); 
  
// change to while?
// environment.get_soil_depths() should ask for the ath element to save calling for a new vector each time
// change environment.get_soil_number_of_depths() change to n or soemtyhing
    double rooting_depth = std::min(height, rooting_depth_max);
  const double root_mass_scale = root_mass_carbon_scale * mass_root_;
    // std::vector<double> Q_root;
    // Q_root.reserve(soil_number_of_depths_);

    double prev_q = 1.0;
    for (int a = 0; a < soil_number_of_depths_; ++a) {
      if(prev_q == 0){
        break;
      }
      const double q = Q(soil_depths_[a], rooting_depth, root_depth_shape_eta);

      mass_root_prop_[a] = root_mass_scale * (prev_q - q);
      prev_q = q;
    }

  // Reuse geometry precomputed by environment; avoids rebuilding z midpoints each call.
  leaf.z_soil_mid_ = environment.get_soil_mid_depths();
  leaf.use_precomputed_z_soil_mid_ = true;

  // Optimise the leaf at a given absorbed radiation: rebuilds physiology and
  // solves the root-collar water potential, leaving the leaf.* outputs
  // (profit_, transpiration_, soil_consumption_, opt_psi_stem_, ...) set. Only
  // the radiation argument varies between calls; every other input is
  // depth-independent and already computed above.
  auto optimise_at = [&](double radiation) {
    leaf.set_physiology(area_leaf_, mass_root_prop_, rho, a_bio, radiation, psi_soil, soil_depths_, leaf_specific_conductance_max, environment.get_atm_vpd(), environment.get_ca(), sapwood_volume_per_leaf_area, environment.get_leaf_temp(), environment.get_atm_o2_kpa(), environment.get_atm_kpa());
    leaf.find_root_collar_psi();
  };

  // Convert canopy openness (0-1) into absorbed radiation: PPFD attenuated by
  // the self-shading coefficient k_I. The light floor (1e-4) matches
  // compute_average_light_environment().
  const double PPFD = environment.get_PPFD();
  auto radiation_at = [&](double light) {
    return k_I * std::max(light, 0.0001) * PPFD;
  };

  // Aggregate the leaf submodel over the crown according to the shading model.
  // The expensive hydraulic optimisation is the unit of work here, so the model
  // choice is about how many times it runs and on what light:
  //  - crown-centre:  one optimisation at the crown-centre light.
  //  - mean-light:    one optimisation at the leaf-area-weighted mean light
  //                   (TF24's established default).
  //  - deep-crown:    one optimisation per crown-depth quadrature point, with
  //                   every leaf output integrated to a leaf-area-weighted mean.
  if (shading_model_ == ShadingModel::CrownCentre) {
    optimise_at(radiation_at(environment.get_environment_at_height(height * eta_c)));
  } else if (shading_model_ == ShadingModel::MeanLight) {
    // Leaf-area-weighted mean canopy openness = integral of (light * q) over the
    // crown (q integrates to one). radiation_at then applies k_I * PPFD, exactly
    // reproducing TF24's established average_radiation.
    auto f = [&](double x) -> double {
      return compute_average_light_environment(x, height, environment);
    };
    optimise_at(radiation_at(function_integrator.integrate(f, 0.0, height)));
  } else { // DeepCrown
    const std::vector<double> nodes =
      function_integrator.integrate_vector_x(0.0, height);
    const size_t nn = nodes.size();
    std::vector<double> profit_y(nn), trans_y(nn), eup_y(nn), psi_y(nn),
      root_psi_y(nn), gco2_y(nn);
    std::vector<std::vector<double>> soil_y(
      soil_number_of_depths_, std::vector<double>(nn));
    for (size_t i = 0; i < nn; ++i) {
      const double qi = q(nodes[i], height);
      optimise_at(radiation_at(environment.get_environment_at_height(nodes[i])));
      profit_y[i]   = leaf.profit_ * qi;
      trans_y[i]    = leaf.transpiration_ * qi;
      eup_y[i]      = leaf.E_up_ * qi;
      psi_y[i]      = leaf.opt_psi_stem_ * qi;
      root_psi_y[i] = leaf.root_collar_psi_ * qi;
      gco2_y[i]     = leaf.stom_cond_CO2_ * qi;
      for (int a = 0; a < soil_number_of_depths_; ++a) {
        soil_y[a][i] = leaf.soil_consumption_[a] * qi;
      }
    }
    // Integrate each leaf output to its leaf-area-weighted crown mean (q
    // integrates to one over the crown). soil_consumption_ feeds the patch
    // water balance, so it must be the depth-integrated total; the rest are
    // diagnostics reported through compute_rates.
    leaf.profit_          = function_integrator.integrate_vector(profit_y, 0.0, height);
    leaf.transpiration_   = function_integrator.integrate_vector(trans_y, 0.0, height);
    leaf.E_up_            = function_integrator.integrate_vector(eup_y, 0.0, height);
    leaf.opt_psi_stem_    = function_integrator.integrate_vector(psi_y, 0.0, height);
    leaf.root_collar_psi_ = function_integrator.integrate_vector(root_psi_y, 0.0, height);
    leaf.stom_cond_CO2_   = function_integrator.integrate_vector(gco2_y, 0.0, height);
    for (int a = 0; a < soil_number_of_depths_; ++a) {
      leaf.soil_consumption_[a] =
        function_integrator.integrate_vector(soil_y[a], 0.0, height);
    }
  }


  //TODO: one point constant ratio and integral width for daylength
  // convert assimilation per leaf area per second (umol m^-2 s^-1) to canopy-level total yearly assimilation (mol yr^-1)
  // converts to canopy area, then years, then mols
  const double assimilation_ = leaf.profit_ * area_leaf_* 60*60*12*365/1e6;
  // const double assimilation_ = assimilation(environment, height, area_leaf_);
  const double respiration_ =
    respiration(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
  const double turnover_ =
    turnover(mass_leaf_, mass_bark_, mass_sapwood_, mass_root_);
  return net_mass_production_dt_A(assimilation_, respiration_, turnover_);
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
  if (util::is_finite(cumulative_mortality)) {
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
    net_mass_production_dt(environment, height_0, area_leaf_0, 1.0 / height_0);
  if (net_mass_production_dt_ > 0) {
    const double tmp = a_d0 * area_leaf_0 / net_mass_production_dt_;
    return 1.0 / (tmp * tmp + 1.0) * decay_over_time;
  } else {
    return 0.0;
  }
}

double TF24_Strategy::compute_competition(double z, double height) const {
  return k_I * area_leaf(height) * Q(z, height, eta);
}

// Ratio-first hot-path overload (see header): receives the cached
// competition_effect (= area_leaf(height)) and height_inverse (= 1/height), so the
// per-call area_leaf() evaluation and z/height division are hoisted out of the
// inner competition loop. Reproduces k_I * area_leaf(height) * Q(z, height, eta).
double TF24_Strategy::compute_competition(double z, double area_leaf_,
                                          double height_inverse) const {
  const double u = z * height_inverse;  // z / height
  if (u > 1.0) {
    return 0.0;
  }
  const double tmp = 1.0 - pow(u, eta);
  return k_I * area_leaf_ * tmp * tmp;
}

// [eqn  9] Probability density of leaf area at height `z`
double TF24_Strategy::q(double z, double height) const {
  const double tmp = pow(z / height, eta);
  return 2 * eta * (1 - tmp) * tmp / z;
}

// [eqn 10] ... Fraction of leaf area above height 'z' for an
//              individual of height 'height'
double TF24_Strategy::Q(double z, double height, double eta_x) const {
  if (z > height) {
    return 0.0;
  }
  const double tmp = 1.0-pow(z / height, eta_x);
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

void TF24_Strategy::prepare_strategy() {

  // Set up the function_integrator
  function_integrator = quadrature::QK(
      // Gauss-Kronrod quadrature integeration rule (see qkrules)
      control.function_integration_rule);

  // Resolve the crown shading model once. The empty Control default maps to
  // TF24's own default (mean-light, its established behaviour); PPA is an
  // FF16-only stepped-light model and is rejected here.
  shading_model_ =
    shading_model_from_string(control.shading_model, ShadingModel::MeanLight);
  // PPA and the flat-top-box variants are FF16-only (they reshape the FF16
  // competition / light profile, which TF24 does not use).
  if (shading_model_ == ShadingModel::PPA ||
      shading_model_ == ShadingModel::FlatTopBox ||
      shading_model_ == ShadingModel::FlatTopSoftBox) {
    throw std::invalid_argument(
      "shading_model '" + control.shading_model +
      "' is not supported for the TF24 strategy");
  }

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
    leaf = Leaf(vcmax_25,  c,  b, psi_crit, root_c, root_b, root_psi_crit, beta2, jmax_25, hk_s, a, curv_fact_elec_trans,curv_fact_colim, control.GSS_tol_abs,
           control.vulnerability_curve_ncontrol,
           control.ci_abs_tol,
           control.ci_niter,g1_TF24, beta_R_H, beta_R_V);
   // set phsyiology possibly here?        
}

TF24_Strategy::ptr make_strategy_ptr(TF24_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<TF24_Strategy>(s);
}
}
