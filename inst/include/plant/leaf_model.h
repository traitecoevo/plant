// -*-c++-*-
#ifndef PLANT_PLANT_LEAF_MODEL_H_
#define PLANT_PLANT_LEAF_MODEL_H_

// TODO: replace with constants
// #define umol_per_mol_to_mol_per_mol 0.000...
// #define umol_per_mol_to_Pa ...
// #define kg_to_mol_h2o ...
// #define kPa_to_Pa ...

#include <plant/models/tf24_environment.h>
#include <plant/qag.h>
#include <plant/uniroot.h>
#include <plant/optimize.h>

namespace plant {

// double check best namespace for constants (private vs global)
// converts vcmax to jmax 25 unitless (Sperry et el. (2017))
static const double vcmax_25_to_jmax_25 = 1.67;

// kJ mol ^-1
static const double vcmax_ha = 60000;
// kJ mol ^-1
static const double vcmax_H_d = 200000;
// kJ mol ^-1
static const double vcmax_d_S = 650;

// kJ mol ^-1
static const double jmax_ha = 30000;
// kJ mol ^-1
static const double jmax_H_d = 200000;
// kJ mol ^-1
static const double jmax_d_S = 650;

// umol ^ -1 mol ^ 1
static const double gamma_25 = 42.75;
// dimensionless
static const double gamma_c = 19.02;
// kJ mol ^-1
static const double gamma_ha = 37.83e3;

// umol mol ^-1
static const double kc_25 = 404.9 ;
// dimensionless
static const double kc_c = 38.05;
// kJ mol ^-1
static const double kc_ha = 79.43e3;

// umol mol ^-1
static const double ko_25 = 278400 ;
// dimensionless
static const double ko_c = 20.30;
// kJ mol ^-1
static const double ko_ha = 36.38e3;

// Pa umol ^ -1 mol ^ 1 
static const double umol_per_mol_to_Pa = 0.1013;

// mol H2o kg ^-1
static const double kg_to_mol_h2o = 55.4939;
// mol mol ^-1 / (umol mol ^-1)
static const double umol_to_mol = 1e-6;
// Pa kPa^-1
static const double kPa_to_Pa = 1000.0;

// universal gas constant J mol^-1 K^-1
static const double R = 8.314;

//convert deg C to deg K
static const double C_to_K = 273.15;

//H20:CO2 stomatal diffusion ratio
static const double H2O_CO2_stom_diff_ratio = 1.67;

const double gravity_head = 9.8e-3; // MPa / m

// number of intergration steps
const double n = 5;

class Leaf {
public:
  //anonymous Leaf function as in canopy.h
  Leaf();
  
  Leaf(double vcmax_25, 
       double c, 
       double b, 
       double psi_crit,
       double root_c,
       double root_b,
       double root_psi_crit,
       double beta2, 
       double jmax_25, 
       double hk_s,
       double a, 
       double curv_fact_elec_trans, 
       double curv_fact_colim,
       double GSS_tol_abs,
       double vulnerability_curve_ncontrol,
       double ci_abs_tol,
       double ci_niter,
      double g1_TF24,
    double beta_R_H,
    double beta_R_V); 
        
  quadrature::QAG integrator;
  interpolator::Interpolator transpiration_from_psi;
  interpolator::Interpolator psi_from_transpiration;
  // pre-computed root vulnerability curve (same role as transpiration_from_psi for xylem)
  interpolator::Interpolator root_vuln_from_psi;
  // cumulative integral of the root vulnerability curve, G(m) = int_0^m f_r(s) ds,
  // indexed by magnitude m = -psi. Lets E_from_Soil_to_Root_Collar obtain the
  // mean conductivity over a potential interval from 2 evals instead of (n+1)
  // (same pre-integrated-curve trick as transpiration_from_psi for the xylem).
  interpolator::Interpolator root_vuln_integral_from_psi;

  // psi_from_E

  double vcmax_25;
  double c;
  double b;
  double psi_crit;  // derived from b and c
  double root_c;
  double root_b;
  double root_psi_crit;
  double beta2;
  double jmax_25;
  double hk_s;
  double a;
  double curv_fact_elec_trans; // unitless - obtained from Smith and Keenan (2020)
  double curv_fact_colim;
  double GSS_tol_abs;
  double vulnerability_curve_ncontrol;
  double ci_abs_tol;
  double ci_niter;
  double g1_TF24;
  double beta_R_H;
  double beta_R_V;
  double soil_number_of_depths_;
  int max_soil_layer;

  double ci_;
  double stom_cond_CO2_;
  double assim_colimited_;
  double transpiration_;
  std::vector<double> soil_consumption_;
  double profit_;
  double psi_stem;
  double lambda_;
  double lambda_analytical_;
  double hydraulic_cost_;
  
  double electron_transport_;
  double gamma_;
  double ko_;
  double kc_;
  double km_;
  double R_d_;
  double leaf_specific_conductance_max_;
  double sapwood_volume_per_leaf_area_;
  double k_s_;
  double root_mass_;
  std::vector<double> c_r_V_;
  std::vector<double> c_r_H_;
  double area_leaf_;
  double rho_;
  double vcmax_;
  double jmax_;
  double lma_; //kg m^-2
  double a_bio_;
  
  std::vector<double> psi_soil_;
  std::vector<double> psi_soil_inverted_;
  // Per-layer cache of root_vuln_integral_from_psi.eval(-psi_soil_inverted_[i]).
  // psi_soil_inverted_ is fixed for the whole find_root_collar_psi solve, so the
  // soil-side endpoint of the cumulative-integral lookup in
  // E_from_Soil_to_Root_Collar is constant across every (re)evaluation of the
  // nested root-finders. Precomputing it once per solve (alongside the
  // P_x_r-side eval, hoisted out of the layer loop) collapses ~2 spline evals
  // per layer to ~1 per call. Rebuilt in find_root_collar_psi.
  std::vector<double> root_vuln_integral_soil_;
  std::vector<double> soil_depth_;
  std::vector<double> z_soil_mid_;
  bool use_precomputed_z_soil_mid_;
  double dz_;
  std::vector<double> r_R_H_min;
        // vertical root resistance
    std::vector<double> r_R_V;

            // cumulative vertical sum of root resistance
    std::vector<double> r_R_V_sum;
    

  double leaf_temp_;
  double PPFD_;
  double atm_vpd_;
  double atm_o2_kpa_;
  double atm_kpa_;
  double ca_;
  double root_collar_psi_;
  double assim_max_;

  
  double opt_psi_stem_;
  double opt_ci_;
  double count;
  double E_up_;

  // 1-entry memo for transpiration(). Within a single root-collar/profit solve
  // leaf_specific_conductance_max_ and the transpiration spline are fixed, so
  // supply-side transpiration depends only on (psi_stem, psi_upstream). That
  // pair is queried repeatedly with identical values per profit evaluation
  // (psi_stem_to_ci -> stom_cond_CO2 -> transpiration, then transpiration again
  // for transpiration_). Caching the last result avoids the redundant spline
  // lookups; it is invalidated in set_physiology when conductance/soil change.
  bool   transpiration_cached_ = false;
  double transpiration_cache_psi_stem_ = 0.0;
  double transpiration_cache_psi_upstream_ = 0.0;
  double transpiration_cache_value_ = 0.0;
  std::vector<double> f_r;
  // TODO: move into environment?

  // TODO: atm_vpd - now set in set_physiology although ideally should be moved to enviroment
  double atm_vpd = 2.0; //kPa
  double ca = 40.0; // Pa
  double atm_kpa = 101.3; //kPa
  //partial pressure o2 (kPa)
  double atm_o2_kpa = 21;
  //leaf temperature (deg C)
  double leaf_temp = 25;
  // density of water 


  // this might end up hard-coded
  void initialize_integrator(int integration_rule = 21,
                             double integration_tol = 1e-3) {

    integrator = quadrature::QAG(integration_rule,
                                 1, // fixed integration
                                 integration_tol, integration_tol);
  }
  
  // set-up functions
  void set_physiology(double area_leaf, const std::vector<double>& mass_root_prop, double rho, double a_bio, double PPFD, const std::vector<double>& psi_soil, const std::vector<double>& soil_depth, double leaf_specific_conductance_max, double atm_vpd, double ca, double sapwood_volume_per_leaf_area, double leaf_temp, double atm_o2_kpa, double atm_kpa);
  void setup_transpiration(double resolution);
  void setup_root_vulnerability(double resolution);
  // Shared builder for the knot grid {0, step, .., <= psi_max} and the
  // cumulative vulnerability integral G(m) = int_0^m exp(-(s/b)^c) ds, seeded
  // from its gamma closed form. Used by both setup_* functions (see #468).
  void build_cumulative_vulnerability_integral(double b, double c,
                                               double resolution,
                                               std::vector<double>& x,
                                               std::vector<double>& y_integral);
  void setup_clean_leaf();
  // std::vector<double> root_collar_psi(std::vector<double> soil_moist_);

  void E_from_Soil_to_Root_Collar(double P_x_r, const std::vector<double>& psi_soil);
  void find_root_collar_psi();
  // Shut-down operating point used by the find_root_collar_psi early-exits: stem
  // held at psi_crit (no transpiration), paying only respiration + hydraulic
  // cost. Only root_collar_psi_ differs between the cases, so it is the argument.
  void set_shutdown_state(double root_collar);
  double find_root_psi(double wettest_soil_layer, const std::vector<double>& psi_soil, int find_root_crit);
  double find_psi_stem_from_psi_root(double psi_root, const std::vector<double>& psi_soil);
  double E_column(double x, const std::vector<double>& psi_soil, double psi_leaf);
  double E_column_zero(double x, const std::vector<double>& psi_soil);
  
  double arrh_curve(double Ea, double ref_value, double leaf_temp) const;
  double peak_arrh_curve(double Ea, double ref_value, double leaf_temp, double H_d, double d_S) const;

  // transpiration functions

  // proportion of conductivity in xylem at a given water potential (return: unitless)
  double proportion_of_conductivity(double psi) const;

  // supply-side transpiration for a given water potential gradient between leaves and soil, 
  // references setup_transpiraiton for values (return: kg h20 s^-1 m^-2 LA)
  // should be renamed to reflect supply-side
  double transpiration(double psi_stem, double psi_upstream);
  // supply-side transpiration for a given water potential gradient between leaves and soil, integrated internally (return: kg h20 s^-1 m^-2 LA)
  // should be renamed to reflect supply-side
  double transpiration_full_integration(double psi_stem, double psi_upstream);                    
  // stomatal conductance rate of c02 (return: mol CO2 m^-2 s^-1)
  double stom_cond_CO2(double psi_stem, double psi_upstream); // define as a constant
  // converts transpiration in kg h20 s^-1 m^-2 LA to psi_stem (return: -MPa)
  double transpiration_to_psi_stem(double transpiration_, double psi_upstream);
  
  // assimilation functions

  //
  double assim_rubisco_limited(double ci_);
  double electron_transport();
  double assim_electron_limited(double ci_);
  double assim_colimited(double ci_);
  double assim_minus_stom_cond_CO2(double x, double psi_stem, double psi_upstream);
  double psi_stem_to_ci(double psi_stem, double psi_upstream);
  void set_leaf_states_rates_from_psi_stem(double psi_stem, double psi_upstream);


// leaf economics functions
  double hydraulic_cost_Sperry(double psi_stem, double psi_upstream);
  double hydraulic_cost_TF(double psi_stem);

  double profit_psi_stem_Sperry(double psi_stem, double psi_upstream);
  double profit_psi_stem_TF(double psi_stem, double psi_upstream);

// optimiser functions
  void optimise_psi_stem_Sperry();
  void optimise_psi_stem_TF();

};
} // namespace plant
#endif
