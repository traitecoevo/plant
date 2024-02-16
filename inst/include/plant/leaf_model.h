// -*-c++-*-
#ifndef PLANT_PLANT_LEAF_MODEL_H_
#define PLANT_PLANT_LEAF_MODEL_H_

// TODO: replace with constants
// #define umol_per_mol_to_mol_per_mol 0.000...
// #define umol_per_mol_to_Pa ...
// #define kg_to_mol_h2o ...
// #define kPa_to_Pa ...

#include <plant/models/ff16_environment.h>
#include <plant/qag.h>
#include <plant/uniroot.h>

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

class Leaf {
public:
  //anonymous Leaf function as in canopy.h
  Leaf();
  
  Leaf(double vcmax_25, 
       double c, 
       double b, 
       double psi_crit,
       double beta1,
       double beta2, 
       double jmax_25, 
       double hk_s,
       double a, 
       double curv_fact_elec_trans, 
       double curv_fact_colim,
       double GSS_tol_abs,
       double vulnerability_curve_ncontrol,
       double ci_abs_tol,
       double ci_niter); 
        
  quadrature::QAG integrator;
  interpolator::Interpolator transpiration_from_psi;
  interpolator::Interpolator psi_from_transpiration;

  // psi_from_E

  double vcmax_25;
  double c;
  double b;
  double psi_crit;  // derived from b and c
  double beta1;
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

  double ci_;
  double stom_cond_CO2_;
  double assim_colimited_;
  double transpiration_;
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
  double rho_;
  double vcmax_;
  double jmax_;
  double lma_; //kg m^-2
  double a_bio_;
  
  double psi_soil_;
  double leaf_temp_;
  double PPFD_;
  double atm_vpd_;
  double atm_o2_kpa_;
  double atm_kpa_;
  double ca_;
  
  double opt_psi_stem_;
  double opt_ci_;
  double count;

  // TODO: move into environment?

  // TODO: atm_vpd - now set in set_physiology although ideally should be moved to enviroment
  double atm_vpd = 2.0; //kPa
  double ca = 40.0; // Pa
  double atm_kpa = 101.3; //kPa
  //partial pressure o2 (kPa)
  double atm_o2_kpa = 21;
  //leaf temperature (deg C)
  double leaf_temp = 25;

  // this might end up hard-coded
  void initialize_integrator(int integration_rule = 21,
                             double integration_tol = 1e-3) {

    integrator = quadrature::QAG(integration_rule,
                                 1, // fixed integration
                                 integration_tol, integration_tol);
  }
  
  // set-up functions
  void set_physiology(double rho, double a_bio, double PPFD, double psi_soil, double leaf_specific_conductance_max, double atm_vpd, double ca, double sapwood_volume_per_leaf_area, double leaf_temp, double atm_o2_kpa, double atm_kpa);
  void setup_transpiration(double resolution);
  void setup_clean_leaf();

  double arrh_curve(double Ea, double ref_value, double leaf_temp) const;
  double peak_arrh_curve(double Ea, double ref_value, double leaf_temp, double H_d, double d_S) const;

  // transpiration functions

  // proportion of conductivity in xylem at a given water potential (return: unitless)
  double proportion_of_conductivity(double psi) const;

  // supply-side transpiration for a given water potential gradient between leaves and soil, 
  // references setup_transpiraiton for values (return: kg h20 s^-1 m^-2 LA)
  // should be renamed to reflect supply-side
  double transpiration(double psi_stem);
  // supply-side transpiration for a given water potential gradient between leaves and soil, integrated internally (return: kg h20 s^-1 m^-2 LA)
  // should be renamed to reflect supply-side
  double transpiration_full_integration(double psi_stem);                    
  // stomatal conductance rate of c02 (return: mol CO2 m^-2 s^-1)
  double stom_cond_CO2(double psi_stem); // define as a constant
  // converts transpiration in kg h20 s^-1 m^-2 LA to psi_stem (return: -MPa)
  double transpiration_to_psi_stem(double transpiration_);
  
  // assimilation functions

  //
  double assim_rubisco_limited(double ci_);
  double electron_transport();
  double assim_electron_limited(double ci_);
  double assim_colimited(double ci_);
  double assim_minus_stom_cond_CO2(double x, double psi_stem);
  double psi_stem_to_ci(double psi_stem);
  void set_leaf_states_rates_from_psi_stem(double psi_stem);


// leaf economics functions
  double hydraulic_cost_Sperry(double psi_stem);
  double hydraulic_cost_TF(double psi_stem);

  double profit_psi_stem_Sperry(double psi_stem);
  double profit_Sperry_ci(double ci_);
  double profit_psi_stem_TF(double psi_stem);

// optimiser functions
  void optimise_psi_stem_Sperry();
  void optimise_ci_Sperry(double ci_guess);
  void optimise_psi_stem_TF();

};
} // namespace plant
#endif
