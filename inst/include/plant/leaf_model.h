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
// unitless
static const double vcmax_25_to_jmax_25 = 1.67;
// unitless - obtained from Sperry et al. (2017)
static const double curv_fact = 0.9;
// unitless
static const double a = 0.3;
// umol ^ -1 mol ^ 1
static const double gamma_25 = 42.75;
// Pa umol ^ -1 mol ^ 1 
static const double umol_per_mol_to_Pa = 0.1013;
// umol mol ^-1
static const double kc_25 = 404.9 ;
//partial pressure o2 (kPa)
static const double atm_o2_kpa = 21 ;
//multiplicative converter from kPa to Pa (Pa kPa^-1)
static const double kPa_2_Pa = 1000;
//umol mol ^-1
static const double ko_25 = 278400;

// Pa
static const double km_25 = (kc_25*umol_per_mol_to_Pa)*(1 + (atm_o2_kpa*kPa_2_Pa)/(ko_25*umol_per_mol_to_Pa));
// mol H2o kg ^-1
static const double kg_to_mol_h2o = 55.4939;
// mol mol ^-1 / (umol mol ^-1)
static const double umol_to_mol = 1e-6;
// Pa kPa^-1
static const double kPa_to_Pa = 1000.0;


class Leaf {
public:
  Leaf(double vcmax        = 100, // umol m^-2 s^-1
       double p_50         = 1.0, // (- MPa)
       double c            = 2.04, //unitless
       double b            = 2.0, // MPa
       double psi_crit     = 3.42,  // derived from b and c (- MPa)
       double huber_value  = 1.57e-4, // m^2 sapwood area m^-2 leaf area
       double K_s = 2, // kg m^-1 s^-1 MPa ^-1 Liu et al. 2010 
       double epsilon_leaf = 0.001); 

  quadrature::QAG integrator;
  interpolator::Interpolator E_from_psi;
  interpolator::Interpolator psi_from_E;

  // psi_from_E

  double vcmax;
  double p_50;
  double c;
  double b;
  double psi_crit;  // derived from b and c
  double huber_value;
  double K_s;

  //actually a control paramaeter and needs to be moved
  double epsilon_leaf;


  double ci;
  double j_;
  double g_c;
  double A_lim_;
  double E;
  double profit;
  double psi_stem;
  double lambda_;
  double lambda_analytical_;
  double PPFD_;
  double atm_vpd_;
  double ca_;
  double k_l_max_;
  double psi_soil_;
  double opt_psi_stem;
  double opt_ci;
  double method;
  double count;
  double GSS_count;


  // TODO: move into environment?

  // TODO: atm_vpd - now set in set_physiology although ideally should be moved to enviroment
  double atm_vpd = 2.0; //kPa
  double ca = 40.0; // Pa
  double atm_kpa = 101.3; //kPa

  // this might end up hard-coded
  void initialize_integrator(int integration_rule = 21,
                             double integration_tol = 1e-3) {

    integrator = quadrature::QAG(integration_rule,
                                 1, // fixed integration
                                 integration_tol, integration_tol);
  }
  
  void set_physiology(double PPFD, double psi_soil, double k_l_max, double atm_vpd, double ca);
  void setup_E_supply(double resolution);


  double calc_cond_vuln(double psi) const;
  double calc_E_supply(double psi_stem);
  double calc_E_supply_full_integration(double psi_stem);                    

  double calc_g_c(double psi_stem); // define as a constant
  double calc_A_c(double ci_);
  double calc_j();
  double calc_A_j(double ci_);
  double A_lim(double ci_);
  double A_lim_analytical(double c_i);
  
  double diff_ci(double x, double psi_stem);
  double convert_psi_stem_to_ci_analytical(double psi_stem);
  void get_leaf_states_rates_from_psi_stem_analytical(double psi_stem);
  double convert_psi_stem_to_ci(double psi_stem);
  void get_leaf_states_rates_from_psi_stem(double psi_stem);
  double convert_E_from_ci_to_psi_stem(double E_ci);

  double calc_hydraulic_cost_Sperry(double psi_stem);

  double profit_psi_stem_Sperry(double psi_stem);
  double profit_psi_stem_Sperry_analytical(double psi_stem);
  double calc_profit_Sperry_ci(double c_i);
  double calc_profit_Sperry_ci_analytical(double c_i);

  void optimise_psi_stem_Sperry();
  void optimise_psi_stem_Sperry_analytical();
  void optimise_psi_stem_Sperry_Newton(double psi_guess);
  void optimise_psi_stem_Sperry_Newton_analytical(double psi_guess);
  void optimise_ci_Sperry(double ci_guess);
  void optimise_ci_Sperry_analytical(double max_ci);
  void optimise_ci_Sperry_Newton(double ci_guess);
  void optimise_ci_Sperry_Newton_analytical(double ci_guess);

  
};
} // namespace plant
#endif
