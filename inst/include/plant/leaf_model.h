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
// unitless
static const double curv_fact = 0.9;
// unitless
static const double a = 0.3;
// Pa
static const double gamma_25 = 42.75;
// Pa umol ^ -1 mol ^ 1 
static const double umol_per_mol_to_Pa = 0.1013;
// Pa
static const double km_25 = 71.56;
// mol H2o kg ^-1
static const double kg_to_mol_h2o = 55.4939;
// mol mol ^-1 / (umol mol ^-1)
static const double umol_per_mol_to_mol_per_mol = 1e-6;
// Pa kPa^-1
static const double kPa_to_Pa = 1000.0;


class Leaf {
public:
  Leaf(double vcmax        = 100, // umol m^-2 s^-1
       double p_50         = 1.0, // MPa
       double c            = 2.04, //unitless
       double b            = 2.0, // MPa
       double psi_crit     = 3.42,  // derived from b and c
       double beta         = 15000.0, // umol m^-2 s^-1
       double beta_2       = 1.0, //unitless
       double huber_value  = 1.57e-4, // m^2 sapwood area m^-2 leaf area
       double K_s = 2); // kg m^-1 s^-1 MPa ^-1 Liu et al. 2010 

  quadrature::QAG integrator;
  interpolator::Interpolator E_curve;

  double vcmax;
  double p_50;
  double c;
  double b;
  double psi_crit;  // derived from b and c
  double beta;
  double beta_2;
  double huber_value;
  double K_s;

  double ci;
  double g_c;
  double A_lim;
  double E;
  double psi;
  double profit;

  // TODO: move into environment?
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

  double p_50_to_K_s() const;

  // double calc_k_l_max(double K_s, double h_v, double h) const;

  double calc_vul_b() const;
  double calc_cond_vuln(double psi) const;

  double calc_E_supply(double k_l_max, double psi_soil,
                       double psi_stem);

  double setup_E_supply(double resolution);

  double calc_g_c(double psi_soil, double psi_stem, double k_l_max); // define as a constant

  double calc_A_c(double ci_);
  double calc_A_j(double PPFD, double ci_);
  double calc_A_lim(double PPFD, double ci_);

  double diff_ci(double PPFD, double x, double psi_soil, double psi_stem, double k_l_max);

  double calc_assim_gross(double PPFD, double psi_soil, double psi_stem, double k_l_max);

  double calc_hydraulic_cost_Sperry(double psi_soil, double psi_stem, double k_l_max);

  double calc_profit_Sperry(double PPFD, double psi_soil, double psi_stem, double k_l_max);

  double calc_hydraulic_cost_Bartlett(double psi_soil, double psi_stem,
                                      double k_l_max);

  double calc_profit_Bartlett(double PPFD, double psi_soil,
                              double psi_stem, double k_l_max);

  double optimise_psi_stem_Bartlett(double PPFD, double psi_soil, double k_l_max);

  double optimise_psi_stem_Sperry(double PPFD, double psi_soil, double k_l_max);

};

} // namespace plant
#endif
