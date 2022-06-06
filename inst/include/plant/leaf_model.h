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
static const double vcmax_25_to_jmax_25 = 1.67;
static const double curv_fact = 0.9;
static const double a = 0.3;
static const double gamma_25 = 42.75;
static const double umol_per_mol_to_Pa = 0.1013;
static const double km_25 = 71.56;
static const double kg_to_mol_h2o = 55.4939;
static const double umol_per_mol_to_mol_per_mol = 1e-6;
static const double kPa_to_Pa = 1000.0;

class Leaf {
public:
  Leaf(double vcmax        = 100,
       double p_50         = 1.0,
       double c            = 2.04,
       double b            = 2.0,
       double psi_crit     = 3.42,  // derived from b and c
       double beta         = 15000.0,
       double beta_2       = 1.0,
       double huber_value  = 1.57e-4,
       double K_s = 2);

  quadrature::QAG integrator;

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
  double atm_vpd = 2.0;
  double ca = 40.0;
  double atm_kpa = 101.3;

  // this might end up hard-coded
  void initialize_integrator(int integration_rule = 21,
                             double integration_tol = 1e-6) {

    integrator = quadrature::QAG(integration_rule,
                                 1, // fixed integration
                                 integration_tol, integration_tol);
  }

  double p_50_to_K_s() const;

  // double calc_k_l_max(double K_s, double h_v, double h) const;

  double calc_vul_b() const;
  double calc_cond_vuln(double psi, double k_l_max) const;

  double calc_E_supply(double k_l_max, double psi_soil,
                       double psi_stem);

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
