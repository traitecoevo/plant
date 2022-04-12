// -*-c++-*-
#ifndef PLANT_PLANT_LEAF_MODEL_H_
#define PLANT_PLANT_LEAF_MODEL_H_

// TODO: replace with constants
// #define umol_per_mol_2_mol_per_mol 0.000...
// #define umol_per_mol_2_Pa ...
// #define kg_2_mol_h2o ...
// #define kPa_2_Pa ...

#include <plant/models/ff16_environment.h>
#include <plant/qag.h>
#include <plant/uniroot.h>

namespace plant {

class Leaf {
public:
  Leaf();

  quadrature::QAG integrator;

  std::vector<double> ci;
  double gc[3];
  double A_lim[3];
  double E[3];
  double psi[3];
  double profit[3];


  // this might end up hard-coded
  void initialize_integrator(int integration_rule = 21, double integration_tol = 1e-6) {

    integrator = quadrature::QAG(integration_rule,
                                 1, // fixed integration
                                 integration_tol, integration_tol);
  }

  double calc_soil_water_pot(const FF16_Environment &environment,
                             double psi_aep, double b_CH);
  double calc_k_l_max(double K_s, double h_v, double h) const;

  double calc_vul_b(double p_50, double c) const;
  double calc_cond_vuln(double psi, double k_l_max, double b, double c) const;

  double calc_E_supply(double k_l_max, double b, double c, double psi_soil,
                       double psi_stem);

  double calc_g_c(const FF16_Environment &environment, double psi_aep,
                  double b_CH, double psi_stem, double k_l_max, double p_50,
                  double c, double b, double atm_kpa, const double kg_2_mol_h2o,
                  double atm_vpd); // define as a constant

  double calc_A_c(double c_i, double vcmax, double gamma_25,
                  double umol_per_mol_2_Pa, double km_25);
  double calc_A_j(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                  double curv_fact, double a, double gamma_25,
                  double umol_per_mol_2_Pa, double c_i);
  double calc_A_lim(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                    double curv_fact, double a, double gamma_25,
                    double umol_per_mol_2_Pa, double c_i, double km_25);

  double diff_ci(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                     double curv_fact, double a, double gamma_25,
                     double umol_per_mol_2_Pa, double x, double km_25,
                     const FF16_Environment &environment, double psi_aep,
                     double b_CH, double psi_stem, double k_l_max, double p_50,
                     double c, double b, const double kg_2_mol_h2o,
                     double umol_per_mol_2_mol_per_mol, double atm_vpd,
                     double ca, double atm_kpa, double kPa_2_Pa);
    
  double solve_A_lim(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                         double curv_fact, double a, double gamma_25,
                         double umol_per_mol_2_Pa, double km_25,
                         const FF16_Environment &environment, double psi_aep,
                         double b_CH, double psi_stem, double k_l_max,
                         double p_50, double c, double b,
                         const double kg_2_mol_h2o,
                         double umol_per_mol_2_mol_per_mol, double atm_vpd,
                         double ca, double atm_kpa, double kPa_2_Pa);
  double calc_hydraulic_cost(const FF16_Environment &environment, double psi_stem, double k_l_max, double b, double c, double psi_aep, double b_CH);
  double calc_lambda(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                         double curv_fact, double a, double gamma_25,
                         double umol_per_mol_2_Pa, double km_25,
                         const FF16_Environment &environment, double psi_aep,
                         double b_CH, double psi_crit, double k_l_max,
                         double p_50, double c, double b,
                         const double kg_2_mol_h2o,
                         double umol_per_mol_2_mol_per_mol, double atm_vpd,
                         double ca, double atm_kpa, double kPa_2_Pa);
  double calc_assim_gross(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                         double curv_fact, double a, double gamma_25,
                         double umol_per_mol_2_Pa, double km_25,
                         const FF16_Environment &environment, double psi_aep,
                         double b_CH, double psi_stem, double k_l_max,
                         double p_50, double c, double b,
                         const double kg_2_mol_h2o,
                         double umol_per_mol_2_mol_per_mol, double atm_vpd,
                         double ca, double atm_kpa, double kPa_2_Pa);
double calc_profit(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                         double curv_fact, double a, double gamma_25,
                         double umol_per_mol_2_Pa, double km_25,
                         const FF16_Environment &environment, double psi_aep,
                         double b_CH, double psi_stem, double k_l_max,
                         double p_50, double c, double b,
                         const double kg_2_mol_h2o,
                         double umol_per_mol_2_mol_per_mol, double atm_vpd,
                         double ca, double atm_kpa, double kPa_2_Pa, double psi_crit);                                                 
double optimise_profit(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                         double curv_fact, double a, double gamma_25,
                         double umol_per_mol_2_Pa, double km_25,
                         const FF16_Environment &environment, double psi_aep,
                         double b_CH, double k_l_max,
                         double p_50, double c, double b,
                         const double kg_2_mol_h2o,
                         double umol_per_mol_2_mol_per_mol, double atm_vpd,
                         double ca, double atm_kpa, double kPa_2_Pa, double psi_crit); 
};

} // namespace plant
#endif
