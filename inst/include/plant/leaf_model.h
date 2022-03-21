// -*-c++-*-
#ifndef PLANT_PLANT_LEAF_MODEL_H_
#define PLANT_PLANT_LEAF_MODEL_H_

#include <Rcpp.h>
#include <vector>
#include <plant/util.h>
#include <plant/qag.h>
#include <plant/models/ff16_environment.h>
#include <plant/leaf_model.h>

namespace plant {

class Leaf {
public:
  Leaf();

  void initialize(double a1, double a2, double e,
                  bool adaptive_integration=false,
                  int integration_rule=21,
                  int iterations=1000,
                  double integration_tol=1e-6) {

    // strategy parameters
    a_p1 = a1;
    a_p2 = a2;
    eta = e;

    // set up integrator
    if(!adaptive_integration) {
      iterations = 1;
    }

    integrator = quadrature::QAG(integration_rule,
                                 iterations,
                                 integration_tol,
                                 integration_tol);
  }

  quadrature::QAG integrator;

  double calc_soil_water_pot(const FF16_Environment &environment, double psi_aep, double b_CH);
  double calc_k_l_max(double K_s, double h_v, double h) const;
  double calc_vul_b(double p_50, double c) const;
  double calc_cond_vuln(double k_l_max, double b, double c, double psi_stem) const;
  double calc_g_c(double psi_stem, double psi_soil, double atm_kpa, double atm_vpd, double psi_aep, double b_CH, const double kg_2_mol_h2o); // define as a constant
  double calc_A_c(double c_i, double vcmax, double gamma_25, double umol_per_mol_2_Pa, double km_25);
  double calc_A_j(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double c_i);
  double calc_A_lim(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double c_i);
  double diff_ci(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double c_i,
                 double psi_stem, double psi_soil, double atm_vpd, double psi_aep, double b_CH, double kg_2_mol_h2o, double ca, double x, double atm_kpa, double kPa_2_Pa);
  double solve_A_lim(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double c_i);
};

}
#endif
