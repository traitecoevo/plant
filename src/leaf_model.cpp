#include <plant/leaf_model.h>

namespace plant {
Leaf::Leaf() {}

// TODO: vectorise for soil layers
double Leaf::calc_soil_water_pot(const FF16_Environment &environment,
                                 double psi_aep, double b_CH) {
  double rel_theta =
      environment.get_soil_water_state()[0]; // this hack assumes one soil layer
  return psi_aep * pow(rel_theta, -b_CH);
}

double Leaf::calc_k_l_max(double K_s, double h_v, double h) const {
  return K_s * h_v / h;
}

double Leaf::calc_vul_b(double p_50, double c) const {
  return p_50/pow(-log(1 - 50.0/100.0),1/c);
}

// integrates
double Leaf::calc_cond_vuln(double psi, double k_l_max, double b,
                            double c) const {
  return k_l_max * exp(-pow((psi / b), c));
}

// replace f with some other function
double Leaf::calc_E_supply(double k_l_max, double b, double c, double psi_soil,
                           double psi_stem) {
  std::function<double(double)> f;
  f = [&](double psi) -> double { return calc_cond_vuln(psi, k_l_max, b, c); };

  return integrator.integrate(f, psi_soil, psi_stem);
}

double Leaf::calc_g_c(const FF16_Environment &environment, double psi_aep,
                      double b_CH, double psi_stem, double k_l_max, double p_50,
                      double c, double b, double atm_kpa,
                      const double kg_2_mol_h2o, double atm_vpd) {
  double psi_soil = calc_soil_water_pot(environment, psi_aep, b_CH);
  double E_supply = calc_E_supply(k_l_max, c, b, psi_soil, psi_stem);

  return atm_kpa * E_supply * kg_2_mol_h2o / atm_vpd / 1.6;
}

double Leaf::calc_A_c(double c_i, double vcmax, double gamma_25,
                      double umol_per_mol_2_Pa, double km_25) {
  return (vcmax * (c_i - gamma_25 * umol_per_mol_2_Pa)) / (c_i + km_25);
}

double Leaf::calc_A_j(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                      double curv_fact, double a, double gamma_25,
                      double umol_per_mol_2_Pa, double c_i) {

  double jmax = vcmax * vcmax_25_2_jmax_25;
  double j = (a * PPFD + jmax -
              sqrt(pow(a * PPFD + jmax, 2) - 4 * curv_fact * a * PPFD * jmax)) /
             (2 * curv_fact);  // check brackets are correct

  return j / 4 *
         ((c_i - gamma_25 * umol_per_mol_2_Pa) /
          (c_i + 2 * gamma_25 * umol_per_mol_2_Pa));
}

double Leaf::calc_A_lim(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                        double curv_fact, double a, double gamma_25,
                        double umol_per_mol_2_Pa, double c_i, double km_25) {

  double A_c = calc_A_c(c_i, vcmax, gamma_25, umol_per_mol_2_Pa, km_25);
  double A_j = calc_A_j(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25,
                        umol_per_mol_2_Pa, c_i);

  return (A_c + A_j - sqrt(pow(A_c + A_j, 2) - 4 * 0.98 * A_c * A_j)) / (2 * 0.98);
}

double Leaf::diff_ci(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                     double curv_fact, double a, double gamma_25,
                     double umol_per_mol_2_Pa, double x, double km_25,
                     const FF16_Environment &environment, double psi_aep,
                     double b_CH, double psi_stem, double k_l_max, double p_50,
                     double c, double b, const double kg_2_mol_h2o,
                     double umol_per_mol_2_mol_per_mol, double atm_vpd,
                     double ca, double atm_kpa, double kPa_2_Pa) {

  double A_lim_ = calc_A_lim(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a,
                             gamma_25, umol_per_mol_2_Pa, x, km_25);

  double g_c_ = calc_g_c(environment, psi_aep, b_CH, psi_stem, k_l_max, p_50, c,
                         b, atm_kpa, kg_2_mol_h2o, atm_vpd);

  return A_lim_ * umol_per_mol_2_mol_per_mol -
         (g_c_ * (ca - x) / (atm_kpa * kPa_2_Pa));
}

// need to fill in tol and max_iteratiosn
double Leaf::solve_A_lim(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                         double curv_fact, double a, double gamma_25,
                         double umol_per_mol_2_Pa, double km_25,
                         const FF16_Environment &environment, double psi_aep,
                         double b_CH, double psi_stem, double k_l_max,
                         double p_50, double c, double b,
                         const double kg_2_mol_h2o,
                         double umol_per_mol_2_mol_per_mol, double atm_vpd,
                         double ca, double atm_kpa, double kPa_2_Pa, int i) {

  // not clear what x is here
  auto target = [&](double x) mutable -> double {
    return diff_ci(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25,
                   umol_per_mol_2_Pa, x, km_25, environment, psi_aep, b_CH,
                   psi_stem, k_l_max, p_50, c, b, kg_2_mol_h2o,
                   umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa,
                   kPa_2_Pa);
  };

  // tol and iterations copied from control defaults (for now)
  ci[i] = util::uniroot(target, 0, 40, 1e-6, 1000);

  A_lim[i] = calc_A_lim(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25,
                    umol_per_mol_2_Pa, ci[i], km_25);

  g_c[i] = calc_g_c(environment, psi_aep, b_CH, psi_stem, k_l_max, p_50, c,
                         b, atm_kpa, kg_2_mol_h2o, atm_vpd);
                         

  return A_lim[i];
}

double Leaf::calc_hydraulic_cost(const FF16_Environment &environment, double psi_stem, double k_l_max, double b, double c, double psi_aep, double b_CH) {
  double psi_soil_ = calc_soil_water_pot(environment, psi_aep, b_CH);
  double k_l_soil_ = calc_cond_vuln(psi_soil_, k_l_max, b, c);
  double k_l_stem_ = calc_cond_vuln(psi_stem, k_l_max, b, c);

  return k_l_soil_ - k_l_stem_;
}

double Leaf::calc_lambda(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double km_25, const FF16_Environment &environment, double psi_aep, double b_CH, double psi_crit, double k_l_max, double p_50, double c, double b, double kg_2_mol_h2o, double umol_per_mol_2_mol_per_mol, double atm_vpd, double ca, double atm_kpa, double kPa_2_Pa) {
  double psi_soil_ = calc_soil_water_pot(environment, psi_aep, b_CH);
  double assim_gross_max_ = solve_A_lim(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, environment, psi_aep, b_CH, psi_crit, k_l_max, p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa);
  double max_diff_cond_ = calc_cond_vuln(psi_soil_, k_l_max, b, c) - k_l_max * 0.05;

  return assim_gross_max_ / max_diff_cond_;
}

double Leaf::calc_assim_gross(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double km_25, const FF16_Environment &environment, double psi_aep, double b_CH, double psi_stem, double k_l_max, double p_50, double c, double b, double kg_2_mol_h2o, double umol_per_mol_2_mol_per_mol, double atm_vpd, double ca, double atm_kpa, double kPa_2_Pa) {
  return solve_A_lim(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, environment, psi_aep, b_CH, psi_stem, k_l_max, p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa);
}

double Leaf::calc_profit(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double km_25, const FF16_Environment &environment, double psi_aep, double b_CH, double psi_stem, double k_l_max, double p_50, double c, double b, double kg_2_mol_h2o, double umol_per_mol_2_mol_per_mol, double atm_vpd, double ca, double atm_kpa, double kPa_2_Pa, double psi_crit, int i) {
  double psi_soil_ = calc_soil_water_pot(environment, psi_aep, b_CH);
  
  double benefit_ = calc_assim_gross(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, environment, psi_aep, b_CH, psi_stem, k_l_max, p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, i);
  double lambda_ = calc_lambda(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, environment, psi_aep, b_CH, psi_crit, k_l_max, p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa);
  double cost_ = calc_hydraulic_cost(environment, psi_stem, k_l_max, b, c, psi_aep, b_CH);
  
  return benefit_ - lambda_*cost_;
}

double Leaf::optimise_profit(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double km_25, 
  const FF16_Environment &environment, double psi_aep, double b_CH, double k_l_max, double p_50, double c, double b, double kg_2_mol_h2o, double umol_per_mol_2_mol_per_mol, 
  double atm_vpd, double ca, double atm_kpa, double kPa_2_Pa, double psi_crit) {
  
  double psi_soil_ = calc_soil_water_pot(environment, psi_aep, b_CH);
  double psi_curr = psi_soil_;

  if (psi_soil_ >= psi_crit) {
    return(psi_curr);
  } else {
   
  double delta = 0.01;
  double profit_curr = calc_profit(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, environment, psi_aep, b_CH, psi_curr, k_l_max, p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit, 1) ;
  double delta_crit = 1e-5;
  double psi_next = psi_soil_;
  double profit_next = psi_soil_;
  
  while (delta > delta_crit) {
    psi_next = psi_curr + delta;
    double psi_next2 = psi_curr + 2*delta;

    profit[2]
    profit_next = calc_profit(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, environment, psi_aep, b_CH, psi_next, k_l_max, 
    p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit, 2) ;

    double profit_next2 = calc_profit(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, environment, psi_aep, b_CH, psi_next2, k_l_max, 
    p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit, 3) ;
  
  
    if (profit_next2 > profit_curr){
      profit_curr = profit_next;
      profit_next = profit_next2;

      psi_curr = psi_next;
      psi_curr = psi_next;
      
      psi[1] = psi[2];
      A_lim[1] = A_lim[2];
      g_c[1] = g_c[2];
      ci[1] = c1[2];

    
    } else {
      delta = delta/2;
    }
  }
  return(profit_curr);
  }
  }
} // namespace plant
