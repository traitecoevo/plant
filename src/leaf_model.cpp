#include <plant/leaf_model.h>

namespace plant {
Leaf::Leaf():
ci(NA_REAL),
g_c(NA_REAL),
A_lim(NA_REAL),
E(NA_REAL),
psi(NA_REAL),
profit(NA_REAL)
{}

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

double Leaf::calc_g_c(double psi_soil, double psi_stem, double k_l_max, double p_50,
                      double c, double b, double atm_kpa,
                      const double kg_2_mol_h2o, double atm_vpd) {
  double E_supply = calc_E_supply(k_l_max, b, c, psi_soil, psi_stem);

  return atm_kpa * E_supply * kg_2_mol_h2o / atm_vpd / 1.6;
}

double Leaf::calc_A_c(double ci_, double vcmax, double gamma_25,
                      double umol_per_mol_2_Pa, double km_25) {
  return (vcmax * (ci_ - gamma_25 * umol_per_mol_2_Pa)) / (ci_ + km_25);
}

double Leaf::calc_A_j(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                      double curv_fact, double a, double gamma_25,
                      double umol_per_mol_2_Pa, double ci_) {

  double jmax = vcmax * vcmax_25_2_jmax_25;
  double j = (a * PPFD + jmax -
              sqrt(pow(a * PPFD + jmax, 2) - 4 * curv_fact * a * PPFD * jmax)) /
             (2 * curv_fact);  // check brackets are correct

  return j / 4 *
         ((ci_ - gamma_25 * umol_per_mol_2_Pa) /
          (ci_ + 2 * gamma_25 * umol_per_mol_2_Pa));
}

double Leaf::calc_A_lim(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                        double curv_fact, double a, double gamma_25,
                        double umol_per_mol_2_Pa, double ci_, double km_25) {

  double A_c = calc_A_c(ci_, vcmax, gamma_25, umol_per_mol_2_Pa, km_25);
  double A_j = calc_A_j(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25,
                        umol_per_mol_2_Pa, ci_);
  double R_d = vcmax*0.015;                      

  return (A_c + A_j - sqrt(pow(A_c + A_j, 2) - 4 * 0.98 * A_c * A_j)) / (2 * 0.98) - R_d;
}

double Leaf::diff_ci(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                     double curv_fact, double a, double gamma_25,
                     double umol_per_mol_2_Pa, double x, double km_25,
                     double psi_soil, double psi_stem, double k_l_max, double p_50,
                     double c, double b, const double kg_2_mol_h2o,
                     double umol_per_mol_2_mol_per_mol, double atm_vpd,
                     double ca, double atm_kpa, double kPa_2_Pa) {

  double A_lim_ = calc_A_lim(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a,
                             gamma_25, umol_per_mol_2_Pa, x, km_25);

  double g_c_ = calc_g_c(psi_soil, psi_stem, k_l_max, p_50, c,
                         b, atm_kpa, kg_2_mol_h2o, atm_vpd);

  return A_lim_ * umol_per_mol_2_mol_per_mol -
         (g_c_ * (ca - x) / (atm_kpa * kPa_2_Pa));
}

// need to fill in tol and max_iteratiosn
double Leaf::calc_assim_gross(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                         double curv_fact, double a, double gamma_25,
                         double umol_per_mol_2_Pa, double km_25,
                         double psi_soil, double psi_stem, double k_l_max,
                         double p_50, double c, double b,
                         const double kg_2_mol_h2o,
                         double umol_per_mol_2_mol_per_mol, double atm_vpd,
                         double ca, double atm_kpa, double kPa_2_Pa) {

  // not clear what x is here
  auto target = [&](double x) mutable -> double {
    return diff_ci(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25,
                   umol_per_mol_2_Pa, x, km_25, psi_soil,
                   psi_stem, k_l_max, p_50, c, b, kg_2_mol_h2o,
                   umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa,
                   kPa_2_Pa);
  };

  // tol and iterations copied from control defaults (for now)
  ci = util::uniroot(target, 0, ca, 1e-8, 1000);

  A_lim = calc_A_lim(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25,
                    umol_per_mol_2_Pa, ci, km_25);

  g_c = calc_g_c(psi_soil, psi_stem, k_l_max, p_50, c,
                         b, atm_kpa, kg_2_mol_h2o, atm_vpd);

  E = g_c * 1.6 * atm_vpd / kg_2_mol_h2o / atm_kpa;

  psi = psi_stem;                         

  return A_lim;
}

//Sperry et al. 2017; Sabot et al. 2020 implementation

double Leaf::calc_hydraulic_cost(double psi_soil, double psi_stem, double k_l_max, double b, double c) {
  double k_l_soil_ = calc_cond_vuln(psi_soil, k_l_max, b, c);
  double k_l_stem_ = calc_cond_vuln(psi_stem, k_l_max, b, c);

  return k_l_soil_ - k_l_stem_;
}

double Leaf::calc_profit(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double km_25, 
double psi_soil, double psi_stem, double k_l_max, double p_50, double c, double b, double kg_2_mol_h2o, 
double umol_per_mol_2_mol_per_mol, double atm_vpd, double ca, double atm_kpa, double kPa_2_Pa, double psi_crit) {
  
  double lambda_ = calc_assim_gross(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, psi_soil, psi_crit, k_l_max, p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa)/
    (calc_cond_vuln(psi_soil, k_l_max, b, c) - calc_cond_vuln(psi_crit, k_l_max, b, c)) ;
  double benefit_ = calc_assim_gross(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, psi_soil, psi_stem, k_l_max, p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa);
  double cost_ = calc_hydraulic_cost(psi_soil, psi_stem, k_l_max, b, c);
  
  return benefit_ - lambda_*cost_;
}

double Leaf::optimise_profit_gss(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double km_25, 
  double psi_soil, double k_l_max, double p_50, double c, double b, double kg_2_mol_h2o, double umol_per_mol_2_mol_per_mol, 
  double atm_vpd, double ca, double atm_kpa, double kPa_2_Pa, double psi_crit) {
  
  double gr = (sqrt(5) + 1)/2;

  profit = calc_profit(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, psi_soil, psi_soil, k_l_max, 
  p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit);

  double bound_a = psi_soil; 
  double bound_b = psi_crit;

  double opt_psi_stem = psi_soil;

  if (psi_soil < psi_crit) {
    double delta_crit = 1e-5;
    
    double bound_c = bound_b - (bound_b - bound_a)/gr;
    double bound_d = bound_a + (bound_b - bound_a)/gr;
  
    while (abs(bound_b - bound_a) > delta_crit) {
    
      if(calc_profit(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, psi_soil, bound_c, k_l_max, 
      p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit) >
      calc_profit(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, psi_soil, bound_d, k_l_max, 
      p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit)) {
        bound_b = bound_d;
      } else {
        bound_a = bound_c;
      }

     bound_c = bound_b - (bound_b - bound_a)/gr;
     bound_d = bound_a + (bound_b - bound_a)/gr;


    }
    opt_psi_stem  = ((bound_b + bound_a)/2);
    
  profit = calc_profit(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25,psi_soil, opt_psi_stem, k_l_max, 
  p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit) ;

  }
  return(profit);
}


//Bartlett et al. implementation


double Leaf::calc_hydraulic_cost_bartlett(double psi_soil, double psi_stem, double k_l_max, double b, double c, double beta, double beta_2, double huber_value, double height) {
  double k_l_soil_ = calc_cond_vuln(psi_soil, k_l_max, b, c);
  double k_l_stem_ = calc_cond_vuln(psi_stem, k_l_max, b, c);

  return beta*huber_value*height*pow((1-k_l_stem_/k_l_soil_), beta_2);
}

double Leaf::calc_profit_bartlett(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double km_25, 
double psi_soil, double psi_stem, double k_l_max, double p_50, double c, double b, double kg_2_mol_h2o, 
double umol_per_mol_2_mol_per_mol, double atm_vpd, double ca, double atm_kpa, double kPa_2_Pa, double psi_crit, double beta, double beta_2, double huber_value, double height) {
  
  double benefit_ = calc_assim_gross(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, psi_soil, psi_stem, k_l_max, p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa);
  double cost_ = calc_hydraulic_cost_bartlett(psi_soil, psi_stem, k_l_max, b, c, beta, beta_2, huber_value, height);
  
  return benefit_ - cost_;
}

double Leaf::optimise_profit_gss_bartlett(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double km_25, 
  double psi_soil, double k_l_max, double p_50, double c, double b, double kg_2_mol_h2o, double umol_per_mol_2_mol_per_mol, 
  double atm_vpd, double ca, double atm_kpa, double kPa_2_Pa, double psi_crit, double beta, double beta_2, double huber_value, double height) {
  
  double gr = (sqrt(5) + 1)/2;

  profit = calc_profit_bartlett(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, psi_soil, psi_soil, k_l_max, 
  p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit, beta, beta_2, huber_value, height);

  double bound_a = psi_soil; 
  double bound_b = psi_crit;

  double opt_psi_stem = psi_soil;

  if (psi_soil < psi_crit) {
    double delta_crit = 1e-5;
    
    double bound_c = bound_b - (bound_b - bound_a)/gr;
    double bound_d = bound_a + (bound_b - bound_a)/gr;
  
    while (abs(bound_b - bound_a) > delta_crit) {
    
      if(calc_profit_bartlett(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, psi_soil, bound_c, k_l_max, 
      p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit, beta, beta_2, huber_value, height) >
      calc_profit_bartlett(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, psi_soil, bound_d, k_l_max, 
      p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit, beta, beta_2, huber_value, height)) {
        bound_b = bound_d;
      } else {
        bound_a = bound_c;
      }

     bound_c = bound_b - (bound_b - bound_a)/gr;
     bound_d = bound_a + (bound_b - bound_a)/gr;


    }
    opt_psi_stem  = ((bound_b + bound_a)/2);
    
  profit = calc_profit_bartlett(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25,psi_soil, opt_psi_stem, k_l_max, 
  p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit, beta, beta_2, huber_value, height) ;

  }
  return(profit);
}





}
// namespace plant
