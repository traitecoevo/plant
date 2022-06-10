#include <plant/leaf_model.h>

namespace plant {
Leaf::Leaf(double vcmax, double p_50, double c, double b,
           double psi_crit, // derived from b and c
           double beta, double beta_2, double huber_value, double K_s)
    : vcmax(vcmax), p_50(p_50), c(c), b(b), psi_crit(psi_crit), beta(beta),
      beta_2(beta_2), huber_value(huber_value), K_s(K_s), ci(NA_REAL), g_c(NA_REAL),
      A_lim(NA_REAL), E(NA_REAL), psi(NA_REAL), profit(NA_REAL) {}

double Leaf::p_50_to_K_s() const {
  return std::pow(1.638999 * (p_50 / 2), -1.38);
}

// double Leaf::calc_k_l_max(double K_s, double h_v, double h) const {
//   return K_s * h_v / h;
// }

double Leaf::calc_vul_b() const {
  return p_50 / pow(-log(1 - 50.0 / 100.0), 1 / c);
}

// integrates
double Leaf::calc_cond_vuln(double psi, double k_l_max) const {
  return k_l_max * exp(-pow((psi / b), c));
}

// replace f with some other function
double Leaf::calc_E_supply(double k_l_max, double psi_soil,
                           double psi_stem) {
  std::function<double(double)> f;
  f = [&](double psi) -> double { return calc_cond_vuln(psi, k_l_max); };

  return integrator.integrate(f, psi_soil, psi_stem);
}

double Leaf::calc_g_c(double psi_soil, double psi_stem, double k_l_max) {
  double E_supply = calc_E_supply(k_l_max, psi_soil, psi_stem);

  return atm_kpa * E_supply * kg_to_mol_h2o / atm_vpd / 1.6;
}

double Leaf::calc_A_c(double ci_) {
  return (vcmax * (ci_ - gamma_25 * umol_per_mol_to_Pa)) / (ci_ + km_25);
}

double Leaf::calc_A_j(double PPFD, double ci_) {

  double jmax = vcmax * vcmax_25_to_jmax_25;
  double j = (a * PPFD + jmax -
              sqrt(pow(a * PPFD + jmax, 2) - 4 * curv_fact * a * PPFD * jmax)) /
             (2 * curv_fact); // check brackets are correct

  return j / 4 *
         ((ci_ - gamma_25 * umol_per_mol_to_Pa) /
          (ci_ + 2 * gamma_25 * umol_per_mol_to_Pa));
}

double Leaf::calc_A_lim(double PPFD, double ci_) {

  double A_c = calc_A_c(ci_);
  double A_j = calc_A_j(PPFD, ci_);
  double R_d = vcmax * 0.015;

  return (A_c + A_j - sqrt(pow(A_c + A_j, 2) - 4 * 0.98 * A_c * A_j)) /
             (2 * 0.98) -
         R_d;
}

double Leaf::diff_ci(double PPFD, double x, double psi_soil, double psi_stem, double k_l_max) {

  double A_lim_ = calc_A_lim(PPFD, x);

  double g_c_ = calc_g_c(psi_soil, psi_stem, k_l_max);

  return A_lim_ * umol_per_mol_to_mol_per_mol -
         (g_c_ * (ca - x) / (atm_kpa * kPa_to_Pa));
}

// need to fill in tol and max_iteratiosn
double Leaf::calc_assim_gross(
    double PPFD, double psi_soil, double psi_stem, double k_l_max) {


  // not clear what x is here
  auto target = [&](double x) mutable -> double {
    return diff_ci(PPFD, x, psi_soil, psi_stem, k_l_max);
  };

  // tol and iterations copied from control defaults (for now)
  ci = util::uniroot(target, 0, ca, 1e-8, 1000);

  A_lim = calc_A_lim(PPFD, ci);

  g_c = calc_g_c(psi_soil, psi_stem, k_l_max);

  E = g_c * 1.6 * atm_vpd / kg_to_mol_h2o / atm_kpa;

  psi = psi_stem;

  return A_lim;
  
}

// Sperry et al. 2017; Sabot et al. 2020 implementation

double Leaf::calc_hydraulic_cost_Sperry(double psi_soil, double psi_stem,
                                 double k_l_max) {
  double k_l_soil_ = calc_cond_vuln(psi_soil, k_l_max);
  double k_l_stem_ = calc_cond_vuln(psi_stem, k_l_max);

  return k_l_soil_ - k_l_stem_;
}

double Leaf::calc_profit_Sperry(double PPFD, double psi_soil, double psi_stem, double k_l_max) {

  double lambda_ =
      calc_assim_gross(PPFD, psi_soil, psi_crit, k_l_max) /
      (calc_cond_vuln(psi_soil, k_l_max) -
       calc_cond_vuln(psi_crit, k_l_max));
  double benefit_ =
      calc_assim_gross(PPFD, psi_soil, psi_stem, k_l_max);
  double cost_ = calc_hydraulic_cost_Sperry(psi_soil, psi_stem, k_l_max);

  return benefit_ - lambda_ * cost_;
}

// Bartlett et al. implementation

double Leaf::calc_hydraulic_cost_Bartlett(double psi_soil, double psi_stem,
                                          double k_l_max) {

  double k_l_soil_ = calc_cond_vuln(psi_soil, k_l_max);
  double k_l_stem_ = calc_cond_vuln(psi_stem, k_l_max);
  double height_ = (K_s * huber_value)/k_l_max;

  return beta * huber_value * height_ * pow((1 - k_l_stem_ / k_l_soil_), beta_2);
}

double Leaf::calc_profit_Bartlett(double PPFD, double psi_soil, double psi_stem,
                                  double k_l_max) {                                  

  double benefit_ =
      calc_assim_gross(PPFD, psi_soil, psi_stem, k_l_max);

  double cost_ =
      calc_hydraulic_cost_Bartlett(psi_soil, psi_stem, k_l_max);

  return benefit_ - cost_;
}

// need docs on Golden Section Search and reference to Bartlett.
double Leaf::optimise_psi_stem_Bartlett(double PPFD, double psi_soil, double k_l_max) {

  double gr = (sqrt(5) + 1) / 2;
  double opt_psi_stem = psi_soil;

  // optimise for stem water potential
  if (psi_soil < psi_crit) {
    double bound_a = psi_soil;
    double bound_b = psi_crit;

    double delta_crit = 1e-3;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;

    while (abs(bound_b - bound_a) > delta_crit) {

      double profit_at_c =
          calc_profit_Bartlett(PPFD, psi_soil, bound_c, k_l_max);

      double profit_at_d =
          calc_profit_Bartlett(PPFD, psi_soil, bound_d, k_l_max);

      if (profit_at_c > profit_at_d) {
        bound_b = bound_d;
      } else {
        bound_a = bound_c;
      }

      bound_c = bound_b - (bound_b - bound_a) / gr;
      bound_d = bound_a + (bound_b - bound_a) / gr;
    }

    opt_psi_stem = ((bound_b + bound_a) / 2);
  }

  return opt_psi_stem;
}


// need docs on Golden Section Search and reference to Bartlett.
double Leaf::optimise_psi_stem_Sperry(double PPFD, double psi_soil, double k_l_max) {

  double gr = (sqrt(5) + 1) / 2;
  double opt_psi_stem = psi_soil;

  // optimise for stem water potential
  if (psi_soil < psi_crit) {
    double bound_a = psi_soil;
    double bound_b = psi_crit;

    double delta_crit = 1e-5;

    double bound_c = bound_b - (bound_b - bound_a) / gr;
    double bound_d = bound_a + (bound_b - bound_a) / gr;

    while (abs(bound_b - bound_a) > delta_crit) {
        std::cout << "bound_b = " << bound_b << ", bound_a = " << bound_a << "\n";
      double profit_at_c =
          calc_profit_Sperry(PPFD, psi_soil, bound_c, k_l_max);

      double profit_at_d =
          calc_profit_Sperry(PPFD, psi_soil, bound_d, k_l_max);

      if (profit_at_c > profit_at_d) {
        bound_b = bound_d;
      } else {
        bound_a = bound_c;
      }

      bound_c = bound_b - (bound_b - bound_a) / gr;
      bound_d = bound_a + (bound_b - bound_a) / gr;
    }

    opt_psi_stem = ((bound_b + bound_a) / 2);
  }

  return opt_psi_stem;
}



} // namespace plant
