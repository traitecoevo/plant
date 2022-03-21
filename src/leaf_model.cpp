#include <plant/leaf_model.h>

namespace plant {
Leaf::Leaf(){}


// TODO: vectorise for soil layers
double Leaf::calc_soil_water_pot(const FF16_Environment &environment, double psi_aep, double b_CH) {
    double rel_theta = environment.get_soil_water_state()[0]; // this hack assumes one soil layer
    return psi_aep*pow(rel_theta, b_CH);
}

double Leaf::calc_k_l_max(double K_s, double h_v, double h) const {
  return K_s*h_v/h;
}

double Leaf::calc_vul_b(double p_50, double c) const {
    return p_50 / pow(-log(1-50/100), 1/c);
}

// integrates
double Leaf::calc_cond_vuln(double psi, double k_l_max, double b, double c) const {
    return k_l_max * exp(-pow((psi / b), c));
}

// replace f with some other function
double calc_E_supply(double k_l_max, double b, double c, double psi_soil, double psi_stem){
    std::function<double(double)> f;
    f = [&] (double psi) -> double {
         return calc_cond_vuln(psi, k_l_max, b, c);
       };

  integrator.integrate(f, psi_soil, psi_stem)
}

//,
double Leaf::calc_g_c(double psi_stem, double psi_soil, double atm_kpa, double atm_vpd, double psi_aep, double b_CH, const double kg_2_mol_h2o){
    double = psi_soil = calc_soil_water_pot(FF16_Environment &environment, psi_aep, b_CH);
    double E_supply = calc_E_supply(psi_stem, psi_soil, k_l_max, p_50, c, b);

    return atm_kpa * E_supply * kg_2_mol_h2o / atm_vpd / 1.6;
}

double Leaf::calc_A_c(double c_i, double vcmax, double gamma_25, double umol_per_mol_2_Pa, double km_25){
    return (vcmax * (c_i - gamma_25*umol_per_mol_2_Pa))/(c_i + km_25);
}

double Leaf::calc_A_j(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double c_i){
    double jmax = vcmax*vcmax_25_2_jmax_25;
    double j = (a*PPFD + jmax - sqrt((a*PPFD+jmax)^2-4*curv_fact*PPFD*jmax))/(2*curv_fact);
    return j / 4 * ((c_i - gamma_25*umol_per_mol_2_Pa)/(c_i + 2*gamma_25*umol_per_mol_2_Pa));
}

double Leaf::calc_A_lim(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double c_i){
     double A_j = calc_A_j(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, c_i);
     double A_c = calc_A_c(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, c_i);

     return (A_c + A_j  - sqrt(pow(A_c + A_j, 2) - 4 * 0.98 * A_c * A_j)/(2 * 0.98));
}

double Leaf::diff_ci(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double c_i,
                     double psi_stem, double psi_soil, double atm_vpd, double psi_aep, double b_CH, double kg_2_mol_h2o, double ca, double x, double atm_kpa, double kPa_2_Pa ){
     double A_lim_ = calc_A_lim(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, c_i);
     double g_c_ = calc_g_c(psi_stem, psi_soil, atm_vpd, psi_aep, b_CH, kg_2_mol_h2o)*(ca - x)/(atm_kpa*kPa_2_Pa)

     return A_lim_ * umol_per_mol_2_mol_per_mol - (g_c_ * (ca - x) / (atm_kpa * kPa_2_Pa));

      // calc_A_lim(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, c_i) *umol_per_mol_2_mol_per_mol - (calc_g_c(psi_stem, psi_soil, atm_vpd, psi_aep, b_CH, kg_2_mol_h2o)*(ca - x)/(atm_kpa*kPa_2_Pa));
}

// need to fill in tol and max_iteratiosn
double Leaf::solve_A_lim(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double c_i){
     double diff_ci_ = diff_ci(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, c_i);
     double ci_root_ = util::uniroot(diff_ci, 0, 40, 1e-8, 1000); // tol and iterations copied from control defaults (for now)
     return calc_A_lim(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, double c_i = ci_root_);
}

// void Leaf::solve_A_lim(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double c_i){
//      double ci_root_ = util::uniroot(diff_ci, 0, 40, tol, max_iterations);
//      return calc_A_lim(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, double c_i = ci_root_);
// }
}
