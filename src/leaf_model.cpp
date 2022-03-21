#include <plant/leaf_model.h>
#include <Rcpp.h>
#include <plant/util_post_rcpp.h>

namespace plant {
namespace leaf {

void leaf_model::calc_soil_water_pot(const FF16_Environment &environment, double psi_aep, double b_CH) {
    std::vector<double> rel_theta = environment.get_soil_water_state();
    return psi_aep*pow(rel_theta, b_CH)

void leaf_model::calc_k_l_max(double K_s, double h_v, double h) const {
  return K_s*h_v/h;
}

void leaf_model::calc_vul_b(double p_50, double c) const {
    return p_50/
    pow(-log(1-50/100)), 1/c);
}

void leaf_model::calc_cond_vuln(double K_s, double h_v, double h, double p_50, double c) const {
    return calc_k_l_max(K_s, h_v, h)*exp(-(x/calc_vul_b(p_50, c)));
}

void leaf_model::calc_g_c(double psi_stem, double psi_soil, double atm_vpd, double psi_aep, double b_CH, const kg_2_mol_h20){
    return atm_kpa*calc_E_supply(psi_stem, psi_soil = calc_soil_water_pot(FF16_Environment &environment, psi_aep, b_CH), K_s, h_v, h, p_50, c)*kg_2_mol_h20/atm_vpd
    /1.6;
}

void leaf_model::calc_A_c(double c_i, double vcmax, double gamma_25, double umol_per_mol_2_Pa, double km_25){
    return (vcmax * (c_i - gamma_25*umol_per_mol_2_Pa))/(c_i + km_25);
}

void leaf_model::calc_A_j(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double c_i){
    jmax = vcmax*vcmax_25_2_jmax_25
    j = (a*PPFD + jmax - sqrt((a*PPFD+jmax)^2-4*curv_fact*PPFD*jmax))/(2*curv_fact)
    return j/4 * ((c_i - gamma_25*umol_per_mol_2_Pa)/(c_i + 2*gamma_25*umol_per_mol_2_Pa));
}

void leaf_model::calc_A_lim(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double c_i){
     double A_j = calc_A_j(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, c_i)
     double A_c = calc_A_c(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, c_i)

     return (A_c + A_j  - sqrt((A_c + A_j)^2-4*0.98*A_c*A_j))/(2*0.98);
}

void leaf_model::diff_ci(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double c_i){
     return calc_A_lim(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, c_i)*umol_per_mol_2_mol_per_mol - (calc_g_c(psi_stem, psi_soil, atm_vpd, psi_aep, b_CH, kg_2_mol_h20)*(ca - x)/(atm_kpa*kPa_2_Pa));
}

// need to fill in tol and max_iteratiosn
void leaf_model::solve_A_lim(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double c_i){
     double ci_root_ = util::uniroot(diff_ci, 0, 40, tol, max_iterations);
     return calc_A_lim(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, double c_i = ci_root_)
}

void leaf_model::solve_A_lim(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double c_i){
     double ci_root_ = util::uniroot(diff_ci, 0, 40, tol, max_iterations);
     return calc_A_lim(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, double c_i = ci_root_)
}