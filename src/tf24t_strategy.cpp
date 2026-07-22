#include <plant/models/tf24t_strategy.h>

namespace plant {

// The base TF24_Strategy constructor sets collect_all_auxiliary and runs the
// base refresh_indices(); we re-run our own so the appended state slots are
// registered, and set the reported name.
TF24t_Strategy::TF24t_Strategy() {
  name = "TF24t";
  refresh_indices();
}

// Build on the base index maps, then register the two appended acclimation
// states after TF24's states (indices TF24 state_size() and +1).
void TF24t_Strategy::refresh_indices() {
  TF24_Strategy::refresh_indices();
  const int idx = static_cast<int>(TF24_Strategy::state_size());
  state_index["acclim_topt"] = idx;
  state_index["acclim_tcrit"] = idx + 1;
  state_idx_acclim_topt = idx;
  state_idx_acclim_tcrit = idx + 1;
}

// Fail fast on a misconfigured acclimation kinetic (a negative/non-finite gain
// or decay would poison the state rate); shared guard for compute_rates and
// set_initial_states.
namespace {
void check_kinetics(double alpha, double beta, double s) {
  if (!util::is_finite(alpha) || alpha < 0.0 ||
      !util::is_finite(beta) || beta < 0.0 ||
      !util::is_finite(s) || s <= 0.0) {
    util::stop("TF24t: acclimation kinetics must be finite with alpha,beta>=0 and softplus_s>0");
  }
}
}

// Enable the leaf thermal-damage layer and copy the thermal traits into the
// leaf, after the base sets up the leaf (and the PM energy balance). Mirrors the
// base's use_energy_balance_ / d_ handling.
void TF24t_Strategy::prepare_strategy() {
  TF24_Strategy::prepare_strategy();
  // Phase 4 (midday evaluation wiring, #566): the Lumry-Eyring damage factor N is
  // a single quasi-steady evaluation at the *midday operating-point* leaf
  // temperature, and that operating point only exists on the Penman-Monteith
  // energy-balance path (set_leaf_states_rates_from_psi_stem recomputes the temp
  // params -- and N -- at Tleaf = f(E) per candidate psi; base off-PM path leaves
  // N at the prescribed-Tair baseline with no cooling/avoidance response). The
  // model hierarchy is TF24 subset TF24+PM subset TF24t (+PM+ATLS), so force PM on
  // here regardless of pars.use_energy_balance -- thermal damage is not physically
  // meaningful without a real midday Tleaf. This mirrors the unconditional
  // use_thermal_damage_ = true below.
  leaf.use_energy_balance_ = true;
  leaf.use_thermal_damage_ = true;
  leaf.topt_offset_ = topt_offset;
  leaf.tcrit_0_ = tcrit_0;
  leaf.k_d1_0_ = k_d1_0;
  leaf.k_r1_0_ = k_r1_0;
  leaf.m_switch_ = m_switch;
  leaf.m_rep_ = m_rep;
  leaf.t_rep_cut_ = t_rep_cut;
  leaf.dTcrit_max_ = dTcrit_max;
  leaf.dTopt_max_ = dTopt_max;
  leaf.K_A_ = K_A;
  // Leaf-side thermal-cost coefficients (Phase 3): acclimation maintenance,
  // repair standing maintenance, repair activity. Added to leaf.R_d_.
  leaf.c_acclim_maint_ = c_acclim_maint;
  leaf.c_repair_maint_ = c_repair_maint;
  leaf.c_repair_flux_  = c_repair_flux;
  // Start the leaf's acclimation inputs at 0 (compute_rates overwrites them with
  // the live states each step; this only defines what establishment_probability
  // reads before any rates run).
  leaf.A_opt_ = 0.0;
  leaf.A_crit_ = 0.0;
}

// Smooth (C-infinity) positive part 0.5*(x + sqrt(x^2 + eps^2)) -> max(0,x),
// matching the form TF24 uses for the positive part of net production. Used for
// the induction cost so the ODE-integrated carbon balance is smooth in the
// acclimation build rate.
namespace {
double positive_part(double x, double eps) {
  return 0.5 * (x + std::sqrt(x * x + eps * eps));
}
}

// Whole-plant carbon balance = base TF24 net production, minus the thermal costs
// that are not leaf-local respiration (those already sit in leaf.R_d_):
//  - tolerance construction cost: a one-off build cost per leaf, amortised
//    through leaf turnover, charged for constitutive tolerance bought above the
//    intrinsic optimum (topt_offset > 0, tcrit_0 > tcrit_ref).
//  - acclimation induction cost: proportional to the smoothed positive build
//    rate of each acclimation state (paying to synthesise the proteins).
// leaf.A_opt_/A_crit_ hold the current acclimation states (compute_rates set
// them before the inherited compute_rates called this); at establishment they
// read 0, a negligible build-from-zero snapshot at default temperatures.
double TF24t_Strategy::net_mass_production_dt(const TF24_Environment& environment,
                                             double height, double area_leaf_,
                                             double height_inverse) {
  const double base =
    TF24_Strategy::net_mass_production_dt(environment, height, area_leaf_,
                                          height_inverse);

  const double turnover_leaf_ = turnover_leaf(mass_leaf(area_leaf_));
  const double construction =
    (c_build_topt * std::max(0.0, topt_offset) +
     c_build_tcrit * std::max(0.0, tcrit_0 - tcrit_ref)) * turnover_leaf_;

  const double T = environment.get_leaf_temp();
  const double dA_opt =
    alpha_opt * Leaf::softplus_(T - t_accl_opt, softplus_s) - beta_opt * leaf.A_opt_;
  const double dA_crit =
    alpha_crit * Leaf::softplus_(T - t_accl_crit, softplus_s) - beta_crit * leaf.A_crit_;
  const double induction =
    c_accl_induct * (positive_part(dA_opt, induct_eps) +
                     positive_part(dA_crit, induct_eps));

  return base - construction - induction;
}

// Feed the current acclimation states into the leaf (so this step's damage
// factor and T_opt shift see them), reuse TF24's rates for the shared states,
// then set the two acclimation-state rates:
//   dA/dt = alpha * softplus(T_regime - T_accl) - beta * A.
// The forcing temperature is the (midday) environmental leaf-temp driver.
void TF24t_Strategy::compute_rates(const TF24_Environment& environment,
                                   Internals& vars) {
  check_kinetics(alpha_opt, beta_opt, softplus_s);
  check_kinetics(alpha_crit, beta_crit, softplus_s);

  const double A_opt = vars.state(state_idx_acclim_topt);
  const double A_crit = vars.state(state_idx_acclim_tcrit);
  leaf.A_opt_ = std::max(A_opt, 0.0);
  leaf.A_crit_ = std::max(A_crit, 0.0);

  TF24_Strategy::compute_rates(environment, vars);

  const double T = environment.get_leaf_temp();
  vars.set_rate(state_idx_acclim_topt,
                alpha_opt * Leaf::softplus_(T - t_accl_opt, softplus_s) - beta_opt * A_opt);
  vars.set_rate(state_idx_acclim_tcrit,
                alpha_crit * Leaf::softplus_(T - t_accl_crit, softplus_s) - beta_crit * A_crit);
}

// Seed the acclimation states at their environmental equilibrium
// (A_eq = (alpha/beta) * softplus(T - T_accl)) so a recruit is born already
// acclimated to its birth environment -- no cold-start transient. Seed the
// shared TF24 states (notably the NSC storage pool, #517) first.
void TF24t_Strategy::set_initial_states(const TF24_Environment& environment,
                                        Internals& vars) {
  TF24_Strategy::set_initial_states(environment, vars);
  check_kinetics(alpha_opt, beta_opt, softplus_s);
  check_kinetics(alpha_crit, beta_crit, softplus_s);

  const double T = environment.get_leaf_temp();
  const double A_opt_eq =
    beta_opt > 0.0 ? alpha_opt / beta_opt * Leaf::softplus_(T - t_accl_opt, softplus_s) : 0.0;
  const double A_crit_eq =
    beta_crit > 0.0 ? alpha_crit / beta_crit * Leaf::softplus_(T - t_accl_crit, softplus_s) : 0.0;
  vars.set_state(state_idx_acclim_topt, A_opt_eq);
  vars.set_state(state_idx_acclim_tcrit, A_crit_eq);
}

TF24t_Strategy::ptr make_strategy_ptr(TF24t_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<TF24t_Strategy>(s);
}

}
