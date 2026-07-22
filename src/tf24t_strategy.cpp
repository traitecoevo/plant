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
