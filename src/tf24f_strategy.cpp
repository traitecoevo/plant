#include <plant/models/tf24f_strategy.h>

namespace plant {

// The base TF24_Strategy constructor sets collect_all_auxiliary and the name
// and runs the *base* refresh_indices(); we re-run our own refresh_indices()
// here so the appended state slot is registered, and set the reported name.
TF24f_Strategy::TF24f_Strategy() {
  name = "TF24f";
  refresh_indices();
}

// Build on the base index maps, then register the extra tracked state. The new
// slot is appended after TF24's states, so its index is TF24's state_size().
void TF24f_Strategy::refresh_indices() {
  TF24_Strategy::refresh_indices();
  const int idx = static_cast<int>(TF24_Strategy::state_size());
  state_index["opt_root_psi_state"] = idx;
  state_idx_opt_root_psi_state = idx;
}

// Reuse TF24's rates for the five shared states; the tracked-state value is fed
// to solve_leaf() (called from the reused net_mass_production_dt) via
// tracked_root_psi_, and the resulting profit gradient comes back in
// dprofit_dpsi_, which becomes the tracked state's rate (gradient ascent).
void TF24f_Strategy::compute_rates(const TF24_Environment& environment,
                                   Internals& vars) {
  tracked_root_psi_ = vars.state(state_idx_opt_root_psi_state);
  TF24_Strategy::compute_rates(environment, vars);
  vars.set_rate(state_idx_opt_root_psi_state, k_acclim * dprofit_dpsi_);
}

// Track instead of optimise: evaluate the leaf at the tracked collar psi and
// finite-difference the profit gradient. evaluate_root_collar_psi clamps to the
// feasible interval, so we perturb from the *clamped* operating value (read back
// as -root_collar_psi_); this keeps the gradient meaningful even when the
// tracked state sits outside the feasible interval (e.g. an uninitialised state
// at 0), pulling it back inside. The final evaluate leaves the leaf outputs at
// the operating point that compute_rates' aux reads expect.
void TF24f_Strategy::solve_leaf() {
  const double h = psi_fd_step;
  const double p0 = leaf.evaluate_root_collar_psi(tracked_root_psi_);
  const double used = -leaf.root_collar_psi_;
  const double p1 = leaf.evaluate_root_collar_psi(used + h);
  dprofit_dpsi_ = (p1 - p0) / h;
  leaf.evaluate_root_collar_psi(used);
}

TF24f_Strategy::ptr make_strategy_ptr(TF24f_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<TF24f_Strategy>(s);
}

}
