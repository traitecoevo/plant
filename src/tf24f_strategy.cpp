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
  // k_acclim is a user-settable gain; a negative value would silently turn the
  // gradient ascent into descent (away from the optimum) and a non-finite value
  // would poison the state rate, so fail fast on misconfiguration.
  if (!util::is_finite(k_acclim) || k_acclim < 0.0) {
    util::stop("TF24f: k_acclim must be finite and >= 0 (got " +
               util::to_string(k_acclim) + ")");
  }
  tracked_root_psi_ = vars.state(state_idx_opt_root_psi_state);
  TF24_Strategy::compute_rates(environment, vars);
  vars.set_rate(state_idx_opt_root_psi_state, k_acclim * dprofit_dpsi_);
}

// Track instead of optimise: evaluate the leaf at the tracked collar psi and
// supply d(profit)/d(psi) for the gradient-ascent rate. evaluate_root_collar_psi
// clamps to the feasible interval, so we read back the *clamped* operating value
// (`used` = -root_collar_psi_) and form the gradient about it; this keeps the
// gradient meaningful even when the tracked state has drifted outside the
// interval (e.g. an uninitialised state at 0), pulling it back inside, and the
// final evaluate leaves the leaf outputs at the operating point that
// compute_rates' aux reads expect. Two gradient methods are available
// (use_ad_gradient): the exact AD/IFT gradient (default, #527) or a centred
// finite difference (#526); see the branches below.
void TF24f_Strategy::solve_leaf() {
  if (initializing_) {
    // Birth initialisation: run the full optimiser so set_initial_states can
    // read the optimum collar psi.
    leaf.find_root_collar_psi();
    return;
  }
  // Establish the operating point (and the clamped collar psi `used`) and leave
  // the leaf outputs there for compute_rates' aux reads.
  leaf.evaluate_root_collar_psi(tracked_root_psi_);
  const double used = -leaf.root_collar_psi_;

  if (use_ad_gradient) {
    // Exact gradient (default, #527): forward-mode AD over the analytic algebra
    // + IFT at the ci root-find + analytic spline derivatives for the transport.
    // No O(h) bias and no finite-difference step to tune.
    dprofit_dpsi_ = leaf.dprofit_droot_collar_psi(used);
    leaf.evaluate_root_collar_psi(used);  // restore operating-point outputs
  } else {
    // Centred finite-difference fallback (#526), perturbing about the clamped
    // operating value `used`. A one-sided difference biases the fixed point to
    // psi* - h/2 (O(h)); the centred form cancels that term (O(h^2)) for one
    // extra leaf evaluation. Near a boundary the clamp inside
    // evaluate_root_collar_psi collapses one arm, degrading it gracefully to a
    // one-sided difference that still points back into the interior.
    const double h = psi_fd_step;
    if (!util::is_finite(h) || h <= 0.0) {
      util::stop("TF24f: psi_fd_step must be finite and > 0 (got " +
                 util::to_string(h) + ")");
    }
    const double p_plus  = leaf.evaluate_root_collar_psi(used + h);
    const double p_minus = leaf.evaluate_root_collar_psi(used - h);
    dprofit_dpsi_ = (p_plus - p_minus) / (2.0 * h);
    leaf.evaluate_root_collar_psi(used);
  }
}

// Seed the tracked state at its optimum: run the base optimiser once (via the
// initializing_ flag, which makes solve_leaf optimise rather than track) and
// store the resulting collar psi as the initial state. leaf.root_collar_psi_ is
// the signed (negative) potential; the tracked state is the positive magnitude.
void TF24f_Strategy::set_initial_states(const TF24_Environment& environment,
                                        Internals& vars) {
  initializing_ = true;
  net_mass_production_dt(environment, vars);
  initializing_ = false;
  vars.set_state(state_idx_opt_root_psi_state, -leaf.root_collar_psi_);
}

TF24f_Strategy::ptr make_strategy_ptr(TF24f_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<TF24f_Strategy>(s);
}

}
