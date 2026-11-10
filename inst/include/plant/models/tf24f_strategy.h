// -*-c++-*-
#ifndef PLANT_PLANT_TF24F_STRATEGY_H_
#define PLANT_PLANT_TF24F_STRATEGY_H_

#include <plant/models/tf24_strategy.h>

namespace plant {

// TF24f ("f" for fast / forecasting): a variant of TF24 that will let the
// optimal leaf hydraulic state chase its optimum via an extra ODE state
// (gradient-ascent), instead of re-solving the nested leaf optimisation /
// root-collar root-find from scratch at every step (issue #525). It inherits
// TF24_Strategy and reuses TF24_Environment + TF24_Pars; only the leaf-solve
// hook and the extra state are overridden.
//
// Adds one extra ODE state, opt_root_psi_state, holding the tracked optimal
// root-collar water potential. Instead of TF24's per-step golden-section
// optimisation, solve_leaf() evaluates the leaf at the tracked state and the
// state relaxes toward its optimum by gradient ascent on carbon profit
// (dpsi/dt = k_acclim * d(profit)/d(psi)); the state is initialised at its optimum at
// birth. The states shared with TF24 are computed by the inherited
// compute_rates, so TF24f does NOT reproduce TF24 bit-for-bit -- it tracks the
// optimum with a (gain-dependent) lag and is faster.
template <typename S = double>
class TF24f_Strategy : public TF24_Strategy<S> {
public:
  using value_type = S;
  using base_type = TF24_Strategy<S>;

  typedef std::shared_ptr<TF24f_Strategy> ptr;
  TF24f_Strategy();

  // Scientific version — compound, because TF24f is a fast *approximation* of
  // TF24 and inherits its equations/parameters. The version is reported as
  // "<TF24 version>.<approximation_revision>" (e.g. "2.1"), so:
  //   * the major component auto-tracks TF24_Strategy::scientific_version, so a
  //     TF24 scientific change also invalidates TF24f (the safe direction);
  //   * bump `approximation_revision` for changes specific to the fast
  //     approximation itself.
  // See plant::model_version() / model_id() and src/strategy_version.cpp. Bump
  // rules match the other models (output-changing science only, not refactors).
  static constexpr int approximation_revision = 1;

  // TF24's states + the tracked root-collar psi (appended last so the
  // inherited state indices are unchanged). state_size()/state_names() are
  // static and resolved on the concrete type by Individual<TF24f, ...>.
  static size_t state_size() { return base_type::state_size() + 1; }
  static std::vector<std::string> state_names() {
    std::vector<std::string> ret = base_type::state_names();
    ret.push_back("opt_root_psi_state");
    return ret;
  }

  // Re-register indices including the new state slot. Base refresh_indices()
  // calls the *static* (base) state_names(), so it would otherwise miss the
  // appended state; we call it then add the extra slot.
  void refresh_indices();

  void compute_rates(const TF24_Environment<S>& environment, Internals<S>& vars);

  // Override the leaf solve: instead of optimising the root-collar psi, evaluate
  // the leaf at the tracked state and finite-difference the profit gradient
  // (left in dprofit_dpsi_ for compute_rates to turn into the state's rate).
  void solve_leaf();

  // Initialise the tracked state at its optimum for a newly introduced individual, so
  // gradient ascent starts at the optimum (no birth transient / no climb from 0,
  // which otherwise lets shaded recruits escape suppression). Runs the base
  // optimiser once via the initializing_ flag below.
  void set_initial_states(const TF24_Environment<S>& environment, Internals<S>& vars);

  // Acclimation gain k in  dpsi/dt = k * d(profit)/d(psi). Exposed to R so the
  // stiffness / accuracy-vs-speed k-sweep (#525) can be driven without a rebuild;
  // large k recovers the quasi-steady-state (TF24) optimum.
  S k_acclim = 1.0;

  // Cached slot for the tracked state, resolved in refresh_indices().
  int state_idx_opt_root_psi_state = -1;

  // Channel between solve_leaf() (writes) and compute_rates() (reads). Default is
  // finite so solve_leaf is safe before the first compute_rates (e.g. during
  // establishment_probability), where the tracked state has not been read yet.
  S tracked_root_psi_ = 0.0;
  S dprofit_dpsi_ = 0.0;
  // When true, solve_leaf() runs the base optimiser (Leaf::optimise on the seated
  // curve) rather than the tracked evaluation; used by set_initial_states to read
  // the optimum.
  bool initializing_ = false;
};

template <typename S>
typename TF24f_Strategy<S>::ptr make_strategy_ptr(TF24f_Strategy<S> s);

// The base TF24_Strategy constructor sets collect_all_auxiliary and the name
// and runs the *base* refresh_indices(); we re-run our own refresh_indices()
// here so the appended state slot is registered, and set the reported name.
template <typename S>
TF24f_Strategy<S>::TF24f_Strategy() {
  this->name = "TF24f";
  this->refresh_indices();
}

// Build on the base index maps, then register the extra tracked state. The new
// slot is appended after TF24's states, so its index is TF24's state_size().
template <typename S>
void TF24f_Strategy<S>::refresh_indices() {
  base_type::refresh_indices();
  const int idx = static_cast<int>(base_type::state_size());
  this->state_index["opt_root_psi_state"] = idx;
  state_idx_opt_root_psi_state = idx;
  // Again on the concrete model: the base checked its own names against its
  // own count, and the slot appended above is outside both.
  check_state_layout(this->state_index, state_size(), "TF24f");
}

// Reuse TF24's rates for the shared states; the tracked-state value is fed
// to solve_leaf() (called from the reused net_mass_production_dt) via
// tracked_root_psi_, and the resulting profit gradient comes back in
// dprofit_dpsi_, which becomes the tracked state's rate (gradient ascent).
template <typename S>
void TF24f_Strategy<S>::compute_rates(const TF24_Environment<S>& environment,
                                   Internals<S>& vars) {
  // k_acclim is a user-settable gain; a negative value would silently turn the
  // gradient ascent into descent (away from the optimum) and a non-finite value
  // would poison the state rate, so fail fast on misconfiguration.
  if (!util::is_finite(k_acclim) || k_acclim < 0.0) {
    util::stop("TF24f: k_acclim must be finite and >= 0 (got " +
               util::to_string(k_acclim) + ")");
  }
  tracked_root_psi_ = vars.state(state_idx_opt_root_psi_state);
  base_type::compute_rates(environment, vars);
  vars.set_rate(state_idx_opt_root_psi_state, k_acclim * dprofit_dpsi_);
}

// Track instead of optimise: evaluate the leaf at the tracked collar psi and
// supply d(profit)/d(psi) for the gradient-ascent rate. The tracked state is
// clamped into the feasible interval here (phylloptim refuses a collar outside
// it), and the gradient is formed about the clamped value `used`. This keeps the
// gradient meaningful even when the tracked state has drifted outside the
// interval (e.g. an uninitialised state at 0, or a wet bound that a drying soil
// moves past the state during shutdown), pulling it back inside, and the final
// evaluate leaves the leaf outputs at the operating point that compute_rates'
// aux reads expect. The gradient is exact: forward-mode AD over the leaf's
// algebra, with the implicit-function theorem at the ci root-find.
// A function-try-block, so every exit of the body is covered: phylloptim's
// infeasible_error becomes odelia's DomainError, a rejected step rather than a
// dead run, for the reason TF24_Strategy::solve_leaf gives. Only that error is
// translated; anything else still propagates.
template <typename S>
void TF24f_Strategy<S>::solve_leaf() try {
  if (initializing_) {
    // Birth initialisation: run the full optimiser so set_initial_states can
    // read the optimum collar psi.
    // The curve-typed form, for the reason TF24_Strategy::solve_leaf gives.
    this->leaf.template find_root_collar_psi_for<
        phylloptim::Leaf::CostCurve::TF24_floor>();
    return;
  }
  // Forward-mode AD over the analytic algebra, the implicit-function theorem at
  // the ci root-find, and analytic spline derivatives for the transport.
  // Establish the operating point (and the clamped collar psi `used`) and leave
  // the leaf outputs there for compute_rates' aux reads. Where feasibility
  // handling fixes the operating point (shutdown / closed interval) `used` is
  // the collar it chose.
  double bound_a, bound_b;
  const bool open = this->leaf.template prepare_collar_solve<
      phylloptim::Leaf::CostCurve::TF24_floor>(bound_a, bound_b);
  const double used =
      open ? std::min(std::max(tracked_root_psi_, bound_a), bound_b)
           : this->leaf.opt_root_psi_;
  dprofit_dpsi_ = this->leaf.dprofit_droot_collar_psi(used);
  // restore operating-point outputs
  if (open) {
    this->leaf.template profit_at_collar_psi<
        phylloptim::Leaf::CostCurve::TF24_floor>(used, bound_a, bound_b);
  } else {
    this->leaf.template prepare_collar_solve<
        phylloptim::Leaf::CostCurve::TF24_floor>(bound_a, bound_b);
  }
} catch (const phylloptim::util::infeasible_error& e) {
  odelia::util::stop_domain(std::string("leaf solve infeasible: ") + e.what());
}

// Initialise the tracked state at its optimum: run the base optimiser once (via the
// initializing_ flag, which makes solve_leaf optimise rather than track) and
// store the resulting collar psi as the initial state. Both are the positive
// magnitude of the potential, in MPa.
template <typename S>
void TF24f_Strategy<S>::set_initial_states(const TF24_Environment<S>& environment,
                                        Internals<S>& vars) {
  // Initialise the shared TF24 states first (notably the NSC storage pool, #517) --
  // set_initial_states is non-virtual, so without this call TF24f seedlings
  // would be born with empty reserves and die immediately.
  base_type::set_initial_states(environment, vars);
  initializing_ = true;
  this->net_mass_production_dt(environment, vars);
  initializing_ = false;
  vars.set_state(state_idx_opt_root_psi_state, this->leaf.opt_root_psi_);
}

template <typename S>
typename TF24f_Strategy<S>::ptr make_strategy_ptr(TF24f_Strategy<S> s) {
  s.prepare_strategy();
  return std::make_shared<TF24f_Strategy<S> >(s);
}

}

#endif
