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
// (dpsi/dt = k_acclim * d(profit)/d(psi)); the state is seeded at its optimum at
// birth. The five states shared with TF24 are computed by the inherited
// compute_rates, so TF24f does NOT reproduce TF24 bit-for-bit -- it tracks the
// optimum with a (gain-dependent) lag and is faster.
class TF24f_Strategy : public TF24_Strategy {
public:
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

  // TF24's five states + the tracked root-collar psi (appended last so the
  // inherited state indices 0..4 are unchanged). state_size()/state_names() are
  // static and resolved on the concrete type by Individual<TF24f, ...>.
  static size_t state_size() { return TF24_Strategy::state_size() + 1; }
  static std::vector<std::string> state_names() {
    std::vector<std::string> ret = TF24_Strategy::state_names();
    ret.push_back("opt_root_psi_state");
    return ret;
  }

  // Re-register indices including the new state slot. Base refresh_indices()
  // calls the *static* (base) state_names(), so it would otherwise miss the
  // appended state; we call it then add the extra slot.
  void refresh_indices();

  void compute_rates(const TF24_Environment& environment, Internals& vars);

  // Override the leaf solve: instead of optimising the root-collar psi, evaluate
  // the leaf at the tracked state and finite-difference the profit gradient
  // (left in dprofit_dpsi_ for compute_rates to turn into the state's rate).
  //
  // Resolves the seated cost curve once (see the definition) and hands the whole
  // body to solve_leaf_for<K>, so every phylloptim call inside it differentiates
  // and evaluates the SAME model.
  void solve_leaf();
  // The tracked solve for one cost curve. Defined in the .cpp and instantiated
  // only through `solve_leaf`'s `with_curve` dispatch, which is why it needs no
  // explicit instantiation list: the lambda there names every arm.
  template <Leaf::CostCurve K> void solve_leaf_for();

  // Seed the tracked state at its optimum for a newly introduced individual, so
  // gradient ascent starts at the optimum (no birth transient / no climb from 0,
  // which otherwise lets shaded recruits escape suppression). Runs the base
  // optimiser once via the initializing_ flag below.
  void set_initial_states(const TF24_Environment& environment, Internals& vars);

  // Acclimation gain k in  dpsi/dt = k * d(profit)/d(psi). Exposed to R so the
  // stiffness / accuracy-vs-speed k-sweep (#525) can be driven without a rebuild;
  // large k recovers the quasi-steady-state (TF24) optimum.
  double k_acclim = 1.0;
  // Finite-difference step (positive magnitude, MPa) for d(profit)/d(psi); used
  // only when use_ad_gradient is false.
  double psi_fd_step = 1e-3;
  // Use the exact AD/IFT gradient (Leaf::dprofit_droot_collar_psi) instead of the
  // finite difference. Exposed to R for A/B comparison (#527).
  bool use_ad_gradient = true;

  // Cached slot for the tracked state, resolved in refresh_indices().
  int state_idx_opt_root_psi_state = -1;

  // Channel between solve_leaf() (writes) and compute_rates() (reads). Default is
  // finite so solve_leaf is safe before the first compute_rates (e.g. during
  // establishment_probability), where the tracked state has not been read yet.
  double tracked_root_psi_ = 0.0;
  double dprofit_dpsi_ = 0.0;
  // When true, solve_leaf() runs the base optimiser (Leaf::optimise on the seated
  // curve) rather than the tracked evaluation; used by set_initial_states to read
  // the optimum.
  bool initializing_ = false;
};

TF24f_Strategy::ptr make_strategy_ptr(TF24f_Strategy s);

}

#endif
