#include <plant/models/tf24f_strategy.h>
#include <plant/models/tf24_production_kernel.h>
#include <XAD/XAD.hpp>

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

// Exact d(dheight/dt)/d(height) at the tracked collar (#472 scope B / #537 A1).
// Two parts compose: (1) d(profit)/d(height) at the FIXED tracked collar, summed
// over the leaf input channels height moves -- kmax and the crown-centre light
// (both exact analytic via the leaf's forward-AD+IFT sensitivities), plus a weak
// E_up term (zero above the rooting-depth clamp, otherwise a stable central FD of
// the smooth closed-form soil uptake -- NOT the amplifying height_dt FD); and (2)
// forward-mode AD of the demographic growth kernel, seeding height and injecting
// profit with derivative dprofit_dh, so the explicit-height allometry/cascade
// terms compose with the profit coupling. Returns NA_REAL on an infeasible leaf or
// non-smooth light so Node::growth_rate_gradient falls back to its finite difference.
double TF24f_Strategy::growth_rate_gradient_height_ad(double height,
                                                      const TF24_Environment& environment) {
  using AD = xad::fwd<double>::active_type;
  // (1a) Operating point at (height, tracked collar): leaves leaf.* at the point
  // and supplies profit, the clamped collar `used`, kmax and E_up.
  net_mass_production_dt(environment, height, area_leaf(height), 1.0 / height);
  const double profit0 = leaf.profit_;
  if (!util::is_finite(profit0)) return NA_REAL;
  const double used  = -leaf.root_collar_psi_;
  const double kmax0 = leaf.leaf_specific_conductance_max_;
  const double dprofit_dkmax = leaf.dprofit_dkmax(used);
  const double dprofit_dPPFD = leaf.dprofit_dPPFD(used);
  const double dprofit_dEup  = leaf.dprofit_dEup(used);

  // (1b) Channel derivatives d(input)/d(height).
  const double dkmax_dh = -kmax0 / height;                     // kmax = K_s*theta/(h*eta_c)
  const double PPFD_top = environment.get_PPFD();
  const double z = height * eta_c;
  const double light = environment.get_environment_at_height(z);
  // d(light)/d(z): use a central FD of the light VALUES the crown-centre model
  // actually reads (get_environment_at_height), not the interpolator's analytic
  // get_environment_deriv_at_height -- the resident canopy spline's value and its
  // analytic derivative are not mutually consistent (knots at cohort heights), and
  // the model's light response is value-based, so a value-consistent slope is the
  // correct one. This FD is over the smooth light spline (a first derivative), NOT
  // the amplifying height_dt FD it replaces.
  const double dz = 1e-5 * z;
  const double dlight_dz =
      (environment.get_environment_at_height(z + dz) -
       environment.get_environment_at_height(z - dz)) / (2.0 * dz);
  const double dPPFD_dh = (light > 1e-4) ? pars.k_I * PPFD_top * dlight_dz * eta_c : 0.0;

  // E_up channel: zero above the 1.5 m rooting clamp; otherwise a benign central FD
  // of the smooth E_up (the leaf optimiser is never re-run -- evaluation at the
  // tracked collar only). Restore the operating point afterwards so the node's aux
  // reads (already taken) are not disturbed by a stale leaf state.
  double dEup_dh = 0.0;
  {
    const double d = 1e-5 * height;
    net_mass_production_dt(environment, height + d, area_leaf(height + d), 1.0 / (height + d));
    const double eup_p = leaf.E_up_;
    net_mass_production_dt(environment, height - d, area_leaf(height - d), 1.0 / (height - d));
    const double eup_m = leaf.E_up_;
    dEup_dh = (eup_p - eup_m) / (2.0 * d);
    net_mass_production_dt(environment, height, area_leaf(height), 1.0 / height);  // restore
  }

  const double dprofit_dh = dprofit_dkmax * dkmax_dh + dprofit_dPPFD * dPPFD_dh +
                            dprofit_dEup * dEup_dh;

  // (2) Forward-AD the growth kernel: lift the (double) prod-pars to AD constants,
  // seed height, inject profit carrying its total height derivative.
  const TF24ProdPars<double> p0 = prod_pars();
  TF24ProdPars<AD> p;
  p.lma=p0.lma; p.rho=p0.rho; p.theta=p0.theta; p.a_b1=p0.a_b1; p.a_r1=p0.a_r1;
  p.eta_c=p0.eta_c; p.r_l=p0.r_l; p.r_s=p0.r_s; p.r_b=p0.r_b; p.r_r=p0.r_r;
  p.k_l=p0.k_l; p.k_b=p0.k_b; p.k_s=p0.k_s; p.k_r=p0.k_r; p.a_bio=p0.a_bio; p.a_y=p0.a_y;
  p.a_l1=p0.a_l1; p.a_l2=p0.a_l2; p.a_f1=p0.a_f1; p.a_f2=p0.a_f2; p.hmat=p0.hmat;
  p.omega=p0.omega; p.a_f3=p0.a_f3; p.d_I=p0.d_I; p.a_dG1=p0.a_dG1; p.a_dG2=p0.a_dG2;

  AD h = height;       xad::derivative(h) = 1.0;
  AD profit = profit0; xad::derivative(profit) = dprofit_dh;
  AD area_leaf_ad = tf24_area_leaf<AD>(p.a_l1, p.a_l2, h);
  AD net = tf24_net_mass_production<AD>(p, h, area_leaf_ad, profit);
  AD dt = tf24_height_dt_from_net<AD>(p, h, area_leaf_ad, net);
  return xad::derivative(dt);
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
  if (use_ad_gradient) {
    // Exact gradient (default, #527): forward-mode AD over the analytic algebra
    // + IFT at the ci root-find + analytic spline derivatives for the transport.
    // No O(h) bias and no finite-difference step to tune.
    // Feasibility gate (mirrors the #526 FD branch below): when the collar solve
    // is infeasible -- drought shutdown (wettest soil drier than psi_crit) or a
    // collapsed interval -- prepare_collar_solve sets the shutdown operating point
    // (collar = -psi_crit) and there is no informative gradient. Crucially, the AD
    // transport-spline derivative would EXTRAPOLATE at that shutdown psi_crit (the
    // transpiration splines are non-extrapolating), a hard "Extrapolation disabled"
    // abort that bit long / dry patches (e.g. TF24f lifetime >= 8). Return 0 there,
    // exactly as the finite-difference path does; the leaf outputs are already at
    // the shutdown operating point. (#527 robustness / #472 scope B.)
    double bound_a, bound_b;
    if (!leaf.prepare_collar_solve(bound_a, bound_b)) {
      dprofit_dpsi_ = 0.0;
      return;
    }
    // Establish the operating point (and the clamped collar psi `used`) and leave
    // the leaf outputs there for compute_rates' aux reads.
    leaf.evaluate_root_collar_psi(tracked_root_psi_);
    const double used = -leaf.root_collar_psi_;
    dprofit_dpsi_ = leaf.dprofit_droot_collar_psi(used);
    leaf.evaluate_root_collar_psi(used);  // restore operating-point outputs
  } else {
    // Centred finite-difference fallback (#526), perturbing about the clamped
    // operating value `used`. A one-sided difference biases the fixed point to
    // psi* - h/2 (O(h)); the centred form cancels that term (O(h^2)) for one
    // extra leaf evaluation. Near a boundary the clamp collapses one arm,
    // degrading it gracefully to a one-sided difference that still points back
    // into the interior.
    const double h = psi_fd_step;
    if (!util::is_finite(h) || h <= 0.0) {
      util::stop("TF24f: psi_fd_step must be finite and > 0 (got " +
                 util::to_string(h) + ")");
    }
    // Share one prepare_collar_solve across the three profit evals (#530): the
    // plant/environment state is fixed within this solve, so the soil-side caches
    // and feasible interval are identical at `used` and `used ± h`. Re-deriving
    // them per eval (the old four-evaluate_root_collar_psi form) was the ~29%
    // cost over the forward difference.
    double bound_a, bound_b;
    if (!leaf.prepare_collar_solve(bound_a, bound_b)) {
      // Operating point fully determined by feasibility handling (shutdown /
      // assim<0 / collapsed interval); no interior interval to perturb in, so the
      // gradient is zero (matching the old form, where every clamped eval
      // returned the same fixed profit_). Leaf outputs are already at the
      // operating point.
      dprofit_dpsi_ = 0.0;
      return;
    }
    // `used` is the tracked state clamped into the feasible interval -- the same
    // value the old leading evaluate_root_collar_psi(tracked_root_psi_) produced.
    const double used = std::min(std::max(tracked_root_psi_, bound_a), bound_b);
    const double p_plus  = leaf.profit_at_collar_psi(used + h, bound_a, bound_b);
    const double p_minus = leaf.profit_at_collar_psi(used - h, bound_a, bound_b);
    dprofit_dpsi_ = (p_plus - p_minus) / (2.0 * h);
    leaf.profit_at_collar_psi(used, bound_a, bound_b);  // restore operating point
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
