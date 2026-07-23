#include <plant/models/tf24t_strategy.h>

namespace plant {

// The base TF24_Strategy constructor sets collect_all_auxiliary and runs the
// base refresh_indices(); we re-run our own so the appended state slots are
// registered, and set the reported name.
TF24t_Strategy::TF24t_Strategy() {
  name = "TF24t";
  refresh_indices();
}

// Build on the base index maps, then register the three appended thermal states
// after TF24's states (indices TF24 state_size(), +1, +2).
void TF24t_Strategy::refresh_indices() {
  TF24_Strategy::refresh_indices();
  const int idx = static_cast<int>(TF24_Strategy::state_size());
  state_index["acclim_thermostab"] = idx;
  state_index["damage_recoverable"] = idx + 1;
  state_index["damage_permanent"] = idx + 2;
  state_idx_acclim_A = idx;
  state_idx_damage_r = idx + 1;
  state_idx_damage_p = idx + 2;
}

// Fail fast on a misconfigured kinetic (a negative/non-finite gain or rate would
// poison the state rates); shared guard for compute_rates and set_initial_states.
namespace {
void check_kinetics(double alpha, double beta, double s) {
  if (!util::is_finite(alpha) || alpha < 0.0 ||
      !util::is_finite(beta) || beta < 0.0 ||
      !util::is_finite(s) || s <= 0.0) {
    util::stop("TF24t: acclimation kinetics must be finite with alpha,beta>=0 and softplus_s>0");
  }
}
void check_damage_rates(double k_i, double k_rec, double k_mat) {
  if (!util::is_finite(k_i) || k_i < 0.0 ||
      !util::is_finite(k_rec) || k_rec < 0.0 ||
      !util::is_finite(k_mat) || k_mat < 0.0) {
    util::stop("TF24t: damage kinetics k_i,k_rec,k_mat must be finite and >= 0");
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
  leaf.m_rep_ = m_rep;
  leaf.t_rep_cut_ = t_rep_cut;
  leaf.k_i_ = k_i;
  leaf.k_rec_ = k_rec;
  leaf.k_mat_ = k_mat;
  // leaf.k_i_ref_ keeps its fixed default (the no-extra-protection baseline);
  // a strategy that buys protection lowers k_i below it and pays c_protect_maint.
  leaf.dTopt_max_ = dTopt_max;
  leaf.K_A_ = K_A;
  // Leaf-side thermal-cost coefficients: acclimation maintenance, repair/protection
  // standing maintenance, repair activity. Added to leaf.R_d_.
  leaf.c_acclim_maint_  = c_acclim_maint;
  leaf.c_repair_maint_  = c_repair_maint;
  leaf.c_protect_maint_ = c_protect_maint;
  leaf.c_repair_flux_   = c_repair_flux;
  // Start the leaf's thermal-state inputs at 0 (compute_rates overwrites them
  // with the live states each step; this only defines what
  // establishment_probability reads before any rates run).
  leaf.A_ = 0.0;
  leaf.I_r_ = 0.0;
  leaf.I_p_ = 0.0;
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

// Acclimation forcing rate dA/dt on the yearly ODE clock (day^-1 constants
// converted via DAYS_PER_YEAR): dA/dt = DAYS_PER_YEAR*(alpha*sp(T-t_accl) - beta*A).
double TF24t_Strategy::acclim_rate(double T, double A) const {
  return TF24t_Strategy::DAYS_PER_YEAR *
    (alpha * Leaf::softplus_(T - t_accl, softplus_s) - beta * A);
}

// Whole-plant carbon balance = base TF24 net production, minus the thermal costs
// that are not leaf-local respiration (those already sit in leaf.R_d_):
//  - tolerance construction cost: a one-off build cost per leaf, amortised
//    through leaf turnover, charged for the thermostability offset bought above
//    the intrinsic optimum (topt_offset > 0). One axis (merged revision).
//  - acclimation induction cost: proportional to the smoothed positive build
//    rate of the single thermostability state (paying to synthesise proteins).
// leaf.A_ holds the current state (compute_rates set it before the inherited
// compute_rates called this); at establishment it reads 0, a negligible
// build-from-zero snapshot at default temperatures.
double TF24t_Strategy::net_mass_production_dt(const TF24_Environment& environment,
                                             double height, double area_leaf_,
                                             double height_inverse) {
  const double base =
    TF24_Strategy::net_mass_production_dt(environment, height, area_leaf_,
                                          height_inverse);

  const double turnover_leaf_ = turnover_leaf(mass_leaf(area_leaf_));
  const double construction =
    c_build * std::max(0.0, topt_offset) * turnover_leaf_;

  const double T = environment.get_leaf_temp();
  const double dA = acclim_rate(T, leaf.A_);
  const double induction = c_accl_induct * positive_part(dA, induct_eps);

  return base - construction - induction;
}

// Feed the current thermal states into the leaf (so this step's capacity
// discount and d_S shift see them), reuse TF24's rates for the shared states,
// then set the three thermal-state rates. phi_d and the gated restorative rate
// k_rec_eff are read back from the leaf at the operating-point Tleaf that the
// base rates just solved. Rate constants in day^-1 are put on the yearly ODE
// clock via DAYS_PER_YEAR (except g/L, which is already per-year).
void TF24t_Strategy::compute_rates(const TF24_Environment& environment,
                                   Internals& vars) {
  check_kinetics(alpha, beta, softplus_s);
  check_damage_rates(k_i, k_rec, k_mat);

  const double A = vars.state(state_idx_acclim_A);
  const double I_r = std::min(std::max(vars.state(state_idx_damage_r), 0.0), 1.0);
  const double I_p = std::min(std::max(vars.state(state_idx_damage_p), 0.0), 1.0);
  leaf.A_ = std::max(A, 0.0);
  leaf.I_r_ = I_r;
  leaf.I_p_ = I_p;

  TF24_Strategy::compute_rates(environment, vars);

  const double T = environment.get_leaf_temp();

  // Acclimation.
  vars.set_rate(state_idx_acclim_A, acclim_rate(T, A));

  // Recoverable damage I_r: preventive leak from the transiently-unfolded pool
  // phi_d into not-yet-damaged protein (1-I_r-I_p), minus resynthesis and
  // maturation. phi_d and k_rec_eff are the operating-point leaf outputs.
  const double survivable = std::max(0.0, 1.0 - I_r - I_p);
  const double dI_r = DAYS_PER_YEAR *
    (k_i * leaf.phi_d_ * survivable - (leaf.k_rec_eff_ + k_mat) * I_r);
  vars.set_rate(state_idx_damage_r, dI_r);

  // Permanent damage I_p: maturation of I_r, cleared by dilution with clean new
  // leaf. g/L = leaf turnover rate (mature limit) + fractional canopy expansion.
  // area_leaf_dt is recovered from the height rate the base just set.
  const double height = vars.state(HEIGHT_INDEX);
  const double area_leaf_ = area_leaf(height);
  const double area_leaf_dt = vars.rate(HEIGHT_INDEX) / dheight_darea_leaf(area_leaf_);
  const double dilution = pars.k_l + std::max(0.0, area_leaf_dt) / area_leaf_;
  const double dI_p = DAYS_PER_YEAR * k_mat * I_r - dilution * I_p;
  vars.set_rate(state_idx_damage_p, dI_p);
}

// Seed the thermostability state at its environmental equilibrium
// (A_eq = (alpha/beta) * softplus(T - t_accl); the DAYS_PER_YEAR factor cancels
// in the ratio) so a recruit is born already acclimated to its birth
// environment -- no cold-start transient. The lasting-damage pools start empty
// (fresh recruits have undamaged leaves). Seed the shared TF24 states (notably
// the NSC storage pool, #517) first.
void TF24t_Strategy::set_initial_states(const TF24_Environment& environment,
                                        Internals& vars) {
  TF24_Strategy::set_initial_states(environment, vars);
  check_kinetics(alpha, beta, softplus_s);
  check_damage_rates(k_i, k_rec, k_mat);

  const double T = environment.get_leaf_temp();
  const double A_eq =
    beta > 0.0 ? alpha / beta * Leaf::softplus_(T - t_accl, softplus_s) : 0.0;
  vars.set_state(state_idx_acclim_A, A_eq);
  vars.set_state(state_idx_damage_r, 0.0);
  vars.set_state(state_idx_damage_p, 0.0);
}

TF24t_Strategy::ptr make_strategy_ptr(TF24t_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<TF24t_Strategy>(s);
}

}
