// -*-c++-*-
#ifndef PLANT_PLANT_TF24T_STRATEGY_H_
#define PLANT_PLANT_TF24T_STRATEGY_H_

#include <plant/models/tf24_strategy.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <type_traits>

namespace plant {

// TF24t ("t" for thermal): TF24 with phylloptim's lasting leaf thermal-damage
// and acclimation layer (the ATLS model, #566) switched on, plus the four slow
// ODE states that layer reads. It inherits TF24_Strategy and reuses
// TF24_Environment and TF24_Pars.
//
// The leaf half is phylloptim's, behind Leaf::use_thermal_damage_: the shifted
// peaked Arrhenius, the survival discount on Vcmax and Jmax, the deactivated
// fraction phi_d_, the repair rate k_rec_eff_ gated off when hot, and four
// maintenance costs added to R_d_. This class turns the layer on, copies its
// traits into the leaf, and integrates, appended after TF24's states so the
// inherited indices are unchanged:
//
//   dA/dt     = alpha * softplus(T_air - t_accl) - beta * A
//   dI_r/dt   = k_i * phi_d * (1 - I_r - I_p) - (k_rec_eff + k_mat) * I_r
//   dI_p/dt   = k_mat * I_r - (k_l + max(0, dL/dt) / L) * I_p
//   dTbar/dt  = k_mean_leaf_temp * (Tleaf - Tbar)
//
// A is the thermostability (acclimation) state, I_r and I_p the recoverable and
// permanent damaged fractions, L the leaf area and Tbar a running mean of the
// operating-point leaf temperature that feeds back into nothing. The current A,
// I_r and I_p are handed to the leaf before each solve; phi_d and k_rec_eff are
// read back after it, at the operating-point leaf temperature.
//
// ⚠️ UNITS. Every rate constant here and on the leaf is per DAY, and the ODE
// clock is per YEAR, so each is multiplied by DAYS_PER_YEAR -- except leaf
// turnover k_l and the growth dilution in dI_p/dt, which are already per year.
//
// ⚠️ THE ENERGY BALANCE IS FORCED ON, whatever pars.use_energy_balance says:
// phi_d is evaluated at the operating-point leaf temperature, which exists only
// on the Penman-Monteith path. The `leaf_temp` driver is then the AIR
// temperature, and it is what the acclimation forcing reads.
//
// ⚠️ PHI_D DRIVES DAMAGE BUT ENTERS NO OBJECTIVE. The optimiser sees damage
// only through capacities already discounted by it, so a leaf cannot trade
// carbon now against damage later.
//
// ⚠️ NO DERIVATIVE PATH CARRIES THE LAYER. phylloptim refuses its derivative
// surface with the gate on, and rebind_from is deleted below, so no AD route
// reaches this model.
template <typename S = double>
class TF24t_Strategy : public TF24_Strategy<S> {
  // Only ever instantiated at double; see rebind_from below.
  static_assert(std::is_same_v<S, double>,
                "TF24t has no differentiable form: the leaf's thermal-damage "
                "layer carries no derivative");

public:
  using value_type = S;
  using base_type = TF24_Strategy<S>;

  typedef std::shared_ptr<TF24t_Strategy> ptr;
  TF24t_Strategy();

  // Puts the per-day rate constants on the per-year ODE clock.
  static constexpr double DAYS_PER_YEAR = 365.25;

  // Compound scientific version, "<TF24 version>.<thermal_revision>": the major
  // component tracks TF24_Strategy::scientific_version, so a TF24 change also
  // invalidates TF24t (the safe direction). Bump thermal_revision for a change
  // to the thermal layer alone. See src/strategy_version.cpp.
  static constexpr int thermal_revision = 3;

  // TF24's states and the four thermal states, appended last. Static, and
  // resolved on the concrete type by Individual<TF24t, ...>.
  static size_t state_size() { return base_type::state_size() + 4; }
  static std::vector<std::string> state_names() {
    std::vector<std::string> ret = base_type::state_names();
    ret.push_back("acclim_thermostab");
    ret.push_back("damage_recoverable");
    ret.push_back("damage_permanent");
    ret.push_back("mean_leaf_temp");
    return ret;
  }

  // The base refresh_indices() reads the base (static) state_names(), so it
  // misses the appended states; this registers them after it.
  void refresh_indices();

  // Runs TF24's preparation, refuses the deep-crown shading model, then turns
  // the leaf's energy balance and thermal-damage layer on and copies the
  // thermal traits into it. Hides the non-virtual base method; make_strategy_ptr
  // calls it on the concrete type.
  void prepare_strategy();

  // Hands the thermal states to the leaf, runs TF24's rates for the shared
  // states (which solve the leaf), then sets the four thermal rates.
  void compute_rates(const TF24_Environment<S>& environment, Internals<S>& vars);

  // dA/dt on the per-year clock, at air temperature T (deg C) and state A.
  double acclim_rate(double T, double A) const;

  // TF24's net production less the two whole-plant thermal costs the leaf's R_d_
  // cannot carry: the construction cost of the constitutive thermostability
  // offset, amortised through leaf turnover, and the induction cost of building
  // acclimation, proportional to the smoothed positive part of dA/dt. The leaf
  // maintenance costs are already in R_d_ and so in leaf.profit_; they are not
  // charged again here.
  //
  // Virtual in the base, so this is reached from the inherited compute_rates and
  // from establishment_probability. It reads the thermal states from the leaf,
  // which holds whatever compute_rates last handed it.
  using base_type::net_mass_production_dt;
  S net_mass_production_dt(const TF24_Environment<S>& environment,
                           const S& height, const S& area_leaf_,
                           const S& height_inverse) override;

  // Seeds the shared TF24 states, then A at its equilibrium in the birth air
  // temperature, I_r = I_p = 0, and the running mean at the birth air
  // temperature (the operating point is not solved until the first
  // compute_rates).
  void set_initial_states(const TF24_Environment<S>& environment,
                          Internals<S>& vars);

  // ⚠️ DELETED, NOT INHERITED. TF24_Strategy::rebind_from returns a
  // TF24_Strategy<U>, so an inherited one would rebind a TF24t patch as a TF24
  // patch: no thermal states, no layer, and a derivative of a different model
  // with every number finite. Deleting it makes Patch's at_scalar and odelia's
  // Rebindable read false, so the stand gradient, the sweep and an implicit
  // stepper's Jacobian all refuse at compile time instead.
  template <class U>
  void rebind_from() const = delete;

  // --- Leaf thermal traits, copied into the leaf in prepare_strategy ---------
  // Rate constants per day. Calibration targets, not measured values.
  double topt_offset = 0.0;  // constitutive shift of T_opt and damage onset, deg C
  double m_rep = 0.4;        // slope of the repair shut-off, deg C^-1
  double t_rep_cut = 45.0;   // midpoint of the repair shut-off, deg C
  double k_i = 0.5;          // leak phi_d -> I_r, day^-1
  double k_rec = 0.25;       // resynthesis of I_r when cool, day^-1
  double k_mat = 0.02;       // maturation I_r -> I_p, day^-1; read by the ODE only
  double dTopt_max = 6.0;    // largest acclimated rise in T_opt, deg C
  double K_A = 1.0;          // half-saturation of that rise in A

  // --- Acclimation kinetics --------------------------------------------------
  double alpha = 0.02;       // gain of the acclimation forcing, day^-1 deg C^-1
  double beta = 0.01;        // decay of A, day^-1
  double t_accl = 30.0;      // forcing threshold, deg C
  double softplus_s = 1.0;   // sharpness of the softplus forcing, deg C^-1

  // --- Running mean of the operating-point leaf temperature ------------------
  // Diagnostic only; 1/k is its memory, about 30 days at the default.
  double k_mean_leaf_temp = 0.0333;  // day^-1

  // --- Thermal costs ----------------------------------------------------------
  // The first four are copied into the leaf and added to R_d_ (umol C m^-2 s^-1
  // per unit of what each multiplies; see phylloptim's Leaf). The last three are
  // charged in net_mass_production_dt.
  double c_acclim_maint = 0.02;   // per unit A
  double c_repair_maint = 1e-5;   // per day^-1 of k_rec
  double c_protect_maint = 1e-5;  // per day^-1 of protection, k_i_ref - k_i
  double c_repair_flux = 1e-4;    // per day^-1 of realised repair flux
  double c_build = 0.02;          // construction premium on leaf turnover, per deg C of topt_offset
  double c_accl_induct = 0.01;    // biomass per unit positive dA/dt (yr^-1)
  double induct_eps = 1e-4;       // smoothing width of that positive part

  // Cached slots for the appended states, resolved in refresh_indices().
  int state_idx_acclim_A = -1;
  int state_idx_damage_r = -1;
  int state_idx_damage_p = -1;
  int state_idx_mean_leaf_temp = -1;

private:
  // A negative or non-finite rate would poison the state rates, so a
  // misconfigured kinetic stops the run.
  void check_kinetics() const;
};

template <typename S>
typename TF24t_Strategy<S>::ptr make_strategy_ptr(TF24t_Strategy<S> s);

template <typename S>
TF24t_Strategy<S>::TF24t_Strategy() {
  this->name = "TF24t";
  this->refresh_indices();
}

template <typename S>
void TF24t_Strategy<S>::refresh_indices() {
  base_type::refresh_indices();
  const int idx = static_cast<int>(base_type::state_size());
  this->state_index["acclim_thermostab"] = idx;
  this->state_index["damage_recoverable"] = idx + 1;
  this->state_index["damage_permanent"] = idx + 2;
  this->state_index["mean_leaf_temp"] = idx + 3;
  state_idx_acclim_A = idx;
  state_idx_damage_r = idx + 1;
  state_idx_damage_p = idx + 2;
  state_idx_mean_leaf_temp = idx + 3;
  // Again on the concrete model: the base checked its own names against its
  // own count, and the appended slots are outside both.
  check_state_layout(this->state_index, state_size(), "TF24t");
}

template <typename S>
void TF24t_Strategy<S>::check_kinetics() const {
  if (!util::is_finite(alpha) || alpha < 0.0 ||
      !util::is_finite(beta) || beta < 0.0 ||
      !util::is_finite(softplus_s) || softplus_s <= 0.0 ||
      !util::is_finite(k_mean_leaf_temp) || k_mean_leaf_temp < 0.0) {
    util::stop("TF24t: acclimation kinetics must be finite with "
               "alpha,beta,k_mean_leaf_temp>=0 and softplus_s>0");
  }
  if (!util::is_finite(k_i) || k_i < 0.0 ||
      !util::is_finite(k_rec) || k_rec < 0.0 ||
      !util::is_finite(k_mat) || k_mat < 0.0) {
    util::stop("TF24t: damage kinetics k_i,k_rec,k_mat must be finite and >= 0");
  }
}

template <typename S>
void TF24t_Strategy<S>::prepare_strategy() {
  base_type::prepare_strategy();
  // phi_d_ and k_rec_eff_ are read after the crown loop. Deep-crown shading
  // solves the leaf once per quadrature node and integrates only the outputs
  // TF24 names, so these two would hold the last node's values rather than a
  // crown mean.
  if (this->shading_model_ == ShadingModel::DeepCrown) {
    throw std::invalid_argument(
      "shading_model '" + this->control.shading_model + "' is not supported "
      "for the TF24t strategy: the leaf's phi_d and k_rec_eff are read after "
      "the crown loop, which under deep-crown shading leaves the last "
      "quadrature node's values rather than a crown mean. Use mean-light (the "
      "default) or crown-centre.");
  }
  Leaf& leaf = this->leaf;
  leaf.use_energy_balance_ = true;
  leaf.use_thermal_damage_ = true;
  leaf.topt_offset_ = topt_offset;
  leaf.m_rep_ = m_rep;
  leaf.t_rep_cut_ = t_rep_cut;
  leaf.k_i_ = k_i;
  leaf.k_rec_ = k_rec;
  // k_i_ref_ keeps the leaf's default: the no-extra-protection baseline the
  // protection cost is measured from.
  leaf.dTopt_max_ = dTopt_max;
  leaf.K_A_ = K_A;
  leaf.c_acclim_maint_ = c_acclim_maint;
  leaf.c_repair_maint_ = c_repair_maint;
  leaf.c_protect_maint_ = c_protect_maint;
  leaf.c_repair_flux_ = c_repair_flux;
  // What a solve before the first compute_rates reads; compute_rates overwrites
  // all three with the plant's own states.
  leaf.A_ = 0.0;
  leaf.I_r_ = 0.0;
  leaf.I_p_ = 0.0;
}

template <typename S>
double TF24t_Strategy<S>::acclim_rate(double T, double A) const {
  return DAYS_PER_YEAR *
    (alpha * Leaf::softplus_(T - t_accl, softplus_s) - beta * A);
}

template <typename S>
S TF24t_Strategy<S>::net_mass_production_dt(const TF24_Environment<S>& environment,
                                            const S& height, const S& area_leaf_,
                                            const S& height_inverse) {
  const S base = base_type::net_mass_production_dt(environment, height,
                                                   area_leaf_, height_inverse);

  const S turnover_leaf_ = this->turnover_leaf(this->mass_leaf(area_leaf_));
  const S construction = c_build * std::max(0.0, topt_offset) * turnover_leaf_;

  // Smooth positive part 0.5*(x + sqrt(x^2 + eps^2)), the form TF24 uses for
  // net production, so the integrated carbon balance stays smooth in dA/dt.
  // ⚠️ IT IS NOT ZERO AT dA/dt = 0: it is eps/2 there, so a plant at its
  // acclimation equilibrium still pays c_accl_induct * induct_eps / 2.
  const double dA = acclim_rate(environment.get_leaf_temp(), this->leaf.A_);
  const double dA_pos =
    0.5 * (dA + std::sqrt(dA * dA + induct_eps * induct_eps));
  const double induction = c_accl_induct * dA_pos;

  return base - construction - induction;
}

template <typename S>
void TF24t_Strategy<S>::compute_rates(const TF24_Environment<S>& environment,
                                      Internals<S>& vars) {
  check_kinetics();

  const double A = vars.state(state_idx_acclim_A);
  const double I_r =
    std::min(std::max(static_cast<double>(vars.state(state_idx_damage_r)), 0.0), 1.0);
  const double I_p =
    std::min(std::max(static_cast<double>(vars.state(state_idx_damage_p)), 0.0), 1.0);
  Leaf& leaf = this->leaf;
  leaf.A_ = std::max(A, 0.0);
  leaf.I_r_ = I_r;
  leaf.I_p_ = I_p;

  base_type::compute_rates(environment, vars);

  vars.set_rate(state_idx_acclim_A, acclim_rate(environment.get_leaf_temp(), A));

  // phi_d_ and k_rec_eff_ are the leaf's outputs at the operating point the
  // base rates just solved.
  const double survivable = std::max(0.0, 1.0 - I_r - I_p);
  vars.set_rate(state_idx_damage_r,
                DAYS_PER_YEAR * (k_i * leaf.phi_d_ * survivable -
                                 (leaf.k_rec_eff_ + k_mat) * I_r));

  // Permanent damage is cleared only by dilution with new, undamaged leaf: leaf
  // turnover plus the fractional rate of canopy expansion, both per year. The
  // expansion is recovered from the height rate the base just set.
  const double height = vars.state(HEIGHT_INDEX);
  const double area_leaf_ = this->area_leaf(height);
  const double area_leaf_dt =
    vars.rate(HEIGHT_INDEX) / this->dheight_darea_leaf(area_leaf_);
  const double dilution =
    this->pars.k_l + std::max(0.0, area_leaf_dt) / area_leaf_;
  vars.set_rate(state_idx_damage_p,
                DAYS_PER_YEAR * k_mat * I_r - dilution * I_p);

  // Relaxes toward the leaf's own temperature at the operating point, which is
  // also TF24's `Tleaf` aux.
  vars.set_rate(state_idx_mean_leaf_temp,
                DAYS_PER_YEAR * k_mean_leaf_temp *
                  (leaf.Tleaf_ - vars.state(state_idx_mean_leaf_temp)));
}

// A_eq = (alpha / beta) * softplus(T - t_accl): the DAYS_PER_YEAR factor cancels
// in the ratio, and with beta = 0 there is no equilibrium, so A starts at zero.
template <typename S>
void TF24t_Strategy<S>::set_initial_states(const TF24_Environment<S>& environment,
                                           Internals<S>& vars) {
  // Non-virtual in the base, so the shared states (notably the storage pool)
  // are seeded only if this calls it.
  base_type::set_initial_states(environment, vars);
  check_kinetics();

  const double T = environment.get_leaf_temp();
  const double A_eq =
    beta > 0.0 ? alpha / beta * Leaf::softplus_(T - t_accl, softplus_s) : 0.0;
  vars.set_state(state_idx_acclim_A, A_eq);
  vars.set_state(state_idx_damage_r, 0.0);
  vars.set_state(state_idx_damage_p, 0.0);
  vars.set_state(state_idx_mean_leaf_temp, T);
}

template <typename S>
typename TF24t_Strategy<S>::ptr make_strategy_ptr(TF24t_Strategy<S> s) {
  s.prepare_strategy();
  return std::make_shared<TF24t_Strategy<S> >(s);
}

}

#endif
