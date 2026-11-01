// -*-c++-*-
#ifndef PLANT_PLANT_TF24T_STRATEGY_H_
#define PLANT_PLANT_TF24T_STRATEGY_H_

#include <plant/models/tf24_strategy.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <type_traits>

namespace plant {

// TF24t ("t" for thermal): TF24 with phylloptim's lasting leaf thermal-damage
// and acclimation layer (the ATLS model, #566) switched on, plus the four slow
// ODE states that layer reads. It inherits TF24_Strategy and reuses
// TF24_Environment and TF24_Pars.
//
// THE LAYER IS phylloptim's, behind Leaf::use_thermal_damage_, and so are its
// DYNAMICS: three fractions of leaf area, functional F <=> impaired I -> damaged
// D (F = 1 - I - D), plus an acclimation state A, all on the LEAF temperature.
// Leaf::damage_rates() is the right-hand side, written once there:
//
//   dA/dt     = alpha * softplus(T_leaf - t_accl; s) - beta * A
//   dI/dt     = k_imp(T_leaf) * F - (k_rec_eff(T_leaf) + k_mat) * I
//   dD/dt     = k_mat * I                         - (k_l + max(0, dL/dt)/L) * D
//   dTbar/dt  = k_mean_leaf_temp * (Tleaf - Tbar)
//
// What THIS class owns is only what a whole plant adds: the dilution of damaged
// leaf by turnover and canopy growth (the last term of dD/dt), the per-day to
// per-year clock, the running mean Tbar (feeds back into nothing), and two
// whole-plant costs. It hands A, I and D to the leaf before each solve, so the
// solve scales both capacities by F, and evaluates damage_rates() at the
// operating point's leaf temperature after it.
//
// ⚠️ UNITS. Every rate constant here and on the leaf is per DAY, and the ODE
// clock is per YEAR, so the leaf's rates are multiplied by DAYS_PER_YEAR --
// dilution is not, being per year already.
//
// ⚠️ THE ENERGY BALANCE IS FORCED ON, whatever pars.use_energy_balance says:
// the dynamics read the operating-point leaf temperature, which exists only on
// the Penman-Monteith path. The `leaf_temp` driver is then the AIR temperature;
// nothing here reads it except the equilibrium A at birth, before any solve.
//
// ⚠️ THE DAMAGE RATES ENTER NO OBJECTIVE. The optimiser sees the states only
// through capacities already scaled by F, so a leaf cannot trade carbon now
// against impairment later.
//
// ⚠️ NO DERIVATIVE PATH CARRIES THIS MODEL. rebind is deleted below, so no
// AD route reaches it. (phylloptim's own trait gradient does carry the layer; a
// TF24t stand gradient would need the four states on the tape, which is not
// built.)
template <typename S = double>
class TF24t_Strategy : public TF24_Strategy<S> {
  // Only ever instantiated at double; see rebind below.
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
  static constexpr int thermal_revision = 5;

  // TF24's states and the four thermal states, appended last. Static, and
  // resolved on the concrete type by Individual<TF24t, ...>.
  static size_t state_size() { return base_type::state_size() + 4; }
  static std::vector<std::string> state_names() {
    std::vector<std::string> ret = base_type::state_names();
    ret.push_back("acclim_thermostab");
    ret.push_back("impaired");
    ret.push_back("damaged");
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

  // The leaf's damage rates at its operating-point temperature, per YEAR:
  // {dA/dt, dI/dt, dD/dt}, dD/dt without dilution. Reads the leaf's seated A_,
  // I_ and D_, which compute_rates hands it.
  std::array<double, 3> leaf_damage_rates_per_year() const;

  // TF24's net production less the two whole-plant thermal costs the leaf's R_d_
  // cannot carry: the construction cost of the constitutive thermostability
  // offset, amortised through leaf turnover, and the induction cost of building
  // acclimation, proportional to leaf area and the smoothed positive part of
  // dA/dt. The leaf
  // maintenance costs are already in R_d_ and so in leaf.profit_; they are not
  // charged again here.
  //
  // Virtual in the base, so this is reached from the inherited compute_rates and
  // from establishment_probability. It reads the thermal states from the leaf,
  // which holds whatever compute_rates last handed it.
  using base_type::net_mass_production_dt;
  S net_mass_production_dt(const TF24_Environment<S>& environment,
                           const S& height, const S& area_leaf_,
                           const S& height_inverse,
                           const S& sapwood_departure,
                           const S& leaf_departure,
                           const S& tracked_marginal) override;

  // Seeds the shared TF24 states, then A at its equilibrium at the birth air
  // temperature -- the operating point, and so the leaf temperature, is not
  // solved until the first compute_rates -- I = D = 0, and the running mean at
  // the birth air temperature.
  void set_initial_states(const TF24_Environment<S>& environment,
                          Internals<S>& vars);

  // ⚠️ DELETED, NOT INHERITED. TF24_Strategy::rebind returns a
  // TF24_Strategy<U>, so an inherited one would rebind a TF24t patch as a TF24
  // patch: no thermal states, no layer, and a derivative of a different model
  // with every number finite. Deleting it makes Patch's at_scalar and odelia's
  // Rebindable read false, so the stand gradient, the sweep and an implicit
  // stepper's Jacobian all refuse at compile time instead.
  template <class U>
  void rebind() const = delete;

  // --- Leaf thermal traits, copied into the leaf in prepare_strategy ---------
  // Rate constants per day. Calibration targets, not measured values. The leaf's
  // field each one sets is named beside it; see phylloptim's Leaf.
  // Impairment is Bison & Michaletz's first-order heat inactivation, centred on
  // the isokinetic temperature; repair and maturation are centred at 25 C. The
  // defaults are the leaf's: their 177-species fit for impairment, placeholders
  // for the other two.
  double topt_offset = 0.0;  // constitutive shift of T_opt and of k_i_Tiso, deg C
  double k_i = 276.0;        // impairment rate at the isokinetic temperature, day^-1
  double k_i_Ea = 316.0e3;   // its activation energy, J mol^-1
  double k_i_Tiso = 50.5;    // the isokinetic temperature, deg C
  double k_rec = 0.25;       // repair I -> F at 25 C, day^-1
  double k_rec_Ea = 0.0;     // its activation energy, J mol^-1 (0: constant)
  double k_mat = 0.02;       // maturation I -> D at 25 C, day^-1
  double k_mat_Ea = 0.0;     // its activation energy, J mol^-1 (0: constant)
  double dTopt_max = 6.0;    // largest acclimated rise in T_opt, deg C
  double dTimpair_max = 6.0; // largest acclimated rise in k_i_Tiso, deg C
  double K_A = 1.0;          // half-saturation of both rises in A

  // --- Acclimation kinetics, on the LEAF temperature -------------------------
  double alpha = 0.02;       // gain of the acclimation forcing, day^-1 deg C^-1 (acclim_alpha_)
  double beta = 0.01;        // decay of A, day^-1 (acclim_beta_)
  double t_accl = 30.0;      // forcing threshold, deg C (acclim_T_)
  double softplus_s = 1.0;   // sharpness of the softplus forcing, deg C^-1 (acclim_s_)

  // --- Running mean of the operating-point leaf temperature ------------------
  // Diagnostic only; 1/k is its memory, about 30 days at the default.
  double k_mean_leaf_temp = 0.0333;  // day^-1

  // --- Thermal costs ----------------------------------------------------------
  // The first four are copied into the leaf and added to R_d_ (umol C m^-2 s^-1
  // per unit of what each multiplies; see phylloptim's Leaf). The last three are
  // charged in net_mass_production_dt.
  double c_acclim_maint = 0.02;   // per unit A
  double c_repair_maint = 1e-5;   // per day^-1 of k_rec
  double c_protect_maint = 1e-5;  // per day^-1 of protection, k_i_baseline_ - k_i
  double c_repair_flux = 1e-4;    // per day^-1 of realised repair flux
  double c_build = 0.02;          // construction premium on leaf turnover, per deg C of topt_offset
  double c_accl_induct = 0.01;    // kg per m^2 of leaf per unit positive dA/dt (yr^-1)
  double induct_eps = 1e-4;       // smoothing width of that positive part

  // Cached slots for the appended states, resolved in refresh_indices().
  int state_idx_acclim_A = -1;
  int state_idx_impaired = -1;
  int state_idx_damaged = -1;
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
  this->state_index["impaired"] = idx + 1;
  this->state_index["damaged"] = idx + 2;
  this->state_index["mean_leaf_temp"] = idx + 3;
  state_idx_acclim_A = idx;
  state_idx_impaired = idx + 1;
  state_idx_damaged = idx + 2;
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
  if (!util::is_finite(k_i_Ea) || k_i_Ea < 0.0 ||
      !util::is_finite(k_rec_Ea) || k_rec_Ea < 0.0 ||
      !util::is_finite(k_mat_Ea) || k_mat_Ea < 0.0) {
    util::stop("TF24t: activation energies k_i_Ea,k_rec_Ea,k_mat_Ea must be "
               "finite and >= 0");
  }
  if (!util::is_finite(k_i_Tiso)) {
    util::stop("TF24t: k_i_Tiso must be finite");
  }
}

template <typename S>
void TF24t_Strategy<S>::prepare_strategy() {
  base_type::prepare_strategy();
  // The damage rates are evaluated at the leaf temperature after the crown
  // loop. Deep-crown shading solves the leaf once per quadrature node and
  // integrates only the outputs TF24 names, so they would be the last node's
  // rather than a crown mean.
  if (this->shading_model_ == ShadingModel::DeepCrown) {
    throw std::invalid_argument(
      "shading_model '" + this->control.shading_model + "' is not supported "
      "for the TF24t strategy: the leaf's damage rates are evaluated after the "
      "crown loop, which under deep-crown shading leaves the last quadrature "
      "node's leaf temperature rather than a crown mean. Use mean-light (the "
      "default) or crown-centre.");
  }
  Leaf& leaf = this->leaf;
  leaf.use_energy_balance_ = true;
  leaf.use_thermal_damage_ = true;
  leaf.topt_offset_ = topt_offset;
  leaf.k_i_ = k_i;
  leaf.k_i_Ea_ = k_i_Ea;
  leaf.k_i_Tiso_ = k_i_Tiso;
  leaf.k_rec_ = k_rec;
  leaf.k_rec_Ea_ = k_rec_Ea;
  leaf.k_mat_ = k_mat;
  leaf.k_mat_Ea_ = k_mat_Ea;
  // k_i_baseline_ keeps the leaf's default: the no-extra-protection baseline the
  // protection cost is measured from.
  leaf.dTopt_max_ = dTopt_max;
  leaf.dTimpair_max_ = dTimpair_max;
  leaf.K_A_ = K_A;
  leaf.acclim_alpha_ = alpha;
  leaf.acclim_beta_ = beta;
  leaf.acclim_T_ = t_accl;
  leaf.acclim_s_ = softplus_s;
  leaf.c_acclim_maint_ = c_acclim_maint;
  leaf.c_repair_maint_ = c_repair_maint;
  leaf.c_protect_maint_ = c_protect_maint;
  leaf.c_repair_flux_ = c_repair_flux;
  // What a solve before the first compute_rates reads; compute_rates overwrites
  // all three with the plant's own states.
  leaf.A_ = 0.0;
  leaf.I_ = 0.0;
  leaf.D_ = 0.0;
}

template <typename S>
std::array<double, 3> TF24t_Strategy<S>::leaf_damage_rates_per_year() const {
  const Leaf& leaf = this->leaf;
  const std::vector<double> r =
    leaf.damage_rates(leaf.Tleaf_, leaf.A_, leaf.I_, leaf.D_);
  return {DAYS_PER_YEAR * r[0], DAYS_PER_YEAR * r[1], DAYS_PER_YEAR * r[2]};
}

template <typename S>
S TF24t_Strategy<S>::net_mass_production_dt(const TF24_Environment<S>& environment,
                                            const S& height, const S& area_leaf_,
                                            const S& height_inverse,
                                            const S& sapwood_departure,
                                            const S& leaf_departure,
                                            const S& tracked_marginal) {
  const S base = base_type::net_mass_production_dt(
      environment, height, area_leaf_, height_inverse, sapwood_departure,
      leaf_departure, tracked_marginal);

  const S turnover_leaf_ = this->turnover_leaf(this->mass_leaf(area_leaf_));
  const S construction = c_build * std::max(0.0, topt_offset) * turnover_leaf_;

  // Smooth positive part 0.5*(x + sqrt(x^2 + eps^2)), the form TF24 uses for
  // net production, so the integrated carbon balance stays smooth in dA/dt.
  // ⚠️ IT IS NOT ZERO AT dA/dt = 0: it is eps/2 there, so a plant at its
  // acclimation equilibrium still pays c_accl_induct * induct_eps / 2 per m^2.
  // At the operating point's leaf temperature: base_type's rate evaluation above
  // solved the leaf.
  //
  // ⚠️ PER UNIT LEAF AREA, because A is built per unit leaf. Charged per plant,
  // the cost does not shrink with the plant: a newborn at 28 C air, whose leaf
  // runs 4 C warmer than the air its A was seeded at, has dA/dt = 16 yr^-1 and
  // would pay 0.16 kg yr^-1 against a net production of 9e-5, which starves
  // every seedling and holds the stepper to tiny steps.
  const double dA = leaf_damage_rates_per_year()[0];
  const double dA_pos =
    0.5 * (dA + std::sqrt(dA * dA + induct_eps * induct_eps));
  const double induction = c_accl_induct * dA_pos * area_leaf_;

  return base - construction - induction;
}

template <typename S>
void TF24t_Strategy<S>::compute_rates(const TF24_Environment<S>& environment,
                                      Internals<S>& vars) {
  check_kinetics();

  const double A = vars.state(state_idx_acclim_A);
  const double I =
    std::min(std::max(static_cast<double>(vars.state(state_idx_impaired)), 0.0), 1.0);
  const double D =
    std::min(std::max(static_cast<double>(vars.state(state_idx_damaged)), 0.0), 1.0);
  Leaf& leaf = this->leaf;
  leaf.A_ = A;
  leaf.I_ = I;
  leaf.D_ = D;

  base_type::compute_rates(environment, vars);

  // The leaf's own dynamics at the operating point the base rates just solved.
  // A is handed over as the state, negative included: the leaf reads a negative
  // A as zero for the shift and decays it as given.
  const std::array<double, 3> r = leaf_damage_rates_per_year();
  vars.set_rate(state_idx_acclim_A, r[0]);
  vars.set_rate(state_idx_impaired, r[1]);

  // Damaged leaf is cleared only by dilution with new, undamaged leaf: leaf
  // turnover plus the fractional rate of canopy expansion, both per year. The
  // expansion is recovered from the height rate the base just set.
  const double height = vars.state(HEIGHT_INDEX);
  const double area_leaf_ = this->area_leaf(height);
  const double area_leaf_dt =
    vars.rate(HEIGHT_INDEX) / this->dheight_darea_leaf(area_leaf_);
  const double dilution =
    this->pars.k_l + std::max(0.0, area_leaf_dt) / area_leaf_;
  vars.set_rate(state_idx_damaged, r[2] - dilution * D);

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
  vars.set_state(state_idx_impaired, 0.0);
  vars.set_state(state_idx_damaged, 0.0);
  vars.set_state(state_idx_mean_leaf_temp, T);
}

template <typename S>
typename TF24t_Strategy<S>::ptr make_strategy_ptr(TF24t_Strategy<S> s) {
  s.prepare_strategy();
  return std::make_shared<TF24t_Strategy<S> >(s);
}

}

#endif
