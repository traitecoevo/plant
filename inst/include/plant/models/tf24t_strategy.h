// -*-c++-*-
#ifndef PLANT_PLANT_TF24t_STRATEGY_H_
#define PLANT_PLANT_TF24t_STRATEGY_H_

#include <plant/models/tf24_strategy.h>

namespace plant {

// TF24t ("t" for thermal): a variant of TF24 that adds a mechanistic leaf
// thermal damage / repair / acclimation layer (the ATLS / Lumry-Eyring model,
// issue #566) on top of TF24's Penman-Monteith leaf energy balance (#523). It
// inherits TF24_Strategy and reuses TF24_Environment + TF24_Pars; only the
// acclimation ODE states and the leaf-thermal wiring are added.
//
// The Leaf carries the merged-revision damage machinery (use_thermal_damage_,
// the deactivated fraction phi_d from the unified Medlyn curve, the lasting-
// damage pools I_r_/I_p_, and the tolerance/repair traits); TF24t (a) turns that
// layer on and copies its traits into the leaf in prepare_strategy, and (b) adds
// four slow ODE states -- acclim_thermostab (a single thermostability state A
// shifting d_S, moving T_opt AND the damage onset together), damage_recoverable
// (I_r, resynthesised over ~4 d), damage_permanent (I_p, cleared only by
// new-leaf dilution), and mean_leaf_temp (a diagnostic exponential moving average
// of the PM operating-point leaf temperature, no feedback) -- appended after
// TF24's states so the inherited indices are unchanged. The damage/acclimation
// states relax as
//   dA/dt   = alpha * softplus(T - t_accl) - beta * A
//   dI_r/dt = k_i * phi_d * (1-I_r-I_p) - (k_rec_eff + k_mat) * I_r
//   dI_p/dt = k_mat * I_r - (g/L) * I_p
// (softplus/logistic, so the forcing is C1). The current state values are fed
// into the leaf each step (leaf.A_/leaf.I_r_/leaf.I_p_) before the reused TF24
// rates run, so the capacity discount and d_S shift see them; phi_d (evaluated
// at the operating-point Tleaf) is read back to drive the I_r leak. Day^-1 rate
// constants are converted to the yearly ODE clock via DAYS_PER_YEAR.
class TF24t_Strategy : public TF24_Strategy {
public:
  typedef std::shared_ptr<TF24t_Strategy> ptr;
  TF24t_Strategy();

  // Convert the day^-1 damage/acclimation rate constants to the SCM's yearly
  // integration clock (k_l is /yr, patch lifetime in years).
  static constexpr double DAYS_PER_YEAR = 365.25;

  // Compound scientific version: "<TF24 version>.<thermal_revision>", so a TF24
  // science change also invalidates TF24t (the safe direction); bump
  // thermal_revision for changes specific to the thermal layer. Revision 2 is
  // the merged-revision rewrite (unified curve + lasting I_r/I_p pools);
  // revision 3 adds the diagnostic running-mean leaf-temperature state.
  static constexpr int thermal_revision = 3;

  // TF24's states + the four appended thermal states (indices unchanged for
  // 0..base state_size()-1). Statics resolve on the concrete type in Individual<>.
  static size_t state_size() { return TF24_Strategy::state_size() + 4; }
  static std::vector<std::string> state_names() {
    std::vector<std::string> ret = TF24_Strategy::state_names();
    ret.push_back("acclim_thermostab");
    ret.push_back("damage_recoverable");
    ret.push_back("damage_permanent");
    ret.push_back("mean_leaf_temp");
    return ret;
  }

  // TF24's aux + one appended diagnostic: the PM operating-point leaf temperature
  // (the temperature the running-mean state relaxes toward). aux_names()/aux_size()
  // are non-virtual, so we hide the base versions; Individual<> holds a concrete
  // TF24t_Strategy::ptr, so name lookup binds to these (and the appended aux is
  // sized and copied out correctly).
  std::vector<std::string> aux_names();
  size_t aux_size() { return aux_names().size(); }

  // Base refresh_indices() uses the base (static) state_names(); re-run it then
  // register the appended state and aux slots.
  void refresh_indices();

  // Enable the leaf thermal-damage layer and copy the thermal traits into the
  // leaf (mirrors TF24_Strategy::prepare_strategy's use_energy_balance_ / d_
  // handling). Hides (does not override) the non-virtual base method; called on
  // the concrete type by make_strategy_ptr.
  void prepare_strategy();

  // Reuse TF24's rates for the shared states, feed the thermal states into the
  // leaf, then set the three thermal-state rates.
  void compute_rates(const TF24_Environment& environment, Internals& vars);

  // Acclimation forcing rate on the yearly ODE clock (day^-1 -> /yr):
  //   dA/dt = DAYS_PER_YEAR * (alpha * softplus(T - t_accl) - beta * A).
  double acclim_rate(double T, double A) const;

  // Whole-plant carbon balance with the thermal costs the leaf R_d_ cannot
  // carry (Phase 3, #566): the tolerance *construction* cost (a one-off build
  // cost per leaf, amortised through leaf turnover) and the acclimation
  // *induction* cost (proportional to the smoothed positive build rate of the
  // acclimation states). Base net_mass_production_dt is virtual, so this
  // override is reached even from the inherited compute_rates and from
  // establishment_probability. The leaf maintenance/activity costs live in
  // Leaf::R_d_ and flow through leaf.profit_, so they are NOT re-applied here.
  // The using-declaration keeps the base's (env, Internals&) convenience
  // overload visible (this 4-arg declaration would otherwise hide it).
  using TF24_Strategy::net_mass_production_dt;
  double net_mass_production_dt(const TF24_Environment& environment,
                                double height, double area_leaf_,
                                double height_inverse) override;

  // Seed the acclimation states at their environmental equilibrium (no birth
  // transient), after seeding the shared TF24 states.
  void set_initial_states(const TF24_Environment& environment, Internals& vars);

  // --- Thermal-damage leaf traits (copied into the leaf in prepare_strategy) --
  // Defaults mirror the Leaf defaults; R-settable so trait sweeps need no
  // rebuild. Rate constants are in day^-1 (converted to the yearly ODE clock in
  // compute_rates); they are CALIBRATION TARGETS pending reconciliation to data.
  double topt_offset = 0.0;   // constitutive thermostability (T_opt/onset) shift, deg C
  double m_rep = 0.4;         // repair-collapse switch slope, deg C^-1
  double t_rep_cut = 45.0;    // temperature above which resynthesis (k_rec) shuts off, deg C
  double k_i = 0.5;           // preventive leak phi_d -> I_r, day^-1
  double k_rec = 0.25;        // restorative resynthesis of I_r, day^-1 (3-5 d rebound)
  double k_mat = 0.02;        // maturation I_r -> I_p, day^-1
  double dTopt_max = 6.0;     // max acclimation rise in T_opt, deg C
  double K_A = 1.0;           // half-saturation of the acclimation response

  // --- Acclimation ODE kinetics (used by compute_rates) ----------------------
  // Single merged thermostability state A. Rates in day^-1.
  double alpha = 0.02;        // gain of the acclimation forcing
  double beta = 0.01;         // decay rate of the acclimation state, day^-1
  double t_accl = 30.0;       // forcing threshold for acclimation, deg C
  double softplus_s = 1.0;    // softplus sharpness for the acclimation forcing

  // --- Running-mean leaf temperature (diagnostic ODE state) -------------------
  // Exponential moving average of the PM operating-point leaf temperature
  // (mean_leaf_temp state): d(mean)/dt = DAYS_PER_YEAR * k_mean_leaf_temp *
  // (T_op - mean). Purely diagnostic -- it does not feed back into any other
  // state. Default gives a ~30-day memory (1/k). R-settable.
  double k_mean_leaf_temp = 0.0333;  // running-mean decay rate, day^-1

  // --- Thermal-cost coefficients (#566) --------------------------------------
  // Ship modest nonzero so the trade-offs bite out of the box; R-settable;
  // CALIBRATION TARGETS, not measured values. Leaf-side coefficients are copied
  // into the leaf in prepare_strategy; the construction/induction coefficients
  // are used directly in net_mass_production_dt.
  double c_acclim_maint  = 0.02;  // -> leaf.c_acclim_maint_ (R_d_, per unit A)
  double c_repair_maint  = 1e-5;  // -> leaf.c_repair_maint_ (R_d_, per day^-1 of k_rec)
  double c_protect_maint = 1e-5;  // -> leaf.c_protect_maint_ (R_d_, per day^-1 of protection)
  double c_repair_flux   = 1e-4;  // -> leaf.c_repair_flux_  (R_d_, per day^-1 resynthesis flux)
  double c_build         = 0.02;  // construction premium on leaf turnover, per deg C thermostability offset
  double c_accl_induct   = 0.01;  // induction cost, biomass per unit positive dA/dt
  double induct_eps      = 1e-4;  // smoothing of the positive-part of dA/dt (C-infinity)

  // Cached slots for the appended states, resolved in refresh_indices().
  int state_idx_acclim_A = -1;
  int state_idx_damage_r = -1;
  int state_idx_damage_p = -1;
  int state_idx_mean_leaf_temp = -1;
  // Cached slot for the appended aux (operating-point leaf temperature).
  int aux_idx_operating_leaf_temp = -1;
};

TF24t_Strategy::ptr make_strategy_ptr(TF24t_Strategy s);

}

#endif
