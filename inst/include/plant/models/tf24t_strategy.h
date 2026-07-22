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
// The Leaf already carries the damage machinery (use_thermal_damage_, the
// Lumry-Eyring functional fraction N, and the tolerance/repair traits); TF24t
// (a) turns that layer on and copies its traits into the leaf in
// prepare_strategy, and (b) adds two slow acclimation states -- acclim_topt
// (raising the photosynthetic optimum T_opt) and acclim_tcrit (raising the
// critical temperature T_crit) -- appended after TF24's states so the inherited
// indices are unchanged. Each state relaxes as
//   dA/dt = alpha * softplus(T_regime - T_accl) - beta * A
// (softplus, not max(0,.), so the forcing is C1). The current state values are
// fed into the leaf each step (leaf.A_opt_ / leaf.A_crit_) before the reused
// TF24 rates run, so the damage factor and T_opt shift see them.
class TF24t_Strategy : public TF24_Strategy {
public:
  typedef std::shared_ptr<TF24t_Strategy> ptr;
  TF24t_Strategy();

  // Compound scientific version: "<TF24 version>.<thermal_revision>", so a TF24
  // science change also invalidates TF24t (the safe direction); bump
  // thermal_revision for changes specific to the thermal layer.
  static constexpr int thermal_revision = 1;

  // TF24's states + the two appended acclimation states (indices unchanged for
  // 0..state_size()-1). Statics resolve on the concrete type in Individual<>.
  static size_t state_size() { return TF24_Strategy::state_size() + 2; }
  static std::vector<std::string> state_names() {
    std::vector<std::string> ret = TF24_Strategy::state_names();
    ret.push_back("acclim_topt");
    ret.push_back("acclim_tcrit");
    return ret;
  }

  // Base refresh_indices() uses the base (static) state_names(); re-run it then
  // register the two appended slots.
  void refresh_indices();

  // Enable the leaf thermal-damage layer and copy the thermal traits into the
  // leaf (mirrors TF24_Strategy::prepare_strategy's use_energy_balance_ / d_
  // handling). Hides (does not override) the non-virtual base method; called on
  // the concrete type by make_strategy_ptr.
  void prepare_strategy();

  // Reuse TF24's rates for the shared states, feed the acclimation states into
  // the leaf, then set the two acclimation-state rates.
  void compute_rates(const TF24_Environment& environment, Internals& vars);

  // Seed the acclimation states at their environmental equilibrium (no birth
  // transient), after seeding the shared TF24 states.
  void set_initial_states(const TF24_Environment& environment, Internals& vars);

  // --- Thermal-damage leaf traits (copied into the leaf in prepare_strategy) --
  // Defaults mirror the Leaf/ATLS defaults; R-settable so trait sweeps need no
  // rebuild.
  double topt_offset = 0.0;   // constitutive T_opt shift, deg C
  double tcrit_0 = 38.0;      // baseline critical temperature, deg C
  double k_d1_0 = 864.0;      // unfolding rate scale, day^-1
  double k_r1_0 = 864.0;      // refold/repair rate scale, day^-1 ("repair" axis)
  double m_switch = 1.0;      // damage switch slope, deg C^-1
  double m_rep = 0.4;         // repair-suppression switch slope, deg C^-1
  double t_rep_cut = 45.0;    // repair shut-off temperature, deg C
  double dTcrit_max = 6.0;    // max acclimation rise in T_crit, deg C
  double dTopt_max = 6.0;     // max acclimation rise in T_opt, deg C
  double K_A = 1.0;           // half-saturation of the acclimation response

  // --- Acclimation ODE kinetics (used by compute_rates) ----------------------
  double alpha_opt = 0.02;    // gain of the T_opt acclimation forcing
  double beta_opt = 0.01;     // decay rate of the T_opt acclimation state, day^-1
  double t_accl_opt = 30.0;   // forcing threshold for T_opt acclimation, deg C
  double alpha_crit = 0.02;   // gain of the T_crit acclimation forcing
  double beta_crit = 0.01;    // decay rate of the T_crit acclimation state, day^-1
  double t_accl_crit = 30.0;  // forcing threshold for T_crit acclimation, deg C
  double softplus_s = 1.0;    // softplus sharpness for the acclimation forcing

  // Cached slots for the appended states, resolved in refresh_indices().
  int state_idx_acclim_topt = -1;
  int state_idx_acclim_tcrit = -1;
};

TF24t_Strategy::ptr make_strategy_ptr(TF24t_Strategy s);

}

#endif
