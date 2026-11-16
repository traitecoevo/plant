// -*-c++-*-
#ifndef PLANT_PLANT_TF24F_MODEL_H_
#define PLANT_PLANT_TF24F_MODEL_H_

#include <plant/models/tf24_model.h>

namespace plant {

// TF24f: TF24 with a plastic hydraulic operating point. Where TF24 re-solves
// the leaf's optimal root-collar potential at every rate evaluation, TF24f
// carries the operating point as one more ODE state that climbs the profit
// gradient at a finite gain, so it lags its optimum; a large gain recovers
// TF24. It inherits TF24_Model and reuses TF24_Environment and TF24_Pars; only
// the leaf-solve hook and the extra state are overridden.
//
// ⚠️ THE STATE IS THE COLLAR'S POSITION, NOT ITS POTENTIAL. A collar potential
// is feasible only on an interval (a, b): from the wet end a, where the soil
// delivers just the leaf's residual flux, to the dry end b, where the stem
// carries all it can. Both ends move with the soil, so a state held in MPa can
// be left outside the interval and has to be clamped back; and near b the
// stem's whole tail is crowded into the last micro-MPa of collar, where the
// marginal profit goes like the log of the distance to the end. So the state
// is z, the logit of the collar's relative position u in the interval:
//
//   psi = a + (b - a) u,      u = 1 / (1 + exp(-z)),
//   dz/dt = k_acclim [ (b - a) dP/dpsi - 2 mu sinh(z) ].
//
// That is gradient ascent on profit in u with a mobility that vanishes at
// either end, du/dt = k_acclim u (1 - u) dP/du, plus a log barrier
// mu ln(u (1 - u)). Every z is a feasible collar, so nothing is clamped. The
// collar moves with the ends as the soil dries. Where profit still rises at an
// end, the state settles where the barrier balances the marginal, a short
// step inside that end, instead of running on. And near the dry end the
// marginal goes like the log of the distance to it, which is linear in z, so
// the rate is not stiff there.
//
// The states shared with TF24 are computed by the inherited compute_rates.
template <typename S = double>
class TF24f_Model : public TF24_Model<S> {
public:
  using value_type = S;
  using base_type = TF24_Model<S>;

  typedef std::shared_ptr<TF24f_Model> ptr;
  TF24f_Model();

  // Scientific version — compound, because TF24f is an *approximation* of
  // TF24 and inherits its equations/parameters. The version is reported as
  // "<TF24 version>.<approximation_revision>" (e.g. "2.1"), so:
  //   * the major component auto-tracks TF24_Model::scientific_version, so a
  //     TF24 scientific change also invalidates TF24f (the safe direction);
  //   * bump `approximation_revision` for changes specific to the tracking
  //     itself.
  // See plant::model_version() / model_id() and src/model_exports.cpp. Bump
  // rules match the other models (output-changing science only, not refactors).
  static constexpr int approximation_revision = 3;

  // TF24's states + the collar's position (appended last so the inherited
  // state indices are unchanged). state_size()/state_names() are static and
  // resolved on the concrete type by Individual<TF24f, ...>.
  static size_t state_size() { return base_type::state_size() + 1; }
  static std::vector<std::string> state_names() {
    std::vector<std::string> ret = base_type::state_names();
    ret.push_back("collar_position");
    return ret;
  }

  // Re-register indices including the new state slot. Base refresh_indices()
  // calls the *static* (base) state_names(), so it would otherwise miss the
  // appended state; we call it then add the extra slot.
  void refresh_indices();

  void compute_rates(const TF24_Environment<S>& environment, Internals<S>& vars);

  // Override the leaf solve: instead of optimising the root-collar psi,
  // evaluate the leaf at the collar the tracked position names and take the
  // exact profit gradient there (left in marginal_ for compute_rates to turn
  // into the state's rate).
  void solve_leaf() override;

  // Initialise the position at the optimum for a newly introduced individual,
  // so tracking starts there (no birth transient, which otherwise lets shaded
  // recruits escape suppression). Runs the base optimiser once via the
  // initializing_ flag below.
  void set_initial_states(const TF24_Environment<S>& environment, Internals<S>& vars);

  // The collar potential at position z between the ends a and b. Written from
  // the nearer end, so the distance to that end keeps its precision.
  static double collar_at(double z, double a, double b) {
    return z > 0.0 ? b - (b - a) / (1.0 + std::exp(z))
                   : a + (b - a) / (1.0 + std::exp(-z));
  }

  // The feasible collar interval at the leaf's current inputs, each end
  // solved to collar precision. False where there is none: the leaf is shut,
  // and is left at its shut state.
  //
  // ⚠️ BOTH ENDS ARE POLISHED, because the collar is measured between them.
  // phylloptim places a bracket's end to its bracket tolerance, and where that
  // root-find stops moves in steps as a trait moves; read at that tolerance,
  // the collar would step with it.
  bool collar_interval(double& a, double& b);

  // Acclimation gain k in dz/dt = k * [(b - a) dP/dpsi - 2 mu sinh(z)]: how
  // fast the position climbs the profit gradient. A large k recovers TF24's
  // instantaneous optimum, to within the barrier's hold at an end.
  S k_acclim = 1.0;

  // The log barrier's weight mu, as a fraction of TF24_cost_scale, the scale
  // of the profit it is added to. It decides how close to an end a state
  // settles when profit still rises there: about mu / (dP/dpsi), in MPa. A
  // state held deeper takes longer to leave, so the weight matters most at a
  // low gain, and less than that gain's own lag: on a stand under dry pulses,
  // offspring production differs between this and a tenth of it by 3e-3 at
  // k_acclim = 10 and 5e-5 at 100.
  static constexpr double barrier_fraction = 1e-5;

  // How far inside either end a newborn is seated, as a fraction of the
  // interval: an optimum on an end has no finite logit.
  static constexpr double seed_margin = 1e-4;

  // Cached slot for the tracked state, resolved in refresh_indices().
  int state_idx_collar_position = -1;

  // Channel between compute_rates() and solve_leaf(): the position going in;
  // the marginal profit at its collar and the interval's width coming back.
  // A species' cohorts share one strategy, so a solve outside compute_rates
  // (establishment_probability's own) reads the position the last rate
  // evaluation left, or the default before any: every position is a feasible
  // collar, so that solve is safe, at another cohort's place in the interval.
  // With no interval the width is zero, and the marginal with it.
  S position_ = 0.0;
  S marginal_ = 0.0;
  S width_ = 0.0;
  // When true, solve_leaf() runs the base optimiser rather than the tracked
  // evaluation, and keeps the interval the optimum sits in; used by
  // set_initial_states to seat the newborn.
  bool initializing_ = false;
  double seed_a_ = 0.0;
  double seed_b_ = 0.0;
  bool seed_open_ = false;
};

template <typename S>
typename TF24f_Model<S>::ptr make_strategy_ptr(TF24f_Model<S> s);

// The base TF24_Model constructor sets collect_all_auxiliary and the name
// and runs the *base* refresh_indices(); we re-run our own refresh_indices()
// here so the appended state slot is registered, and set the reported name.
template <typename S>
TF24f_Model<S>::TF24f_Model() {
  this->name = "TF24f";
  this->refresh_indices();
}

// Build on the base index maps, then register the extra tracked state. The new
// slot is appended after TF24's states, so its index is TF24's state_size().
template <typename S>
void TF24f_Model<S>::refresh_indices() {
  base_type::refresh_indices();
  const int idx = static_cast<int>(base_type::state_size());
  this->state_index["collar_position"] = idx;
  state_idx_collar_position = idx;
  // Again on the concrete model: the base checked its own names against its
  // own count, and the slot appended above is outside both.
  check_state_layout(this->state_index, state_size(), "TF24f");
}

// Reuse TF24's rates for the shared states; the position is fed to
// solve_leaf() (called from the reused net_mass_production_dt) via position_,
// and the marginal profit and the interval's width come back for its rate.
template <typename S>
void TF24f_Model<S>::compute_rates(const TF24_Environment<S>& environment,
                                   Internals<S>& vars) {
  // k_acclim is a user-settable gain; a negative value would silently turn the
  // gradient ascent into descent (away from the optimum) and a non-finite value
  // would poison the state rate, so fail fast on misconfiguration.
  const double k = odelia::util::to_passive(k_acclim);
  if (!util::is_finite(k) || k < 0.0) {
    util::stop("TF24f: k_acclim must be finite and >= 0 (got " +
               util::to_string(k) + ")");
  }
  // The barrier is the only thing that stops a position running on where
  // profit rises at an end, and its weight is a fraction of the cost scale.
  if (!(odelia::util::to_passive(this->pars.TF24_cost_scale) > 0.0)) {
    util::stop("TF24f: TF24_cost_scale must be > 0 (got " +
               util::format_double(
                   odelia::util::to_passive(this->pars.TF24_cost_scale)) +
               "): the barrier that holds the tracked position inside its "
               "interval is a fraction of it");
  }
  position_ = vars.state(state_idx_collar_position);
  base_type::compute_rates(environment, vars);
  // Gradient ascent on profit in the relative position, less the barrier's
  // pull from either end. A shut leaf has no interval to climb, so only the
  // barrier acts, drawing the position slowly to the middle.
  using std::sinh;
  const S barrier = barrier_fraction * this->pars.TF24_cost_scale;
  vars.set_rate(state_idx_collar_position,
                k_acclim * (width_ * marginal_ - 2.0 * barrier * sinh(position_)));
}

template <typename S>
bool TF24f_Model<S>::collar_interval(double& a, double& b) {
  Leaf& leaf = this->leaf;
  if (!leaf.template prepare_collar_solve<
          phylloptim::Leaf::CostCurve::TF24_floor>(a, b)) {
    return false;
  }
  a = leaf.polish_bracket_bound(a, leaf.supply_psi_soil(), 0);
  b = leaf.polish_bracket_bound(b, leaf.supply_psi_soil(), 1);
  if (!(b > a)) {
    odelia::util::stop_domain(
        "TF24f: the collar interval closed when its ends were solved to "
        "collar precision: wet end " + util::format_double(a) +
        " MPa, dry end " + util::format_double(b) + " MPa");
  }
  return true;
}

// Track instead of optimise: evaluate the leaf at the collar the position
// names and supply d(profit)/d(psi) there for the position's rate. The final
// evaluation leaves the leaf outputs at that operating point, which
// compute_rates' aux reads expect. The gradient is exact: forward-mode AD over
// the leaf's algebra, the implicit-function theorem at the ci root-find, and
// analytic spline derivatives for the transport.
//
// The leaf solves in double at every scalar, on the position's value.
//
// A function-try-block, so every exit of the body is covered: phylloptim's
// infeasible_error becomes odelia's DomainError, a rejected step rather than a
// dead run, for the reason TF24_Model::solve_leaf gives. Only that error is
// translated; anything else still propagates.
template <typename S>
void TF24f_Model<S>::solve_leaf() try {
  constexpr auto curve = phylloptim::Leaf::CostCurve::TF24_floor;
  Leaf& leaf = this->leaf;
  if (initializing_) {
    // Birth: the interval first, because the optimiser prepares its own and
    // leaves the leaf at the optimum set_initial_states reads. The
    // curve-typed form, for the reason TF24_Model::solve_leaf gives.
    seed_open_ = collar_interval(seed_a_, seed_b_);
    leaf.template find_root_collar_psi_for<curve>();
    return;
  }
  double a, b;
  if (!collar_interval(a, b)) {
    // Shut: the leaf is at its own operating point, whatever the position.
    marginal_ = 0.0;
    width_ = 0.0;
    return;
  }
  const double collar =
      collar_at(odelia::util::to_passive(position_), a, b);
  marginal_ = leaf.dprofit_droot_collar_psi(collar);
  width_ = b - a;
  // The marginal moved the leaf; put its outputs back at the collar.
  leaf.template profit_at_collar_psi<curve>(collar, a, b);
} catch (const phylloptim::util::infeasible_error& e) {
  odelia::util::stop_domain(std::string("leaf solve infeasible: ") + e.what());
}

// Seat a newborn at its optimum: run the base optimiser once (via the
// initializing_ flag, which makes solve_leaf optimise rather than track) and
// store the optimum's position in the interval it was found in.
template <typename S>
void TF24f_Model<S>::set_initial_states(const TF24_Environment<S>& environment,
                                        Internals<S>& vars) {
  // Every leaf solve here is the optimiser's, TF24's own included: a newborn
  // has no position yet, and the strategy's is another cohort's. Put back on
  // every exit, a rejected step's among them, because the strategy is shared.
  struct restore_flag {
    bool& flag;
    ~restore_flag() { flag = false; }
  };
  initializing_ = true;
  const restore_flag restore{initializing_};
  // Initialise the shared TF24 states first (notably the NSC storage pool) --
  // set_initial_states is non-virtual, so without this call TF24f seedlings
  // would be born with empty reserves and die immediately.
  base_type::set_initial_states(environment, vars);
  this->net_mass_production_dt(environment, vars);
  // With no interval the leaf is shut wherever the position is: the middle.
  double z = 0.0;
  if (seed_open_) {
    // An optimum the solve left on a bracket end it did not polish can sit a
    // bracket tolerance outside the polished interval, so the position is
    // held to [0, 1] before the margin seats it inside.
    const double u = std::min(
        std::max((this->leaf.opt_root_psi_ - seed_a_) / (seed_b_ - seed_a_),
                 0.0),
        1.0);
    const double seated = seed_margin + (1.0 - 2.0 * seed_margin) * u;
    z = std::log(seated / (1.0 - seated));
  }
  vars.set_state(state_idx_collar_position, z);
}

template <typename S>
typename TF24f_Model<S>::ptr make_strategy_ptr(TF24f_Model<S> s) {
  s.prepare_strategy();
  return std::make_shared<TF24f_Model<S> >(s);
}

}

#endif
