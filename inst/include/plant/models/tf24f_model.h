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

  // Another TF24f's values at this scalar, which rebind() uses. TF24's
  // constructor copies the leaf and the table; this adds what TF24f holds.
  template <class S1>
  TF24f_Model(tf24_strategy_from_t tag, const TF24f_Model<S1>& src)
      : base_type(tag, src) {
    assign_tracking_from(src);
  }

  // Another TF24f's values, written into this one. Declared here, so it hides
  // TF24's: a TF24f patch is brought to each width through this.
  template <class S1>
  void assign_from(const TF24f_Model<S1>& src) {
    base_type::assign_from(src);
    assign_tracking_from(src);
  }

  // ⚠️ THIS STRATEGY AT ANOTHER SCALAR, AS A TF24f. TF24's rebind returns a
  // TF24_Model, so an inherited one would rebind a TF24f patch as a TF24 patch:
  // no tracked state, and the leaf re-optimised at every evaluation.
  template <class U>
  TF24f_Model<U> rebind() const {
    return TF24f_Model<U>(tf24_strategy_from_t{}, *this);
  }

  // TF24's census, read off the states TF24f shares with it.
  static const auto& census_metrics() {
    return base_type::template census_metrics_for<TF24f_Model<S>>();
  }
  // Each species' net reproduction ratio and offspring production are stand
  // metrics, refereed for TF24f in test-tf24f-stand-gradient.R. Declared rather
  // than inherited, because the flag says this model's channel was checked.
  static constexpr bool differentiable_fitness = true;

  // TF24's columns with the acclimation gain before the birth rate, which stays
  // last as in TF24. k_acclim reaches the position's rate and nothing else.
  size_t ad_column_count() const { return base_type::ad_column_count() + 1; }
  std::vector<S*> ad_parameters() {
    std::vector<S*> ret = base_type::ad_parameters();
    ret.insert(ret.end() - (this->has_birth_rate_column() ? 1 : 0), &k_acclim);
    return ret;
  }
  std::vector<std::string> ad_parameter_names() {
    std::vector<std::string> ret = base_type::ad_parameter_names();
    ret.insert(ret.end() - (this->has_birth_rate_column() ? 1 : 0), "k_acclim");
    return ret;
  }

  // Every active value this strategy holds: TF24's, and the position's
  // working, which carries the state's and the ends' rows into its rate.
  template <class F>
  void for_each_active(F&& f) {
    base_type::for_each_active(f);
    odelia::ode::visit_active(f, k_acclim, position_, marginal_, width_,
                              seed_position_);
  }

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

  // The collar the leaf is billed at on a differentiated pass, with its rows,
  // and the rows of the position's rate (the definition says how). At birth
  // it is TF24's operating point, and the seed's position takes its rows.
  S collar_on_tape(const Leaf::SupplyDraw<S>& draw,
                   const phylloptim::leaf_pars<S>& in,
                   const phylloptim::SupplyInputs<S>& supply, double slope,
                   Leaf::HeldCoords<S>& held) override;

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

  // The step in z over which a differentiated pass differences the marginal
  // for its slope in the position. Near an end a step in z is a step in the
  // log of the distance to it, which is the coordinate the marginal is smooth
  // in there. Against differences of runs pinned to the recorded steps, a step
  // of 1e-2 leaves errors of 1e-4 in the columns whose traits move the ends,
  // and 1e-3 or 1e-4 leaves 1e-6.
  static constexpr double marginal_step = 1e-4;

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
  // Whether the interval's dry end is the energy balance's flux ceiling rather
  // than the stem's end, as the leaf classified it when the interval was solved:
  // the last interval's, the seed's, and the tracked evaluation's.
  bool cold_ = false;
  bool seed_cold_ = false;
  bool bound_cold_ = false;

  // What solve_leaf() found, for collar_on_tape to give its rows: the
  // interval's ends, whether there was one, and whether the marginal there
  // was a real derivative rather than phylloptim's 0.0 for no flow.
  double bound_a_ = 0.0;
  double bound_b_ = 0.0;
  bool open_ = false;
  bool marginal_feasible_ = false;
  // The seed's position with its rows, where a differentiated pass placed it.
  S seed_position_ = 0.0;
  bool seed_on_tape_ = false;

private:
  template <class S1>
  void assign_tracking_from(const TF24f_Model<S1>& src) {
    k_acclim = S(odelia::util::to_passive(src.k_acclim));
    // The working below is rewritten by the evaluation that reads it, so a
    // rebind carries none of it.
    initializing_ = false;
    position_ = S(0.0);
    marginal_ = S(0.0);
    width_ = S(0.0);
    seed_position_ = S(0.0);
    seed_on_tape_ = false;
    seed_open_ = false;
    open_ = false;
    cold_ = false;
    seed_cold_ = false;
    bound_cold_ = false;
    marginal_feasible_ = false;
    // TF24's assign_from ran TF24's refresh_indices, which misses the slot
    // TF24f appends.
    refresh_indices();
  }

  // The two ends of the interval at the caller's scalar, each a root with its
  // rows in the traits and the soil.
  struct ends_on_tape {
    S a;
    S b;
  };
  ends_on_tape interval_on_tape(double a, double b, bool cold,
                                const phylloptim::leaf_pars<S>& in,
                                const phylloptim::SupplyInputs<S>& supply);
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
  cold_ = leaf.dry_bound_kind() == phylloptim::Leaf::OperatingPointKind::BoundaryCold;
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
// The leaf solves in double at every scalar, on the position's value; a
// differentiated pass gives the collar and the marginal their rows in
// collar_on_tape, from what is recorded here.
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
    seed_cold_ = cold_;
    leaf.template find_root_collar_psi_for<curve>();
    return;
  }
  double a, b;
  open_ = collar_interval(a, b);
  bound_cold_ = cold_;
  if (!open_) {
    // Shut: the leaf is at its own operating point, whatever the position.
    marginal_ = 0.0;
    width_ = 0.0;
    return;
  }
  const double collar =
      collar_at(odelia::util::to_passive(position_), a, b);
  bool feasible = false;
  marginal_ = leaf.dprofit_droot_collar_psi(collar, &feasible);
  marginal_feasible_ = feasible;
  bound_a_ = a;
  bound_b_ = b;
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
  seed_on_tape_ = false;
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
  if constexpr (std::is_same_v<S, double>) {
    vars.set_state(state_idx_collar_position, z);
  } else {
    // The same value, with the rows of the optimum and of both ends where a
    // differentiated pass recorded them: at an interior optimum the collar's
    // are the implicit function theorem's. Where the record refused, or the
    // position was held to [0, 1], the value alone.
    vars.set_state(state_idx_collar_position,
                   seed_on_tape_
                       ? S(S(z) + (seed_position_ -
                                   S(odelia::util::to_passive(seed_position_))))
                       : S(z));
  }
}

// The interval's two ends with their rows: the wet end is where uptake meets
// the residual flux, the dry end where the soil delivers what the stem can
// still carry. Each takes a supply draw at its own collar.
//
// ⚠️ WITH THE ENERGY BALANCE ON, the dry end can be the collar whose uptake is
// the flux ceiling instead, a condition phylloptim's surface has no rows for
// (it refuses TF24's own pin there). Refused here rather than differentiated
// as the stem's end; `cold` is the leaf's own classification of that end.
template <typename S>
typename TF24f_Model<S>::ends_on_tape TF24f_Model<S>::interval_on_tape(
    double a, double b, bool cold, const phylloptim::leaf_pars<S>& in,
    const phylloptim::SupplyInputs<S>& supply) {
  constexpr auto curve = phylloptim::Leaf::CostCurve::TF24_floor;
  Leaf& leaf = this->leaf;
  const Leaf::SupplyDraw<S> at_a = leaf.template supply_draw_at<S>(S(a), supply);
  const Leaf::SupplyDraw<S> at_b = leaf.template supply_draw_at<S>(S(b), supply);
  if (cold) {
    this->refuse("TF24f gradient: the dry end of the collar interval, " +
                 util::format_double(b) + " MPa, is where uptake reaches the "
                 "energy balance's flux ceiling, and that bound has no "
                 "derivative rows");
  }
  return {leaf.template wet_bound_at<S>(a, at_a, in),
          leaf.template dry_end_collar_at<curve, S>(b, at_b, in)};
}

// The collar solve_leaf() evaluated the leaf at, with its rows, and the rows
// of the marginal M and of the interval's width, so that the position's rate
// k [(b - a) M - 2 mu sinh z] is differentiated as compute_rates forms it.
//
// The collar is a + (b - a) u(z): its rows are the position's and the two
// ends'. At a shut leaf it is the leaf's own collar, as TF24 places it.
//
// M's rows in the traits, the soil and the light come from the marginal at the
// collar held; its row in the collar is a supplied slope, dM/dz over dpsi/dz,
// with dM/dz a central difference of the double marginal in z at fixed ends.
// The tape's own dM/dcollar is wrong by a factor of about eight (phylloptim's
// profit_curvature_at_collar says why). Where no flow reaches the collar the
// marginal is phylloptim's 0.0, a constant.
template <typename S>
S TF24f_Model<S>::collar_on_tape(const Leaf::SupplyDraw<S>& draw,
                                 const phylloptim::leaf_pars<S>& in,
                                 const phylloptim::SupplyInputs<S>& supply,
                                 double slope, Leaf::HeldCoords<S>& held) {
  using odelia::util::to_passive;
  using std::exp;
  using std::log;
  constexpr auto curve = phylloptim::Leaf::CostCurve::TF24_floor;
  Leaf& leaf = this->leaf;
  if (initializing_) {
    const S collar = base_type::collar_on_tape(draw, in, supply, slope, held);
    if (seed_open_) {
      const double u =
          (to_passive(collar) - seed_a_) / (seed_b_ - seed_a_);
      if (u >= 0.0 && u <= 1.0) {
        const ends_on_tape ends = interval_on_tape(seed_a_, seed_b_, seed_cold_, in, supply);
        const S seated = seed_margin + (1.0 - 2.0 * seed_margin) *
                                           ((collar - ends.a) / (ends.b - ends.a));
        seed_position_ = log(seated / (1.0 - seated));
        seed_on_tape_ = true;
      }
    }
    return collar;
  }
  const double marginal = to_passive(marginal_);
  if (!open_) {
    marginal_ = S(marginal);
    width_ = S(0.0);
    return base_type::collar_on_tape(draw, in, supply, slope, held);
  }
  const double at = draw.at;
  const double z = to_passive(position_);
  const ends_on_tape ends = interval_on_tape(bound_a_, bound_b_, bound_cold_, in, supply);
  const S width = ends.b - ends.a;
  width_ = width;
  // From the nearer end, as collar_at writes it; the value is the collar the
  // leaf was evaluated at.
  const S live = z > 0.0 ? S(ends.b - width / (1.0 + exp(position_)))
                         : S(ends.a + width / (1.0 + exp(-position_)));
  const S collar = S(at) + (live - S(to_passive(live)));
  if (!marginal_feasible_) {
    marginal_ = S(marginal);
    return collar;
  }
  // Placed once at the collar and read by the marginal here and by the billing
  // after, which adds only the collar's channel.
  leaf.template place_held_coords<S>(collar, draw, in, held);
  const S M_held =
      leaf.template dprofit_dcollar_at<curve, S>(S(at), draw, in, held);
  // dM/dz at fixed ends. The marginal moves the leaf, so its outputs are put
  // back at the collar before anything reads them.
  bool up_ok = false, down_ok = false;
  const double up = leaf.dprofit_droot_collar_psi(
      collar_at(z + marginal_step, bound_a_, bound_b_), &up_ok);
  const double down = leaf.dprofit_droot_collar_psi(
      collar_at(z - marginal_step, bound_a_, bound_b_), &down_ok);
  leaf.template profit_at_collar_psi<curve>(at, bound_a_, bound_b_);
  // One side is enough: next to an end with no flow, the far neighbour has no
  // marginal to read.
  double dM_dz = std::numeric_limits<double>::quiet_NaN();
  if (up_ok && down_ok) {
    dM_dz = (up - down) / (2.0 * marginal_step);
  } else if (up_ok) {
    dM_dz = (up - marginal) / marginal_step;
  } else if (down_ok) {
    dM_dz = (marginal - down) / marginal_step;
  }
  const double dpsi_dz =
      (bound_b_ - bound_a_) / (2.0 + std::exp(z) + std::exp(-z));
  const double dM_dcollar = dM_dz / dpsi_dz;
  if (!std::isfinite(dM_dcollar)) {
    this->refuse("TF24f gradient: the marginal profit's slope at the tracked "
                 "collar " + util::format_double(at) + " MPa is " +
                 util::format_double(dM_dcollar) +
                 ", so the position's rate has no row in the collar");
    marginal_ = S(marginal);
    return collar;
  }
  const odelia::input_and_derivative<S> rows[]{{M_held, 1.0},
                                               {collar, dM_dcollar}};
  const odelia::supply_report report =
      odelia::supply_derivatives<S>(marginal, rows, marginal_);
  if (!report.whole) {
    this->refuse("TF24f gradient: the position's rate cannot be recorded at "
                 "the collar " + util::format_double(at) + " MPa: " + report.why);
  }
  return collar;
}

template <typename S>
typename TF24f_Model<S>::ptr make_strategy_ptr(TF24f_Model<S> s) {
  s.prepare_strategy();
  return std::make_shared<TF24f_Model<S> >(s);
}

}

#endif
