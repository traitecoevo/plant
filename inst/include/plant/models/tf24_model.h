// -*-c++-*-
#ifndef PLANT_PLANT_TF24_MODEL_H_
#define PLANT_PLANT_TF24_MODEL_H_

#include <plant/census.h>
#include <plant/model.h>
#include <odelia/ode_util.hpp>
#include <plant/models/tf24_environment.h>
#include <plant/census_gradient.h>
#include <plant/qag.h>
#include <plant/leaf_model.h>
#include <plant/canopy_shape.h>
#include <plant/stem_hydraulics.h>
#include <odelia/implicit_node.hpp>
#include <plant/with_slope.h>
#include <array>
// cstdio/cstdlib for the environment-gated curvature comparison in
// record_leaf_outputs, and nothing else here.
#include <cstdio>
#include <cstdlib>
#include <string_view>
#include <utility>
#include <type_traits>
#include <plant/models/tf24_pars.h>

namespace plant {

// Templated on the scalar S the state, the traits and everything derived from
// them carry; double is production. The embedded Leaf, the Control tolerances
// and the extrinsic drivers stay double.
// A cursor over the operating points one rate evaluation solves for. The points
// themselves are the SOLVER's record, beside the times, the sizes and the states,
// because they are part of what the run did; this walks the list the solver hands
// over for the evaluation about to run.
//
// At most one handle is open, and which one is the constness of what was handed
// over. Nothing here remembers a mode, so nothing can outlive the list it
// described.
// What a pass records about one leaf solve: the collar it placed, and the ARM it
// placed it on.
//
// Two values rather than a class, because that is all a recorded decision is --
// the collar MOVES and is re-derived from, the arm SELECTS and is replayed. What
// else a pass might record is rederivable from the collar alone, which is what
// phylloptim's evaluate_root_collar_psi does.
struct leaf_solved_point {
  double collar = 0.0;
  Leaf::OperatingPointKind kind = Leaf::OperatingPointKind::Unsolved;
};

class leaf_solved_points {
public:
  using points = std::vector<leaf_solved_point>;

  void store_into(points& into) {
    solved = 0;
    storing = &into;
    loading = nullptr;
    storing->clear();
  }
  void load_from(const points& from) {
    solved = 0;
    storing = nullptr;
    loading = &from;
  }
  // One evaluation is one list, so this ends where the evaluation does. A solve
  // outside a rate evaluation then stores nothing and loads nothing.
  void end() {
    solved = 0;
    storing = nullptr;
    loading = nullptr;
  }

  // The operating point the run found for the next solve of this evaluation, or one
  // holding nothing -- which is what a pass filling the record gets, and is the
  // caller's signal to solve for it.
  //
  // ⚠️ RUNNING OFF THE END IS A FAULT, not a fall back to solving. The run made one
  // entry per solve at this evaluation, so a pass that asks for more is a pass
  // whose model disagrees with the record it was handed, and searching from here
  // would tape a discretisation the run never took with every number finite.
  leaf_solved_point load() {
    if (loading == nullptr) {
      return leaf_solved_point();
    }
    if (solved >= loading->size()) {
      util::stop("leaf_solved_points: this rate evaluation asked for operating "
                 "point " + util::to_string(solved + 1) + " where the run it is "
                 "replaying solved for " + util::to_string(loading->size()));
    }
    ++placed;
    return (*loading)[solved++];
  }

  // How many points this record placed. Counted because a record that engages and
  // one that quietly does not are the same green suite: every number a placement
  // produces is the number a search produces, so nothing else can tell them apart.
  size_t placements() const { return placed; }

  void store(const leaf_solved_point& point) {
    if (storing != nullptr) {
      storing->push_back(point);
    }
  }

private:
  // The list this evaluation stores into, and the one it loads from. At most one is
  // ever open.
  points* storing = nullptr;
  const points* loading = nullptr;
  size_t solved = 0;
  size_t placed = 0;
};

// Selects TF24_Model's constructor from another strategy, which rebind() uses.
struct tf24_strategy_from_t {};

template <typename S = double>
class TF24_Model: public Model<TF24_Environment<S>> {
public:
  using value_type = S;

  typedef std::shared_ptr<TF24_Model> ptr;
  TF24_Model();

  // Another strategy's values at this scalar, as assign_from() writes them.
  //
  // ⚠️ THE LEAF IS COPY-CONSTRUCTED, NOT DEFAULTED AND THEN ASSIGNED. A default
  // phylloptim::Leaf computes both vulnerability tables from incomplete gamma
  // functions, about 60 times the cost of copying them, and assign_from() would
  // overwrite them at once. The sweep rebinds the patch on the order of once
  // per recorded step, so each build it skips is one per step.
  template <class S1>
  TF24_Model(tf24_strategy_from_t, const TF24_Model<S1>& src)
      : leaf(src.leaf) {
    assign_from(src);
  }

  // Scientific version. Bump only when equations or default parameters change
  // the simulation output for identical inputs, never for a refactor or a
  // performance, interface or serialisation change: a bump invalidates logpile's
  // cache for this model (model_version(), model_id()). TF24f, TF24t and TF24w
  // carry compound versions that follow this one. NEWS.md records each bump and
  // what it moved, with two gaps: v11, which added R_d_25, moved the root curve
  // to (root_P50, root_c) and vulnerability_curve_ncontrol from 100 to 400; and
  // the half of v9 that took the light field's and the census's trapezium widths
  // from the coordinate the density is carried in.
  static constexpr int scientific_version = 15;

  // calculate the amount of water transpired relativised by leaf area index.

  S evapotranspiration_dt(const S& area_leaf_, int soil_layer);


  // Overrides ----------------------------------------------

  // update this when the length of state_names changes
  static size_t state_size () { return 9; }
  // update this when the length of aux_names changes
  size_t aux_size () { return aux_names().size(); }

  static std::vector<std::string> state_names() {
    return  std::vector<std::string>({
      "height",
      "mortality",
      "fecundity",
      "area_heartwood",
      "mass_heartwood",
      "storage",
      "log_area_leaf_departure",
      "log_area_sapwood_departure",
      "leaf_marginal_return_tracked"
      });
  }

  // TF24 carries reserves that feed back into growth, so its growth rate is not
  // a function of size alone and the height coordinate's compression term is a
  // different derivative rather than a worse approximation of the right one
  // (#590). The two coordinates' offspring production differ by several times
  // here, where FF16's and K93's agree to 5e-5 and 5.6e-7, so this is a
  // correctness default and not a performance one. TF24f, TF24t and TF24w
  // inherit it.
  static constexpr bool density_in_birth_date_default = true;

  // Each species' net reproduction ratio and offspring production are stand
  // metrics (SCM::stand_metrics), refereed in test-tf24-fitness-gradient.R.
  static constexpr bool differentiable_fitness = true;

  // The metrics a census of this model sums, in the order it reports them. Beside
  // the state and aux names because all three say what this model calls its own
  // quantities -- and every slot is read by its cached index, which is how every
  // other reader of these slots reaches them.
  //
  // A table rather than a list built per call, for the reason
  // ad_parameter_fields is one: a census is read from several entry points and
  // each was allocating its own copy of this.
  //
  // ⚠️ THE LEAF AREA IS THE ONE THE PLANT HAS, read from the departure state, not
  // the one its height implies, and the sapwood is likewise at the Huber value it
  // carries (#516). At rest both departures are exactly zero and exp(0) is
  // exactly one, so these are the height-only expressions bit for bit.
  static const auto& census_metrics() {
    using strategy = TF24_Model<S>;
    static constexpr std::array<census_metric<strategy>, 3> metrics{{
      {"leaf_area",
       [](const strategy& p, const Internals<S>& vars) -> S {
         return p.area_leaf_carried(vars);
       }},
      {"mass_above_ground",
       [](const strategy& p, const Internals<S>& vars) -> S {
         const S height = vars.state(HEIGHT_INDEX);
         const S area_leaf = p.area_leaf_carried(vars);
         return p.mass_above_ground(
             p.mass_leaf(area_leaf),
             p.mass_bark(p.area_bark(area_leaf), height),
             p.mass_sapwood(
                 p.area_sapwood(
                     area_leaf,
                     vars.state(p.state_idx_log_area_sapwood_departure)),
                 height),
             vars.state(p.state_idx_mass_heartwood));
       }},
      {"area_stem",
       [](const strategy& p, const Internals<S>& vars) -> S {
         const S area_leaf = p.area_leaf_carried(vars);
         return p.area_stem(
             p.area_bark(area_leaf),
             p.area_sapwood(area_leaf,
                            vars.state(p.state_idx_log_area_sapwood_departure)),
             vars.state(p.state_idx_area_heartwood));
       }},
    }};
    return metrics;
  }

  // The leaf area a plant carries: A*(h) exp(phi). update_dependent_aux keeps the
  // same number in the competition_effect aux; this reads the states directly,
  // for a caller holding states whose aux may not have been refreshed.
  S area_leaf_carried(const Internals<S>& vars) const {
    using std::exp;
    return area_leaf(vars.state(HEIGHT_INDEX)) *
      exp(vars.state(state_idx_log_area_leaf_departure));
  }

  // The storage pool. Its rate holds the exact flow inside [0, S_max], so a
  // negative value is a step that overshot rather than a state the model has,
  // and the solver rejects it and retries smaller.
  static std::vector<std::string> non_negative_states() { return {"storage"}; }

  // Virtual so that a subclass can append its own readings (TF24w's per-column
  // fluxes, #681): refresh_indices builds aux_index from this, and aux_size()
  // counts it, so one override reaches both.
  virtual std::vector<std::string> aux_names() {
    std::vector<std::string> ret({
      "competition_effect",
      "height_inverse",
      "net_mass_production_dt",
      "root_mass",
      "opt_psi_leaf",
      "opt_root_psi",
      "transpiration",
      "E_up_",
      "profit",
      // The part of the cost the objective deducted that is a PRICE rather than
      // a realised carbon loss (umol C m^-2 s^-1): `TF24_floor_lambda_o * E` on
      // the seated curve, and exactly zero at the default price. Reported beside
      // `profit` because without it the two cannot be separated after the fact
      // from plant's own output, and the distinction is the whole point of the
      // curve: the carbon the plant actually kept is `profit + shadow_cost`,
      // which is what net_mass_production_dt grows on.
      "shadow_cost",
      "stom_cond_CO2",
      // Net CO2 assimilation at the optimal operating point, per unit leaf
      // area (umol CO2 m^-2 s^-1). Net, not gross: Leaf::assim_colimited()
      // subtracts day respiration R_d_, so gross = assimilation + R_d_, with
      // R_d_ = R_d_25 * q10^((leaf_temp - 25) / 10) inside phylloptim.
      "assimilation",
      // Leaf temperature at the optimal operating point (deg C) -- an OUTPUT,
      // not the `leaf_temp` driver, and the distinction is the point of
      // reporting it. With pars.use_energy_balance off this equals the driver,
      // so the column is flat at TF24's defaults; with it on the leaf solves
      // its own temperature from its transpiration per operating point. It is
      // the one quantity the Penman-Monteith path exists to produce, so it is
      // reported rather than left to be inferred from an assimilation that
      // cannot be explained without it.
      //
      // Under the deep-crown shading model this is the leaf-area-weighted crown
      // mean, integrated alongside the other leaf outputs; it is NOT the
      // temperature of any single leaf, and a canopy with a hot top and a cool
      // base reports the mean of the two. The depth profile itself is not an
      // aux (a fixed-width slot cannot carry a per-quadrature-node vector).
      "Tleaf",
      // The marginal leaf's carbon balance, which the replacement gate reads.
      // Dimensionless and zero at break-even, so its sign says whether the last
      // leaf is paying for itself.
      "leaf_marginal_return",
      // The marginal return of SAPWOOD area, in growth rather than production
      // (#516). Dimensionless, and zero at the Huber value that maximises the
      // plant's growth rate, so its sign says which way the stem should move.
      // Zero unless sapwood acclimation is on (a_sw > 0), because computing it
      // costs two further leaf solves.
      "sapwood_marginal_return"
    });
    // add the associated computation to compute_rates and compute there
    if (this->collect_all_auxiliary) {
      ret.push_back("area_sapwood");
    }
    return ret;
  }

  using ad_parameter = typename TF24_Pars<S>::ad_parameter;

  // A constant birth rate, at the working scalar so that it is a gradient
  // column; a variable one is a driver read as a double and has no column.
  S birth_rate_ = 1.0;
  bool has_birth_rate_column() const { return !this->is_variable_birth_rate; }
  S birth_rate_at(double t) const {
    return has_birth_rate_column()
             ? birth_rate_
             : S(this->extrinsic_drivers.evaluate("birth_rate", t));
  }

  // How many gradient columns one of these carries, for a caller sizing a batch.
  size_t ad_column_count() const {
    return TF24_Pars<S>::column_count + (has_birth_rate_column() ? 1 : 0);
  }

  // Addresses of the parameters a gradient can be taken with respect to, in the
  // order ad_parameter_names() gives: the table entries whose role has a column,
  // then the birth rate where it is one. Both allocate, so take them once per
  // gradient evaluation and hold them for the run rather than per block; the
  // strategy is shared and the fields do not move. Index against .size().
  std::vector<S*> ad_parameters() {
    std::vector<S*> ret;
    ret.reserve(ad_column_count());
    for (const ad_parameter& p : TF24_Pars<S>::ad_parameter_columns) {
      ret.push_back(&(pars.*p.at));
    }
    if (has_birth_rate_column()) {
      ret.push_back(&birth_rate_);
    }
    return ret;
  }

  // The name each gradient column carries, one per ad_parameters() entry and in
  // that order. Which parameters are left out, and why, is `undifferentiable`.
  std::vector<std::string> ad_parameter_names() {
    std::vector<std::string> ret;
    ret.reserve(ad_column_count());
    for (const ad_parameter& p : TF24_Pars<S>::ad_parameter_columns) {
      ret.push_back(p.name);
    }
    if (has_birth_rate_column()) {
      ret.push_back("birth_rate");
    }
    return ret;
  }

  // The parameters no gradient exists for, each with the sentence that says
  // why, so a caller asking for one is refused in those words rather than told
  // the name is unknown.
  static std::vector<std::pair<std::string, std::string>> undifferentiable_reasons() {
    std::vector<std::pair<std::string, std::string>> ret;
    for (const auto& n : TF24_Pars<S>::undifferentiable) {
      ret.emplace_back(n.name, n.why);
    }
    return ret;
  }


  void refresh_indices();


  // TF24 Methods  ----------------------------------------------

  // Leaf area from height
  S area_leaf(const S& height) const;

  // Leaf mass from leaf area
  S mass_leaf(const S& area_leaf) const;

  // Sapwood area per unit leaf area -- the Huber value, `theta` displaced by
  // the sapwood departure state. ALWAYS form the ratio this way rather than as
  // area_sapwood / area_leaf: the two are algebraically equal but not bitwise,
  // because a multiply followed by a divide does not recover its operand, and
  // that difference is enough to lose the exactness of the resting model.
  S sapwood_per_leaf_area(const S& sapwood_departure) const {
    using std::exp;
    return pars.theta * exp(sapwood_departure);
  }

  // Area and mass of sapwood. The one-argument form is the resting
  // allometry (used where a plant is on its trajectory by definition, e.g. the
  // germination inversion); the two-argument form reads the departure state.
  S area_sapwood(const S& area_leaf) const;
  S area_sapwood(const S& area_leaf, const S& sapwood_departure) const;
  S mass_sapwood(const S& area_sapwood, const S& height) const;

  // Area and mass of bark
  S area_bark(const S& area_leaf) const;
  S mass_bark (const S& area_bark, const S& height) const;

  S area_stem(const S& area_bark, const S& area_sapwood,
                            const S& area_heartwood) const;
  S diameter_stem(const S& area_stem) const;

  // Mass of (fine) roots
  S mass_root(const S& area_leaf) const;

  // Mass of coarse (structural) roots, as a fraction of sapwood mass
  // (see TF24_Pars::a_cr1). Zero unless a_cr1 is set.
  S mass_coarse_root(const S& mass_sapwood) const;

  // Total mass
  S mass_live(const S& mass_leaf, const S& mass_bark,
              const S& mass_sapwood, const S& mass_root,
              const S& mass_coarse_root) const;

  S mass_total(const S& mass_leaf, const S& mass_bark, const S& mass_sapwood,
               const S& mass_heartwood, const S& mass_root,
               const S& mass_coarse_root) const;

  // Above-ground mass = leaf + all stem components (bark + sapwood +
  // heartwood); excludes roots.
  S mass_above_ground(const S& mass_leaf, const S& mass_bark,
                      const S& mass_sapwood, const S& mass_heartwood) const;

  void compute_rates(const TF24_Environment<S>& environment,
                Internals<S>& vars);

  void update_dependent_aux(const int index, Internals<S>& vars);

  // * Mass production
  // Total maintenance respiration
  S respiration(const S& mass_leaf, const S& mass_sapwood,
                const S& mass_bark, const S& mass_root,
                const S& mass_coarse_root) const;

  S respiration_leaf(const S& mass) const;
  S respiration_bark(const S& mass) const;
  S respiration_sapwood(const S& mass) const;
  S respiration_root(const S& mass) const;
  S respiration_coarse_root(const S& mass) const;

  // Total turnover
  // Turnover the plant actually PAYS for. `replacement_leaf` and
  // `replacement_sapwood` are the fractions of the leaf and sapwood flux it
  // chooses to rebuild; the rest is not charged, and the matching tissue is lost
  // through the departure states instead.
  //
  // ⚠️ ONLY A POOL WHOSE SIZE IS A STATE MAY BE GATED. Bark and fine-root
  // turnover stay fully charged because neither has a state to lose, so
  // withholding their cost would save carbon that no tissue paid for -- a free
  // lunch, and a silent one. Coarse roots follow the sapwood's mass, so they do
  // shrink with it, but they stay charged in full: over-charging is the safe
  // side, and the pool is empty at the default a_cr1 = 0.
  S turnover(const S& replacement_leaf, const S& replacement_sapwood,
             const S& mass_leaf, const S& mass_bark,
             const S& mass_sapwood, const S& mass_root,
             const S& mass_coarse_root) const;
  S turnover_leaf(const S& mass) const;
  S turnover_bark(const S& mass) const;
  S turnover_sapwood(const S& mass) const;
  S turnover_root(const S& mass) const;
  S turnover_coarse_root(const S& mass) const;

  // The fraction of leaf and sapwood turnover a plant rebuilds, as a smooth
  // function of the tracked marginal leaf return: ~1 while the marginal leaf
  // pays for itself, falling to 1 - a_pl0 as it stops. a_pl0 = 0 gives exactly 1
  // and so exactly the fixed allometry.
  //
  // Two forms of one logistic, chosen by the sign of the exponent, because a
  // healthy plant's return sits tens of gate widths above the centre: there
  // e^x overflows to inf, which the value survives and its derivative,
  // inf/inf, does not -- and a NaN partial multiplied by a_pl0 = 0 is still NaN.
  S replacement_fraction(const S& marginal_return) const {
    using std::exp;
    const S x = (marginal_return - pars.a_pl1) / plasticity_gate_width;
    if (x > 0.0) {
      const S e = exp(-x);
      return 1.0 - pars.a_pl0 * e / (1.0 + e);
    }
    return 1.0 - pars.a_pl0 / (1.0 + exp(x));
  }

  // How much of the withheld replacement actually happens, tapering to nothing
  // at the canopy floor. Same functional form as the rebuild share, mirrored:
  // rebuilding vanishes at the ceiling (phi = 0) and shedding at the floor
  // (phi = ln a_pl3), so phi is confined to [ln a_pl3, 0] by the SHAPE of the
  // flow rather than by any clamp -- and both bounds are therefore differentiable
  // and cost the solver nothing.
  S shedding_gate(const S& leaf_departure) const {
    using std::exp;
    using std::log;
    const S floor_departure = log(pars.a_pl3);
    return 1.0 - exp((floor_departure - leaf_departure) / pars.a_pl2);
  }

  // Sapwood replaces a scaled share of what its conversion to heartwood
  // removes: the same carbon gate as leaves, times exp(-psi / a_pl2), so a stem
  // carrying more conducting area than its canopy needs declines to renew the
  // excess even when carbon is ample.
  //
  // This is what gives the Huber value a restoring force. Without it psi has
  // only terms that push it around and none that return it, and measured, it
  // drifted to leave plants hydraulically UNDER-built (#516).
  //
  // ⚠️ psi >= 0 IS INVARIANT, and that is the reason this form needs no bound.
  // At psi = 0 the factor is 1, so the two replacement fractions coincide and
  //     dpsi/dt = (1 - replacement) * (k_l - k_s) >= 0,
  // since k_l > k_s. The flow at the boundary points inward, so the "too little
  // sapwood" branch is unreachable, the factor stays within (0, 1], and no
  // over-replacement (a fraction above 1) can arise to cost unbounded carbon.
  // If k_s were ever raised above k_l that argument reverses and this needs
  // revisiting.
  S sapwood_replacement_fraction(const S& replacement_leaf,
                                 const S& sapwood_departure) const {
    using std::exp;
    return replacement_leaf * exp(-sapwood_departure / pars.a_pl2);
  }

  // Carbon cost of a unit of leaf area bought at UNCHANGED height: the leaf,
  // the fine root, and the sapwood and bark that supply it, at the ratio the
  // plant prefers rather than the one it currently has. Buying the whole
  // package is what stops rebuilding a canopy from diluting the Huber value.
  S mass_live_per_area_leaf_at_height(const S& height) const {
    return pars.lma + pars.a_r1 +
      pars.theta * height * eta_c * pars.rho * (1.0 + pars.a_b1);
  }

  // Net production
  S net_mass_production_dt_A(const S& assimilation, const S& respiration,
                             const S& turnover) const;

  // Takes the three #516 states beside the sizes: the sapwood departure moves
  // the conductance handed to the leaf and the sapwood carried, and the leaf
  // departure and the tracked marginal return set what turnover is charged. A
  // plant on its allometry with no history passes zeros.
  virtual S net_mass_production_dt(const TF24_Environment<S>& environment,
                                const S& height, const S& area_leaf_,
                                const S& height_inverse,
                                const S& sapwood_departure,
                                const S& leaf_departure,
                                const S& tracked_marginal);

  // The soil water potentials the leaf's supply solve reads, one per layer.
  // TF24 reads the patch's shared column; TF24w (#681) overrides this to hand
  // over the cohort's own root-zone column, computed from its own water states
  // before the base compute_rates runs. Virtual for the reason solve_leaf is:
  // net_mass_production_dt is reused unchanged by the subclass.
  virtual const std::vector<S>& soil_water_potential_for_leaf(
      const TF24_Environment<S>& environment) {
    return environment.get_soil_water_potential_state();
  }

  // Resolve the leaf operating point on the already-set-up `leaf` (i.e. after
  // leaf.set_physiology(...)). Base TF24 optimises the root-collar psi via
  // golden-section search; the TF24f variant overrides this to make the optimum
  // chase a tracked ODE state (#525). Called per crown light point from
  // net_mass_production_dt, so it must be virtual to dispatch to the override
  // when net_mass_production_dt is reused unchanged by the subclass.
  virtual void solve_leaf();

  // Read how the leaf's outputs respond to what it was given, and record each
  // output that re-enters the active chain carrying that response. The leaf is
  // already driven and solved; nothing here re-supplies it.
  void record_leaf_outputs(const S& PPFD_int, const std::vector<S>& psi_soil,
                           const S& conductance_max, double NIR_per_PPFD_int);
  // Two steps of record_leaf_outputs, and the marginal its curvature guard
  // reports (tf24_leaf.h says what each does).
  phylloptim::leaf_pars<S> leaf_pars_on_tape(const S& PPFD_int,
                                             const S& conductance_max,
                                             double NIR_per_PPFD_int) const;
  double marginal_at_collar() const;
  void refuse_unusable_curvature(double slope);

  // The steps of net_mass_production_dt, in the order it takes them, and
  // seat_leaf, which two of them use (tf24_leaf.h says what each computes).
  S update_leaf_specific_conductance_max(const S& height,
                                         const S& sapwood_departure);
  void build_root_network(const TF24_Environment<S>& environment,
                          const S& height);
  void seat_leaf(const TF24_Environment<S>& environment, double PPFD,
                 const std::vector<double>& psi, double conductance);
  S solve_leaf_over_crown(const TF24_Environment<S>& environment,
                          const S& height, const std::vector<S>& psi_soil,
                          const S& leaf_specific_conductance_max,
                          double PPFD_above, double NIR_per_PPFD_int);
  void measure_dprofit_dkmax(const TF24_Environment<S>& environment,
                             const S& PPFD_int_used,
                             const std::vector<S>& psi_soil,
                             const S& leaf_specific_conductance_max);
  void update_marginal_leaf_return(const S& assimilation_, const S& mass_leaf_,
                                   const S& mass_bark_, const S& mass_root_);

  // Two steps of compute_rates (tf24_rates.h says what each computes).
  void write_leaf_readings(Internals<S>& vars, const S& area_leaf_);
  S sapwood_marginal_return(const S& area_leaf_, const S& area_leaf_star,
                            const S& sapwood_departure, const S& mass_sapwood_,
                            const S& withheld_sapwood, const S& growth_flux,
                            const S& net_mass_production_dt_,
                            const S& darea_leaf_dmass_live_);
  // Model-agnostic entry point used by Individual<TF24> (#266): reads the
  // height state and the cached aux slots itself, so the generic Individual
  // does not need to know TF24's state/aux layout.
  S net_mass_production_dt(const TF24_Environment<S>& environment,
                                const Internals<S>& vars) {
    return net_mass_production_dt(
        environment, vars.state(HEIGHT_INDEX),
        vars.aux(aux_idx_competition_effect),
        vars.aux(aux_idx_height_inverse),
        vars.state(state_idx_log_area_sapwood_departure),
        vars.state(state_idx_log_area_leaf_departure),
        vars.state(state_idx_leaf_marginal_tracked));
  }

  // Fraction of production allocated to reproduction
  virtual S fraction_allocation_reproduction(const S& height) const;

  // The height above which this strategy allocates at least 1 - epsilon of its
  // maximum, a_f1, to reproduction: where a continuous patch's absorbing top
  // layer sits (#519). Growth's share of surplus there is 1 - (1 - epsilon) a_f1,
  // within epsilon a_f1 of its floor -- and epsilon itself only where a_f1 is one.
  // The inverse of fraction_allocation_reproduction()'s logistic, and NaN where
  // the logistic does not rise with height (a_f2 <= 0) or hmat is not positive.
  double asymptotic_height(double epsilon) const {
    using odelia::util::to_passive;
    const double hmat = to_passive(pars.hmat);
    const double a_f2 = to_passive(pars.a_f2);
    if (!(hmat > 0.0) || !(a_f2 > 0.0)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return hmat * (1.0 + std::log((1.0 - epsilon) / epsilon) / a_f2);
  }
  S fraction_allocation_growth(const S& height) const;
  // Rate of offspring production
  S fecundity_dt(const S& net_mass_production_dt,
                 const S& fraction_allocation_reproduction) const;

  // Leaf area bought per unit live mass. Takes the sapwood departure because
  // the stem that comes with a new leaf is priced at the Huber value the plant
  // ACTUALLY carries, not at pars.theta. Without that a plant holding extra
  // conducting area gets its hydraulic benefit free, which is most of why net
  // production appeared to peak at 3.3x the pipe-model ratio.
  S darea_leaf_dmass_live(const S& area_leaf,
                          const S& sapwood_departure) const;

  // change in height per change in leaf area
  S dheight_darea_leaf(const S& area_leaf) const;
  // Mass of leaf needed for new unit area leaf, d m_l / d a_l
  S dmass_leaf_darea_leaf(const S& area_leaf) const;
  // Mass of stem needed for new unit area leaf, d m_s / d a_l
  S dmass_sapwood_darea_leaf(const S& area_leaf,
                             const S& sapwood_departure) const;
  // Mass of bark needed for new unit area leaf, d m_b / d a_l
  // Bark is pinned to LEAF area, not to sapwood, so it keeps pars.theta and
  // does not move with the departure. Sharing the sapwood expression would give
  // it the Huber value's excursion, which it has no business having.
  S dmass_bark_darea_leaf(const S& area_leaf) const;
  // Mass of root needed for new unit area leaf, d m_r / d a_l
  S dmass_root_darea_leaf(const S& area_leaf) const;
  // Mass of coarse root needed for new unit area leaf, d m_cr / d a_l.
  // == a_cr1 * dmass_sapwood_darea_leaf, so it carries the same height scaling,
  // and the same Huber value: coarse roots continue the sapwood the plant has.
  S dmass_coarse_root_darea_leaf(const S& area_leaf,
                                 const S& sapwood_departure) const;
  // Growth rate of heartwood area: it depends on the sapwood already there,
  // not on leaf-area growth.
  S area_heartwood_dt(const S& area_leaf, const S& sapwood_departure) const;
  S mass_heartwood_dt(const S& mass_sapwood) const;

  S mass_live_given_height(const S& height) const;
  S height_given_mass_leaf(const S& mass_leaf_) const;


  S mortality_dt(const S& relative_reserves, const S& cumulative_mortality) const;
  S mortality_growth_independent_dt()const ;
  // Storage-dependent growth mortality (#517): rises smoothly as relative
  // reserves r = S/S_max deplete, bounded in [a_dG1*e^-a_dG2, a_dG1].
  S mortality_storage_dependent_dt(const S& relative_reserves) const;
  // NSC storage capacity S_max = a_st1 * mass_sapwood [kg NSC], at the
  // sapwood the plant carries.
  S storage_capacity(const S& area_leaf, const S& height,
                     const S& sapwood_departure) const;
  // Initialise the storage state for a newly germinated individual (#517).
  void set_initial_states(const TF24_Environment<S>& environment, Internals<S>& vars);
  // Survival of seedlings during establishment, from the carbon a
  // seedling produces at birth size. This form works that carbon out.
  S establishment_probability(const TF24_Environment<S>& environment);
  // The same, for a newborn whose rates have just been computed. A newborn is
  // already at birth size, so compute_rates has left that carbon in aux and the
  // leaf need not be solved there twice.
  S establishment_probability(const TF24_Environment<S>& environment,
                              const Internals<S>& vars) {
    return establishment_probability(environment,
                                     vars.aux(aux_idx_net_mass_production_dt));
  }
  // The equation the two above share.
  S establishment_probability(const TF24_Environment<S>& environment,
                              const S& net_mass_production_dt_);

  // * Competitive environment
  // Total projected leaf area above height `z` for a given plant,
  // from the plant's leaf area and 1/height. The overload below reads both from
  // the cached aux slots, where competition_effect is the leaf area the plant
  // carries (area_leaf_carried), so nothing is recomputed per call.
  S compute_competition(const S& z, const S& area_leaf_,
                        const S& height_inverse) const;
  // Model-agnostic entry point used by Individual<TF24> (#266): reads the
  // cached competition_effect and height_inverse aux slots itself.
  S compute_competition(const S& z, const Internals<S>& vars) const {
    return compute_competition(z, vars.aux(aux_idx_competition_effect),
                               vars.aux(aux_idx_height_inverse));
  }

  // The competition contribution and its vertical derivative from one pass, so
  // u^eta is evaluated once. `value` is bit-for-bit the one
  // compute_competition() returns: both read the shading model's own profile,
  // which is what keeps them equal under a flat-top one.
  with_slope<S> compute_competition_and_slope(const S& z, const Internals<S>& vars) const {
    const S& area_leaf_ = vars.aux(aux_idx_competition_effect);
    const S height_inverse = vars.aux(aux_idx_height_inverse);
    const S scale = pars.k_I * area_leaf_;
    const std::pair<S, S> Qq =
      canopy_shape.Q_and_q(z * height_inverse, z, height_inverse);
    // Negated on the way out: q is -dQ/dz and everything above here carries the
    // signed slope.
    return {scale * Qq.first, -(scale * Qq.second)};
  }


  // The fraction of root mass below soil depth `z`, for a plant rooted to
  // `rooting_depth` with shape exponent `eta_x` (pars.root_depth_shape_eta). The
  // canopy's own cumulative form is CanopyShape::Q, at pars.eta.
  S Q(const S& z, const S& rooting_depth, const S& eta_x) const;

  // The inverse of dheight_darea_leaf, so the allometry has one source.
  S darea_leaf_dheight(const S& area_leaf) const {
    return 1.0 / dheight_darea_leaf(area_leaf);
  }

  // The aim is to find a plant height that gives the correct seed mass.
  double height_seed(void) const;

  // The seed's height and leaf area at the current scalar.
  //
  // Preparation solves the height in plain arithmetic and cannot run at an active
  // scalar, so on a differentiated path the height is declared by the residual that
  // defines it -- live mass equals seed mass -- and the leaf area is derived from
  // the height it returns. The two come from one call because the leaf area's own
  // partials in the allometric constants and its chain through the height are the
  // same channel: taking either with the other fixed mixes them.
  struct SeedGeometry { S height; S area_leaf; };
  SeedGeometry seed_geometry() const {
    if constexpr (std::is_same_v<S, double>) {
      return {height_0, area_leaf_0};
    } else {
      // ⚠️ THE SLOPE IS READ, NOT TAKEN HERE. Taking it needs the strategy at a
      // tangent, and a rebind CONSTRUCTS one -- whose Leaf member builds two
      // vulnerability curves from an incomplete gamma per knot, which
      // assign_from then overwrites with a copy. This runs per rate evaluation,
      // so prepare_strategy takes the tangent once and this reads it.
      const S h = odelia::implicit_value<S>(
        height_0, dmass_dheight_0,
        [this](const S& y) -> S {
          return mass_live_given_height(y) - pars.omega;
        });
      return {h, area_leaf(h)};
    }
  }

  // Set constants within TF24_Model
  void prepare_strategy();
  // Two steps of prepare_strategy.
  void refuse_unsupported_traits() const;
  void construct_leaf();

  // The same strategy at scalar U.

  // Another strategy's values, written into this one. prepare_strategy() is
  // refused at an active scalar, so what it produced is carried rather than
  // rebuilt: the leaf, and the canopy shape below.
  //
  // A rebind hands back a fresh strategy, so anything it leaves default has to be
  // set here rather than left alone -- an assignment writes into one that already
  // exists, and a member it does not write keeps a previous state's working.
  template <class S1>
  void assign_from(const TF24_Model<S1>& src) {
    // Qualified: these are the base's, and an unqualified name is not looked up
    // in a dependent base.
    this->birth_rate_x = src.birth_rate_x;
    this->birth_rate_y = src.birth_rate_y;
    this->is_variable_birth_rate = src.is_variable_birth_rate;
    this->collect_all_auxiliary = src.collect_all_auxiliary;
    this->size_0 = src.size_0;
    this->control = src.control;
    this->name = src.name;
    this->extrinsic_drivers = src.extrinsic_drivers;

    TF24_Pars<S1> from_pars = src.pars;
    std::vector<S1*> from = from_pars.field_ptrs();
    std::vector<S*> to = pars.field_ptrs();
    // Both lists come from the same table at two scalars, so this says the
    // rebind still describes the same struct.
    util::check_length(from.size(), to.size());
    util::check_length(to.size(), TF24_Pars<S>::field_count);
    for (size_t i = 0; i < from.size(); ++i) {
      *to[i] = S(odelia::util::to_passive(*from[i]));
    }

    shading_model_ = src.shading_model_;
    eta_c = S(odelia::util::to_passive(src.eta_c));
    birth_rate_ = S(odelia::util::to_passive(src.birth_rate_));
    canopy_shape.initialise(pars.eta, shading_model_);
    height_0 = src.height_0;
    dmass_dheight_0 = src.dmass_dheight_0;
    area_leaf_0 = S(odelia::util::to_passive(src.area_leaf_0));
    leaf = src.leaf;
    storage_gate_width = src.storage_gate_width;
    storage_prod_eps = src.storage_prod_eps;
    plasticity_gate_width = src.plasticity_gate_width;
    beta_R_H = src.beta_R_H;
    beta_R_V = src.beta_R_V;
    function_integrator = src.function_integrator;
    crown = src.crown;
    // Shared, not copied: a clamp on this path fires inside a per-unit copy that
    // is then discarded, so the count has to land in storage the run still owns.
    // The missing-row flag is shared for the same reason.
    clamps.differentiated = src.clamps.differentiated;
    curvature_margin = src.curvature_margin;
    recorded_refusal = src.recorded_refusal;
    leaf_points = src.leaf_points;

    // Sized, not copied: these hold one right-hand side's working, and are
    // rewritten before they are read.
    root_carbon_per_leaf_area_.assign(src.root_carbon_per_leaf_area_.size(),
                                      S(0.0));

    // The index maps and the slot numbers are a function of the names, so they
    // are derived rather than carried.
    refresh_indices();
  }

  // This strategy at another scalar, carrying no derivative -- assign_from takes
  // every parameter through to_passive. It is what an implicit node asks a
  // parameter set for when it needs a slope with the parameters fixed.
  template <class U>
  TF24_Model<U> rebind() const {
    return TF24_Model<U>(tf24_strategy_from_t{}, *this);
  }

  // Birth height of a (germinated) seed. Model-agnostic accessor used by
  // the templated Individual; here height_0 is derived in prepare_strategy().
  double initial_height() const { return height_0; }

  // Crown shading model, resolved once from control.shading_model in
  // prepare_strategy(). TF24 supports deep-crown, mean-light (its default)
  // and crown-centre; PPA is not available for TF24.
  ShadingModel shading_model_ = ShadingModel::MeanLight;

  // Biological (user-settable) parameters; TF24_Pars is in tf24_pars.h.
  TF24_Pars<S> pars;

  // The exponent of the stem path integral, beta = 2*D_c + theta_c. The factor 2
  // on D_c is the packing limit: under a conserved lumen fraction, widening is
  // paid for by proportionally fewer conduits, so sapwood-specific conductivity
  // scales as D^2 and not the D^4 of Hagen-Poiseuille. See
  // plant/stem_hydraulics.h.
  //
  // What it implies: leaf-specific resistance grows as H^(1-beta) rather than
  // linearly with height. beta = 0 is the linear case; the default beta = 0.4
  // gives H^0.6; beta >= 1 saturates, so resistance approaches a finite limit no
  // matter how tall the plant grows. Larger beta therefore means a weaker height
  // penalty on carbon gain. Named for the path rather than for either parameter,
  // since it belongs to neither.
  //
  // DO NOT cache this in a member. rebind() runs before the directions are
  // set on ad_parameters(), so a copy taken at rebind carries the passive
  // value and D_c's column comes back an exact zero -- a number that reads as an
  // answer. Read here, it carries whatever tape identity pars.D_c has at the
  // call. prepare_strategy() validates (beta, L_tip) off this same expression,
  // so nothing reaches stem_hydraulics::effective_path_length that the guards
  // there have not seen.
  S stem_path_exponent() const { return S(2.0) * pars.D_c + pars.theta_c; }

  // Derived / precomputed in prepare_strategy() (NOT user-set) -------------
  S eta_c     = NA_REAL; // crown shape factor, precomputed from pars.eta
  CanopyShape<S> canopy_shape;
  // Height and leaf area of a (germinated) seed
  double height_0  = NA_REAL;
  // dM/dh at the seed height above, which is what closes the implicit function
  // theorem on it in seed_geometry. Derived here rather than there because the
  // tangent that takes it needs a whole strategy rebound, and seed_geometry runs
  // on the rate-evaluation path. NA_REAL until prepare_strategy runs, so a
  // strategy that skipped it stops in implicit_value's own dF/dy guard rather
  // than returning a number.
  double dmass_dheight_0 = NA_REAL;
  S area_leaf_0;

  // Embedded leaf hydraulic/photosynthesis sub-model, built in prepare_strategy()
  Leaf leaf;

  // Width of the smooth reserve gate G(r) on growth (#517); small relative to
  // [0,1] so the switch about the growth threshold a_st2 is fairly sharp but
  // differentiable. storage_prod_eps smooths the positive-part of net production
  // (replacing the old hard net>0 growth cutoff) for AD-readiness.
  double storage_gate_width = 0.1;
  double storage_prod_eps   = 1e-4;
  // Width of the replacement gate, in units of the marginal leaf return (#516).
  double plasticity_gate_width = 0.02;

  // The clamp sites this strategy reaches; the list itself is shared with the
  // environment, which reaches the others.
  clamp_counter clamps;

  // Count a clamp against the path it fired on. Which path this is is a property
  // of the scalar, so it is decided at compile time.
  void note_clamp(int site) {
    if constexpr (std::is_same_v<S, double>) {
      ++clamps.forward[site];
    } else {
      ++(*clamps.differentiated)[site];
    }
  }

  // How close the interior derivation came to dividing by nothing. The guard
  // below refuses on a floor, and a floor nothing approaches and a floor nothing
  // reaches report the same green -- so the margin is carried out beside the
  // count rather than left to be assumed. Shared for the same reason the counts
  // are: the copy that measures it is discarded.
  std::shared_ptr<double> curvature_margin = std::make_shared<double>(-1.0);
  void note_curvature(double value) {
    const double m = std::abs(value);
    if (*curvature_margin < 0.0 || m < *curvature_margin) {
      *curvature_margin = m;
    }
  }

  // Below this the collar's own response is amplification rather than an answer.
  // It lives on Control because it moves which states answer, so two gradients
  // taken at different values are gradients of different functions and
  // stand_gradient() has to be able to say so.
  //
  // The units are the profit's, so a reparameterisation that rescales profit
  // needs it re-measured -- which is what curvature_margin is for.
  double curvature_floor() const { return this->control.gradient_curvature_floor; }

  // Where a row a census metric needs does not exist. The values still do, so the
  // recording carries them on and leaves the missing rows off the tape.
  //
  // Recorded rather than thrown, and a recording is the reason: a throw from here
  // unwinds through a live tape and the whole sweep to deliver what a poll after
  // it delivers, and there is nothing for it to save, because a metric is a sum
  // and a sum has no defined value with an undefined term.
  //
  // Shared for the reason the clamp counts are: the strategy that records it is a
  // per-unit copy the sweep discards, so a refusal on the copy is one nothing can
  // read. Cleared where the call that reads it starts, so it latches for one
  // gradient call rather than for one operating point.
  std::shared_ptr<refusal> recorded_refusal = std::make_shared<refusal>();

  // How many leaf solves landed in each operating-point kind, indexed by the
  // enum. Diagnostic and deliberately outside rebind: a block copies the
  // strategy per unit and discards it, so carrying the tally across would count
  // the sweep's copies as well as the run.
  std::vector<size_t> operating_point_counts =
    std::vector<size_t>(Leaf::operating_point_kind_count, 0);

  // The root-architecture model's two constants. The leaf takes resistances and
  // this strategy owns the model that produces them: the 1/3 : 2/3 root split
  // and the dz^2 vertical scaling are not gas exchange.
  double beta_R_H = 3.4e2;
  double beta_R_V = 9.4e3;

  // Cached aux/state indices, resolved once in refresh_indices(), so the hot
  // compute_rates path does not do a std::map<string,int>::at (string compare)
  // lookup per ODE derivs evaluation per individual (profile hot spot).
  int aux_idx_competition_effect = -1;
  int aux_idx_height_inverse = -1;
  int aux_idx_net_mass_production_dt = -1;
  int aux_idx_root_mass = -1;
  int aux_idx_opt_psi_leaf = -1;
  int aux_idx_opt_root_psi = -1;
  int aux_idx_transpiration = -1;
  int aux_idx_E_up = -1;
  int aux_idx_profit = -1;
  int aux_idx_shadow_cost = -1;
  int aux_idx_stom_cond_CO2 = -1;
  int aux_idx_assimilation = -1;
  int aux_idx_Tleaf = -1;
  int aux_idx_leaf_marginal_return = -1;
  int aux_idx_sapwood_marginal_return = -1;
  int aux_idx_area_sapwood = -1;       // only present when collect_all_auxiliary
  int state_idx_area_heartwood = -1;
  int state_idx_mass_heartwood = -1;
  int state_idx_storage        = -1;
  // ln(area_leaf / area_leaf(height)): how far the canopy sits from the leaf
  // area the height allometry prefers. Signed, so not a non_negative_state.
  int state_idx_log_area_leaf_departure = -1;
  // ln(area_sapwood / (theta * area_leaf)): how far the conducting stem sits
  // from the leaf area it supports, i.e. the departure of the Huber value from
  // its preferred one. Also signed.
  int state_idx_log_area_sapwood_departure = -1;
  // Running mean of the marginal leaf's carbon balance. The gate reads THIS,
  // not the instantaneous value, so a plant does not shed and re-flush with
  // every passing dry spell.
  int state_idx_leaf_marginal_tracked = -1;

  // For integrating functions with using Gauss-Kronrod quadrature
  quadrature::QK function_integrator;
  // The quadrature over this strategy's unit crown, which the mean light reads
  // the light field on.
  crown_rule crown;

  // Reusable per-layer root-carbon buffer, refilled (not reallocated) each
  // net_mass_production_dt call to avoid a heap allocation per derivs eval.
  // Carries S: per-layer root carbon is one of the state directions the leaf's
  // supplied Jacobian has rows for, so it is a live gradient channel.
  std::vector<S> root_carbon_per_leaf_area_;

  // The same numbers with the derivative stripped, because the architecture
  // model and the leaf both take double. Held rather than made per call for the
  // reason the buffer above is.
  std::vector<double> root_carbon_value_;

  // Filled in place by the architecture model each call and taken by the leaf as
  // a const reference, so it must not be moved from: that is what keeps its
  // capacity across calls.
  phylloptim::RootNetwork root_network_;

  // The soil potentials with the derivative stripped, for the reason
  // root_carbon_value_ is.
  std::vector<double> psi_soil_value_;

  // The first reason is kept: a later one is a consequence of the same
  // degeneracy.
  void refuse(const std::string& why) {
    if (!recorded_refusal->happened()) {
      recorded_refusal->reason = why;
    }
  }

  // The leaf's own sites, one site at a time. Reading them as a folded vector
  // allocates, and a differentiated leaf evaluation reads them twice.
  using leaf_clamp_tally =
      std::array<std::size_t, phylloptim::CLAMP_SITE_COUNT>;
  static_assert(CLAMP_LEAF_FIRST + phylloptim::CLAMP_SITE_COUNT ==
                    CLAMP_SITE_COUNT,
                "the leaf's sites are the last of plant's, so this tally maps "
                "onto them by adding CLAMP_LEAF_FIRST");

  leaf_clamp_tally leaf_clamps() const {
    leaf_clamp_tally out{};
    for (std::size_t s = 0; s < out.size(); ++s) {
      out[s] = leaf.clamp_count(static_cast<int>(s));
    }
    return out;
  }

  // What the leaf clamped between two readings of its own tally. Everything
  // between re-supplies and re-solves it many times, so anything clamped there
  // is the gradient's and the forward run's share is the total less this.
  void note_leaf_clamps(const leaf_clamp_tally& before) {
    const leaf_clamp_tally after = leaf_clamps();
    for (std::size_t s = 0; s < after.size(); ++s) {
      if (after[s] > before[s]) {
        (*clamps.differentiated)[CLAMP_LEAF_FIRST + s] += after[s] - before[s];
      }
    }
  }

  // What this species' leaf solves found, kept against the rate evaluation that
  // found them so a pass re-running the model over the same states places the point
  // instead of searching for it again. The patch hands the address down.
  //
  // Shared, not copied, for the reason the clamp counts are: the pass that reads the
  // record runs on a rebound strategy.
  std::shared_ptr<leaf_solved_points> leaf_points =
      std::make_shared<leaf_solved_points>();
  // What one of this species' rate evaluations solves for.
  using solved_values = leaf_solved_points::points;
  void store_solved(solved_values& into) { leaf_points->store_into(into); }
  void load_solved(const solved_values& from) { leaf_points->load_from(from); }
  void end_solved() { leaf_points->end(); }
  size_t leaf_placements() const { return leaf_points->placements(); }


  // The leaf's two outputs on the active chain, carrying its supplied Jacobian.
  // Written by net_mass_production_dt before compute_rates reads either. The
  // carbon is the leaf's, the objective plus the price it deducted for its water
  // (`shadow_cost`), as phylloptim bills it (net_mass_production_dt says why).
  S leaf_carbon_;
  std::vector<S> leaf_soil_consumption_;

  // Flexible allometry's working (#516), written by net_mass_production_dt and
  // read by compute_rates in the same evaluation.
  //
  // The marginal leaf's carbon balance; dimensionless, zero at break-even.
  S marginal_leaf_return_ = 0.0;
  // The withheld shares the budget actually charged, kept so compute_rates
  // cannot form them differently from the budget that used them. The carbon not
  // spent and the tissue lost are two halves of one decision, so they must come
  // from one number.
  S withheld_leaf_ = 0.0;
  S withheld_sapwood_ = 0.0;
  // Respiration plus turnover (kg yr^-1), the strictly-positive scale the
  // sapwood controller normalises by.
  S maintenance_flux_ = 0.0;
  // The leaf-specific max hydraulic conductance actually handed to the leaf
  // (kg m^-2 LA s^-1 MPa^-1), stored so the sapwood controller differentiates
  // the value the model used rather than a recomputed copy of it.
  S leaf_specific_conductance_max_ = 0.0;
  // d(profit)/d(kmax) at the leaf's own operating point, umol m^-2 s^-1 per unit
  // conductance. A finite difference through the leaf, so a double; zero unless
  // sapwood acclimation is on.
  double dprofit_dkmax_ = 0.0;

  // Every active value this strategy holds: the whole parameter table, the
  // constant birth rate, the leaf outputs, and the per-layer root carbon. `leaf`
  // is not among them -- it solves in double, which is what keeps the leaf off
  // the tape.
  // ⚠️ A LEAF OUTPUT LEFT OFF THIS LIST CARRIES NO ROWS, SILENTLY.
  template <class F>
  void for_each_active(F&& f) {
    for (const ad_parameter& field : pars.ad_parameter_table()) {
      f(pars.*field.at);
    }
    f(birth_rate_);
    odelia::ode::visit_active(f, leaf_carbon_, leaf_soil_consumption_,
                              root_carbon_per_leaf_area_,
                              marginal_leaf_return_, withheld_leaf_,
                              withheld_sapwood_, maintenance_flux_,
                              leaf_specific_conductance_max_);
  }
};

template <typename S>
typename TF24_Model<S>::ptr make_strategy_ptr(TF24_Model<S> s);

// Named here for clarity; promotion to user-tunable traits (RcppR6) is a
// deliberate follow-up (see vignettes/models/code_review_leaf_tf24.qmd #9).
// rescales total fine-root mass into the per-layer carbon units expected by the
// root hydraulic network in Leaf::set_physiology.
inline const double root_mass_carbon_scale = 83.26 * 0.5;

// TODO: Document consistent argument order: l, b, s, h, r
// TODO: Document ordering of different types of variables (size
// before physiology, before compound things?)
// TODO: Consider moving to activating as an initialisation list?
template <typename S>
TF24_Model<S>::TF24_Model() {
  this->collect_all_auxiliary = false;
  // build the string state/aux name to index map
  refresh_indices();
  this->name = "TF24";
}

template <typename S>
void TF24_Model<S>::refresh_indices () {
    // Create and fill the name to state index maps
  this->state_index = std::map<std::string,int>();
  this->aux_index   = std::map<std::string,int>();
  std::vector<std::string> aux_names_vec = aux_names();
  std::vector<std::string> state_names_vec = state_names();
  for (size_t i = 0; i < state_names_vec.size(); i++) {
    this->state_index[state_names_vec[i]] = i;
  }
  for (size_t i = 0; i < aux_names_vec.size(); i++) {
    this->aux_index[aux_names_vec[i]] = i;
  }

  // Cache integer indices for the keys used in the hot compute_rates path, so it
  // does no std::map<string,int> lookup per derivs evaluation.
  aux_idx_competition_effect    = this->aux_index.at("competition_effect");
  aux_idx_height_inverse        = this->aux_index.at("height_inverse");
  aux_idx_net_mass_production_dt = this->aux_index.at("net_mass_production_dt");
  aux_idx_root_mass             = this->aux_index.at("root_mass");
  aux_idx_opt_psi_leaf          = this->aux_index.at("opt_psi_leaf");
  aux_idx_opt_root_psi          = this->aux_index.at("opt_root_psi");
  aux_idx_transpiration         = this->aux_index.at("transpiration");
  aux_idx_E_up                  = this->aux_index.at("E_up_");
  aux_idx_profit                = this->aux_index.at("profit");
  aux_idx_shadow_cost           = this->aux_index.at("shadow_cost");
  aux_idx_stom_cond_CO2         = this->aux_index.at("stom_cond_CO2");
  aux_idx_assimilation          = this->aux_index.at("assimilation");
  aux_idx_Tleaf                 = this->aux_index.at("Tleaf");
  aux_idx_leaf_marginal_return  = this->aux_index.at("leaf_marginal_return");
  aux_idx_sapwood_marginal_return =
    this->aux_index.at("sapwood_marginal_return");
  // area_sapwood is only registered when collect_all_auxiliary is set.
  aux_idx_area_sapwood = this->aux_index.count("area_sapwood") ? this->aux_index.at("area_sapwood") : -1;
  state_idx_area_heartwood      = this->state_index.at("area_heartwood");
  state_idx_mass_heartwood      = this->state_index.at("mass_heartwood");
  state_idx_storage             = this->state_index.at("storage");
  state_idx_log_area_leaf_departure =
    this->state_index.at("log_area_leaf_departure");
  state_idx_log_area_sapwood_departure =
    this->state_index.at("log_area_sapwood_departure");
  state_idx_leaf_marginal_tracked =
    this->state_index.at("leaf_marginal_return_tracked");
  check_state_layout(this->state_index, state_size(), "TF24");
}

template <typename S>
void TF24_Model<S>::prepare_strategy() {

  // Set up the function_integrator
  function_integrator = quadrature::QK(
      // Gauss-Kronrod quadrature integeration rule (see qkrules)
      this->control.function_integration_rule);

  // Resolve the crown shading model once. The empty Control default maps to
  // TF24's own default (mean-light, its established behaviour); PPA is an
  // FF16-only stepped-light model and is rejected here.
  shading_model_ =
    shading_model_from_string(this->control.shading_model, ShadingModel::MeanLight);
  // PPA and the flat-top-box variants are FF16-only (they reshape the FF16
  // competition / light profile, which TF24 does not use).
  if (shading_model_ == ShadingModel::PPA ||
      shading_model_ == ShadingModel::FlatTopBox ||
      shading_model_ == ShadingModel::FlatTopSoftBox) {
    throw std::invalid_argument(
      "shading_model '" + this->control.shading_model +
      "' is not supported for the TF24 strategy");
  }

  // ⚠️ Sapwood acclimation needs d(profit)/d(kmax), which is measured by
  // re-solving the leaf at ONE light. Under deep-crown the profit is an integral
  // over crown positions instead, so that single-point difference is not its
  // derivative. Refuse rather than skip: skipping leaves the controller with
  // only its (negative) cost term, so the stem would shrink without bound and
  // the run would look plausible the whole way down.
  if (odelia::util::to_passive(pars.a_sw) > 0.0 &&
      shading_model_ == ShadingModel::DeepCrown) {
    throw std::invalid_argument(
      "TF24 sapwood acclimation (a_sw > 0) is not supported with the "
      "deep-crown shading model: the sapwood sensitivity is measured at a "
      "single light and deep-crown integrates over the crown");
  }

  canopy_shape.initialise(pars.eta, shading_model_);
  crown = canopy_shape.unit_crown(function_integrator);

  eta_c = CanopyShape<S>::eta_c(pars.eta);
  // NOTE: Also pre-computing, though less trivial
  height_0 = height_seed();
  area_leaf_0 = area_leaf(height_0);
  // The residual seed_geometry closes is five composed allometries, so its slope
  // is not one term to write down: it is a tangent through the same expression,
  // taken here where the choice is visible and the scalar is already double.
  {
    using tangent = odelia::ode::tangent_scalar<double>;
    const TF24_Model<tangent> at_tangent = rebind<tangent>();
    tangent probe = height_0;
    odelia::ode::set_direction(probe, 1.0);
    dmass_dheight_0 = odelia::ode::derivative_along(
        at_tangent.mass_live_given_height(probe));
  }

  refuse_unsupported_traits();

  this->set_birth_rate_driver();
  if (!this->is_variable_birth_rate) {
    birth_rate_ = this->birth_rate_y[0];
  }
  construct_leaf();
}

// The trait values the model cannot run, refused once the seed size is known,
// since L_tip is measured against it.
template <typename S>
void TF24_Model<S>::refuse_unsupported_traits() const {
  // theta_c is declared but NOT YET USABLE. It profiles theta along the flow
  // path, and theta is not a hydraulics-only trait: it also sets area_sapwood,
  // area_bark, their growth rates, mass_sapwood (hence construction cost,
  // respiration, turnover and NSC capacity) and the hard-coded
  // dmass_sapwood_darea_leaf derivative. Those all still read a flat pars.theta.
  //
  // Applying the profile to the hydraulic term alone would give a plant whose
  // stem conducts as though theta varied while it is built and respired as
  // though theta were constant -- two different plants sharing one trait. There
  // is no staged version of this worth having, so it is refused rather than
  // half-applied. Lift the guard in the same change that profiles theta
  // everywhere.
  //
  // Read through to_passive because prepare_strategy also runs on a rebound
  // strategy carrying an active scalar, where the comparison is of the value
  // either way and an active operand only records a statement to reach it.
  if (odelia::util::to_passive(pars.theta_c) != 0.0) {
    throw std::invalid_argument(
      "theta_c is not implemented yet: theta also sets sapwood and bark area, "
      "construction cost, respiration and storage capacity, and those still "
      "use a constant theta. A hydraulics-only theta profile would be "
      "physically inconsistent, so it is refused rather than half-applied. Use "
      "D_c to vary the height dependence of resistance.");
  }

  // A negative offset can zero the storage rate's divisor or flip its sign.
  if (!(odelia::util::to_passive(pars.storage_relaxation_offset) >= 0.0)) {
    throw std::invalid_argument("storage_relaxation_offset must be >= 0");
  }

  // a_pl0 is the largest fraction of turnover a plant can leave unreplaced.
  // Above 1 replacement goes negative, destroying tissue faster than it dies;
  // below 0 a struggling plant replaces more than dies. Written to refuse NaN.
  const double a_pl0_value = odelia::util::to_passive(pars.a_pl0);
  if (!(a_pl0_value >= 0.0 && a_pl0_value <= 1.0)) {
    throw std::invalid_argument(
      "a_pl0 must be in [0, 1]: it is the largest fraction of leaf and sapwood "
      "turnover left unreplaced, and a_pl0 = " + std::to_string(a_pl0_value) +
      " would make replacement negative or greater than turnover");
  }

  if (odelia::util::to_passive(stem_path_exponent()) != 0.0) {
    // L_tip is the anchor of both profiles, so it cannot be zero once either is
    // active: k_s(L) = K_s*(L/L_tip)^(2*D_c) diverges everywhere as L_tip -> 0,
    // giving zero resistance. That is a degenerate configuration to reject, not
    // a numerical edge case to tolerate.
    //
    // Written to reject zero and NaN as well as negatives: `L_tip < 0.0` alone
    // would let a zero through, which is the case this guard exists for.
    const double L_tip_value = odelia::util::to_passive(pars.L_tip);
    if (L_tip_value <= 0.0 || std::isnan(L_tip_value)) {
      throw std::invalid_argument(
        "L_tip must be > 0 when D_c or theta_c is non-zero: the within-plant "
        "profiles are defined relative to the terminal segment, and L_tip -> 0 "
        "sends sapwood-specific conductivity to infinity everywhere");
    }
    // A plant cannot be shorter than one terminal segment. Worth catching here
    // rather than downstream: height_0 is solved from seed mass, so a user
    // sweeping omega down at a fixed L_tip will eventually cross this, and the
    // symptom is a negative path length, hence a negative conductance, hence an
    // unattributable NaN twenty frames inside the leaf solver.
    if (!(L_tip_value < odelia::util::to_passive(height_0 * eta_c))) {
      throw std::invalid_argument(
        "L_tip must be shorter than the birth-size flow path (height_0*eta_c): "
        "a plant cannot be smaller than one terminal segment");
    }
  }
}

// The leaf this strategy solves, built from its traits and seated on the
// layered supply and the TF24_floor cost curve.
template <typename S>
void TF24_Model<S>::construct_leaf() {
  if constexpr (std::is_same_v<S, double>) {
    // (P50, c) per curve, and phylloptim derives b and psi_crit itself. Handing
    // over the derived pair instead would state the curve twice and let the two
    // disagree silently, which is why they are not in the trait set.
    leaf = Leaf(pars.vcmax_25, pars.stem_c, pars.stem_P50,
                pars.root_c, pars.root_P50,
                pars.TF24_beta2, pars.jmax_25, pars.phi_J,
                pars.curv_fact_elec_trans, pars.curv_fact_colim,
                this->control.GSS_tol_abs, this->control.vulnerability_curve_ncontrol,
                this->control.ci_abs_tol, this->control.ci_niter,
                pars.TF24_cost_scale);
    // phylloptim's Leaf defaults to the single-potential supply on the leaf
    // route; this strategy builds a root network, so it asks for the layered
    // supply, which also seats the collar route.
    leaf.set_supply_multilayer();
    // Not a constructor argument, and written before any physiology is set, so
    // the temperature block derives R_d_ from it on the first call.
    leaf.R_d_25 = pars.R_d_25;
    // Not a constructor argument either. The leaf's own default is zero, so at
    // the default this writes the value the leaf already held.
    leaf.g_min = pars.g_min;
    // Penman-Monteith leaf energy balance (#523): enable per pars (default off,
    // backward-compatible) and pass the leaf-dimension trait. Wind speed is a
    // per-timestep driver, set from the environment before each set_physiology.
    //
    // ⚠️ BOTH ARE SETTABLE AND REACH NOTHING IF THIS IS DROPPED. Moving the model
    // body out of src/tf24_strategy.cpp and into this header lost the pair, and
    // the symptom was Tleaf reporting the air temperature exactly on every run
    // with the balance switched on -- an answer, not an error.
    leaf.use_energy_balance_ = (odelia::util::to_passive(pars.use_energy_balance) != 0.0);
    leaf.d_ = odelia::util::to_passive(pars.d);
    leaf.alpha_PAR_ = odelia::util::to_passive(pars.alpha_PAR);
    leaf.alpha_NIR_ = odelia::util::to_passive(pars.alpha_NIR);
    // ⚠️ SEATED ON TF24_floor, WHICH IS TF24 AT lambda_o = 0 -- and exactly, not
    // approximately: each term is TF24's own expression, so zeroing the price
    // adds an exact zero to TF24's exact value. The default price IS zero, so no
    // run moves; what changes is that the price now reaches the model at all.
    // Without this the parameter is settable and inert, which is the silent
    // failure the whole (P50, c) reparameterisation was about avoiding.
    leaf.TF24_floor_lambda_o = pars.TF24_floor_lambda_o;
    leaf.set_model(phylloptim::Leaf::CostCurve::TF24_floor, true);
  } else {
    static_assert(std::is_same_v<S, double>,
                  "Leaf carries double; an active strategy must supply the leaf's "
                  "local Jacobian across this boundary, not template Leaf.");
  }
}

template <typename S>
typename TF24_Model<S>::ptr make_strategy_ptr(TF24_Model<S> s) {
  s.prepare_strategy();
  return std::make_shared<TF24_Model<S> >(s);
}

}

// The member definitions, by what they compute.
#include <plant/models/tf24_allometry.h>
#include <plant/models/tf24_leaf.h>
#include <plant/models/tf24_rates.h>

#endif
