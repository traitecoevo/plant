// -*-c++-*-
#ifndef PLANT_PLANT_FF16_STRATEGY_H_
#define PLANT_PLANT_FF16_STRATEGY_H_

#include <plant/census.h>
#include <plant/strategy.h>
#include <plant/models/ff16_environment.h>
#include <plant/qag.h>
#include <plant/canopy_shape.h>
#include <plant/models/ff16_production_kernel.h>
#include <plant/with_slope.h>
#include <odelia/ode_util.hpp>
#include <odelia/implicit_node.hpp>
#include <array>
#include <cmath>
#include <string_view>
#include <type_traits>
#include <utility>

namespace plant {

// Biological (user-settable) parameters for the FF16 strategy. Held as a value
// member `pars` on FF16_Strategy and exposed to R as a nested RcppR6 list class
// (so R access is `s$pars$lma`). Derived/precomputed quantities (eta_c,
// height_0, canopy_shape, ...) are NOT here -- they are outputs of
// prepare_strategy() and stay as plain members on the strategy.
//
// Templated on the scalar S so that a gradient can be taken with respect to
// every member: ad_parameter_fields below lists them, and the sizeof assert
// after the struct makes that list total.
template <typename S = double>
struct FF16_Pars {
  using value_type = S;

  // * Core traits
  S lma       = 0.1978791;  // Leaf mass per area [kg / m2]
  S rho       = 608.0;      // Wood density [kg/m3]
  S hmat      = 16.5958691; // Height at maturation [m]
  S omega     = 3.8e-5;     // Seed mass [kg]
  // * Individual allometry
  // Canopy shape parameter
  S eta       = 12.0; // [dimensionless]
  // Sapwood area per leaf area
  // Ratio sapwood area area to leaf area
  S theta     = 1.0/4669; // [dimensionless]
  // Height - leaf mass scaling
  S a_l1        = 5.44; // height with 1m2 leaf [m]
  S a_l2        = 0.306; // dimensionless scaling of height with leaf area
  // Root mass per leaf area
  S a_r1        = 0.07;  //[kg / m]
  // Ratio of bark area : sapwood area
  S a_b1         = 0.17; // [dimensionless]

  // * Production
  // Ratio of leaf dark respiration to leaf mass [mol CO2 / yr  / kg]
  // =  [mol CO2 / m2 / yr]  |  (39.27 = 2100 * 0.00187)  | narea * photosynthesis_per_nitrogen
  //    / [kg(leaf) / m2 ]   |    / (0.1978791)           | lma
  // Hard coded in value of lma here so that this value doesn't change
  // if that trait changes above.
  S r_l    = 39.27 / 0.1978791;
  // Root respiration per mass [mol CO2 / yr / kg]
  S r_r    = 217.0;
  // Sapwood respiration per stem mass  [mol CO2 / yr / kg]
  // = respiration per volume [mol CO2 / m3 / yr]
  // /  wood density [kg/m3]
  S r_s    = 4012.0 / 608.0;
  // Bark respiration per stem mass
  // assumed to be twice rate of sapwood
  // (NOTE that there is a re-parametrisation here relative to the paper
  // -- r_b is defined (new) as 2*r_s, whereas the paper assumes a
  // fixed multiplication by 2)
  // Spelled out rather than as 2.0 * r_s: r_b and r_s are two independent
  // fields (the hyperparameter function sets each), and the literal keeps the
  // default readable in a struct whose members are a template scalar.
  S r_b    = 2.0 * (4012.0 / 608.0);
  // Carbon conversion parameter
  S a_y    = 0.7;
  // Constant converting assimilated CO2 to dry mass [kg / mol]
  // (12E-3 / 0.49)
  S a_bio  = 2.45e-2;
  // Leaf turnover [/yr]
  S k_l    =  0.4565855;
  // Bark turnover [/yr]
  S k_b    = 0.2;
  // Sapwood turnover [/yr]
  S k_s           = 0.2;
  // Root turnover [/yr]
  S k_r    = 1.0;
  // Parameters of the hyperbola for annual LRC
  S a_p1   = 151.177775377968; // [mol CO2 / yr / m2]
  S a_p2   = 0.204716166503633; // [dimensionless]

  // * Seed production
  // Accessory cost of reproduction
  S a_f3  = 3.0 *  3.8e-5; // [kg per seed]
  // Maximum allocation to reproduction
  S a_f1   = 1.0; //[dimensionless]
  // Size range across which individuals mature
  S a_f2   = 50; // [dimensionless]

  // * Mortality parameters
  // Probability of survival during dispersal
  S S_D   = 0.25; // [dimensionless]
  // Parameter for seedling survival
  S a_d0    = 0.1; //[kg / yr / m2]
  // Baseline for intrinsic mortality
  S d_I    = 0.01; // [ / yr]
  // Baseline rate for growth-related mortality
  S a_dG1    = 5.5; // [ / yr]
  // Risk coefficient for dry mass production (per area)
  S a_dG2    = 20.0;// [yr m2 / kg ]

  // Germination
  S recruitment_decay = 0.0;

  // * Light capture parameters
  S k_I = 0.5;

  // One registered parameter: where it lives and what its column is called.
  struct ad_parameter {
    // A member pointer rather than an address in one instance, which is what
    // lets the table below be constexpr.
    S FF16_Pars::* at;
    const char* name;
  };

  // The parameters no gradient exists for, each with the sentence that says why.
  // Each is excluded by name rather than given a silent zero, so asking for one
  // errors in these words. Nothing else is excluded: every other member reaches
  // a rate, the seed's residual or the fitness metrics' own reading -- S_D
  // multiplies the net reproduction ratio directly, and a_f1, a_f2 and a_f3
  // reach fecundity_dt, which the fitness metrics read and the census does not,
  // so their census columns are zero on any trajectory rather than absent.
  struct no_gradient {
    const char* name;
    const char* why;
  };
  static constexpr std::array undifferentiable{
      no_gradient{"eta",
                  "the canopy profile's u^eta reaches u = 0, where the recorded "
                  "derivative u^eta*log(u) is 0*(-inf) and the guard returns a "
                  "constant zero instead, so a row here would be a silently "
                  "wrong zero rather than a NaN; and at the specialised eta "
                  "values the profile is a multiplication chain eta does not "
                  "appear in at all"}};

  // Every member above, in declaration order. field_ptrs() and both
  // ad_parameter_* projections read this, so a member reaches all of them or
  // none, and the sizeof assert below makes the list total.
#define PLANT_FF16_AD_PARAMETER(x) ad_parameter{&FF16_Pars::x, #x}
  static constexpr std::array ad_parameter_fields{
      PLANT_FF16_AD_PARAMETER(lma),
      PLANT_FF16_AD_PARAMETER(rho),
      PLANT_FF16_AD_PARAMETER(hmat),
      PLANT_FF16_AD_PARAMETER(omega),
      PLANT_FF16_AD_PARAMETER(eta),
      PLANT_FF16_AD_PARAMETER(theta),
      PLANT_FF16_AD_PARAMETER(a_l1),
      PLANT_FF16_AD_PARAMETER(a_l2),
      PLANT_FF16_AD_PARAMETER(a_r1),
      PLANT_FF16_AD_PARAMETER(a_b1),
      PLANT_FF16_AD_PARAMETER(r_l),
      PLANT_FF16_AD_PARAMETER(r_r),
      PLANT_FF16_AD_PARAMETER(r_s),
      PLANT_FF16_AD_PARAMETER(r_b),
      PLANT_FF16_AD_PARAMETER(a_y),
      PLANT_FF16_AD_PARAMETER(a_bio),
      PLANT_FF16_AD_PARAMETER(k_l),
      PLANT_FF16_AD_PARAMETER(k_b),
      PLANT_FF16_AD_PARAMETER(k_s),
      PLANT_FF16_AD_PARAMETER(k_r),
      PLANT_FF16_AD_PARAMETER(a_p1),
      PLANT_FF16_AD_PARAMETER(a_p2),
      PLANT_FF16_AD_PARAMETER(a_f3),
      PLANT_FF16_AD_PARAMETER(a_f1),
      PLANT_FF16_AD_PARAMETER(a_f2),
      PLANT_FF16_AD_PARAMETER(S_D),
      PLANT_FF16_AD_PARAMETER(a_d0),
      PLANT_FF16_AD_PARAMETER(d_I),
      PLANT_FF16_AD_PARAMETER(a_dG1),
      PLANT_FF16_AD_PARAMETER(a_dG2),
      PLANT_FF16_AD_PARAMETER(recruitment_decay),
      PLANT_FF16_AD_PARAMETER(k_I)
  };
#undef PLANT_FF16_AD_PARAMETER

  // Every entry of the table, which is the whole parameter set: a rebind carries
  // it across a scalar change. ad_parameters() is the subset with a column.
  std::vector<S*> field_ptrs() {
    std::vector<S*> ret;
    ret.reserve(ad_parameter_fields.size());
    for (const ad_parameter& f : ad_parameter_fields) {
      ret.push_back(&(this->*f.at));
    }
    return ret;
  }

  static constexpr size_t field_count = ad_parameter_fields.size();

  // A parameter carries a column unless the list above says no gradient exists
  // for it.
  static constexpr bool has_column(std::string_view name) {
    for (const no_gradient& n : undifferentiable) {
      if (name == n.name) {
        return false;
      }
    }
    return true;
  }

  static constexpr size_t column_count = [] {
    size_t n = 0;
    for (const ad_parameter& f : ad_parameter_fields) {
      if (has_column(f.name)) ++n;
    }
    return n;
  }();

  // A name matching no member excludes nothing, so the two counts stop adding
  // up: a parameter renamed or removed is a compile error here rather than a
  // list entry that has quietly stopped excluding anything.
  static_assert(column_count + undifferentiable.size() == field_count,
                "undifferentiable names a parameter ad_parameter_fields does "
                "not carry, or names one twice");

  static constexpr std::array<ad_parameter, column_count> ad_parameter_columns =
      [] {
        std::array<ad_parameter, column_count> ret{};
        std::size_t n = 0;
        for (const ad_parameter& f : ad_parameter_fields) {
          if (has_column(f.name)) {
            ret[n++] = f;
          }
        }
        return ret;
      }();

  static constexpr const auto& ad_parameter_table() { return ad_parameter_fields; }
};

// Every member of FF16_Pars is an S, so one added without extending
// ad_parameter_fields changes this size and is refused rather than dropped.
static_assert(sizeof(FF16_Pars<double>) ==
              FF16_Pars<double>::field_count * sizeof(double),
              "FF16_Pars has a member ad_parameter_fields does not list");

// Templated on the scalar S the state, the traits and everything derived from
// them carry; double is production, and the whole model is in this header so
// that one body serves both. The double path is bit-identical to the model as
// it stood in src/ff16_strategy.cpp: every expression is transcribed in its
// operand order.
template <typename S = double>
class FF16_Strategy: public Strategy<FF16_Environment<S>> {
public:
  using value_type = S;

  typedef std::shared_ptr<FF16_Strategy> ptr;
  FF16_Strategy();

  // Scientific version. Bump ONLY when equations or default parameters change
  // the simulation output for identical inputs. Do NOT bump for refactors,
  // performance, interface, or serialisation changes. Bumping invalidates
  // logpile's cache for this model (see plant::model_version() / model_id()).
  // v2: the reductions over the size distribution take their trapezium widths
  // from the coordinate the density is carried in, so a run on the birth-date
  // coordinate integrates over birth dates rather than over heights. Measured on
  // the full-lifetime deep-crown anchor, offspring production moves
  // 16.884586 -> 17.172004, i.e. **+1.70%**; on the height coordinate it is
  // unchanged to the last bit.
  static constexpr int scientific_version = 2;

  // Fixed integer slots for the hot ODE rate path, used instead of
  // state_index.at("...") / aux_index.at("...") string-map lookups (those map
  // lookups showed up in profiling, see #466). These MUST stay in sync with
  // the order of state_names() and aux_names() below: *_INDEX is the position
  // of that name in state_names(), *_AUX_INDEX the position in aux_names(). If
  // you add/reorder a name there, update these constants (refresh_indices()
  // still validates the named maps used by the R-facing paths).
  static constexpr int AREA_HEARTWOOD_INDEX = 3;
  static constexpr int MASS_HEARTWOOD_INDEX = 4;
  static constexpr int COMPETITION_EFFECT_AUX_INDEX = 0;
  static constexpr int HEIGHT_INVERSE_AUX_INDEX = 1;
  static constexpr int NET_MASS_PRODUCTION_DT_AUX_INDEX = 2;
  static constexpr int AREA_SAPWOOD_AUX_INDEX = 3;

  // Overrides ----------------------------------------------

  // update this when the length of state_names changes
  static size_t state_size () { return 5; }
  // update this when the length of aux_names changes
  size_t aux_size () { return aux_names().size(); }

  static std::vector<std::string> state_names() {
    return  std::vector<std::string>({
      "height",
      "mortality",
      "fecundity",
      "area_heartwood",
      "mass_heartwood"
      });
  }

  std::vector<std::string> aux_names() {
    std::vector<std::string> ret({
      "competition_effect",
      "height_inverse",
      "net_mass_production_dt"
    });
    // add the associated computation to compute_rates and compute there
    if (this->collect_all_auxiliary) {
      ret.push_back("area_sapwood");
    }
    return ret;
  }

  // Each species' net reproduction ratio and offspring production are offered
  // as differentiable stand metrics beside the census (SCM::stand_metrics).
  // Declared here because the fecundity channel is refereed for this model:
  // test-ff16-fitness-gradient.R checks every route into it against a pinned
  // central difference.
  static constexpr bool differentiable_fitness = true;

  // The metrics a census of this model sums, in the order it reports them: the
  // same three TF24 declares, read off this model's own slots.
  static const auto& census_metrics() {
    using strategy = FF16_Strategy<S>;
    static constexpr std::array<census_metric<strategy>, 3> metrics{{
      {"leaf_area",
       [](const strategy& p, const Internals<S>& vars) -> S {
         return p.area_leaf(vars.state(HEIGHT_INDEX));
       }},
      {"mass_above_ground",
       [](const strategy& p, const Internals<S>& vars) -> S {
         const S height = vars.state(HEIGHT_INDEX);
         const S area_leaf = p.area_leaf(height);
         return p.mass_above_ground(
             p.mass_leaf(area_leaf),
             p.mass_bark(p.area_bark(area_leaf), height),
             p.mass_sapwood(p.area_sapwood(area_leaf), height),
             vars.state(MASS_HEARTWOOD_INDEX));
       }},
      {"area_stem",
       [](const strategy& p, const Internals<S>& vars) -> S {
         const S height = vars.state(HEIGHT_INDEX);
         const S area_leaf = p.area_leaf(height);
         return p.area_stem(p.area_bark(area_leaf), p.area_sapwood(area_leaf),
                            vars.state(AREA_HEARTWOOD_INDEX));
       }},
    }};
    return metrics;
  }

  using ad_parameter = typename FF16_Pars<S>::ad_parameter;

  // The birth rate a constant-birth-rate species is run at, held at the
  // working scalar so that it is a gradient column: it reaches the inflow
  // condition (the density every cohort is born at) and the offspring
  // integral's scalars. A variable birth rate is a driver read as a double and
  // carries no column, so that species' gradient has one column fewer rather
  // than a zero that reads as an answer. Set by prepare_strategy() from
  // birth_rate_y, carried by assign_from().
  S birth_rate_ = 1.0;
  bool has_birth_rate_column() const { return !this->is_variable_birth_rate; }
  S birth_rate_at(double t) const {
    return has_birth_rate_column()
             ? birth_rate_
             : S(this->extrinsic_drivers.evaluate("birth_rate", t));
  }

  // How many gradient columns one of these carries, for a caller sizing a batch.
  size_t ad_column_count() const {
    return FF16_Pars<S>::column_count + (has_birth_rate_column() ? 1 : 0);
  }

  // Addresses of the parameters a gradient can be taken with respect to, in the
  // order ad_parameter_names() gives: the parameter table's columns, then the
  // birth rate where it is one. The strategy is shared and the fields do not
  // move, so take them once per gradient evaluation.
  std::vector<S*> ad_parameters() {
    std::vector<S*> ret;
    ret.reserve(ad_column_count());
    for (const ad_parameter& p : FF16_Pars<S>::ad_parameter_columns) {
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
    for (const ad_parameter& p : FF16_Pars<S>::ad_parameter_columns) {
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
    for (const auto& n : FF16_Pars<S>::undifferentiable) {
      ret.emplace_back(n.name, n.why);
    }
    return ret;
  }

  // Translate generic methods to FF16 strategy leaf area methods

  S competition_effect(const S& height) const {
    return area_leaf(height);
  }

  void refresh_indices();


  // FF16 Methods  ----------------------------------------------

  // [eqn 2] area_leaf (inverse of [eqn 3])
  // Inline (header) so it can inline into the hot competition/assimilation
  // paths that reach it from templated Individual<FF16> code (no LTO build).
  S area_leaf(const S& height) const {
    // Single source: scalar-templated kernel (#472 scope B). Bit-identical to
    // the previous std::pow(height/a_l1, 1/a_l2); validated by the FF16
    // reference-comparison test.
    return ff16_area_leaf(pars.a_l1, pars.a_l2, height);
  }

  // [eqn 1] mass_leaf (inverse of [eqn 2])
  S mass_leaf(const S& area_leaf) const;

  // [eqn 4] area and mass of sapwood
  S area_sapwood(const S& area_leaf) const;
  S mass_sapwood(const S& area_sapwood, const S& height) const;

  // [eqn 5] area and mass of bark
  S area_bark(const S& area_leaf) const;
  S mass_bark (const S& area_bark, const S& height) const;

  S area_stem(const S& area_bark, const S& area_sapwood,
              const S& area_heartwood) const;
  S diameter_stem(const S& area_stem) const;

  // [eqn 7] Mass of (fine) roots
  S mass_root(const S& area_leaf) const;

  // [eqn 8] Total Mass
  S mass_live(const S& mass_leaf, const S& mass_bark,
              const S& mass_sapwood, const S& mass_root) const;

  S mass_total(const S& mass_leaf, const S& mass_bark, const S& mass_sapwood,
               const S& mass_heartwood, const S& mass_root) const;

  // Above-ground mass = leaf + all stem components (bark + sapwood +
  // heartwood); excludes roots.
  S mass_above_ground(const S& mass_leaf, const S& mass_bark,
                      const S& mass_sapwood, const S& mass_heartwood) const;

  void compute_rates(const FF16_Environment<S>& environment,
                     Internals<S>& vars);

  // Inline (header): called per state-set / ODE-state update from templated
  // Individual<FF16> code, so inlining avoids a cross-TU call (no LTO build)
  // and lets the now-inline area_leaf fold in.
  void update_dependent_aux(const int index, Internals<S>& vars) {
    if (index == HEIGHT_INDEX) {
      S height = vars.state(HEIGHT_INDEX);
      vars.set_aux(COMPETITION_EFFECT_AUX_INDEX, area_leaf(height));
      vars.set_aux(HEIGHT_INVERSE_AUX_INDEX, 1.0 / height);
    }
  }

  // * Mass production
  // [eqn 12] Gross annual CO2 assimilation. Thin dispatcher: forwards to the
  // shading-model implementation bound once in prepare_strategy(), so the
  // choice of crown model costs a single predicted indirect call per
  // derivative evaluation and never a string comparison.
  S assimilation(const FF16_Environment<S>& environment, const S& height,
                 const S& area_leaf, const S& height_inverse) {
    return (this->*assimilation_fn)(environment, height, area_leaf,
                                    height_inverse);
  }
  // Shading-model implementations of assimilation (selected in
  // prepare_strategy via assimilation_fn).
  //  - deep crown: integrate photosynthesis over crown depth (Yokozawa q).
  //  - crown top:  single evaluation of the light at the crown centre. Used by
  //                both crown-centre and PPA; they differ only in how the patch
  //                light profile is built (smooth vs stepped), which this read
  //                picks up transparently through the environment.
  S assimilation_deep_crown(const FF16_Environment<S>& environment,
                            const S& height, const S& area_leaf,
                            const S& height_inverse);
  //  - average light: integrate the light over crown depth to a leaf-area-
  //                weighted mean, then a single photosynthesis evaluation.
  S assimilation_average_light(const FF16_Environment<S>& environment,
                               const S& height, const S& area_leaf,
                               const S& height_inverse);
  S assimilation_crown_top(const FF16_Environment<S>& environment,
                           const S& height, const S& area_leaf,
                           const S& height_inverse);

  typedef S (FF16_Strategy::*assimilation_fn_t)(const FF16_Environment<S>&,
                                                const S&, const S&, const S&);
  // Bound once in prepare_strategy(); defaults to the deep-crown integral.
  assimilation_fn_t assimilation_fn = &FF16_Strategy::assimilation_deep_crown;

  // [Appendix S6] Per-leaf photosynthetic rate.
  S assimilation_leaf(const S& x) const;

  // [eqn 13] Total maintenance respiration
  S respiration(const S& mass_leaf, const S& mass_sapwood,
                const S& mass_bark, const S& mass_root) const;

  S respiration_leaf(const S& mass) const;
  S respiration_bark(const S& mass) const;
  S respiration_sapwood(const S& mass) const;
  S respiration_root(const S& mass) const;

  // [eqn 14] Total turnover
  S turnover(const S& mass_leaf, const S& mass_bark,
             const S& mass_sapwood, const S& mass_root) const;
  S turnover_leaf(const S& mass) const;
  S turnover_bark(const S& mass) const;
  S turnover_sapwood(const S& mass) const;
  S turnover_root(const S& mass) const;

  // [eqn 15] Net production
  S net_mass_production_dt_A(const S& assimilation, const S& respiration,
                             const S& turnover) const;

  S net_mass_production_dt(const FF16_Environment<S>& environment,
                           const S& height, const S& area_leaf_);
  S net_mass_production_dt(const FF16_Environment<S>& environment,
                           const S& height, const S& area_leaf_,
                           const S& height_inverse);
  // Worker overload that also reports the sapwood intermediates so callers
  // (compute_rates) can reuse them instead of recomputing. Bit-identical: the
  // out-refs receive exactly the values the body already computed.
  S net_mass_production_dt(const FF16_Environment<S>& environment,
                           const S& height, const S& area_leaf_,
                           const S& height_inverse,
                           S& area_sapwood_, S& mass_sapwood_);
  // Strategy-agnostic entry point used by Individual<FF16> (#266): reads the
  // height state and the cached aux slots itself, so the generic Individual
  // does not need to know FF16's state/aux layout.
  S net_mass_production_dt(const FF16_Environment<S>& environment,
                           const Internals<S>& vars) {
    return net_mass_production_dt(environment, vars.state(HEIGHT_INDEX),
                                  vars.aux(COMPETITION_EFFECT_AUX_INDEX),
                                  vars.aux(HEIGHT_INVERSE_AUX_INDEX));
  }

  // [eqn 16] Fraction of whole plan growth that is leaf
  S fraction_allocation_reproduction(const S& height) const;
  S fraction_allocation_growth(const S& height) const;
  // [eqn 17] Rate of offspring production
  S fecundity_dt(const S& net_mass_production_dt,
                 const S& fraction_allocation_reproduction) const;

  // [eqn 18] Fraction of mass growth that is leaves
  S darea_leaf_dmass_live(const S& area_leaf) const;

  // change in height per change in leaf area
  S dheight_darea_leaf(const S& area_leaf) const;
  // d(leaf area)/d(height), for a caller pulling a leaf-area adjoint back.
  S darea_leaf_dheight(const S& area_leaf) const {
    return 1.0 / dheight_darea_leaf(area_leaf);
  }
  // Mass of leaf needed for new unit area leaf, d m_s / d a_l
  S dmass_leaf_darea_leaf(const S& area_leaf) const;
  // Mass of stem needed for new unit area leaf, d m_s / d a_l
  S dmass_sapwood_darea_leaf(const S& area_leaf) const;
  // Overloads taking a precomputed pow(area_leaf, a_l2): compute_rates needs
  // this term for both the height rate and the live-mass partition, so it
  // evaluates the (non-integer, libm) pow once and shares it -- see issue #361.
  S darea_leaf_dmass_live(const S& area_leaf,
                          const S& area_leaf_pow_a_l2) const;
  S dheight_darea_leaf(const S& area_leaf, const S& area_leaf_pow_a_l2) const;
  S dmass_sapwood_darea_leaf(const S& area_leaf,
                             const S& area_leaf_pow_a_l2) const;
  // Mass of bark needed for new unit area leaf, d m_b / d a_l
  S dmass_bark_darea_leaf(const S& area_leaf) const;
  // Mass of root needed for new unit area leaf, d m_r / d a_l
  S dmass_root_darea_leaf(const S& area_leaf) const;
  // Growth rate of basal diameter_stem per unit stem area
  S ddiameter_stem_darea_stem(const S& area_stem) const;
  // Growth rate of components per unit time:
  S area_sapwood_dt(const S& area_leaf_dt) const;
  S area_heartwood_dt(const S& area_leaf) const;
  S area_bark_dt(const S& area_leaf_dt) const;
  S area_stem_dt(const S& area_leaf, const S& area_leaf_dt) const;
  S diameter_stem_dt(const S& area_stem, const S& area_stem_dt) const;
  S mass_root_dt(const S& area_leaf, const S& area_leaf_dt) const;
  S mass_live_dt(const S& fraction_allocation_reproduction,
                 const S& net_mass_production_dt) const;
  S mass_total_dt(const S& fraction_allocation_reproduction,
                  const S& net_mass_production_dt,
                  const S& mass_heartwood_dt) const;
  S mass_above_ground_dt(const S& area_leaf,
                         const S& fraction_allocation_reproduction,
                         const S& net_mass_production_dt,
                         const S& mass_heartwood_dt,
                         const S& area_leaf_dt) const;

  S mass_heartwood_dt(const S& mass_sapwood) const;

  S mass_live_given_height(const S& height) const;
  S height_given_mass_leaf(const S& mass_leaf_) const;


  S mortality_dt(const S& productivity_area, const S& cumulative_mortality) const;
  S mortality_growth_independent_dt() const;
  S mortality_growth_dependent_dt(const S& productivity_area) const;

  // The seed's height, written at birth. On double this is the height the
  // Individual was constructed with; on an active scalar it is declared by its
  // own residual so that it carries the seed mass's and the allometry's
  // derivative -- see seed_geometry().
  void set_initial_states(const FF16_Environment<S>& environment,
                          Internals<S>& vars);

  // [eqn 20] Survival of seedlings during establishment, from the carbon a
  // seedling produces at birth size. This form works that carbon out.
  S establishment_probability(const FF16_Environment<S>& environment);
  // The same, for a newborn whose rates have just been computed. A newborn is
  // already at birth size, so compute_rates has left that carbon in aux and the
  // production need not be computed there twice.
  S establishment_probability(const FF16_Environment<S>& environment,
                              const Internals<S>& vars) {
    return establishment_probability(environment,
                                     vars.aux(NET_MASS_PRODUCTION_DT_AUX_INDEX));
  }
  // The equation the two above share. The seed's leaf area is derived once by
  // whichever caller has it: on an active scalar seed_geometry() is an
  // implicit node per call, and the one-argument form had made two.
  S establishment_probability(const FF16_Environment<S>& environment,
                              const S& net_mass_production_dt_);
  S establishment_probability(const FF16_Environment<S>& environment,
                              const S& net_mass_production_dt_,
                              const S& seed_area_leaf);

  // * Competitive environment
  // [eqn 11] total projected leaf area above height above height `z` for given plant
  // Inline (header) so the per-node hot competition path called from
  // Individual<FF16>::compute_competition can inline these tiny helpers
  // instead of paying a cross-TU call each iteration (no LTO build).
  S compute_competition(const S& z, const S& height) const {
    return compute_competition(z, area_leaf(height), 1.0 / height);
  }
  S compute_competition(const S& z, const S& area_leaf_,
                        const S& height_inverse) const {
    return compute_competition_by_ratio(z * height_inverse, area_leaf_);
  }
  S compute_competition_by_ratio(const S& z_over_height,
                                 const S& area_leaf_) const {
    return pars.k_I * area_leaf_ * canopy_shape.leaf_area_above(z_over_height);
  }
  // Strategy-agnostic entry point used by Individual<FF16> (#266): reads the
  // cached competition_effect (= area_leaf) and height_inverse aux slots
  // itself. Inline (header) to keep the per-node hot competition path free of
  // a cross-TU call (no LTO build).
  S compute_competition(const S& z, const Internals<S>& vars) const {
    return compute_competition(z, vars.aux(COMPETITION_EFFECT_AUX_INDEX),
                               vars.aux(HEIGHT_INVERSE_AUX_INDEX));
  }

  // The competition contribution and its vertical derivative from one pass, so
  // u^eta is evaluated once. The first entry is bit-for-bit the one
  // compute_competition() returns.
  with_slope<S>
  compute_competition_and_slope(const S& z, const Internals<S>& vars) const {
    const S area_leaf_ = vars.aux(COMPETITION_EFFECT_AUX_INDEX);
    const S height_inverse = vars.aux(HEIGHT_INVERSE_AUX_INDEX);
    const S scale = pars.k_I * area_leaf_;
    const std::pair<S, S> Qq =
      canopy_shape.Q_and_q(z * height_inverse, z, height_inverse);
    return {scale * Qq.first, -(scale * Qq.second)};
  }

  // [      ] Inverse of Q: height above which fraction 'x' of leaf found
  S Qp(const S& x, const S& height) const;

  // The aim is to find a plant height that gives the correct seed mass.
  double height_seed(void) const;

  // The seed's height and leaf area at the current scalar.
  //
  // Preparation solves the height in plain arithmetic and cannot run at an active
  // scalar, so on a differentiated path the height is declared by the residual that
  // defines it -- live mass equals seed mass -- and the leaf area is derived from
  // the height it returns. The two come from one call because the leaf area's own
  // partials in the allometric constants and its chain through the height are the
  // same channel: taking either against the other held fixed mixes them.
  struct SeedGeometry { S height; S area_leaf; };
  SeedGeometry seed_geometry() const {
    if constexpr (std::is_same_v<S, double>) {
      return {height_0, area_leaf_0};
    } else {
      // The slope is read, not taken here: prepare_strategy takes it once by a
      // tangent through the same residual, and this runs per rate evaluation.
      const S h = odelia::implicit_value<S>(
        height_0, dmass_dheight_0,
        [this](const S& y) -> S {
          return mass_live_given_height(y) - pars.omega;
        });
      return {h, area_leaf(h)};
    }
  }

  // Set constants within FF16_Strategy
  void prepare_strategy();

  // Another strategy's values, written into this one. prepare_strategy() is
  // refused at an active scalar, so what it produced is carried rather than
  // rebuilt. A rebind hands back a fresh strategy, so anything it leaves default
  // has to be set here rather than left alone.
  template <class S1>
  void assign_from(const FF16_Strategy<S1>& src) {
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

    FF16_Pars<S1> from_pars = src.pars;
    std::vector<S1*> from = from_pars.field_ptrs();
    std::vector<S*> to = pars.field_ptrs();
    // Both lists come from the same table at two scalars, so this says the
    // rebind still describes the same struct.
    util::check_length(from.size(), to.size());
    util::check_length(to.size(), FF16_Pars<S>::field_count);
    for (size_t i = 0; i < from.size(); ++i) {
      *to[i] = S(odelia::util::to_passive(*from[i]));
    }

    // eta carries no column, so the shape factor is carried passive, as the
    // canopy shape's own eta is.
    eta_c = S(odelia::util::to_passive(src.eta_c));
    birth_rate_ = S(odelia::util::to_passive(src.birth_rate_));
    bind_shading_model(src.shading_model_);
    height_0 = src.height_0;
    dmass_dheight_0 = src.dmass_dheight_0;
    area_leaf_0 = S(odelia::util::to_passive(src.area_leaf_0));
    function_integrator = src.function_integrator;

    // The index maps and the slot numbers are a function of the names, so they
    // are derived rather than carried.
    refresh_indices();
  }

  // This strategy at another scalar, carrying no derivative -- assign_from takes
  // every parameter through to_passive.
  template <class U>
  FF16_Strategy<U> rebind_from() const {
    FF16_Strategy<U> out;
    out.assign_from(*this);
    return out;
  }

  // Birth height of a (germinated) seed. Strategy-agnostic accessor used by
  // the templated Individual; here height_0 is derived in prepare_strategy().
  double initial_height() const { return height_0; }

  // Every active value this strategy holds: the whole parameter table. Nothing
  // else on it carries a derivative between rate evaluations.
  template <class F>
  void for_each_active(F&& f) {
    for (const ad_parameter& field : pars.ad_parameter_table()) {
      f(pars.*field.at);
    }
    f(birth_rate_);
  }

  // Biological (user-settable) parameters; see FF16_Pars above.
  FF16_Pars<S> pars;

  // Crown shading model, resolved once from control.shading_model in
  // prepare_strategy() and carried across a rebind.
  ShadingModel shading_model_ = ShadingModel::DeepCrown;

  // Derived / precomputed in prepare_strategy() (NOT user-set) -------------
  // Crown shape factor, precomputed from pars.eta
  S eta_c     = NA_REAL; // [dimensionless]
  CanopyShape<S> canopy_shape;
  // Height and leaf area of a (germinated) seed
  double height_0  = NA_REAL;
  // dM/dh at the seed height above, which is what closes the implicit function
  // theorem on it in seed_geometry. NA_REAL until prepare_strategy runs, so a
  // strategy that skipped it stops in implicit_value's own dF/dy guard rather
  // than returning a number.
  double dmass_dheight_0 = NA_REAL;
  S area_leaf_0;

  // For integrating functions with using Gauss-Kronrod quadrature
  quadrature::QK function_integrator;

private:
  // Bind the competition profile and the assimilation implementation to one
  // shading model. Called by prepare_strategy() from the Control string and by
  // assign_from() from the source's resolved model, so the two cannot drift.
  void bind_shading_model(ShadingModel shading_model) {
    shading_model_ = shading_model;
    // canopy_shape also selects the competition contribution: smooth Q for every
    // model except flat-top-box, which casts a step (see CanopyShape).
    canopy_shape.initialise(pars.eta, shading_model);
    switch (shading_model) {
    case ShadingModel::DeepCrown:
      assimilation_fn = &FF16_Strategy::assimilation_deep_crown;
      break;
    case ShadingModel::MeanLight:
      assimilation_fn = &FF16_Strategy::assimilation_average_light;
      break;
    case ShadingModel::CrownCentre:
    case ShadingModel::FlatTopBox:
    case ShadingModel::FlatTopSoftBox:
    case ShadingModel::PPA:
      // All evaluate assimilation at the crown centre. They differ in the light
      // profile they read: crown-centre from the smooth profile, flat-top-box / -soft-
      // box from a profile built with (hard / smoothed) box competition, PPA from
      // a stepped profile.
      assimilation_fn = &FF16_Strategy::assimilation_crown_top;
      break;
    }
  }
};

template <typename S>
typename FF16_Strategy<S>::ptr make_strategy_ptr(FF16_Strategy<S> s);

// ---- Definitions ---------------------------------------------------------

template <typename S>
FF16_Strategy<S>::FF16_Strategy() {
  this->collect_all_auxiliary = false;
  // build the string state/aux name to index map
  refresh_indices();
  this->name = "FF16";
}

template <typename S>
void FF16_Strategy<S>::refresh_indices () {
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
  check_state_layout(this->state_index, state_size(), "FF16");
}

// area_leaf() is defined inline in the class (hot path).

// [eqn 1] mass_leaf (inverse of [eqn 2])
template <typename S>
S FF16_Strategy<S>::mass_leaf(const S& area_leaf) const {
  return area_leaf * pars.lma;
}

// [eqn 4] area and mass of sapwood
template <typename S>
S FF16_Strategy<S>::area_sapwood(const S& area_leaf) const {
  return area_leaf * pars.theta;
}

template <typename S>
S FF16_Strategy<S>::mass_sapwood(const S& area_sapwood, const S& height) const {
  return area_sapwood * height * eta_c * pars.rho;
}

// [eqn 5] area and mass of bark
template <typename S>
S FF16_Strategy<S>::area_bark(const S& area_leaf) const {
  return pars.a_b1 * area_leaf * pars.theta;
}

template <typename S>
S FF16_Strategy<S>::mass_bark(const S& area_bark, const S& height) const {
  return area_bark * height * eta_c * pars.rho;
}

template <typename S>
S FF16_Strategy<S>::area_stem(const S& area_bark, const S& area_sapwood,
                              const S& area_heartwood) const {
  return area_bark + area_sapwood + area_heartwood;
}

template <typename S>
S FF16_Strategy<S>::diameter_stem(const S& area_stem) const {
  using std::sqrt;
  return sqrt(4 * area_stem / M_PI);
}

// [eqn 7] Mass of (fine) roots
template <typename S>
S FF16_Strategy<S>::mass_root(const S& area_leaf) const {
  return pars.a_r1 * area_leaf;
}

// [eqn 8] Total mass
template <typename S>
S FF16_Strategy<S>::mass_live(const S& mass_leaf, const S& mass_bark,
                              const S& mass_sapwood, const S& mass_root) const {
  return mass_leaf + mass_sapwood + mass_bark + mass_root;
}

template <typename S>
S FF16_Strategy<S>::mass_total(const S& mass_leaf, const S& mass_bark,
                               const S& mass_sapwood, const S& mass_heartwood,
                               const S& mass_root) const {
  return mass_leaf + mass_bark + mass_sapwood +  mass_heartwood + mass_root;
}

template <typename S>
S FF16_Strategy<S>::mass_above_ground(const S& mass_leaf, const S& mass_bark,
                                      const S& mass_sapwood, const S& mass_heartwood) const {
  return mass_leaf + mass_bark + mass_sapwood + mass_heartwood;
}

// update_dependent_aux() is defined inline in the class (hot path).

// one-shot update of the scm variables
// i.e. setting rates of ode vars from the state and updating aux vars
template <typename S>
void FF16_Strategy<S>::compute_rates(const FF16_Environment<S>& environment,
                                     Internals<S>& vars) {
  using std::pow;

  S height = vars.state(HEIGHT_INDEX);
  S area_leaf_ = vars.aux(COMPETITION_EFFECT_AUX_INDEX);
  S height_inverse = vars.aux(HEIGHT_INVERSE_AUX_INDEX);

  // Reuse the sapwood intermediates the worker already computes (for
  // respiration/turnover) rather than recomputing them below; bit-identical.
  S area_sapwood_, mass_sapwood_;
  const S net_mass_production_dt_ =
    net_mass_production_dt(environment, height, area_leaf_, height_inverse,
                           area_sapwood_, mass_sapwood_);

  // store the aux sate
  vars.set_aux(NET_MASS_PRODUCTION_DT_AUX_INDEX, net_mass_production_dt_);

  if (net_mass_production_dt_ > 0) {

    const S fraction_allocation_reproduction_ = fraction_allocation_reproduction(height);
    // dheight_darea_leaf and the sapwood/bark terms in darea_leaf_dmass_live all
    // share pow(area_leaf, pars.a_l2); evaluate this libm pow once and reuse it for
    // both rates rather than paying it twice per node per step (issue #361).
    const S area_leaf_pow_a_l2 = pow(area_leaf_, pars.a_l2);
    const S darea_leaf_dmass_live_ = darea_leaf_dmass_live(area_leaf_, area_leaf_pow_a_l2);
    const S fraction_allocation_growth_ = fraction_allocation_growth(height);
    const S area_leaf_dt = net_mass_production_dt_ * fraction_allocation_growth_ * darea_leaf_dmass_live_;

    vars.set_rate(HEIGHT_INDEX, dheight_darea_leaf(area_leaf_, area_leaf_pow_a_l2) * area_leaf_dt);
    vars.set_rate(FECUNDITY_INDEX,
      fecundity_dt(net_mass_production_dt_, fraction_allocation_reproduction_));

    vars.set_rate(AREA_HEARTWOOD_INDEX, area_heartwood_dt(area_leaf_));
    vars.set_rate(MASS_HEARTWOOD_INDEX, mass_heartwood_dt(mass_sapwood_));

    if (this->collect_all_auxiliary) {
      vars.set_aux(AREA_SAPWOOD_AUX_INDEX, area_sapwood_);
    }
  } else {
    vars.set_rate(HEIGHT_INDEX, S(0.0));
    vars.set_rate(FECUNDITY_INDEX, S(0.0));
    vars.set_rate(AREA_HEARTWOOD_INDEX, S(0.0));
    vars.set_rate(MASS_HEARTWOOD_INDEX, S(0.0));
  }
  // [eqn 21] - Instantaneous mortality rate
  vars.set_rate(MORTALITY_INDEX,
      mortality_dt(net_mass_production_dt_ / area_leaf_, vars.state(MORTALITY_INDEX)));
}

// [eqn 12] Gross annual CO2 assimilation -- deep-crown model.
// Integrate photosynthesis over crown depth: for a given height in the crown,
// take photosynthesis at that depth multiplied by the amount of leaf there.
template <typename S>
S FF16_Strategy<S>::assimilation_deep_crown(const FF16_Environment<S>& environment,
                                            const S& height,
                                            const S& area_leaf,
                                            const S& height_inverse) {

  S A = 0.0;

  // Define an anonymous function to integrate.
  // Keep the lambda's own closure type (do not wrap in std::function) so the
  // templated QK::integrate inlines the integrand at each quadrature point
  // instead of making a type-erased indirect call.
  // Hoist the light-spline upper bound (canopy top) out of the integrand: it
  // is invariant across the quadrature, so fetch it once and pass it into the
  // capped get_environment_at_height() overload rather than re-reading
  // spline.max() at every quadrature point.
  const double canopy_top = environment.max_environment_height();
  auto f = [&](const S& z) -> S {
    return assimilation_leaf(environment.get_environment_at_height(z, canopy_top)) *
      canopy_shape.q(z * height_inverse, z);
  };

  // Integrate over crown depth using using Gauss-Kronrod quadrature.
  // The number of points used in the integration is determined by the control parameter
  // function_integration_rule. Rules defined in qk_rules.cpp
  A = function_integrator.integrate(f, S(0.0), height);

  return area_leaf * A;
}

// [eqn 12] Gross annual CO2 assimilation -- mean-light model.
// Integrate the *light* over crown depth, weighted by the leaf-area density q
// (which integrates to one over the crown), to get the leaf-area-weighted mean
// light the crown experiences, then take a single photosynthesis evaluation of
// that mean. This sits between deep-crown (which integrates the concave
// photosynthetic rate itself) and crown-centre (a single point evaluation): it
// captures the mean light exactly but ignores the curvature of photosynthesis
// across the within-crown light distribution.
template <typename S>
S FF16_Strategy<S>::assimilation_average_light(const FF16_Environment<S>& environment,
                                               const S& height,
                                               const S& area_leaf,
                                               const S& height_inverse) {
  const double canopy_top = environment.max_environment_height();
  auto f = [&](const S& z) -> S {
    return environment.get_environment_at_height(z, canopy_top) *
      canopy_shape.q(z * height_inverse, z);
  };
  const S mean_light = function_integrator.integrate(f, S(0.0), height);
  return area_leaf * assimilation_leaf(mean_light);
}

// [eqn 12] Gross annual CO2 assimilation -- crown-top model (crown-centre and PPA).
// Leaf area is treated as a thin layer at the crown centre, so a single
// evaluation of the light environment there replaces the crown-depth integral.
// The crown-centre and PPA models share this code: under crown-centre the environment
// returns the smooth light profile, under PPA it returns the stepped (layered)
// profile, so the only difference between them lives in the environment build.
// (height_inverse is unused but kept for a common dispatch signature.)
template <typename S>
S FF16_Strategy<S>::assimilation_crown_top(const FF16_Environment<S>& environment,
                                           const S& height,
                                           const S& area_leaf,
                                           const S& /* height_inverse */) {
  const S E = environment.get_environment_at_height(height * eta_c);
  return area_leaf * assimilation_leaf(E);
}

// Photosynthetic rate per leaf area
// `x` is openness, ranging from 0 to 1.
template <typename S>
S FF16_Strategy<S>::assimilation_leaf(const S& x) const {
  // Single source: scalar-templated kernel (#472 scope B, Milestone A).
  return ff16_assimilation_leaf(pars.a_p1, pars.a_p2, x);
}

// [eqn 13] Total maintenance respiration
// NOTE: In contrast with Falster ref model, we do not normalise by pars.a_y*pars.a_bio.
template <typename S>
S FF16_Strategy<S>::respiration(const S& mass_leaf, const S& mass_sapwood,
                                const S& mass_bark, const S& mass_root) const {
  // Single source: scalar-templated kernel (#472 scope B, Milestone A).
  return ff16_respiration(mass_leaf, mass_sapwood, mass_bark, mass_root,
                          pars.r_l, pars.r_s, pars.r_b, pars.r_r);
}

template <typename S>
S FF16_Strategy<S>::respiration_leaf(const S& mass) const {
  return pars.r_l * mass;
}

template <typename S>
S FF16_Strategy<S>::respiration_bark(const S& mass) const {
  return pars.r_b * mass;
}

template <typename S>
S FF16_Strategy<S>::respiration_sapwood(const S& mass) const {
  return pars.r_s * mass;
}

template <typename S>
S FF16_Strategy<S>::respiration_root(const S& mass) const {
  return pars.r_r * mass;
}

// [eqn 14] Total turnover
template <typename S>
S FF16_Strategy<S>::turnover(const S& mass_leaf, const S& mass_bark,
                             const S& mass_sapwood, const S& mass_root) const {
   // Single source: scalar-templated kernel (#472 scope B, Milestone A).
   return ff16_turnover(mass_leaf, mass_bark, mass_sapwood, mass_root,
                        pars.k_l, pars.k_b, pars.k_s, pars.k_r);
}

template <typename S>
S FF16_Strategy<S>::turnover_leaf(const S& mass) const {
  return pars.k_l * mass;
}

template <typename S>
S FF16_Strategy<S>::turnover_bark(const S& mass) const {
  return pars.k_b * mass;
}

template <typename S>
S FF16_Strategy<S>::turnover_sapwood(const S& mass) const {
  return pars.k_s * mass;
}

template <typename S>
S FF16_Strategy<S>::turnover_root(const S& mass) const {
  return pars.k_r * mass;
}

// [eqn 15] Net production
//
// NOTE: Translation of variable names from the Falster 2011.  Everything
// before the minus sign is SCM's N, our `net_mass_production_dt` is SCM's P.
template <typename S>
S FF16_Strategy<S>::net_mass_production_dt_A(const S& assimilation, const S& respiration,
                                             const S& turnover) const {
  // Single source: scalar-templated kernel (#472 scope B, Milestone A).
  return ff16_net_production_A(pars.a_bio, pars.a_y, assimilation, respiration, turnover);
}

// One shot calculation of net_mass_production_dt
// Used by establishment_probability() and compute_rates().
template <typename S>
S FF16_Strategy<S>::net_mass_production_dt(const FF16_Environment<S>& environment,
                                           const S& height, const S& area_leaf_) {
  return net_mass_production_dt(environment, height, area_leaf_,
                                1.0 / height);
}

template <typename S>
S FF16_Strategy<S>::net_mass_production_dt(const FF16_Environment<S>& environment,
                                           const S& height, const S& area_leaf_,
                                           const S& height_inverse) {
  S area_sapwood_, mass_sapwood_;
  return net_mass_production_dt(environment, height, area_leaf_, height_inverse,
                                area_sapwood_, mass_sapwood_);
}

template <typename S>
S FF16_Strategy<S>::net_mass_production_dt(const FF16_Environment<S>& environment,
                                           const S& height, const S& area_leaf_,
                                           const S& height_inverse,
                                           S& area_sapwood_, S& mass_sapwood_) {
  const S mass_leaf_    = mass_leaf(area_leaf_);
  area_sapwood_ = area_sapwood(area_leaf_);
  mass_sapwood_ = mass_sapwood(area_sapwood_, height);
  const S area_bark_    = area_bark(area_leaf_);
  const S mass_bark_    = mass_bark(area_bark_, height);
  const S mass_root_    = mass_root(area_leaf_);
  const S assimilation_ =
    assimilation(environment, height, area_leaf_, height_inverse);
  const S respiration_ =
    respiration(mass_leaf_, mass_sapwood_, mass_bark_, mass_root_);
  const S turnover_ =
    turnover(mass_leaf_, mass_bark_, mass_sapwood_, mass_root_);
  return net_mass_production_dt_A(assimilation_, respiration_, turnover_);
}

// [eqn 16] Fraction of production allocated to reproduction
template <typename S>
S FF16_Strategy<S>::fraction_allocation_reproduction(const S& height) const {
  using std::exp;
  return pars.a_f1 / (1.0 + exp(pars.a_f2 * (1.0 - height / pars.hmat)));
}

// Fraction of production allocated to growth
template <typename S>
S FF16_Strategy<S>::fraction_allocation_growth(const S& height) const {
  return 1.0 - fraction_allocation_reproduction(height);
}

// [eqn 17] Rate of offspring production
template <typename S>
S FF16_Strategy<S>::fecundity_dt(const S& net_mass_production_dt,
                                 const S& fraction_allocation_reproduction) const {
  return net_mass_production_dt * fraction_allocation_reproduction /
    (pars.omega + pars.a_f3);
}

template <typename S>
S FF16_Strategy<S>::darea_leaf_dmass_live(const S& area_leaf) const {
  using std::pow;
  return darea_leaf_dmass_live(area_leaf, pow(area_leaf, pars.a_l2));
}

template <typename S>
S FF16_Strategy<S>::darea_leaf_dmass_live(const S& area_leaf,
                                          const S& area_leaf_pow_a_l2) const {
  // dmass_bark_darea_leaf(area_leaf) == pars.a_b1 * dmass_sapwood_darea_leaf(area_leaf),
  // so compute the shared pow(area_leaf, pars.a_l2) term once rather than twice.
  const S dmass_sapwood_darea_leaf_ =
    dmass_sapwood_darea_leaf(area_leaf, area_leaf_pow_a_l2);
  return 1.0/(  dmass_leaf_darea_leaf(area_leaf)
              + dmass_sapwood_darea_leaf_
              + pars.a_b1 * dmass_sapwood_darea_leaf_
              + dmass_root_darea_leaf(area_leaf));
}

template <typename S>
S FF16_Strategy<S>::dheight_darea_leaf(const S& area_leaf) const {
  using std::pow;
  return pars.a_l1 * pars.a_l2 * pow(area_leaf, pars.a_l2 - 1);
}

template <typename S>
S FF16_Strategy<S>::dheight_darea_leaf(const S& area_leaf,
                                       const S& area_leaf_pow_a_l2) const {
  // pow(area_leaf, pars.a_l2 - 1) == pow(area_leaf, pars.a_l2) / area_leaf, so reuse the
  // shared power instead of a second libm pow (a reciprocal-multiply reorder,
  // so not bit-identical to the single-argument form).
  return pars.a_l1 * pars.a_l2 * area_leaf_pow_a_l2 / area_leaf;
}

// Mass of leaf needed for new unit area leaf, d m_s / d a_l
template <typename S>
S FF16_Strategy<S>::dmass_leaf_darea_leaf(const S& /* area_leaf */) const {
  return pars.lma;
}

// Mass of stem needed for new unit area leaf, d m_s / d a_l
template <typename S>
S FF16_Strategy<S>::dmass_sapwood_darea_leaf(const S& area_leaf) const {
  using std::pow;
  return dmass_sapwood_darea_leaf(area_leaf, pow(area_leaf, pars.a_l2));
}

template <typename S>
S FF16_Strategy<S>::dmass_sapwood_darea_leaf(const S& /* area_leaf */,
                                             const S& area_leaf_pow_a_l2) const {
  return pars.rho * eta_c * pars.a_l1 * pars.theta * (pars.a_l2 + 1.0) * area_leaf_pow_a_l2;
}

// Mass of bark needed for new unit area leaf, d m_b / d a_l
template <typename S>
S FF16_Strategy<S>::dmass_bark_darea_leaf(const S& area_leaf) const {
  return pars.a_b1 * dmass_sapwood_darea_leaf(area_leaf);
}

// Mass of root needed for new unit area leaf, d m_r / d a_l
template <typename S>
S FF16_Strategy<S>::dmass_root_darea_leaf(const S& /* area_leaf */) const {
  return pars.a_r1;
}

// Growth rate of basal diameter_stem per unit time
template <typename S>
S FF16_Strategy<S>::ddiameter_stem_darea_stem(const S& area_stem) const {
  using std::pow;
  return pow(M_PI * area_stem, -0.5);
}

// Growth rate of sapwood area at base per unit time
template <typename S>
S FF16_Strategy<S>::area_sapwood_dt(const S& area_leaf_dt) const {
  return area_leaf_dt * pars.theta;
}

// Note, unlike others, heartwood growth does not depend on leaf area growth, but
// rather existing sapwood
template <typename S>
S FF16_Strategy<S>::area_heartwood_dt(const S& area_leaf) const {
  return pars.k_s * area_sapwood(area_leaf);
}

// Growth rate of bark area at base per unit time
template <typename S>
S FF16_Strategy<S>::area_bark_dt(const S& area_leaf_dt) const {
  return pars.a_b1 * area_leaf_dt * pars.theta;
}

// Growth rate of stem basal area per unit time
template <typename S>
S FF16_Strategy<S>::area_stem_dt(const S& area_leaf,
                                 const S& area_leaf_dt) const {
  return area_sapwood_dt(area_leaf_dt) +
    area_bark_dt(area_leaf_dt) +
    area_heartwood_dt(area_leaf);
}

// Growth rate of basal diameter_stem per unit time
template <typename S>
S FF16_Strategy<S>::diameter_stem_dt(const S& area_stem, const S& area_stem_dt) const {
  return ddiameter_stem_darea_stem(area_stem) * area_stem_dt;
}

// Growth rate of root mass per unit time
template <typename S>
S FF16_Strategy<S>::mass_root_dt(const S& area_leaf,
                                 const S& area_leaf_dt) const {
  return area_leaf_dt * dmass_root_darea_leaf(area_leaf);
}

template <typename S>
S FF16_Strategy<S>::mass_live_dt(const S& fraction_allocation_reproduction,
                                 const S& net_mass_production_dt) const {
  return (1 - fraction_allocation_reproduction) * net_mass_production_dt;
}

template <typename S>
S FF16_Strategy<S>::mass_total_dt(const S& fraction_allocation_reproduction,
                                  const S& net_mass_production_dt,
                                  const S& mass_heartwood_dt) const {
  return mass_live_dt(fraction_allocation_reproduction, net_mass_production_dt) +
    mass_heartwood_dt;
}

// TODO(#480): Do we not track root mass change?
template <typename S>
S FF16_Strategy<S>::mass_above_ground_dt(const S& area_leaf,
                                         const S& fraction_allocation_reproduction,
                                         const S& net_mass_production_dt,
                                         const S& mass_heartwood_dt,
                                         const S& area_leaf_dt) const {
  const S mass_root_dt =
    area_leaf_dt * dmass_root_darea_leaf(area_leaf);
  return mass_total_dt(fraction_allocation_reproduction, net_mass_production_dt,
                        mass_heartwood_dt) - mass_root_dt;
}

template <typename S>
S FF16_Strategy<S>::mass_heartwood_dt(const S& mass_sapwood) const {
  return turnover_sapwood(mass_sapwood);
}


template <typename S>
S FF16_Strategy<S>::mass_live_given_height(const S& height) const {
  S area_leaf_ = area_leaf(height);
  return mass_leaf(area_leaf_) +
         mass_bark(area_bark(area_leaf_), height) +
         mass_sapwood(area_sapwood(area_leaf_), height) +
         mass_root(area_leaf_);
}

template <typename S>
S FF16_Strategy<S>::height_given_mass_leaf(const S& mass_leaf) const {
  using std::pow;
  return pars.a_l1 * pow(mass_leaf / pars.lma, pars.a_l2);
}

template <typename S>
S FF16_Strategy<S>::mortality_dt(const S& productivity_area,
                                 const S& cumulative_mortality) const {

  // NOTE: When plants are extremely inviable, the rate of change in
  // mortality can be Inf, because net production is negative, leaf
  // area is small and so we get exp(big number).  However, most of
  // the time that happens we should get infinite mortality variable
  // levels and the rate of change won't matter.  It is possible that
  // we will need to trim this to some large finite value, but for
  // now, just checking that the actual mortality rate is finite.
  // ⚠️ THE CEILING IS TESTED BESIDE FINITENESS, and the pair must not be
  // reduced back to one test. compute_initial_conditions holds an unestablished
  // recruit's hazard at establishment_failure_hazard rather than at the +Inf
  // -log(0) gives, so this is the test that keeps such a cohort's rate parked --
  // which is what makes the finite hazard bit-identical to the infinite one.
  if (util::is_finite(cumulative_mortality) &&
      cumulative_mortality < establishment_failure_hazard) {
    return
      mortality_growth_independent_dt() +
      mortality_growth_dependent_dt(productivity_area);
 } else {
    // Mortality probability is 1, so the rate calculations have nothing left to
    // describe and the state does not move again.
    return 0.0;
  }
}

template <typename S>
S FF16_Strategy<S>::mortality_growth_independent_dt() const {
  return pars.d_I;
}

template <typename S>
S FF16_Strategy<S>::mortality_growth_dependent_dt(const S& productivity_area) const {
  using std::exp;
  return pars.a_dG1 * exp(-pars.a_dG2 * productivity_area);
}

// The seed's height is written here, not inherited from the height an
// Individual was constructed with. Construction runs in plain arithmetic and an
// active strategy receives its results rather than deriving them, so a height
// that arrives that way carries no trait derivative and every rate the newborn
// goes on to have inherits the loss. Declaring it by its own residual at this
// point puts it where the parameters are already differentiable inputs, and the
// value written is the same one either way.
template <typename S>
void FF16_Strategy<S>::set_initial_states(const FF16_Environment<S>& environment,
                                          Internals<S>& vars) {
  (void)environment;
  const SeedGeometry seed = seed_geometry();
  vars.set_state(HEIGHT_INDEX, seed.height);
  // A height written into vars leaves the slots it determines holding what the
  // constructor's plain height derived, which is the same number carrying no
  // derivative. Re-deriving them costs one allometry and is the only thing that
  // puts the seed's row into the leaf area every rate at birth size is scaled by.
  update_dependent_aux(HEIGHT_INDEX, vars);
}

// [eqn 20] Survival of seedlings during establishment
template <typename S>
S FF16_Strategy<S>::establishment_probability(const FF16_Environment<S>& environment) {
  const SeedGeometry seed = seed_geometry();
  return establishment_probability(
    environment,
    net_mass_production_dt(environment, seed.height, seed.area_leaf,
                           1.0 / seed.height),
    seed.area_leaf);
}

template <typename S>
S FF16_Strategy<S>::establishment_probability(const FF16_Environment<S>& environment,
                                              const S& net_mass_production_dt_) {
  return establishment_probability(environment, net_mass_production_dt_,
                                   seed_geometry().area_leaf);
}

// Every form above ends here. The carbon is birth-size carbon either way,
// whatever height the caller's plant happens to be at.
template <typename S>
S FF16_Strategy<S>::establishment_probability(const FF16_Environment<S>& environment,
                                              const S& net_mass_production_dt_,
                                              const S& seed_area_leaf) {
  using std::exp;

  S decay_over_time = exp(-pars.recruitment_decay * environment.time);

  if (net_mass_production_dt_ > 0) {
    const S tmp = pars.a_d0 * seed_area_leaf / net_mass_production_dt_;
    return 1.0 / (tmp * tmp + 1.0) * decay_over_time;
  } else {
    return 0.0;
  }
}

// compute_competition() overloads and compute_competition_by_ratio() are
// defined inline in the class (per-node hot path).

// (inverse of [eqn 10]; return the height above which fraction 'x' of
// the leaf mass would be found).
template <typename S>
S FF16_Strategy<S>::Qp(const S& x, const S& height) const { // x in [0,1], unchecked.
  return canopy_shape.Qp(x, height);
}

// The aim is to find a plant height that gives the correct seed mass.
template <typename S>
double FF16_Strategy<S>::height_seed(void) const {

  // Note, these are not entirely correct bounds. Ideally we would use height
  // given *total* mass, not leaf mass, but that is difficult to calculate.
  // Using "height given leaf mass" will expand upper bound, but that's ok
  // most of time. Only issue is that could break with obscure parameter
  // values for LMA or height-leaf area scaling. Could instead use some
  // absolute maximum height for new seedling, e.g. 1m?
  const S
    h0 = height_given_mass_leaf(std::numeric_limits<double>::min()),
    h1 = height_given_mass_leaf(pars.omega);

  const double tol = this->control.offspring_production_tol;
  const size_t max_iterations = this->control.offspring_production_iterations;

  auto target = [&] (double x) mutable -> S {
    return mass_live_given_height(x) - pars.omega;
  };

  if constexpr (std::is_same_v<S, double>) {
    return util::uniroot(target, h0, h1, tol, max_iterations);
  } else {
    // Bisection is affine in its bracket and blind to the residual's values, so
    // recording the search would return d(h1)/d(trait) rather than the height at
    // which mass_live equals the seed mass. Declare it by the residual instead.
    static_assert(std::is_same_v<S, double>,
                  "height_seed() finds its root by iteration; an active scalar must "
                  "declare it through implicit_value on the residual.");
  }
}

template <typename S>
void FF16_Strategy<S>::prepare_strategy() {

  // Set up the function_integrator
  function_integrator = quadrature::QK(
      // Gauss-Kronrod quadrature integeration rule (see qkrules)
      this->control.function_integration_rule);

  // Resolve the crown shading model once (string -> enum), then bind both hot
  // paths to it: canopy_shape handles competition (leaf_area_above), and
  // assimilation_fn selects the matching assimilation implementation. After
  // this, neither path compares the model string per call.
  bind_shading_model(
    shading_model_from_string(this->control.shading_model, ShadingModel::DeepCrown));

  eta_c = CanopyShape<S>::eta_c(pars.eta);
  // NOTE: Also pre-computing, though less trivial
  height_0 = height_seed();
  area_leaf_0 = area_leaf(height_0);
  // The residual seed_geometry closes is four composed allometries, so its slope
  // is not one term to write down: it is a tangent through the same expression,
  // taken here where the choice is visible and the scalar is already double.
  {
    using tangent = odelia::ode::tangent_scalar<double>;
    const FF16_Strategy<tangent> at_tangent = rebind_from<tangent>();
    tangent probe = height_0;
    odelia::ode::seed_direction(probe, 1.0);
    dmass_dheight_0 = odelia::ode::derivative_along(
        at_tangent.mass_live_given_height(probe));
  }

  if (this->is_variable_birth_rate) {
    this->extrinsic_drivers.set_variable("birth_rate", this->birth_rate_x,
                                         this->birth_rate_y,
                                         odelia::drivers::Slopes::monotone);
  } else {
    this->extrinsic_drivers.set_constant("birth_rate", this->birth_rate_y[0]);
    birth_rate_ = this->birth_rate_y[0];
  }
}

template <typename S>
typename FF16_Strategy<S>::ptr make_strategy_ptr(FF16_Strategy<S> s) {
  s.prepare_strategy();
  return std::make_shared<FF16_Strategy<S> >(s);
}

}

#endif
