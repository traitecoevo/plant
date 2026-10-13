// Built from  inst/include/plant/models/ff16_strategy.h on Fri Jul 24 10:23:19 2020 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_K93_STRATEGY_H_
#define PLANT_PLANT_K93_STRATEGY_H_

#include <plant/strategy.h>
#include <plant/models/k93_environment.h>
#include <plant/canopy_shape.h>
#include <plant/with_slope.h>

namespace plant {

// Biological (user-settable) parameters for the K93 strategy. Held as a value
// member `pars` on K93_Strategy and exposed to R as a nested RcppR6 list class
// (access as `s$pars$b_0`). Derived quantities (canopy_shape) stay as plain
// members on the strategy.
struct K93_Pars {
  // Initial seedling size (dbh cm)
  double height_0 = 2.0;
  // * Growth
  double b_0 = 0.059;    // Growth intercept (yr-1)
  double b_1 = 0.012;    // Growth asymptote (yr-1.(ln cm)-1)
  double b_2 = 0.00041;  // Growth suppression rate (m2.cm-2.yr-1)
  // * Mortality
  double c_0 = 0.008;    // Intercept (yr-1)
  double c_1 = 0.00044;  // Suppression rate (m2.cm-2.yr-1)
  // * Reproduction
  double d_0 = 0.00073;  // Recruitment rate (cm2.yr-1)
  double d_1 = 0.044;    // Reduction from suppression (m2.cm-2.yr-1)
  // Probability of survival during dispersal (required by scm.h)
  double S_D = 1.0;
  // Smoothing / canopy shape parameter
  double eta = 12;
  // Light capture parameter
  double k_I = 0.01;
};

class K93_Strategy: public Strategy<K93_Environment> {
public:
  using value_type = double;

  typedef std::shared_ptr<K93_Strategy> ptr;
  K93_Strategy();

  // Scientific version. Bump ONLY when equations or default parameters change
  // the simulation output for identical inputs. Do NOT bump for refactors,
  // performance, interface, or serialisation changes. Bumping invalidates
  // logpile's cache for this model (see plant::model_version() / model_id()).
  // v2: as FF16 v2 -- the reductions integrate over the coordinate the density
  // is carried in, so a birth-date run's output moves and a height run's does
  // not.
  static constexpr int scientific_version = 2;

  // Direct aux indices for the hot path, avoiding aux_index.at("...") string-map
  // lookups (these showed up in profiling; see #466). MUST stay in sync with the
  // order of aux_names() below. refresh_indices() still fills the named maps used
  // by the R-facing paths.
  static constexpr int COMPETITION_EFFECT_AUX_INDEX = 0;
  static constexpr int HEIGHT_INVERSE_AUX_INDEX = 1;

  // update this when the length of state_names changes
  static size_t state_size () { return 3; }
  // update this when the length of aux_names changes
  size_t aux_size () { return aux_names().size(); }

  static std::vector<std::string> state_names() {
    return  std::vector<std::string>({
      "height",
      "mortality",
      "fecundity"
      });
  }

  std::vector<std::string> aux_names() {
    return std::vector<std::string>({"competition_effect", "height_inverse"});
  }

  void compute_rates(const K93_Environment& environment, Internals<double>& vars);

  void refresh_indices();

  double establishment_probability(const K93_Environment& environment);
  // The rates-carrying entry point Individual calls for every strategy. K93's
  // establishment reads no carbon budget, so the rates are unused here.
  double establishment_probability(const K93_Environment& environment,
                                   const Internals<double>&) {
    return establishment_probability(environment);
  }
  double net_mass_production_dt(const K93_Environment& environment,
                                double size, double cumulative_basal_area);
  double net_mass_production_dt(const K93_Environment& environment,
                                double size, double cumulative_basal_area,
                                double height_inverse);
  // Strategy-agnostic entry point used by Individual<K93> (#266). K93 has no
  // carbon budget, so the worker ignores these arguments and returns NA; the
  // wrapper exists to keep Individual's interface uniform across strategies.
  double net_mass_production_dt(const K93_Environment& environment,
                                const Internals<double>& vars) {
    return net_mass_production_dt(environment, vars.state(HEIGHT_INDEX),
                                  vars.aux(COMPETITION_EFFECT_AUX_INDEX),
                                  vars.aux(HEIGHT_INVERSE_AUX_INDEX));
  }

  double compute_competition(double z, double size) const;
  // Hot path (called per node from Species::compute_competition): defined inline
  // here so the no-LTO build folds it into the loop instead of paying a cross-TU
  // call (mirrors FF16_Strategy; see the profile-plant skill).
  double compute_competition(double z, double whole_plant_competition,
                             double height_inverse) const {
    return compute_competition_by_ratio(z * height_inverse, whole_plant_competition);
  }
  double compute_competition_by_ratio(double z_over_size,
                                      double whole_plant_competition) const {
    // Competition only felt if plant bigger than target size z.
    return whole_plant_competition * canopy_shape.leaf_area_above(z_over_size);
  }
  // Strategy-agnostic entry point used by Individual<K93> (#266): reads the
  // cached competition_effect and height_inverse aux slots itself.
  double compute_competition(double z, const Internals<double>& vars) const {
    return compute_competition(z, vars.aux(COMPETITION_EFFECT_AUX_INDEX),
                               vars.aux(HEIGHT_INVERSE_AUX_INDEX));
  }

  // The competition contribution and its vertical derivative from one pass, so
  // u^eta is evaluated once. The first entry is bit-for-bit the one
  // compute_competition() returns: both read the shading model's own profile,
  // which is what keeps them equal under a flat-top one.
  with_slope<double>
  compute_competition_and_slope(double z, const Internals<double>& vars) const {
    const double whole_plant_competition = vars.aux(COMPETITION_EFFECT_AUX_INDEX);
    const double height_inverse = vars.aux(HEIGHT_INVERSE_AUX_INDEX);
    const std::pair<double, double> Qq =
      canopy_shape.Q_and_q(z * height_inverse, z, height_inverse);
    return {whole_plant_competition * Qq.first,
            -(whole_plant_competition * Qq.second)};
  }

  void update_dependent_aux(const int index, Internals<double>& vars);


  // K93 Methods  ----------------------------------------------
  double size_to_basal_area(double size) const;

  // Growth rate of a plant per unit time:
  double size_dt(double size, double cumulative_basal_area) const;

  // Rate of offspring production
  double fecundity_dt(double size, double cumulative_basal_area) const;

  // Rate of mortality over time
  double mortality_dt(double cumulative_basal_area,
                      double cumulative_mortality) const;

  // Set constants within K93_Strategy
  void prepare_strategy();

  // Birth height (dbh, cm) of a seedling. Strategy-agnostic accessor used by
  // the templated Individual; for K93 this is the user-set pars.height_0.
  double initial_height() const { return pars.height_0; }

  // Biological (user-settable) parameters; see K93_Pars above.
  K93_Pars pars;

  // Derived / precomputed in prepare_strategy() (NOT user-set)
  CanopyShape<double> canopy_shape;
};

// Definitions ----------------------------------------------------------
// Header-only, as FF16_Strategy is: the model is read in one place and the
// per-node calls can inline into the SCM loop.

inline K93_Strategy::K93_Strategy() {
   // Empirical parameter defaults (Table 1) live in the K93_Pars struct
   // member initialisers.

   // build the string state/aux name to index map
   refresh_indices();
   name = "K93";
}

// Signatures fixed in plant.h
inline void K93_Strategy::update_dependent_aux(const int index, Internals<double>& vars) {
  if (index == HEIGHT_INDEX) {
    double height = vars.state(HEIGHT_INDEX);
    vars.set_aux(COMPETITION_EFFECT_AUX_INDEX,
                 compute_competition_by_ratio(
                   0.0, pars.k_I * size_to_basal_area(height)));
    vars.set_aux(HEIGHT_INVERSE_AUX_INDEX, 1.0 / height);
  }
}

inline double K93_Strategy::compute_competition(double z, double size) const {
  return compute_competition(z, pars.k_I * size_to_basal_area(size), 1.0 / size);
}

inline double K93_Strategy::establishment_probability(const K93_Environment& environment){
  (void) environment;
  // K93 (Kohyama 1993) has no carbon-budget establishment filter:
  // establishment is deterministic, so every dispersed seed establishes.
  // Demographic filtering instead happens via growth suppression (size_dt is
  // clamped to >= 0 under competition) and the g > 0 density guard in
  // Node::compute_initial_conditions. This value is consumed by the SCM
  // (initial mortality / density) and the stochastic germination test.
  // TODO(#480): may want to make this dependent on achieving positive growth rate.
  return 1.0;
}

inline double K93_Strategy::net_mass_production_dt(const K93_Environment& environment,
                                            double height, double area_leaf_) {
  (void) environment;
  (void) height;
  (void) area_leaf_;
  // K93 models growth directly (size_dt); it has no carbon mass-production
  // budget, so net mass production is undefined for this strategy. The value
  // is never consumed by compute_rates -- it is reachable only via the R
  // diagnostic accessor and resource_compensation_point(). Return NA to make
  // "not applicable" explicit rather than returning a misleading number.
  return NA_REAL;
}

inline double K93_Strategy::net_mass_production_dt(const K93_Environment& environment,
                                            double height, double area_leaf_,
                                            double height_inverse) {
  (void) height_inverse;
  return net_mass_production_dt(environment, height, area_leaf_);
}

inline void K93_Strategy::refresh_indices () {
  // Create and fill the name to state index maps
  state_index = std::map<std::string,int>();
  aux_index   = std::map<std::string,int>();
  std::vector<std::string> aux_names_vec = aux_names();
  std::vector<std::string> state_names_vec = state_names();
  for (size_t i = 0; i < state_names_vec.size(); i++) {
    state_index[state_names_vec[i]] = i;
  }
  for (size_t i = 0; i < aux_names_vec.size(); i++) {
    aux_index[aux_names_vec[i]] = i;
  }
  check_state_layout(state_index, state_size(), "K93");
}

// i.e. setting rates of ode vars from the state and updating aux vars
inline void K93_Strategy::compute_rates(const K93_Environment& environment, Internals<double>& vars) {

  double height = vars.state(HEIGHT_INDEX);

  // suppression integral mapped [0, 1] using adaptive spline
  // back transform to basal area and add suppression from self
  double competition = environment.get_environment_at_height(height);

  double cumulative_basal_area = -log(competition) / pars.k_I;

  if (!util::is_finite(cumulative_basal_area))
  {
    util::stop("Environmental interpolator has gone out of bounds, try lowering the extinction coefficient k_I");
  }

  vars.set_rate(HEIGHT_INDEX,
    size_dt(height, cumulative_basal_area));

  vars.set_rate(FECUNDITY_INDEX,
    fecundity_dt(height, cumulative_basal_area));

  vars.set_rate(MORTALITY_INDEX,
    mortality_dt(cumulative_basal_area, vars.state(MORTALITY_INDEX)));
}

inline double K93_Strategy::size_to_basal_area(double size) const {
  return M_PI / 4 * pow(size, 2);
}

// [eqn 10] Growth
inline double K93_Strategy::size_dt(double size,
                             double cumulative_basal_area) const {

  double growth = size * (pars.b_0 - pars.b_1 * log(size) - pars.b_2 * cumulative_basal_area);

  if(growth < 0.0) {
    growth = 0.0;
  }

  return growth;
}

// [eqn 12] Reproduction
inline double K93_Strategy::fecundity_dt(double size,
                                  double cumulative_basal_area) const {
  double basal_area = size_to_basal_area(size);
  return pars.d_0 * basal_area * exp(-pars.d_1 * cumulative_basal_area);
}

// [eqn 11] Mortality
inline double K93_Strategy::mortality_dt(double cumulative_basal_area,
                                  double cumulative_mortality) const {
  // Mortality probability is 1 past either test, so the rate calculations have
  // nothing left to describe and the state does not move again.
  //
  // ⚠️ THE CEILING IS TESTED BESIDE FINITENESS, and the pair must not be reduced
  // back to one test. compute_initial_conditions holds an unestablished
  // recruit's hazard at establishment_failure_hazard rather than at the +Inf
  // -log(0) gives, so this is the test that keeps such a cohort's rate parked.
  if (util::is_finite(cumulative_mortality) &&
      cumulative_mortality < establishment_failure_hazard) {
    double mu = -pars.c_0 + pars.c_1 * cumulative_basal_area;
    return (mu > 0)? mu:0.0;
 } else {
    return 0.0;
  }
}

// useful for pre-computing expensive objects
inline void K93_Strategy::prepare_strategy() {
  canopy_shape.initialise(pars.eta);

  if (is_variable_birth_rate) {
    extrinsic_drivers.set_variable("birth_rate", birth_rate_x, birth_rate_y,
                                   odelia::drivers::Slopes::monotone);
  } else {
    extrinsic_drivers.set_constant("birth_rate", birth_rate_y[0]);
  }
}

inline K93_Strategy::ptr make_strategy_ptr(K93_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<K93_Strategy>(s);
}

}

#endif
