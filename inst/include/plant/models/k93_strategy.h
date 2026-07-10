// Built from  inst/include/plant/models/ff16_strategy.h on Fri Jul 24 10:23:19 2020 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_K93_STRATEGY_H_
#define PLANT_PLANT_K93_STRATEGY_H_

#include <plant/strategy.h>
#include <plant/models/k93_environment.h>
#include <plant/canopy_shape.h>

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
  typedef std::shared_ptr<K93_Strategy> ptr;
  K93_Strategy();

  // Scientific version. Bump ONLY when equations or default parameters change
  // the simulation output for identical inputs. Do NOT bump for refactors,
  // performance, interface, or serialisation changes. Bumping invalidates
  // logpile's cache for this model (see plant::model_version() / model_id()).
  static constexpr int scientific_version = 1;

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

  void compute_rates(const K93_Environment& environment, Internals& vars);

  void refresh_indices();

  double establishment_probability(const K93_Environment& environment);
  double net_mass_production_dt(const K93_Environment& environment,
                                double size, double cumulative_basal_area);
  double net_mass_production_dt(const K93_Environment& environment,
                                double size, double cumulative_basal_area,
                                double height_inverse);
  // Strategy-agnostic entry point used by Individual<K93> (#266). K93 has no
  // carbon budget, so the worker ignores these arguments and returns NA; the
  // wrapper exists to keep Individual's interface uniform across strategies.
  double net_mass_production_dt(const K93_Environment& environment,
                                const Internals& vars) {
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
    return whole_plant_competition * canopy_shape.Q(z_over_size);
  }
  // Strategy-agnostic entry point used by Individual<K93> (#266): reads the
  // cached competition_effect and height_inverse aux slots itself.
  double compute_competition(double z, const Internals& vars) const {
    return compute_competition(z, vars.aux(COMPETITION_EFFECT_AUX_INDEX),
                               vars.aux(HEIGHT_INVERSE_AUX_INDEX));
  }

  void update_dependent_aux(const int index, Internals& vars);


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
  CanopyShape canopy_shape;
};

K93_Strategy::ptr make_strategy_ptr(K93_Strategy s);

}

#endif
