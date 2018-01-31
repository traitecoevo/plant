// -*-c++-*-
#ifndef PLANT_PLANT_PLANT_PLUS_H_
#define PLANT_PLANT_PLANT_PLUS_H_

#include <vector>
#include <plant/ff16_strategy.h>
#include <plant/plant_plus_internals.h>
#include <plant/ode_interface.h>
#include <plant/plant.h>
#include <plant/environment.h>

namespace plant {

template <typename T>
class PlantPlus {
public:
  typedef PlantPlus_internals internals;
  typedef T  strategy_type;
  typedef typename strategy_type::ptr strategy_type_ptr;
  PlantPlus(strategy_type_ptr s)
    : strategy(s) {
    set_height(strategy->height_0);
  }

  PlantPlus(Plant<T> p)
    : strategy(make_strategy_ptr(p.r_get_strategy())) {
    std::vector<double> state(ode_size());
    p.ode_state(state.begin());
    set_ode_state(state.begin());
  }

  PlantPlus(Plant<T> p, const Environment& environment)
    : PlantPlus(p) {
    compute_vars_phys(environment);
    compute_vars_growth();
  }

  Plant<T> to_plant() const {
    Plant<T> p(strategy);
    std::vector<double> state(ode_size());
    ode_state(state.begin());
    p.set_ode_state(state.begin());
    return p;
  }

  // * Individual size
  // [eqn 1-8] Update size variables to a new leaf mass.

  // * Individual size
  double height() const {return vars.height;}
  double height_dt() const {return vars.height_dt;}
  void set_height(double x) {
    // TODO: in the original version of the model, all of plant size
    // was driven only by height, so we only needed to check this
    // here.  But now we have two size variables (heartwood mass being
    // the other) we need to be more careful. For the new plant type,
    // replace this height function with a new function that takes
    // both size variables as arguments and update all the size
    // variables at that point.
    vars.height = x;
    compute_vars_size();
  }

  double mortality() const {return vars.mortality;}
  double mortality_dt() const {return vars.mortality_dt;}
  void set_mortality(double x) {vars.mortality = x;}

  double fecundity() const {return vars.fecundity;}
  double fecundity_dt() const {return vars.fecundity_dt;}
  void set_fecundity(double x) {vars.fecundity = x;}

  double area_heartwood() const {return vars.area_heartwood;}
  double area_heartwood_dt() const {return vars.area_heartwood_dt;}
  void set_area_heartwood(double x) {
    vars.area_heartwood = x;
    compute_vars_size();
  }

  double mass_heartwood() const {return vars.mass_heartwood;}
  double mass_heartwood_dt() const {return vars.mass_heartwood_dt;}
  void set_mass_heartwood(double x) {
    vars.mass_heartwood = x;
    compute_vars_size();
  }

  // * Competitive environment
  double area_leaf() const {return vars.area_leaf;}

  // [      ] Leaf area (not fraction) above height `z`
  double area_leaf_above(double z) const {
    return strategy->area_leaf_above(z, vars.height, vars.area_leaf);
  }

  // * Mass production
  // [eqn 12-19,21] Update physiological variables
  void compute_vars_phys(const Environment& environment, bool
                         reuse_intervals=false);
  // TODO: This is temporary -- it should be called by
  // compute_vars_phys, but I don't want that to always happen, so do
  // it manually for now.
  void compute_vars_growth();

  // * Births and deaths
  // [eqn 20] Survival of seedlings during germination
  double germination_probability(const Environment& environment) {
    return strategy->germination_probability(environment);
  }

  // * ODE interface
  //   The ode dimensions are:
  //   1. Height
  //   2. Mortality
  //   3. Fecundity
  //   4. Heartwood area
  //   5. Heartwood mass
  static size_t       ode_size() {return 5;}
  ode::const_iterator set_ode_state(ode::const_iterator it) {
    set_height(*it++);
    set_mortality(*it++);
    set_fecundity(*it++);
    set_area_heartwood(*it++);
    set_mass_heartwood(*it++);
    return it;
  }
  ode::iterator ode_state(ode::iterator it) const {
    *it++ = height();
    *it++ = mortality();
    *it++ = fecundity();
    *it++ = area_heartwood();
    *it++ = mass_heartwood();
    return it;
  }
  ode::iterator ode_rates(ode::iterator it) const {
    *it++ = height_dt();
    *it++ = mortality_dt();
    *it++ = fecundity_dt();
    *it++ = area_heartwood_dt();
    *it++ = mass_heartwood_dt();
    return it;
  }
  // Optional, but useful
  static std::vector<std::string> ode_names() {
    return std::vector<std::string>({"height", "mortality", "fecundity",
          "area_heartwood", "mass_heartwood"});
  }

  // * R interface
  strategy_type r_get_strategy() const {return *strategy.get();}
  internals r_internals() const {return vars;}
  const Control& control() const {return strategy->control;}

  // * Used by tools:
  double net_mass_production_dt() const {return vars.net_mass_production_dt;}

  std::string strategy_name() const {return strategy->name;}

private:
  // * Individual size
  // [eqn 1-8] Update size variables to a new leaf mass.
  void compute_vars_size();

  strategy_type_ptr strategy;
  internals vars;
};

// * Individual size
// [eqn 1-8] Update size variables to a new leaf mass.
template <typename T>
void PlantPlus<T>::compute_vars_size() {
  vars.area_leaf = strategy->area_leaf(vars.height);
  vars.mass_leaf = strategy->mass_leaf(vars.area_leaf);
  vars.area_sapwood =  strategy->area_sapwood(vars.area_leaf);
  vars.mass_sapwood =  strategy->mass_sapwood(vars.area_sapwood, vars.height);
  vars.area_bark =  strategy->area_bark(vars.area_leaf);
  vars.mass_bark = strategy->mass_bark(vars.area_bark, vars.height);
  vars.mass_root = strategy->mass_root(vars.area_leaf);
  vars.mass_live = strategy->mass_live(vars.mass_leaf, vars.mass_sapwood,
                                       vars.mass_bark , vars.mass_root);
  vars.area_stem = strategy->area_stem(vars.area_bark, vars.area_sapwood,
                                         vars.area_heartwood);
  vars.mass_total = strategy->mass_total(vars.mass_leaf, vars.mass_bark, vars.mass_sapwood,
                                         vars.mass_heartwood, vars.mass_root);
  vars.mass_above_ground = strategy->mass_above_ground(vars.mass_leaf, vars.mass_bark,
                                                       vars.mass_sapwood, vars.mass_root);
  vars.diameter_stem = strategy->diameter_stem(vars.area_stem);
}

// * Mass production
// [eqn 12-19,21] Update physiological variables given the current
// light environment (and given the current set of size variables).
template <typename T>
void PlantPlus<T>::compute_vars_phys(const Environment& environment,
                                     bool reuse_intervals) {
  // [eqn 12] Gross annual CO2 assimilation
  vars.assimilation = strategy->assimilation(environment, vars.height,
                                             vars.area_leaf, reuse_intervals);

  // [eqn 13] Total maintenance respiration
  vars.respiration = strategy->respiration(vars.mass_leaf,
                                           vars.mass_sapwood,
                                           vars.mass_bark,
                                           vars.mass_root);

  // [eqn 14] Total turnover
  vars.turnover = strategy->turnover(vars.mass_leaf, vars.mass_sapwood,
                                     vars.mass_bark, vars.mass_root);

  // [eqn 15] Net production
  vars.net_mass_production_dt =
    strategy->net_mass_production_dt_A(vars.assimilation,
                                       vars.respiration,
                                       vars.turnover);

  if (vars.net_mass_production_dt > 0) {
    // [eqn 16] - Fraction of whole plant growth that is leaf
    vars.fraction_allocation_reproduction =
      strategy->fraction_allocation_reproduction(vars.height);

    // [eqn 17] - Rate of offspring production
    //
    // NOTE: In prototype from Falster 2010, was multiplied by S_D (survival during
    // dispersal), but we do not here.
    // NOTE: This is also a hyper-parametrisation and should move into
    // the initialisation function.
    vars.fecundity_dt =
      strategy->fecundity_dt(vars.net_mass_production_dt,
                              vars.fraction_allocation_reproduction);

    // [eqn 19] - Growth rate in leaf height
    // different to Falster 2010, which was growth rate in leaf mass
    vars.darea_leaf_dmass_live =
      strategy->darea_leaf_dmass_live(vars.area_leaf);
    vars.fraction_allocation_growth = strategy->fraction_allocation_growth(vars.height);

    vars.area_leaf_dt =
      vars.net_mass_production_dt * vars.fraction_allocation_growth *
      vars.darea_leaf_dmass_live;
   vars.height_dt =
      strategy->dheight_darea_leaf(vars.area_leaf) *
      vars.area_leaf_dt;
   vars.area_heartwood_dt =
     strategy->area_heartwood_dt(vars.area_leaf);
   vars.mass_heartwood_dt =
      strategy->mass_heartwood_dt(vars.mass_sapwood);
  } else {
    vars.fraction_allocation_reproduction = 0.0;
    vars.fraction_allocation_growth       = 0.0;
    vars.fecundity_dt        = 0.0;
    vars.area_leaf_dt = 0.0;
    vars.darea_leaf_dmass_live = 0.0;
    vars.height_dt    = 0.0;
    vars.area_heartwood_dt   = 0.0;
    vars.mass_heartwood_dt   = 0.0;
  }

  // [eqn 21] - Instantaneous mortality rate
  vars.mortality_dt =
      strategy->mortality_dt(vars.net_mass_production_dt / vars.area_leaf,
                             vars.mortality);
}

template <typename T>
void PlantPlus<T>::compute_vars_growth() {
  const strategy_type *s = strategy.get(); // for brevity.
  const double area_leaf = vars.area_leaf,
    area_leaf_dt = vars.area_leaf_dt;

  const double area_stem = vars.area_bark +
                           vars.area_sapwood +
                           vars.area_heartwood;

  // Changes with leaf area:
  vars.dheight_darea_leaf       = s->dheight_darea_leaf(area_leaf);
  vars.dmass_sapwood_darea_leaf = s->dmass_sapwood_darea_leaf(area_leaf);
  vars.dmass_bark_darea_leaf    = s->dmass_bark_darea_leaf(area_leaf);
  vars.dmass_root_darea_leaf    = s->dmass_root_darea_leaf(area_leaf);

  // Changes over time:
  vars.area_sapwood_dt    = s->area_sapwood_dt(area_leaf_dt);
  vars.area_bark_dt       = s->area_bark_dt(area_leaf_dt);
  vars.area_stem_dt       = s->area_stem_dt(area_leaf,
                                              area_leaf_dt);
  vars.diameter_stem_dt   = s->diameter_stem_dt(area_stem,
                                                  vars.area_stem_dt);
  vars.mass_root_dt       = s->mass_root_dt(area_leaf,
                                             area_leaf_dt);
  vars.mass_live_dt       = s->mass_live_dt(vars.fraction_allocation_reproduction,
                                              vars.net_mass_production_dt);
  vars.mass_total_dt      = s->mass_total_dt(vars.fraction_allocation_reproduction,
                                              vars.net_mass_production_dt,
                                              vars.mass_heartwood_dt);
  vars.mass_above_ground_dt =
    s->mass_above_ground_dt(area_leaf,
                             vars.fraction_allocation_reproduction,
                             vars.net_mass_production_dt,
                             vars.mass_heartwood_dt,
                             area_leaf_dt);

  // Odd one out:
  vars.ddiameter_stem_darea_stem = s->ddiameter_stem_darea_stem(area_stem);
}

template <typename T>
PlantPlus<T> make_plant_plus(T s) {
  return PlantPlus<T>(make_strategy_ptr(s));
}

}

#endif
