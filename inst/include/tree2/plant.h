// -*-c++-*-
#ifndef TREE_PLANT_H_
#define TREE_PLANT_H_

#include <memory> // std::shared_ptr
#include <vector>
#include <tree2/strategy.h>
#include <tree2/qag_internals.h> // quadrature::intervals_type
#include <tree2/environment.h>
#include <tree2/ode_interface.h>
#include <RcppCommon.h>

namespace tree2 {

class Plant {
public:
  typedef Strategy                       strategy_type;
  typedef std::shared_ptr<strategy_type> strategy_ptr_type;
  Plant(strategy_ptr_type s);

  // * Individual size
  // [eqn 1-8] Update size variables to a new leaf mass.

  // * Individual size
  double height() const;
  double height_rate() const;
  void set_height(double height_);

  double mortality() const;
  double mortality_rate() const;
  void set_mortality(double x);

  double fecundity() const;
  double fecundity_rate() const;
  void set_fecundity(double x);

  double heartwood_area() const;
  double heartwood_area_rate() const;
  void set_heartwood_area(double x);

  double heartwood_mass() const;
  double heartwood_mass_rate() const;
  void set_heartwood_mass(double x);

  // Mortality functions
  // Declared as static method so it can be accessed and queried externally
  static double mortality_growth_independent(double d0);
  static double mortality_growth_dependent(double d2, double d3,
                                           double productivity);

  // These are derived from mortality() -- see design.md.
  double mortality_probability() const;
  double survival_probability() const;

  // * Competitive environment
  double leaf_area() const;
  // [      ] Leaf area (not fraction) above height `z`
  double leaf_area_above(double z) const;

  // * Mass production
  // [eqn 12-19,21] Update physiological variables
  void compute_vars_phys(const Environment& environment, bool
			 reuse_intervals=false);

  // * Births and deaths
  // [eqn 20] Survival of seedlings during germination
  double germination_probability() const;
  double germination_probability(const Environment& environment);

  // * ODE interface
  static size_t ode_size() {return 5;}
  ode::const_iterator set_ode_state(ode::const_iterator it);
  ode::iterator       ode_state(ode::iterator it) const;
  ode::iterator       ode_rates(ode::iterator it) const;
  // Optional, but useful
  static std::vector<std::string> ode_names();

  // * Set constants within Strategy
  static void prepare_strategy(strategy_ptr_type s);
  static double height_seed(strategy_ptr_type s);

  // * R interface
  strategy_type r_get_strategy() const;
  SEXP r_get_vars_size() const;
  SEXP r_get_vars_phys() const;
  SEXP r_get_vars_growth() const;

  // * Used by tools:
  double net_production() const {return vars.net_production;}

  Plant r_copy() const;

  const Control& control() const;

private:
  // * Individual size
  // [eqn 1-8] Update size variables to a new leaf mass.
  void compute_vars_size(double height_);

  // * Competitive environment
  // [eqn  9] Probability density of leaf area at height `z`
  double q(double z) const;
  // [eqn 10] Fraction of leaf area above height `z`
  double Q(double z) const;
  // [      ] Inverse of Q: height above which fraction 'x' of leaf found
  double Qp(double x) const;

  // * Mass production
  // [eqn 12] Gross annual CO2 assimilation
  double assimilation(const Environment& environment, bool reuse_intervals);
  // Used internally, corresponding to the inner term in [eqn 12]
  double compute_assimilation_x(double x, const Environment& environment) const;
  double compute_assimilation_h(double h, const Environment& environment) const;
  double compute_assimilation_p(double p, const Environment& environment) const;
  // [Appendix S6] Per-leaf photosynthetic rate.
  double assimilation_leaf(double x) const;

  // TODO: Move into Strategy by computing this statelessly?  Requires
  // moving the size functions over (@dfalster).
  double live_mass_given_height(double h);

  // To simplify my life, I'm making a small internal-only class that
  // contains some implementation details here.

  // TODO: I'll organise exporting this as a RcppR6 list class I
  // think...
  class internals {
  public:
    internals();
    // * Individual size
    // Mass of leaves.  This is the core independent variable
    double leaf_mass;      // [eqn 1]
    // Other size variables that follow directly from `leaf_mass`:
    double leaf_area;      // [eqn 2]
    double height;         // [eqn 3]
    double sapwood_mass;   // [eqn 4]
    double bark_mass;      // [eqn 5]
    double heartwood_mass; // [eqn 6]
    double heartwood_area;
    double root_mass;      // [eqn 7] (fine roots)
    double live_mass;      // [eqn 8]
    // * Mass production
    double assimilation;   // [eqn 12] Gross annual CO2 assimilation
    double respiration;    // [eqn 13] Total maintenance respiration
    double turnover;       // [eqn 14] Total turnover
    double net_production; // [eqn 15] Net production
    double reproduction_fraction; // [eqn 16]
    double fecundity_rate; // [eqn 17] Rate of offspring production
    double leaf_fraction;  // [eqn 18] Fraction of mass growth that is leaves
    double leaf_mass_growth_rate; // [eqn 19] Growth rate in leaf mass
    double height_growth_rate;    // [doc/details.md]
    double heartwood_area_rate;
    double heartwood_mass_rate;
    // * Mortality
    double mortality_rate; // [eqn 21]
    // * Variables
    double mortality;
    double fecundity;
  };

  strategy_ptr_type strategy;
  internals vars;

  // The ode dimensions are:
  // 1. Height
  // 2. Mortality
  // 3. Fecundity
  // 4. Heartwood area
  // 5. Heartwood mass
  static const int ode_dimension = 5;
};

Plant::strategy_ptr_type make_strategy_ptr(Plant::strategy_type s);
Plant make_plant(Plant::strategy_type s);

}

#endif
