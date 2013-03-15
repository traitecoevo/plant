// -*-c++-*-
#ifndef TREE_PLANT_
#define TREE_PLANT_

#include <Rcpp.h>

#include "strategy.h"
#include "spline.h"

namespace model {
    
class Plant {
public:
  Plant(Strategy  s);
  Plant(Strategy *s);

  // Copy constructor, assigment and destructor (rule of three)
  Plant(const Plant &other);
  Plant& operator=(const Plant &rhs);
  ~Plant();

  double get_mass_leaf() const;
  void set_mass_leaf(double x);

  // static const int size = 3;

  Rcpp::NumericVector r_get_vars_size() const;
  Rcpp::NumericVector r_get_vars_phys() const;
  double r_compute_assimilation(spline::Spline env);
  double r_compute_assimilation_x(double x, spline::Spline env);
  Rcpp::List r_get_parameters() const;

  // * Competitive environment
  // [eqn  9] Probability density of leaf area at height `z`
  double q(double z) const;
  // [eqn 10] Fraction of leaf area above height `z`
  double Q(double z) const;
  // [      ] Inverse of Q: height above which fraction 'x' of leaf found
  double Qp(double x) const;
  // [      ] Leaf area (not fraction) above height `z`
  double leaf_area_above(double z) const;

  // * Mass production
  // [Appendix S6] Per-leaf photosynthetic rate.
  double assimilation_leaf(double x) const;

private:
  // * Individual size
  // [eqn 1-8] Update size variables to a new leaf mass.
  void compute_vars_size(double mass_leaf_);

  // * Mass production
  // [eqn 12] Gross annual CO2 assimilation
  double compute_assimilation(spline::Spline *env);
  double compute_assimilation_x(double x, spline::Spline *env);
  // [eqn 13] Total maintenance respiration
  double compute_respiration() const;
  // [eqn 14] Total turnover
  double compute_turnover() const;
  // [eqn 16] Fraction of whole plan growth that is leaf
  double compute_reproduction_fraction() const;
  // [eqn 18] Fraction of mass growth that is leaves
  double compute_leaf_fraction() const;

  // [eqn 12-19,21] Update physiological variables
  void compute_vars_phys(spline::Spline *env);

  bool standalone;
  Strategy *strategy;

  // * Individual size
  // Mass of leaves.  This is the core independent variable
  double mass_leaf;      // [eqn 1]
  // Other size variables that follow directly from `mass_leaf`:
  double leaf_area;      // [eqn 2]
  double height;         // [eqn 3]
  double mass_sapwood;   // [eqn 4]
  double mass_bark;      // [eqn 5]
  double mass_heartwood; // [eqn 6]
  double mass_root;      // [eqn 7] (fine roots)
  double mass_total;     // [eqn 8]

  // * Mass production
  double assimilation;   // [eqn 12] Gross annual CO2 assimilation
  double respiration;    // [eqn 13] Total maintenance respiration
  double turnover;       // [eqn 14] Total turnover
  double net_production; // [eqn 15] Net production
  double reproduction_fraction; // [eqn 16]
  double fecundity_rate; // [eqn 17] Rate of offspring production
  double leaf_fraction;  // [eqn 18] Fraction of mass growth that is leaves
  double growth_rate;    // [eqn 19] Growth rate in leaf mass

  // * Mortality
  double mortality_rate; // [eqn 21]

  // State variables resulting from integration of the corresponding
  // *_rate variables (fecundity_rate and mortality_rate).
  double fecundity, mortality;
};

}

#endif
