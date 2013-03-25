// -*-c++-*-
#ifndef TREE_PLANT_
#define TREE_PLANT_

#include <Rcpp.h>

#include "ode_solver.h"
#include "integrator.h"
#include "strategy.h"
#include "spline.h"
#include "functor.h"

namespace model {
    
class Plant {
public:
  Plant(Strategy  s);
  Plant(Strategy *s);

  // Copy constructor, assigment and destructor (rule of three)
  Plant(const Plant &other);
  Plant& operator=(const Plant &rhs);
  ~Plant();

  static void prepare_strategy(Strategy *s);

  double get_mass_leaf() const;
  double get_height() const;

  // * Individual size
  // [eqn 1-8] Update size variables to a new leaf mass.
  void set_mass_leaf(double mass_leaf_);

  // * Mass production
  // [eqn 12-19,21] Update physiological variables
  void compute_vars_phys(spline::Spline *env);

  Rcpp::NumericVector r_get_vars_size() const;
  Rcpp::NumericVector r_get_vars_phys() const;
  double r_compute_assimilation(spline::Spline env) const;
  double r_compute_assimilation_x(double x, spline::Spline env) const;
  Rcpp::List r_get_parameters() const;
  void r_compute_vars_phys(spline::Spline env);

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

  // ODE interface
  static const size_t ode_size = 3;
  ode::iter_const set_values(ode::iter_const it, bool &changed);
  ode::iter       get_values(ode::iter it) const;
  ode::iter       get_rates(ode::iter it)  const;

private:
  // * Individual size
  // [eqn 1-8] Update size variables to a new leaf mass.
  void compute_vars_size(double mass_leaf_);

  // * Mass production
  // [eqn 12] Gross annual CO2 assimilation
  double compute_assimilation(spline::Spline *env) const;
  // Used internally, corresponding to the inner term in [eqn 12]
  double compute_assimilation_x(double x, spline::Spline *env) const;
  // [eqn 13] Total maintenance respiration
  double compute_respiration() const;
  // [eqn 14] Total turnover
  double compute_turnover() const;
  // [eqn 16] Fraction of whole plan growth that is leaf
  double compute_reproduction_fraction() const;
  // [eqn 18] Fraction of mass growth that is leaves
  double compute_leaf_fraction() const;

  // Update a number of constants within the model.  This is a work in
  // progress.
  static double mass_leaf_seed(Strategy *s);
  double compute_mass_total(double m_);

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

// To prepare for the integration in `compute_assimilation` we need to
// convert the function `compute_assimilation_x(double, util::Spline*)
// to take just a double as an argument.  Boost has the ability to
// bind arguments which would be nice here, but we're avoiding
// depending on that for the time being.
// 
// This binds the second argument, assuming a const method.  All a bit
// of a hack, but it does seem to work correctly.
template <class T, class T2, double (T::*target)(double, T2) const>
class FunctorBind1 : public util::DFunctor {
public:
  FunctorBind1(const T *obj, T2 arg2) : obj(obj), arg2(arg2) {}
  virtual double operator()(double x) {
    return (obj->*target)(x, arg2);
  }
  
private:
  const T* obj;
  T2 arg2;
};

}

RCPP_EXPOSED_CLASS(model::Plant)

#endif
