// -*-c++-*-
#ifndef TREE_PLANT_H_
#define TREE_PLANT_H_

#include <Rcpp.h>

#include "ode_target.h"
#include "integrator.h"
#include "strategy.h"
#include "spline.h"
#include "functor.h"

namespace model {

class PlantBase : public ode::OdeTarget {
public:
  virtual ~PlantBase() {};
  // * "Key" variables.
  virtual double mass_leaf() const = 0;
  virtual void set_mass_leaf(double mass_leaf_) = 0;
  virtual double mass_leaf_rate() const = 0;
  virtual double mortality() const = 0;
  virtual void set_mortality(double x) = 0;
  virtual double mortality_rate() const = 0;
  virtual double fecundity() const = 0;
  virtual void set_fecundity(double x) = 0;
  virtual double fecundity_rate() const = 0;
  // * The rest
  virtual double get_height() const = 0;
  virtual double mortality_probability() const = 0;
  virtual double survival_probability() const = 0;
  virtual double leaf_area_above(double z) const = 0;
  virtual void compute_vars_phys(spline::Spline *env) = 0;
  virtual double germination_probability(spline::Spline *env) = 0;
  virtual int offspring() = 0;
  virtual bool died() = 0;
  // * R interface
  virtual Strategy r_get_strategy() const = 0;
  virtual Rcpp::NumericVector r_get_vars_size() const = 0;
  virtual Rcpp::NumericVector r_get_vars_phys() const = 0;
  virtual void r_compute_vars_phys(spline::Spline env) = 0;
  virtual double r_germination_probability(spline::Spline env) = 0;
  virtual bool r_died() = 0;
};

class Plant : public PlantBase {
public:
  Plant(Strategy  s);
  Plant(Strategy *s);

  // Equivalence operator
  bool operator==(const Plant &rhs);

  // * Individual size
  // [eqn 1-8] Update size variables to a new leaf mass.
  void set_mass_leaf(double mass_leaf_);

  // Accessors to the key variables / rates.
  double mass_leaf() const;
  // void set_mass_leaf(); -- see above.
  double mass_leaf_rate() const; // this is "growth rate"

  double mortality() const;
  void set_mortality(double x);
  double mortality_rate() const;

  double fecundity() const;
  void set_fecundity(double x);
  double fecundity_rate() const;

  // These are derived from mortality() -- see design.md.
  double mortality_probability() const;
  double survival_probability() const;

  // Also need access to height.
  double get_height() const;

  // * Competitive environment
  // [      ] Leaf area (not fraction) above height `z`
  double leaf_area_above(double z) const;

  // * Mass production
  // [eqn 12-19,21] Update physiological variables
  void compute_vars_phys(spline::Spline *env);

  // * Births and deaths
  int offspring();
  bool died();
  // [eqn 20] Survival of seedlings during germination
  double germination_probability(spline::Spline *env);

  // * Access the Control parameter.
  const Control& control() const;

  // * ODE interface
  size_t ode_size() const;
  ode::iter_const ode_values_set(ode::iter_const it);
  ode::iter       ode_values(ode::iter it) const;
  ode::iter       ode_rates(ode::iter it)  const;

  // * Set constants within Strategy
  static void prepare_strategy(Strategy *s);

  // * R interface
  Strategy r_get_strategy() const;
  Rcpp::NumericVector r_get_vars_size() const;
  Rcpp::NumericVector r_get_vars_phys() const;
  void r_compute_vars_phys(spline::Spline env);
  double r_germination_probability(spline::Spline env);
  bool r_died();

private:
  // * Individual size
  // [eqn 1-8] Update size variables to a new leaf mass.
  void compute_vars_size(double mass_leaf_);

  // * Competitive environment
  // [eqn  9] Probability density of leaf area at height `z`
  double q(double z) const;
  // [eqn 10] Fraction of leaf area above height `z`
  double Q(double z) const;
  // [      ] Inverse of Q: height above which fraction 'x' of leaf found
  double Qp(double x) const;

  // * Mass production
  // [eqn 12] Gross annual CO2 assimilation
  double compute_assimilation(spline::Spline *env) const;
  // Used internally, corresponding to the inner term in [eqn 12]
  double compute_assimilation_x(double x, spline::Spline *env) const;
  // [Appendix S6] Per-leaf photosynthetic rate.
  double assimilation_leaf(double x) const;
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

  // To simplify my life, I'm making a small internal-only class that
  // contains some implementation details here.  This is going to
  // simplify the assignment operator / copy constructor etc.

  class internals {
  public:
    internals();
    bool operator==(const Plant::internals &rhs);
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
    // * Variables
    double mortality;
    double fecundity;
  };

  Strategy::ptr strategy;
  internals vars;

  static const int ode_dimension = 3;
};

namespace test {
bool test_plant(Strategy s, bool copy, bool ptr);
}


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
