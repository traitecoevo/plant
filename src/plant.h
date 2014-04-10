// -*-c++-*-
#ifndef TREE_PLANT_H_
#define TREE_PLANT_H_

#include <Rcpp.h>

#include "environment.h"
#include "ode_target.h"
#include "strategy.h"
#include "interpolator.h"
#include "functor.h"

namespace model {

class PlantBase : public ode::OdeTarget {
public:
  typedef std::vector<double> state;

  virtual ~PlantBase();
  // * "Key" variables.
  virtual double height() const = 0;
  virtual void set_height(double height_) = 0;
  virtual double height_rate() const = 0;
  virtual double mortality() const = 0;
  virtual void set_mortality(double x) = 0;
  virtual double mortality_rate() const = 0;
  virtual double fecundity() const = 0;
  virtual void set_fecundity(double x) = 0;
  virtual double fecundity_rate() const = 0;
  // * The rest
  virtual double mortality_probability() const = 0;
  virtual double survival_probability() const = 0;
  virtual double leaf_area() const = 0;
  virtual double leaf_area_above(double z) const = 0;
  virtual void compute_vars_phys(const Environment& environment) = 0;
  virtual double germination_probability(const Environment& environment) = 0;
  virtual int offspring() = 0;
  virtual bool died() = 0;
  virtual double
  assimilation_given_height(double h,
			    const Environment &environment) = 0;
  // * R interface
  virtual Strategy r_get_strategy() const = 0;
  virtual Rcpp::NumericVector r_get_vars_size() const = 0;
  virtual Rcpp::NumericVector r_get_vars_phys() const = 0;
  virtual bool r_died() = 0;
  // Going to need this in a few places, from the look of it.

  virtual size_t state_size() const = 0;
  virtual state::iterator get_state(state::iterator it) const = 0;
  virtual state::const_iterator set_state(state::const_iterator it) = 0;
};

class Plant : public PlantBase {
public:
  Plant(Strategy  s);
  Plant(Strategy *s);

  // Equivalence operator
  bool operator==(const Plant &rhs) const;

  // * Individual size
  // [eqn 1-8] Update size variables to a new leaf mass.

  // * Individual size
  double height() const;
  void set_height(double height_);
  double height_rate() const;

  double mortality() const;
  void set_mortality(double x);
  double mortality_rate() const;

  double fecundity() const;
  void set_fecundity(double x);
  double fecundity_rate() const;

  double heartwood_area() const;
  void set_area_heartwood(double x);
  double dheartwood_area_dt() const;

  double mass_heartwood() const;
  void set_mass_heartwood(double x);
  double sapwood_turnover() const;   // Sapwood turnover



  // These are derived from mortality() -- see design.md.
  double mortality_probability() const;
  double survival_probability() const;

  // * Competitive environment
  double leaf_area() const;
  // [      ] Leaf area (not fraction) above height `z`
  double leaf_area_above(double z) const;

  // * Mass production
  // [eqn 12-19,21] Update physiological variables
  void compute_vars_phys(const Environment& environment);

  // * Births and deaths
  int offspring();
  bool died();
  // [eqn 20] Survival of seedlings during germination
  double germination_probability(const Environment& environment);

  // * ODE interface
  size_t ode_size() const;
  ode::iterator_const set_ode_values(double time, ode::iterator_const it);
  ode::iterator       ode_values(ode::iterator it) const;
  ode::iterator       ode_rates(ode::iterator it)  const;

  // * Set constants within Strategy
  static void prepare_strategy(Strategy *s);
  static void compute_assimilation_fn(Strategy *s,
				      double hmin, double hmax,
				      const Environment &environment);
  static void rescale_assimilation_fn(Strategy *s,
				      double hmin, double hmax,
				      const Environment &environment);

  double assimilation_given_height(double h,
				   const Environment &environment);

  // * R interface
  Strategy r_get_strategy() const;
  Rcpp::NumericVector r_get_vars_size() const;
  Rcpp::NumericVector r_get_vars_phys() const;
  Rcpp::NumericVector r_get_vars_growth_decomp() const;

  double r_germination_probability(interpolator::Interpolator env);
  bool r_died();
  Plant r_copy() const;

  size_t state_size() const;
  state::iterator       get_state(state::iterator       it) const;
  state::const_iterator set_state(state::const_iterator it);

  // * Access the Control parameter (who needs this?)
protected:
  const Control& control() const;
  integration::intervals_type get_last_integration_intervals() const;
  void set_integration_intervals(integration::intervals_type x);

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
  double assimilation(const Environment& environment);
  double compute_assimilation(const Environment& environment);
  // Used internally, corresponding to the inner term in [eqn 12]
  double compute_assimilation_x(double x, const Environment& environment) const;
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
  // change in height per change in leaf area
  double dheight_dleaf_area() const;
  // Mass of stem needed for new unit mass leaf, d m_s / d m_l
  double dmass_sapwood_dmass_leaf() const;
  // Mass of bark needed for new unit mass leaf, d m_b / d m_l
  double dmass_bark_dmass_leaf() const;
  // Mass of root needed for new unit mass leaf, d m_r / d m_l
  double dmass_root_dmass_leaf() const;
  // Growth rate of leaf area per unit time
  double dleaf_area_dt() const;
  // Growth rate of spawood area at base per unit time
  double dsapwood_area_dt() const;
  // Growth rate of bark area at base per unit time
  double dbark_area_dt() const;
  // Growth rate of stem basal per unit time
  double dbasal_area_dt() const;
  // Growth rate of basal dimater per unit basal area
  double dbasal_diam_dbasal_area() const;
  // Growth rate of basal dimaterper unit time
  double dbasal_diam_dt() const;

  // Sapwood area
  double sapwood_area() const;
  // Bark area
  double bark_area() const;
  // basal area
  double basal_area() const;

  // Update a number of constants within the model.  This is a work in
  // progress.
  static double height_seed(Strategy *s);
  double mass_live_given_height(double h);
  double height_given_mass_leaf(double mass_leaf_) const;

  // To simplify my life, I'm making a small internal-only class that
  // contains some implementation details here.

  class internals {
  public:
    internals();
    bool operator==(const Plant::internals &rhs) const;
    // * Individual size
    // Mass of leaves.  This is the core independent variable
    double mass_leaf;      // [eqn 1]
    // Other size variables that follow directly from `mass_leaf`:
    double leaf_area;      // [eqn 2]
    double height;         // [eqn 3]
    double mass_sapwood;   // [eqn 4]
    double mass_bark;      // [eqn 5]
    double mass_heartwood; // [eqn 6]
    double area_heartwood;
    double mass_root;      // [eqn 7] (fine roots)
    double mass_live;      // [eqn 8]
    // * Mass production
    double assimilation;   // [eqn 12] Gross annual CO2 assimilation
    double respiration;    // [eqn 13] Total maintenance respiration
    double turnover;       // [eqn 14] Total turnover
    double net_production; // [eqn 15] Net production
    double reproduction_fraction; // [eqn 16]
    double fecundity_rate; // [eqn 17] Rate of offspring production
    double leaf_fraction;  // [eqn 18] Fraction of mass growth that is leaves
    double mass_leaf_growth_rate; // [eqn 19] Growth rate in leaf mass
    double height_growth_rate;    // [doc/details.md]
    // * Mortality
    double mortality_rate; // [eqn 21]
    // * Variables
    double mortality;
    double fecundity;
  };

  Strategy::ptr strategy;
  internals vars;
  integration::intervals_type integration_intervals;

  static const int ode_dimension = 5;
};

namespace test {
bool test_plant(Strategy s, bool copy, bool ptr);
interpolator::Interpolator
compute_assimilation_fn(Strategy s, double hmin, double hmax,
			const Environment &environment);
}

}

RCPP_EXPOSED_CLASS_NODECL(model::Plant)

#endif
