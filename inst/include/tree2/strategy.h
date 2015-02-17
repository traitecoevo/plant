// -*-c++-*-
#ifndef TREE2_STRATEGY_H_
#define TREE2_STRATEGY_H_

#include <tree2/control.h>
#include <tree2/environment.h>
#include <tree2/qag_internals.h> // quadrature::intervals_type
#include <tree2/environment.h>

namespace tree2 {

struct Strategy {
public:
  Strategy();

  // * Size

  // [eqn 2] leaf_area (inverse of [eqn 3])
  double leaf_area(double height) const;

  // [eqn 1] leaf_mass (inverse of [eqn 2])
  double leaf_mass(double leaf_area) const;

  // [eqn 4] area and mass of sapwood
  double sapwood_area(double leaf_area) const;
  double sapwood_mass(double sapwood_area, double height) const;

  // [eqn 5] area and mass of bark
  double bark_area(double leaf_area) const;
  double bark_mass (double bark_area, double height) const;

  double basal_area(double bark_area, double sapwood_area,
                            double heartwood_area) const;
  double diameter(double basal_area) const;

  // [eqn 7] Mass of (fine) roots
  double root_mass(double leaf_area) const;

  // [eqn 8] Total Mass
  double live_mass(double leaf_mass, double bark_mass,
                   double sapwood_mass, double root_mass) const;

  double total_mass(double leaf_mass, double bark_mass, double sapwood_mass,
                    double heartwood_mass, double root_mass) const;

  double above_ground_mass(double leaf_mass, double bark_mass,
                           double sapwood_mass, double root_mass) const;

   // * Mass production
  // [eqn 12] Gross annual CO2 assimilation
  double assimilation(const Environment& environment, double height,
                      double leaf_area, bool reuse_intervals);
  // Used internally, corresponding to the inner term in [eqn 12]
  double compute_assimilation_x(double x, double height,
                                const Environment& environment) const;
  double compute_assimilation_h(double h, double height,
                                const Environment& environment) const;
  double compute_assimilation_p(double p, double height,
                                const Environment& environment) const;
  // [Appendix S6] Per-leaf photosynthetic rate.
  double assimilation_leaf(double x) const;



  // [eqn 13] Total maintenance respiration
  double respiration(double leaf_mass, double sapwood_mass,
                     double bark_mass, double root_mass) const;

  double respiration_leaf(double mass) const;
  double respiration_bark(double mass) const;
  double respiration_sapwood(double mass) const;
  double respiration_root(double mass) const;

  // [eqn 14] Total turnover
  double turnover(double leaf_mass, double bark_mass,
                  double sapwood_mass, double root_mass) const;
  double turnover_leaf(double mass) const;
  double turnover_bark(double mass) const;
  double turnover_sapwood(double mass) const;
  double turnover_root(double mass) const;

  // [eqn 15] Net production
  double net_production(double assimilation, double respiration,
                        double turnover);
  // [eqn 16] Fraction of whole plan growth that is leaf
  double reproduction_fraction(double height) const;
  double growth_fraction(double height) const;
  // [eqn 17] Rate of offspring production
  double dfecundity_dt(double net_production,
                       double reproduction_fraction) const;

  // [eqn 18] Fraction of mass growth that is leaves
  double leaf_area_deployment_mass(double leaf_area) const;

  // change in height per change in leaf area
  double dheight_dleaf_area(double leaf_area) const;
  // Mass of stem needed for new unit area leaf, d m_s / d a_l
  double dsapwood_mass_dleaf_area(double leaf_area) const;
  // Mass of bark needed for new unit area leaf, d m_b / d a_l
  double dbark_mass_dleaf_area(double leaf_area) const;
  // Mass of root needed for new unit area leaf, d m_r / d a_l
  double droot_mass_dleaf_area(double leaf_area) const;
  // Growth rate of basal diameter per unit basal area
  double dbasal_diam_dbasal_area(double bark_area, double sapwood_area,
                            double heartwood_area) const;
  // Growth rate of components per unit time:
  double dleaf_area_dt(double leaf_area_growth_rate) const;
  double dsapwood_area_dt(double leaf_area_growth_rate) const;
  double dheartwood_area_dt(double leaf_area) const;
  double dbark_area_dt(double leaf_area_growth_rate) const;
  double dbasal_area_dt(double leaf_area, double leaf_area_growth_rate) const;
  double dbasal_diam_dt(double leaf_area, double bark_area, double sapwood_area,
                        double heartwood_area, double leaf_area_growth_rate) const;
  double droot_mass_dt(double leaf_area,
                       double leaf_area_growth_rate) const;
  double dlive_mass_dt(double reproduction_fraction,
                       double net_production) const;
  double dtotal_mass_dt(double reproduction_fraction,
                        double net_production,
                        double dheartwood_mass_dt) const;
  double dabove_ground_mass_dt(double leaf_area,
                               double reproduction_fraction,
                               double net_production,
                               double dheartwood_mass_dt,
                               double leaf_area_growth_rate) const;

  double dheartwood_mass_dt(double sapwood_mass) const;

  double live_mass_given_height(double height) const;
  double height_given_leaf_mass(double leaf_mass_) const;


  double mortality_dt(double productivity_area) const;
  double mortality_growth_independent_dt(double d0)const ;
  double mortality_growth_dependent_dt(double d2, double d3,
                                       double productivity_area
                                       ) const;
  // [eqn 20] Survival of seedlings during germination
  double germination_probability(double leaf_area,
                                 double net_production) const;

  // * Competitive environment
  // [eqn 11] total leaf area above height above height `z` for given plant
  double leaf_area_above(double z, double height, double leaf_area) const;
  // [eqn  9] Probability density of leaf area at height `z`
  double q(double z, double height) const;
  // [eqn 10] Fraction of leaf area above height `z`
  double Q(double z, double height) const;
  // [      ] Inverse of Q: height above which fraction 'x' of leaf found
  double Qp(double x, double height) const;

  // The aim is to find a plant height that gives the correct seed mass.
  double height_seed(void) const;

  // Set constants within Strategy
  void prepare_strategy();

  // Every Strategy needs a set of Control objects -- these govern
  // things to do with how numerical calculations are performed,
  // rather than the biological control that this class has.
  Control control;

  // Previously there was an "integrator" here.  I'm going to stick
  // that into Control or Environment instead.

  // * Core traits
  double lma, rho, hmat, s, n_area;
  // * Individual allometry
  // Canopy shape parameters
  double eta, eta_c;
  // Leaf area per sapwood area
  double theta;
  // Empirical constants for scaling relationships
  double a1, B1, a3;
  // Bark area per sapwood area
  double b;
  // * Production
  // Respiration constants
  double c_Rs, c_Rb, c_Rr, c_Rl;
  // Yield = carbon fixed in tissue per carbon assimilated;
  double Y;
  // Conversion factor
  double c_bio;
  // Leaf, bark sapwood, and root turnover rates
  double k_l, k_b, k_s, k_r;
  // Leaf productivity parameters  - only used when no N reallocation
  double c_p1, c_p2;

  // * Seed production
  // Accessory cost of reproduction, per seed
  double c_acc;
  // Proportion production alloctaed to reproduction
  double c_r1;
  // Size range across which individuals mature
  double c_r2;

  // * Mortality
  // Parameter for seedling mortality
  double c_s0;
  // Baseline structural mortality rate
  double c_d0;
 // Baseline for growth mortality rate
  double c_d2;
  // Coefficient for dry mass production in mortality function
  double c_d3;

  // Height of a (germinated) seed
  double height_0;
};

}

#endif
