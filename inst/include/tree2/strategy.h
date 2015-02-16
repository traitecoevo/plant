// -*-c++-*-
#ifndef TREE2_STRATEGY_H_
#define TREE2_STRATEGY_H_

#include <tree2/control.h>

namespace tree2 {

struct Strategy {
public:
  Strategy();

  // * Size
  // [eqn 4] Sapwood area
  double sapwood_area(double leaf_area) const;

  // * Mass production

  // [eqn 13] Total maintenance respiration
  double respiration(double leaf_area, double sapwood_mass,
                     double bark_mass, double root_mass) const;
  // [eqn 14] Total turnover
  double turnover(double leaf_mass, double bark_mass,
                  double sapwood_mass, double root_mass) const;
  // [eqn 15] Net production
  double net_production(double assimilation, double respiration,
                        double turnover);
  // [eqn 16] Fraction of whole plan growth that is leaf
  double reproduction_fraction(double height) const;
  // [eqn 17] Rate of offspring production
  double dfecundity_dt(double net_production,
                       double reproduction_fraction) const;

  // [eqn 18] Fraction of mass growth that is leaves
  double leaf_mass_fraction(double leaf_area) const;
  double leaf_area_fraction(double leaf_area) const;

  // change in height per change in leaf area
  double dheight_dleaf_area(double leaf_area) const;
  // Mass of stem needed for new unit mass leaf, d m_s / d m_l
  double dsapwood_mass_dleaf_mass(double leaf_area) const;
  // Mass of bark needed for new unit mass leaf, d m_b / d m_l
  double dbark_mass_dleaf_mass(double leaf_area) const;
  // Mass of root needed for new unit mass leaf, d m_r / d m_l
  double droot_mass_dleaf_mass(double leaf_area) const;
  // Growth rate of basal diameter per unit basal area
  double dbasal_diam_dbasal_area(double leaf_area,
                                 double heartwood_area) const;
  // Growth rate of components per unit time:
  double dleaf_area_dt(double leaf_mass_growth_rate) const;
  double dsapwood_area_dt(double leaf_mass_growth_rate) const;
  double dheartwood_area_dt(double leaf_area) const;
  double dbark_area_dt(double leaf_mass_growth_rate) const;
  double dbasal_area_dt(double leaf_area, double leaf_mass_growth_rate) const;
  double dbasal_diam_dt(double leaf_area, double heartwood_area,
                        double leaf_mass_growth_rate) const;
  double droot_mass_dt(double leaf_area,
                       double leaf_mass_growth_rate) const;
  double dlive_mass_dt(double reproduction_fraction,
                       double net_production) const;
  double dtotal_mass_dt(double reproduction_fraction,
                        double net_production,
                        double dheartwood_mass_dt) const;
  double dabove_ground_mass_dt(double leaf_area,
                               double reproduction_fraction,
                               double net_production,
                               double dheartwood_mass_dt,
                               double leaf_mass_growth_rate) const;

  // Bark area
  double bark_area(double leaf_area) const;
  // basal area
  double basal_area(double leaf_area, double heartwood_area) const;
  double dheartwood_mass_dt(double sapwood_mass) const;
  double sapwood_turnover(double sapwood_mass) const;

  // double live_mass_given_height(double h); [TODO]
  double height_given_leaf_mass(double leaf_mass_) const;

  quadrature::QAG& integrator() {return control.integrator;}
  // Every Strategy needs a set of Control objects -- these govern
  // things to do with how numerical calculations are performed,
  // rather than the biological control that this class has.
  Control control;

  // Previously there was an "integrator" here.  I'm going to stick
  // that into Control or Environment instead.

  // * Core traits
  double lma, rho, hmat, s, n_area;

  // * Deafult values for core traits
  double lma_0, rho_0, hmat_0, s_0, n_area_0;

  // * Individual allometry
  // Canopy shape parameters
  double eta, eta_c;
  // Leaf area per sapwood area
  double theta;
  // Empirical constants for scaling relationships
  double a1, B1, a3, k_l0, B4, k_s0, B5;
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
  // Accessory cost of reproduction - multiplication factor
  double c_acc;
  // Scaling of seed accessory costs with seed mass
  double B7;
  // Proportion production alloctaed to reproduction
  double c_r1;
  // Size range across which individuals mature
  double c_r2;

  // * Mortality
  // Parameter for seedling mortality
  double c_s0;
  // Baseline structural mortality rate
  double c_d0;
  // Coeffcieint for wood density in mortality function
  double c_d1;
  // Coeffcieint for height in mortality function
  double B6;
  // Baseline for growth mortality rate
  double c_d2;
  // Coefficient for dry mass production in mortality function
  double c_d3;

  // Height of a (germinated) seed
  double height_0;
};

}

#endif
