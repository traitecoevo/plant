// -*-c++-*-
#ifndef TREE2_STRATEGY_H_
#define TREE2_STRATEGY_H_

#include <tree2/control.h>

namespace tree2 {

struct Strategy {
public:
  Strategy();
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
