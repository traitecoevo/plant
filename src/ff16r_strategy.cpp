// Built from  src/ff16_strategy.cpp on Fri Jul  3 08:14:35 2020 using the scaffolder, from the strategy:  FF16
#include <plant/models/ff16r_strategy.h>

namespace plant {

// TODO: Document consistent argument order: l, b, s, h, r
// TODO: Document ordering of different types of variables (size
// before physiology, before compound things?)
// TODO: Consider moving to activating as an initialisation list?
FF16r_Strategy::FF16r_Strategy() {
  // * Core traits - default values
  lma       = 0.1978791;  // Leaf mass per area [kg / m2]
  rho       = 608.0;      // Wood density [kg/m3]
  hmat      = 16.5958691; // Height at maturation [m]
  omega     = 3.8e-5;     // Seed mass [kg]

  // * Individual allometry
  // Canopy shape parameter (extra calculation here later)
  eta       = 12.0; // [dimensionless]
  // Ratio sapwood area area to leaf area
  theta     = 1.0/4669; // [dimensionless]
  // Height - leaf mass scaling
  a_l1        = 5.44; // height with 1m2 leaf [m]
  a_l2        = 0.306; // dimensionless scaling of height with leaf area
  // Root mass per leaf area
  a_r1        = 0.07;  //[kg / m]
  // Ratio of bark area : sapwood area
  a_b1         = 0.17; // [dimensionless]

  // * Production
  // Ratio of leaf dark respiration to leaf mass [mol CO2 / yr  / kg]
  // =  [mol CO2 / m2 / yr]  |  (39.27 = 2100 * 0.00187)  | narea * photosynthesis_per_nitrogen
  //    / [kg(leaf) / m2 ]   |    / (0.1978791)           | lma
  // Hard coded in value of lma here so that this value doesn't change
  // if that trait changes above.
  r_l   = 39.27 / 0.1978791;
  // Root respiration per mass [mol CO2 / yr / kg]
  r_r   = 217.0;
  // Sapwood respiration per stem mass  [mol CO2 / yr / kg]
  // = respiration per volume [mol CO2 / m3 / yr]
  // /  wood density [kg/m3]
  r_s   = 4012.0 / 608.0;
  // Bark respiration per stem mass
  // assumed to be twice rate of sapwood
  // (NOTE that there is a re-parametrisation here relative to the paper
  // -- r_b is defined (new) as 2*r_s, whereas the paper assumes a
  // fixed multiplication by 2)
  r_b   = 2.0 * r_s;
  // Carbon conversion parameter
  a_y      = 0.7;
  // Constant converting assimilated CO2 to dry mass [kg / mol]
  // (12E-3 / 0.49)
  a_bio  = 2.45e-2;
  // Leaf turnover [/yr]
  k_l=  0.4565855;
 // Bark turnover [/yr]
  k_b    = 0.2;
  // Sapwood turnover [/yr]
  k_s     = 0.2;
  // Root turnover [/yr]
  k_r    = 1.0;
  // Parameters of the hyperbola for annual LRC
  a_p1   = 151.177775377968; // [mol CO2 / yr / m2]
  a_p2   = 0.204716166503633; // [dimensionless]

  // * Seed production
  // Accessory cost of reproduction
  a_f3  = 3.0 *  3.8e-5; // [kg per seed]

  // Maximum allocation to reproduction
  a_f1   = 1.0; //[dimensionless]
  // height above hmat at which allocation to reproduction is half its max value
  a_f2   = 2; // [dimensionless]

// * Mortality parameters
  // Probability of survival during dispersal
  S_D   = 0.25; // [dimensionless]
  // Parameter for seedling survival
  a_d0    = 0.1; //[kg / yr / m2]
  // Baseline for intrinsic mortality
  d_I    = 0.01; // [ / yr]
 // Baseline rate for growth-related mortality
  a_dG1    = 5.5; // [ / yr]
  // Risk coefficient for dry mass production (per area)
  a_dG2    = 20.0;// [yr m2 / kg ]

  // Will get computed properly by prepare_strategy
  height_0 = NA_REAL;
  eta_c    = NA_REAL;

  collect_all_auxillary = false;
  // build the string state/aux name to index map
  refresh_indices();
  name = "FF16r";
}


// [eqn 16] Fraction of production allocated to reproduction
double FF16r_Strategy::fraction_allocation_reproduction(double height) const {
  if(height <= hmat)
    return 0.0;
  else
    return a_f1 * (height - hmat) / (a_f2  + (height - hmat));
}


FF16r_Strategy::ptr make_strategy_ptr(FF16r_Strategy s) {
  s.prepare_strategy();
  return std::make_shared<FF16r_Strategy>(s);
}
}
