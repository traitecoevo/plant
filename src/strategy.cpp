#include "strategy.h"

namespace model {

Strategy::Strategy() {
  reset();
}

void Strategy::reset() {
  // * Core traits
  lma  = NA_REAL;
  rho  = NA_REAL;
  hmat = NA_REAL;
  s    = NA_REAL;

  // * Individual allometry

  // Canopy shape parameter (extra calculation here later)
  eta = 12;
  // ratio leaf area to sapwood area
  theta  = 4669;
  // Height - leaf mass scaling
  a1     = 5.44;
  B1     = 0.306;
  // Leaf area - stem volume scaling
  a2     = 6.67e-5;
  B2      = 1.75;
  // Root - leaf scaling
  a3     = 0.07;
  // Ratio of bark area : sapwood area
  b      = 0.17;

  // * Production
  // TODO: Missing A0?  [But did it actually do anything?]
  // Leaf nitrogen per area (= Plant::v) [kg / m2]
  n_area = 1.87e-3;
  // Ratio of leaf dark respiration to leaf nitrogen mass
  // [mol CO2 / kgN / yr] (6.66e-4 * (365*24*60*60))
  c_Rl   = 2.1e4;
  // Root respiration per mass
  c_Rr   = 217;                 //mol CO2 / kg /yr
  // Sapwood respiration per stem volume [mol CO2 / m3 / yr]
  c_Rs   = 4012;
  // Bark respiration per stem volume [mol CO2 / m3 / yr]
  // (note, new since paper -- see respiration calculation)
  c_Rb   = 2 * c_Rs;
  // Carbon conversion parameter
  Y      = 0.7;
  // Constant converting assimilated CO2 to dry mass [kg / mol]
  // (12E-3 / 0.49)
  c_bio  = 2.45e-2;
  // Leaf turnover - LMA scaling
  a4     = 0.0286;
  B4     = 1.71;
  // Bark turnover
  k_b    = 0.2;
  // Root turnover
  k_r    = 1.0;
  // Parameters of the hyperbola for annual LRC
  c_p1   = 150.36;
  c_p2   = 0.19;

  // * Seed production
  // Accessory cost of reproduction - multiplication factor
  c_acc  = 4.0;
  // Maximum alloction to reproduction
  c_r1   = 1.0;
  // Size range across which individuals mature
  c_r2   = 50;

  // * Mortality parameters
  // Survival during dispersal
  Pi_0    = 0.25;
  // Parameter for seedling survival
  c_s0    = 0.1;
  // Coeffcieint for wood density in mortality function
  c_d1    = 0.0065;
  // Baseline for intrinsic mortality
  // TODO: In the paper this is given as 0.52 only.
  c_d0    = 0.01/exp(-c_d1*608.0);
  // Baseline rate for growth-related mortality
  c_d2    = 5.5;
  // Risk coefficient for dry mass production (per area)
  c_d3    = 20.0;

  assimilation_over_distribution = false;

  compute_constants();
}

void Strategy::compute_constants() {
  eta_c = 1 - 2/(1 + eta) + 1/(1 + 2*eta);
  k_l = a4 * pow(lma, -B4);
}

void Strategy::do_build_lookup() {
  lookup_table["lma"] = &lma;
  lookup_table["hmat"] = &hmat;
  lookup_table["rho"] = &rho;
  lookup_table["s"] = &s;
  lookup_table["eta"] = &eta;
  lookup_table["theta"] = &theta;
  lookup_table["a1"] = &a1;
  lookup_table["B1"] = &B1;
  lookup_table["a2"] = &a2;
  lookup_table["B2"] = &B2;
  lookup_table["a3"] = &a3;
  lookup_table["b"]  = &b;
  lookup_table["a4"] = &a4;
  lookup_table["B4"] = &B4;
  lookup_table["n_area"] = &n_area;
  lookup_table["c_p1"] = &c_p1;
  lookup_table["c_p2"] = &c_p2;
  lookup_table["c_Rl"] = &c_Rl;
  lookup_table["c_Rs"] = &c_Rs;
  lookup_table["c_Rb"] = &c_Rb;
  lookup_table["c_Rr"] = &c_Rr;
  lookup_table["k_b"] = &k_b;
  lookup_table["k_r"] = &k_r;
  lookup_table["Y"] = &Y;
  lookup_table["c_bio"] = &c_bio;
  lookup_table["c_acc"] = &c_acc;
  lookup_table["c_r1"] = &c_r1;
  lookup_table["c_r2"] = &c_r2;
  lookup_table["Pi_0"] = &Pi_0;
  lookup_table["c_s0"] = &c_s0;
  lookup_table["c_d0"] = &c_d0;
  lookup_table["c_d1"] = &c_d1;
  lookup_table["c_d2"] = &c_d2;
  lookup_table["c_d3"] = &c_d3;
  // lookup_table["assimilation_over_distribution"] = &assimilation_over_distribution;
}

}
