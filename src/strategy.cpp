#include <tree2/strategy.h>
#include <RcppCommon.h> // NA_REAL

namespace tree2 {

// TODO: There is some fairly major re-plumbing needed here; we need
// to separate out hyperparameters from the ones the model cares
// about, and possibly define a post-parameter setting hook so that
// other intermediates can be generated.

Strategy::Strategy() {
  // * Core traits - default values
  lma_0  = 0.1978791;  // leaf mas per area
  rho_0  = 608;        // wood density
  hmat_0 = 16.5958691; // Height at maturation
  s_0    = 3.8e-5;  // Seed size
  n_area_0 = 1.87e-3; // Leaf nitrogen per area (= Plant::v) [kg / m2]

  // * To start with set actiual values to default
  // TODO: This needs replacing with hyperparameters?
  lma    = lma_0;
  rho    = rho_0;
  hmat   = hmat_0;
  s      = s_0;
  n_area = n_area_0;

  // * Individual allometry
  // Canopy shape parameter (extra calculation here later)
  eta = 12;
  // ratio leaf area to sapwood area
  theta  = 4669;
  // Height - leaf mass scaling
  a1     = 5.44;
  B1     = 0.306;
  // Root - leaf scaling
  a3     = 0.07;
  // Ratio of bark area : sapwood area
  b      = 0.17;

  // * Production
  // Ratio of leaf dark respiration to leaf nitrogen mass
  // [mol CO2 / kgN / yr] (6.66e-4 * (365*24*60*60))
  c_Rl   = 2.1e4;
  // Root respiration per mass [mol CO2 / kg / yr]
  c_Rr   = 217;
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
  k_l0   =  0.4565855;
  B4     = 1.71;
  // Bark turnover
  k_b    = 0.2;
  // Sapwood turnover
  k_s0    = 0.2;
  B5      = 0;
  // Root turnover
  k_r    = 1.0;
  // Parameters of the hyperbola for annual LRC
  c_p1   = 150.36;
  c_p2   = 0.19;

  // * Seed production
  // Accessory cost of reproduction, kg for seed
  // with mass s_0
  // TODO: Deal with hyperparameters here?
  c_acc  = 3.0 * s_0;
  // Scaling of seed accessory costs with seed mass
  B7     = 1.0;

  // Maximum alloction to reproduction
  c_r1   = 1.0;
  // Size range across which individuals mature
  c_r2   = 50;

  // * Mortality parameters
  // Parameter for seedling survival
  c_s0    = 0.1;
  // Baseline for intrinsic mortality
  c_d0    = 0.01;
  // Coefficient for wood density in mortality function
  c_d1    = 0.0;
  // Coefficient for height density in mortality function
  B6    = 0.0;
  // Baseline rate for growth-related mortality
  c_d2    = 5.5;
  // Risk coefficient for dry mass production (per area)
  c_d3    = 20.0;

  // Will get computed properly by Plant.
  height_0    = NA_REAL;

  // NOTE: These need setting at some point
  // TODO: These possibly need setting as hyperparameters?
  eta_c = NA_REAL;
  k_l = NA_REAL;
  k_s = NA_REAL;
}

}
