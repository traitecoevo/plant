// -*-c++-*-
#ifndef TREE_STRATEGY_H_
#define TREE_STRATEGY_H_

#include <Rcpp.h>
#include <string>

#include "control.h"

#include "lookup.h"
#include "util.h"

namespace model {

class Strategy : public util::Lookup {
public:
  typedef util::PtrWrapper<Strategy> ptr;
  Strategy();
  Strategy(Rcpp::List x);

  // This will be accessed by a number of things.
  // TODO: evaluate if this should be public.  Probably better in
  // constructor?
  Control control;

  // All the rest of the class can be accessed only by Plant.
  friend class Plant;

private:
  // * Core traits
  double lma, rho, hmat, s;

  // * Individual allometry
  // Canopy shape parameters
  double eta, eta_c;
  // Leaf area per sapwood area
  double theta;
  // Empirical constants for scaling relationships
  double a1, B1, a2, B2, a3, a4, B4;
  // Bark area per sapwood area
  double b;
  
  // * Production
  // Leaf nitrogen per area
  double n_area;
  // Respiration constants
  double c_Rs, c_Rb, c_Rr, c_Rl;
  // Yield = carbon fixed in tissue per carbon assimilated;
  double Y;
  // Conversion factor
  double c_bio;
  // Leaf, bark and root turnover rates
  double k_l, k_b, k_r;
  // Leaf productivity parameters  - only used when no N reallocation
  double c_p1, c_p2;

  // * Seed production
  // Accessory cost of reproduction - multiplication factor
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
  // Coeffcieint for wood density in mortality function
  double c_d1;
  // Baseline for growth mortality rate
  double c_d2;
  // Coefficient for dry mass production in mortality function
  double c_d3;

  // Leaf mass of a seed
  double mass_leaf_0;

  // These things are really not to be used by anything, but are all
  // harmless (except for reset, actually).  They're used in
  // construction, and in mantaining the lookup table used to enable
  // {get,set}_parameters().
  void do_build_lookup();
  void reset();
  void set_parameters_post_hook();
};

}

RCPP_EXPOSED_CLASS(model::Strategy)

#endif
