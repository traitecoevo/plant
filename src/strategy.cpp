#include "strategy.h"
#include "plant.h" // for Plant::prepare_strategy;

namespace model {

Strategy::Strategy()
  : integrator(integrator_from_control(control)) {
  reset();
  set_parameters_post_hook();
}

Strategy::Strategy(Rcpp::List x)
  : integrator(integrator_from_control(control)) {
  reset();
  set_parameters(x);
}

void Strategy::set_control(Control x) {
  control = x;
  integrator = integrator_from_control(control);
}

const Control& Strategy::get_control() const {
  return control;
}

Control Strategy::r_control() const {
  return control;
}

interpolator::Interpolator Strategy::r_assimilation_fn() const {
  return assimilation_fn;
}
void Strategy::r_set_assimilation_fn(interpolator::Interpolator x) {
  assimilation_fn = x;
}

Strategy Strategy::r_copy() const {
  return *this;
}

void Strategy::reset() {
  // * Core traits - default values
  lma_0  = 0.1978791;  // leaf mas per area
  rho_0  = 608;        // wood density
  hmat_0 = 16.5958691; // Height at maturation
  s_0    = 3.8e-5;  // Seed size
  n_area_0 = 1.87e-3; // Leaf nitrogen per area (= Plant::v) [kg / m2]

  // * To start with set actiual values to default
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
  a4     = 0.0286;
  B4     = 1.71;
  // Bark turnover
  k_b    = 0.2;
  // Sapwood turnover
  k_s    = 0.2;
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
  // Parameter for seedling survival
  c_s0    = 0.1;
  // Baseline for intrinsic mortality
  c_d0    = 0.01;
  // Coefficient for wood density in mortality function
  c_d1    = 0.0065;
  // Baseline rate for growth-related mortality
  c_d2    = 5.5;
  // Risk coefficient for dry mass production (per area)
  c_d3    = 20.0;

  // Will get computed properly by Plant.
  height_0    = NA_REAL;

  eta_c = NA_REAL;
  k_l = NA_REAL;
}

double Strategy::assimilation_fn_lookup(double h) const {
  return assimilation_fn.eval(h);
}

void Strategy::do_build_lookup() {
  lookup_table["lma"]    = &lma;
  lookup_table["hmat"]   = &hmat;
  lookup_table["rho"]    = &rho;
  lookup_table["s"]      = &s;
  lookup_table["eta"]    = &eta;
  lookup_table["theta"]  = &theta;
  lookup_table["a1"]     = &a1;
  lookup_table["B1"]     = &B1;
  lookup_table["a3"]     = &a3;
  lookup_table["b"]      = &b;
  lookup_table["a4"]     = &a4;
  lookup_table["B4"]     = &B4;
  lookup_table["n_area"] = &n_area;
  lookup_table["c_p1"]   = &c_p1;
  lookup_table["c_p2"]   = &c_p2;
  lookup_table["c_Rl"]   = &c_Rl;
  lookup_table["c_Rs"]   = &c_Rs;
  lookup_table["c_Rb"]   = &c_Rb;
  lookup_table["c_Rr"]   = &c_Rr;
  lookup_table["k_b"]    = &k_b;
  lookup_table["k_s"]    = &k_s;
  lookup_table["k_r"]    = &k_r;
  lookup_table["Y"]      = &Y;
  lookup_table["c_bio"]  = &c_bio;
  lookup_table["c_acc"]  = &c_acc;
  lookup_table["c_r1"]   = &c_r1;
  lookup_table["c_r2"]   = &c_r2;
  lookup_table["c_s0"]   = &c_s0;
  lookup_table["c_d0"]   = &c_d0;
  lookup_table["c_d1"]   = &c_d1;
  lookup_table["c_d2"]   = &c_d2;
  lookup_table["c_d3"]   = &c_d3;

  // It might be nice to flag some variables as readonly, that way we
  // could return height_0, etc_c, and k_l.
}

void Strategy::set_parameters_post_hook() {
  Plant::prepare_strategy(this);
}

// Check that all parameters are non-negative (zeros are OK, even
// though they will cause some problems).
bool Strategy::validate_parameters(Rcpp::List x) const {
  const std::vector<std::string> names = x.names();
  for (int i = 0; i < x.size(); ++i) {
    double tmp = Rcpp::as<double>(x[i]);
    if (tmp < 0)
      Rcpp::stop("Parameter " + names[static_cast<size_t>(i)] +
		 " must be non-negative");
  }
  return true;
}

// Static method:
integration::QAG Strategy::integrator_from_control(const Control& control) {
  if (control.plant_assimilation_adaptive) {
    integration::QAG ret(control.plant_assimilation_rule,
			 control.plant_assimilation_iterations,
			 control.plant_assimilation_tol,
			 control.plant_assimilation_tol);
    return ret;
  } else {
    integration::QAG ret(control.plant_assimilation_rule);
    return ret;
  }
}

}
