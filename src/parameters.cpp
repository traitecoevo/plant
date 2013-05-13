#include "util.h"   // check_bounds_r

#include "parameters.h"

namespace model {

Parameters::Parameters() {
  reset();
}

void Parameters::reset() {
  // Mean disturbance interval
  mean_disturbance_interval = 30.0;
  // Light extinction coefficient
  c_ext = 0.5;
  // Patch area (m^2?)
  patch_area = 10.0;
  // Survival during dispersal
  Pi_0 = 0.25;
  // Number of patches
  n_patches = 10;

  strategies.clear();
}

size_t Parameters::size() const {
  return strategies.size();
}

void Parameters::add_strategy(Strategy s) {
  strategies.push_back(s);
}

Strategy Parameters::r_at(size_t idx) {
  return strategies.at(util::check_bounds_r(idx, size()));
}

Rcpp::List Parameters::r_get_strategies() {
  Rcpp::List ret;
  for ( std::vector<Strategy>::iterator it = strategies.begin();
	it != strategies.end(); it++ )
    ret.push_back(Rcpp::wrap(*it));
  return ret;
}

void Parameters::do_build_lookup() {
  lookup_table["mean_disturbance_interval"] =
    &mean_disturbance_interval;
  lookup_table["c_ext"]      = &c_ext;
  lookup_table["patch_area"] = &patch_area;
  lookup_table["Pi_0"]       = &Pi_0;
  lookup_table["n_patches"]  = &n_patches;
}

}
