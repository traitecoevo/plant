#include "util.h"   // check_bounds

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

  strategies.clear();
}

void Parameters::add_strategy(Strategy s) {
  strategies.push_back(s);
}

Strategy Parameters::r_get_strategy(size_t idx) {
  util::check_bounds(idx, strategies.size());
  return strategies[idx];
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
}

}
