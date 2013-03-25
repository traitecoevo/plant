#include "util.h"   // check_bounds
#include "plant.h"  // prepare_strategy

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

  strategies.clear();
}

// TODO: Should this take a Strategy instead (symmetry with
// get_strategy?).  What is easiest?
void Parameters::add_strategy(Rcpp::List x) {
  Strategy s;
  s.set_parameters(x);
  Plant::prepare_strategy(&s);
  strategies.push_back(s);
}

Strategy Parameters::get_strategy(int idx) {
  util::check_bounds(idx, strategies.size());
  return strategies[idx];
}

Rcpp::List Parameters::get_strategies() {
  Rcpp::List ret;
  for ( std::vector<Strategy>::iterator it = strategies.begin();
	it != strategies.end(); it++ )
    ret.push_back(Rcpp::wrap(*it));
  return ret;
}

// TODO: Not sure that this is the correct name.  
// 
// TODO: Not sure that this is even the right way forward.  Treating
// strategies as essentially read-only might be better.
void Parameters::set_strategy(Rcpp::List x, int idx) {
  util::check_bounds(idx, strategies.size());
  strategies[idx].set_parameters(x);
  Plant::prepare_strategy(&strategies[idx]);
}

void Parameters::do_build_lookup() {
  lookup_table["mean_disturbance_interval"] =
    &mean_disturbance_interval;
  lookup_table["c_ext"]      = &c_ext;
  lookup_table["patch_area"] = &patch_area;
}

}
