#include "util.h"
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
  // TODO: What to do about the new species?  How many to reset down
  // to?  Probably zero.
  strategies.clear();
}

void Parameters::add_strategy(Rcpp::List x) {
  Strategy s;
  s.set_parameters(x);
  strategies.push_back(s);
}

// This is generally fucking up quite badly, returning corrupted
// memory.  Probably need to do a valgrind.
Rcpp::List Parameters::get_strategy(int idx) {
  util::check_bounds(idx, strategies.size());
  return strategies[idx].get_parameters();
}

Rcpp::List Parameters::get_strategies() {
  Rcpp::List ret;
  for ( std::vector<Strategy>::iterator it = strategies.begin();
	it != strategies.end(); it++ )
    ret.push_back(it->get_parameters());
  return ret;
}

void Parameters::set_strategy(Rcpp::List x, int idx) {
  util::check_bounds(idx, strategies.size());
  strategies[idx].set_parameters(x);
}

void Parameters::do_build_lookup() {
  lookup_table["mean_disturbance_interval"] =
    &mean_disturbance_interval;
  lookup_table["c_ext"] = &c_ext;
}

}
