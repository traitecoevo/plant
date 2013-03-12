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
  s.set_params(x);
  strategies.push_back(s);
}

Rcpp::List Parameters::get_params() const {
  Rcpp::List ret;
  ret["mean_disturbance_interval"] = mean_disturbance_interval;
  ret["c_ext"] = c_ext;
  return ret;
}

// This is generally fucking up quite badly, returning corrupted
// memory.  Probably need to do a valgrind.
Rcpp::List Parameters::get_strategy(int idx) const {
  ::check_bounds(idx, strategies.size());
  return strategies[idx].get_params();
}

Rcpp::List Parameters::get_strategies() const {
  Rcpp::List ret;
  for ( std::vector<Strategy>::const_iterator it = strategies.begin();
	it != strategies.end(); it++ )
    ret.push_back(it->get_params());
  return ret;
}

// Not done with a lookup table (cf Strategy) because there are only
// two cases at the moment.  However, if I did do it the same way, I
// either need to move the support into the `utils.cpp` file as
// standalone functions, or as a small class that just wraps around
// map.  That might actually deal with some of the copy issues,
// actually.
void Parameters::set_params(Rcpp::List x) {
  std::vector<std::string> names = x.names();
  for ( int i = 0; i < x.size(); i++ ) {
    std::string key = names[i];
    double *tmp;
    if ( names[i] == "mean_disturbance_interval" )
      tmp = &mean_disturbance_interval;
    else if ( names[i] == "c_ext" )
      tmp = &mean_disturbance_interval;
    else
      ::Rf_error("Key %s not found", key.c_str());
    *tmp = Rcpp::as<double>(x[i]);
  }
}

void Parameters::set_strategy(Rcpp::List x, int idx) {
  ::check_bounds(idx, strategies.size());
  strategies[idx].set_params(x);
}


}
