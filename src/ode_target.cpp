#include <Rcpp.h>
#include "ode_target.h"

namespace ode {

bool OdeTarget::r_ode_values_set(std::vector<double> y) {
  if ( y.size() != ode_size() )
    Rf_error("Incorrect size input (expected %d, recieved %d)",
             ode_size(), y.size());
  
  bool changed = false;  
  ode_values_set(y.begin(), changed);
  return changed;
}
std::vector<double> OdeTarget::r_ode_values() const {
  std::vector<double> ret(ode_size());
  ode_values(ret.begin());
  return ret;
}
std::vector<double> OdeTarget::r_ode_rates() const {
  std::vector<double> ret(ode_size());
  ode_rates(ret.begin());
  return ret;
}

}
