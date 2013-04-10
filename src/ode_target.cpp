#include <Rcpp.h>
#include "ode_target.h"
#include "util.h"

namespace ode {

void OdeTarget::derivs(double time, iter_const y, iter dydt) {
  Rf_error("derivs method not implemented");
}

std::vector<double> OdeTarget::r_derivs(double time, 
					std::vector<double> y) {
  util::check_length(y.size(), ode_size());
  std::vector<double> dydt(y.size());
  derivs(time, y.begin(), dydt.begin());
  return dydt;
}

bool OdeTarget::r_ode_values_set(std::vector<double> y) {
  util::check_length(y.size(), ode_size());
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
