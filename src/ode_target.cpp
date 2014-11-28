#include <Rcpp.h>
#include <tree/ode_target.h>
#include <tree/util.h>

namespace ode {

void OdeTarget::derivs(double time, iterator_const y, iterator dydt) {
  set_ode_values(time, y);
  ode_rates(dydt);
}

std::vector<double> OdeTarget::r_derivs(double time, 
					std::vector<double> y) {
  util::check_length(y.size(), ode_size());
  std::vector<double> dydt(y.size());
  derivs(time, y.begin(), dydt.begin());
  return dydt;
}

void OdeTarget::r_set_ode_values(double time, std::vector<double> y) {
  util::check_length(y.size(), ode_size());
  set_ode_values(time, y.begin());
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
