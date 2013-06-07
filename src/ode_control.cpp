#include "ode_control.h"

#include <Rcpp.h>    // Rf_error
#include <algorithm> // std::max
#include <cmath>     // fabs

namespace ode {

OdeControl::OdeControl()
  : eps_abs(1e-8), // hopefully sane default
    eps_rel(1e-8), // hopefully sane default
    a_y(1.0),      // tune step size based on changes in y...
    a_dydt(0.0),   // ...but not based on changes in dydt
    step_size_shrank_(false) {
}

void OdeControl::set_eps_abs(double x) {
  if ( x < 0 )
    Rf_error("eps_abs is negative");
  eps_abs = x;
}

void OdeControl::set_eps_rel(double x) {
  if ( x < 0 )
    Rf_error("eps_rel is negative");
  eps_rel = x;
}

void OdeControl::set_a_y(double x) {
  if ( x < 0 )
    Rf_error("a_y is negative");
  a_y = x;
}

void OdeControl::set_a_dydt(double x) {
  if ( x < 0 )
    Rf_error("a_dydt is negative");
  a_dydt = x;
}

double OdeControl::adjust_step_size(size_t dim, unsigned int ord, 
				    double step_size,
				    const std::vector<double> &y,
				    const std::vector<double> &yerr,
				    const std::vector<double> &dydt) {
  double rmax = DBL_MIN;
  const double S = 0.9;

  for ( size_t i = 0; i < dim; i++ ) {
    const double D0 = errlevel(y[i], dydt[i], step_size);
    const double r = fabs(yerr[i]) / fabs(D0);
    rmax = std::max(r, rmax);
  }

  if ( rmax > 1.1 ) {
    // decrease step, no more than factor of 5, but a fraction S more
    // than scaling suggests (for better accuracy).
    double r = S / pow(rmax, 1.0 / ord);
    if (r < 0.2)
      r = 0.2;
    step_size *= r;
    step_size_shrank_ = true;
  } else if ( rmax < 0.5 ) {
    // increase step, no more than factor of 5
    double r = S / pow (rmax, 1.0 / (ord + 1.0));
    if ( r > 5.0 )
      r = 5.0;
    if ( r < 1.0 ) // Don't allow any decrease caused by S<1
      r = 1.0;
    step_size *= r;
    step_size_shrank_ = false;
  } else { 
    // otherwise no change to the size
    step_size_shrank_ = false;
  }

  return step_size;
}

double OdeControl::errlevel(double y, double dydt, double h) {
  const double errlev = eps_rel * (a_y    * fabs(y       )  + 
				   a_dydt * fabs(h * dydt)) + eps_abs;
  if ( errlev <= 0.0 )
    Rf_error("errlev <= zero");
  return errlev;
}

bool OdeControl::step_size_shrank() const {
  return step_size_shrank_;
}

}
