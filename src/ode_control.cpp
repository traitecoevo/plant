#include "ode_control.h"

#include <Rcpp.h>    // Rf_error
#include <algorithm> // std::max
#include <cmath>     // fabs

namespace ode {

OdeControl::OdeControl()
  : tol_abs(1e-8), // hopefully sane default
    tol_rel(1e-8), // hopefully sane default
    a_y(1.0),      // tune step size based on changes in y...
    a_dydt(0.0),   // ...but not based on changes in dydt
    step_size_min(1e-8), // hopefully sane default
    step_size_max(1),    // should be inf?
    last_step_size_shrank(false) {
}

OdeControl::OdeControl(double tol_abs, double tol_rel,
		       double a_y, double a_dydt,
		       double step_size_min, double step_size_max)
  : tol_abs(tol_abs),
    tol_rel(tol_rel),
    a_y(a_y),
    a_dydt(a_dydt),
    step_size_min(step_size_min),
    step_size_max(step_size_max),
    last_step_size_shrank(false) {
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
    last_step_size_shrank = true;
    if (step_size < step_size_min) {
      ::Rf_error("Step size became too small");
      step_size = step_size_min; // NOTE: won't get here
    }
  } else if ( rmax < 0.5 ) {
    // increase step, no more than factor of 5
    double r = S / pow (rmax, 1.0 / (ord + 1.0));
    if ( r > 5.0 )
      r = 5.0;
    if ( r < 1.0 ) // Don't allow any decrease caused by S<1
      r = 1.0;
    step_size *= r;
    if (step_size > step_size_max)
      step_size = step_size_max;
    last_step_size_shrank = false;
  } else { 
    // otherwise no change to the size
    last_step_size_shrank = false;
  }

  return step_size;
}

double OdeControl::errlevel(double y, double dydt, double h) const {
  const double errlev = tol_rel * (a_y    * fabs(y       )  +
				   a_dydt * fabs(h * dydt)) + tol_abs;
  if ( errlev <= 0.0 )
    Rf_error("errlev <= zero");
  return errlev;
}

bool OdeControl::step_size_shrank() const {
  return last_step_size_shrank;
}

void OdeControl::do_build_lookup() {
  lookup_table["tol_abs"]       = &tol_abs;
  lookup_table["tol_rel"]       = &tol_rel;
  lookup_table["a_y"]           = &a_y;
  lookup_table["a_dydt"]        = &a_dydt;
  lookup_table["step_size_min"] = &step_size_min;
  lookup_table["step_size_max"] = &step_size_max;
}

}
