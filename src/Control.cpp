#include <Rcpp.h>    // Rf_error
#include <algorithm> // std::max
#include <cmath>     // fabs

#include "Control.h"


Control::Control() {
  eps_abs = 1e-8; // hopefully sane default
  eps_rel = 1e-8; // hopefully sane default
  a_y    = 1.0;   // tune off changes in y
  a_dydt = 0.0;   // but not off changes in dydt
  step_size_shrank_ = false;
}

void Control::set_eps_abs(double x) {
  if ( x < 0 )
    Rf_error("eps_abs is negative");
  eps_abs = x;
}

void Control::set_eps_rel(double x) {
  if ( x < 0 )
    Rf_error("eps_rel is negative");
  eps_rel = x;
}

void Control::set_a_y(double x) {
  if ( x < 0 )
    Rf_error("a_y is negative");
  a_y = x;
}

void Control::set_a_dydt(double x) {
  if ( x < 0 )
    Rf_error("a_dydt is negative");
  a_dydt = x;
}

// This did used to signal decrease, increase or no change, but let's
// find out how that worked before doing that.
double Control::adjust_step_size(int dim, unsigned int ord, 
				 double step_size,
				 const std::vector<double> &y,
				 const std::vector<double> &yerr,
				 const std::vector<double> &dydt) {
  double rmax = DBL_MIN;
  const double S = 0.9;

  for ( int i = 0; i < dim; i++ ) {
    // TODO: This is the same as that computed by errlevel(), but
    // without the check?
    //   const double D0 = errlevel(*y++, *dydt++, h);
    // Could pull that out into an inline function, but I doubt the
    // check costs that much.
    const double D0 =
      eps_rel * (a_y*fabs(y[i]) + a_dydt*fabs(step_size * dydt[i])) +
      eps_abs;
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
    // or std::max(0.2, S / pow(rmax, 1.0 / ord));
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

double Control::errlevel(double y, double dydt, double h) {
  const double errlev = eps_rel * (a_y    * fabs(y       )  + 
				   a_dydt * fabs(h * dydt)) + eps_abs;
  if ( errlev <= 0.0 )
    Rf_error("errlev <= zero");
  return errlev;
}

bool Control::step_size_shrank() const {
  return step_size_shrank_;
}
