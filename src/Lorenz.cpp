#include "Lorenz.h"

Lorenz::Lorenz(double sigma, double R, double b) :
  sigma(sigma), R(R), b(b) {}

void Lorenz::derivs(double time, 
		    std::vector<double>::const_iterator y,
		    std::vector<double>::iterator dydt) {
  const double y0 = *y++;
  const double y1 = *y++;
  const double y2 = *y++;
  
  *dydt++ = sigma * ( y1 - y0 );
  *dydt++ = R * y0 - y1 - y0 * y2;
  *dydt++ = -b * y2 + y0 * y1;
}

// This is something that can move up to a different level, as it's
// just straight up boilerplate.
std::vector<double> Lorenz::r_derivs(double time,
				     std::vector<double> y) {
  if ( y.size() != size() )
    Rf_error("Incorrect parameter vector; expected %d, got %d",
	     size(), y.size());
  std::vector<double> dydt(size());
  derivs(time, y.begin(), dydt.begin());
  return dydt;
}
