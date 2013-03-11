// -*-c++-*-
#ifndef LORENZ_H_
#define LORENZ_H_

#include <vector>
#include <Rcpp.h>

// This is only used as testing code; where do we actually want to put
// it?

class Lorenz {
public:
  Lorenz(double sigma, double R, double b);
  void derivs(double time,
	      std::vector<double>::const_iterator y,
	      std::vector<double>::iterator dydt);
  unsigned int size() const { return 3; }
  
  std::vector<double> r_derivs(double time, 
			       std::vector<double> y);

private:
  const double sigma, R, b;
};


#endif
