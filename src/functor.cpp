#include <Rcpp.h>
#include "functor.h"

namespace util {

double helper_functor(double x, void *data) {
  DFunctor *f = static_cast<DFunctor *>(data);
  return (*f)(x);
}

namespace test {

// Possibly better via modules?
std::vector<double> test_functor(std::vector<double> x, 
				 std::vector<double> pars) {
  if ( (int)pars.size() != 3 )
    Rf_error("Expected parameters of length 3");
  Quadratic obj(pars[0], pars[1], pars[2]);
  Functor<Quadratic, &Quadratic::mytarget> fun(&obj);
  
  std::vector<double> ret;
  for ( std::vector<double>::iterator it = x.begin(); it != x.end(); it++ )
    ret.push_back(fun(*it));
  
  return ret;
}

}

}
