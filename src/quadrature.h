// -*-c++-*-
#ifndef TREE_QUADRATURE_H_
#define TREE_QUADRATURE_H_

#include "functor.h"

namespace integration {

// Gauss-Kronrod Quadrature.
class QK {
public:
  QK(size_t rule);
  double integrate(util::DFunctor *f, double a, double b);

  // These are super simple, so inline:
  double get_last_result()     const {return last_result;    }
  double get_last_error()      const {return last_error;     }
  double get_last_result_abs() const {return last_result_abs;}
  double get_last_result_asc() const {return last_result_asc;}

  // * R interface
  double r_integrate(util::RFunctionWrapper fun, double a, double b);

private:
  void initialise(size_t rule);
  double rescale_error(double err, double result_abs,
		       double result_asc) const;

  double eval(util::DFunctor *f, double x) {return (*f)(x);}

  double last_result;
  double last_result_abs;
  double last_result_asc;
  double last_error;

  size_t n;
  std::vector<double> fv1;
  std::vector<double> fv2;
  std::vector<double> xgk, wg, wgk;
};

}

#endif
