#include "gradient.h"
#include "util.h"

namespace util {

double gradient_fd_forward(DFunctor *f, double x, double dx, double fx) {
  const double x2 = x + dx;
  const double fx2 = (*f)(x2);
  return (fx2 - fx) / dx;
}

double gradient_fd_forward(DFunctor *f, double x, double dx) {
  return gradient_fd_forward(f, x, dx, (*f)(x));
}

double gradient_fd_centre(DFunctor *f, double x, double dx) {
  const double x0 = x - dx/2, x1 = x + dx/2;
  const double fx0 = (*f)(x0), fx1 = (*f)(x1);
  return (fx1 - fx0) / dx;
}

namespace test {
double test_gradient(double x, double dx, bool centre,
		     std::vector<double> pars) {
  util::check_length(pars.size(), 3);
  Quadratic obj(pars[0], pars[1], pars[2]);
  Functor<Quadratic, &Quadratic::mytarget> fun(&obj);

  return centre ? gradient_fd_centre(&fun, x, dx) : 
    gradient_fd_forward(&fun, x, dx);
}
}

}
