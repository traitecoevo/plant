#include <tree/gradient.h>
#include <tree/util.h>

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

double gradient_fd_backward(DFunctor *f, double x, double dx, double fx) {
  return gradient_fd_forward(f, x, -dx, fx);
}

double gradient_fd_backward(DFunctor *f, double x, double dx) {
  return gradient_fd_forward(f, x, -dx);
}

// Based on code in R's numDeriv::grad (actually in grad.default).
//
// First order derivatives are stored in the vector a[r], for r rounds
// of improvement.
//
// We start with deviation from x of d * x, unless x is almost zero
// (determined by being smaller than zero_tol) in which case we use
// d as the absolute deviation.

//------------------------------------------------------------------------
//   Applying Richardson Extrapolation to improve the accuracy of
//   the first and second order derivatives. The algorithm as follows:
//
//   --  For each column of the derivative matrix a,
//	  say, A1, A2, ..., Ar, by Richardson Extrapolation, to calculate a
//	  new sequence of approximations B1, B2, ..., Br used the formula
//
//	     B(i) =( A(i+1)*4^m - A(i) ) / (4^m - 1) ,  i=1,2,...,r-m
//
//		N.B. This formula assumes v=2.
//
//   -- Initially m is taken as 1  and then the process is repeated
//	 restarting with the latest improved values and increasing the
//	 value of m by one each until m equals r-1
//
//-------------------------------------------------------------------------
double gradient_richardson(DFunctor *f, double x, double d, size_t r) {
  const size_t v = 2; // this value is required by scheme (above)
  const double zero_tol = sqrt(DOUBLE_EPS)/7e-7;

  // Initial offset (see above).
  double h = fabs(d * x) + d * (fabs(x) < zero_tol);

  std::vector<double> a;
  for (size_t i = 0; i < r; i++, h /= v)
    a.push_back(((*f)(x + h) - (*f)(x - h))/(2*h));

  for (size_t m = 1; m < r; ++m) {
    const double four_m = pow(4.0, m);
    std::vector<double> a_next;
    for (size_t i = 0; i < r - m; ++i)
      a_next.push_back((a[i+1]*four_m - a[i])/(four_m - 1));
    a = a_next;
  }

  return a.front();
}

namespace test {
double test_gradient(double x, double dx, int type,
		     std::vector<double> pars) {
  util::check_length(pars.size(), 3);
  Quadratic obj(pars[0], pars[1], pars[2]);
  Functor<Quadratic, &Quadratic::mytarget> fun(&obj);
  double ret;

  if (type == -1)
    ret = gradient_fd_backward(&fun, x, dx);
  else if (type == 0)
    ret = gradient_fd_centre(&fun, x, dx);
  else
    ret = gradient_fd_forward(&fun, x, dx);

  return ret;
}

double test_gradient_richardson(double x, double d, size_t r,
				std::vector<double> pars) {
  util::check_length(pars.size(), 3);
  Quadratic obj(pars[0], pars[1], pars[2]);
  Functor<Quadratic, &Quadratic::mytarget> fun(&obj);

  return gradient_richardson(&fun, x, d, r);
}
}

}
