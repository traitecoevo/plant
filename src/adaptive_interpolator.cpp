#include <tree/adaptive_interpolator.h>

namespace interpolator {

AdaptiveInterpolator::AdaptiveInterpolator(double atol_, double rtol_,
					   int nbase_, int max_depth_,
					   bool akima, bool linear)
  : target(NULL),
    atol(atol_),
    rtol(rtol_),
    nbase(nbase_),
    max_depth(max_depth_),
    dx(NA_REAL),
    dxmin(NA_REAL),
    interpolator(akima, linear) {
}

// Evaluate the underlying function (double -> double).
double AdaptiveInterpolator::eval_target(double x) const {
  return (*target)(x);
}

// Adaptively refine a interpolator that spans from a to b so that the
// midpoints of evaluated points have sufficiently low error.
Interpolator AdaptiveInterpolator::construct(util::DFunctor *target_,
					     double a, double b) {
  target = target_;
  check_bounds(a, b);
  dx = (b - a) / (nbase - 1);
  dxmin = dx / pow(2, max_depth);

  xx.clear();
  yy.clear();
  zz.clear();

  // TODO: If we template seq_len better, we can avoid a copy.  But
  // using seq_len here guarantees that the first and last element are
  // *exactly* 'a' and 'b'.
  // TODO: Why is nbase not size_t?
  std::vector<double> tmp = util::seq_len(a, b, static_cast<size_t>(nbase));
  for (size_t i = 0; i < static_cast<size_t>(nbase); i++) {
    const double xi = tmp[i];
    xx.push_back(xi);
    yy.push_back(eval_target(xi));
    zz.push_back(i > 0);
  }

  compute();

  bool flag = true;
  while (flag)
    flag = refine();

  return interpolator;
}

// Given our current set of x/y points, construct the interpolated
// fucntion.  This is the only function that actually modifies the
// interpolator object.
void AdaptiveInterpolator::compute() {
  std::vector<double> x(xx.begin(), xx.end());
  std::vector<double> y(yy.begin(), yy.end());
  interpolator.init(x, y);
}

// Refine the interpolated function by adding points in all intervals
// (xx[i-1], xx[i]) where z[i] is true, and check to see if this new
// point has sufficiently good error that this interval needs no
// further refinement.
bool AdaptiveInterpolator::refine() {
  dx /= 2;

  if (dx < dxmin)
    Rcpp::stop("Interpolated function as refined as currently possible");

  bool flag = false;

  std::list<double>::iterator xi = xx.begin(), yi = yy.begin();
  std::list<bool>::iterator zi = zz.begin();
  for (; xi != xx.end(); ++xi, ++yi, ++zi) {
    if (*zi) {
      const double x_mid = *xi - dx;
      const double y_mid = eval_target(x_mid);
      const double p_mid = interpolator.eval(x_mid);

      // Always insert the new points.
      xx.insert(xi, x_mid);
      yy.insert(yi, y_mid);

      // Flag for refinement based on
      const bool flag_mid = !check_err(y_mid, p_mid);
      // If error was OK/not OK (flag_mid is true/false), say that
      // this interval is OK/not OK...
      *zi = flag_mid;
      // ...and that that the interval implied by the mid point is
      // OK/not OK.
      zz.insert(zi, flag_mid);

      flag = flag || flag_mid;
    }
  }

  // Recompute interpolator to use new points added during refinement.
  compute();

  return flag;
}

void AdaptiveInterpolator::check_bounds(double a, double b) {
  if (a >= b)
    Rcpp::stop("Impossible bounds");
  if (!util::is_finite(a) || !util::is_finite(b))
    Rcpp::stop("Infinite bounds");
}

// Determines if difference between predicted and true values falls
// within error bounds.
bool AdaptiveInterpolator::check_err(double y_true, double y_pred) const {
  const double err_abs = fabs(y_true - y_pred);
  const double err_rel = fabs(1 - y_pred / y_true);
  return err_abs < atol || err_rel < rtol;
}

namespace test {

Interpolator test_adaptive_interpolator(Rcpp::Function fun,
					double a, double b,
					bool akima, bool linear) {
  util::RFunctionWrapper obj(fun);
  // Hopefully sensible defaults:
  const double atol = 1e-6, rtol = 1e-6;
  const int nbase = 17, max_depth = 16;
  AdaptiveInterpolator generator(atol, rtol, nbase, max_depth,
				 akima, linear);
  return generator.construct(&obj, a, b);
}

// Just exposing some control parameters:
Interpolator run_adaptive_interpolator(Rcpp::Function fun,
				       double a, double b,
				       double tol,
				       int nbase, int max_depth) {
  util::RFunctionWrapper obj(fun);
  const double atol = tol, rtol = tol;
  const bool akima = false, linear = false;
  AdaptiveInterpolator generator(atol, rtol, nbase, max_depth,
				 akima, linear);
  return generator.construct(&obj, a, b);
}

}

}
