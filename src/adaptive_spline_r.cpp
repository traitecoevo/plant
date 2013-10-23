#include "adaptive_spline_r.h"

namespace spline {

namespace test {

Spline test_adaptive_spline(SEXP fun, SEXP env,
			    double a, double b,
			    bool akima_spline) {
  util::RFunctionWrapper obj(fun, env);
  // Hopefully sensible defaults:
  const double atol = 1e-6, rtol = 1e-6;
  const int nbase = 17, max_depth = 16;
  AdaptiveSpline generator(atol, rtol, nbase, max_depth, akima_spline);
  Spline spline = generator.construct_spline(&obj, a, b);
  return spline;
}

}

}
