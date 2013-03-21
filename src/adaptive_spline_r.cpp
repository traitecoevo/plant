#include "adaptive_spline_r.h"

namespace spline {

namespace test {

Spline test_adaptive_spline(SEXP fun, SEXP env,
			    double a, double b) {
  RFunctionWrapper obj(fun, env);
  util::Functor<RFunctionWrapper, &RFunctionWrapper::target>
    target(&obj);
  AdaptiveSpline spline;
  spline.set_bounds(a, b);
  spline.set_target(&target);
  spline.construct_spline(); // should be construct(&target, a, b);
  return spline;
}

}

}
