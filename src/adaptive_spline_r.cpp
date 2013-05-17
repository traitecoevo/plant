#include "adaptive_spline_r.h"

namespace spline {

namespace test {

Spline test_adaptive_spline(SEXP fun, SEXP env,
			    double a, double b) {
  RFunctionWrapper obj(fun, env);
  util::Functor<RFunctionWrapper, &RFunctionWrapper::target>
    target(&obj);
  AdaptiveSpline generator(&target);
  Spline spline = generator.construct_spline(a, b);
  return spline;
}

}

}
