
// Used by AdaptiveSpline::r_set_test
class RFunForSpline {
public:
  RFunForSpline(SEXP fun_) : fun(fun_) {}
  double target(double x) {
    SEXP call = PROTECT(Rf_lang2(fun, Rf_ScalarReal(x))); 
    SEXP res = PROTECT(Rf_eval(call, R_GlobalEnv));
    UNPROTECT(2);
    return REAL(res)[0];
  }
private:
  SEXP fun;
};


// High-level R testing function that takes some R function and builds
// an adaptive spline out of it.  Not terribly useful, but exercises
// most of the code.
AdaptiveSpline r_spline_test(double a_, double b_, 
			     double atol_, double rtol_,
			     int nbase_, int max_depth_, 
			     SEXP fun) {
  RFunForSpline *data_ = new RFunForSpline;
  data_->setup(fun);


  AdaptiveSpline s;

  ASFun target_ = helper_spline<RFunForSpline, &RFunForSpline::target>;

  set_bounds(a_, b_);
  set_target(target_, data_);
  set_control(atol_, rtol_, nbase_, max_depth_);
}
