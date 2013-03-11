#include <Rcpp.h>

#include "Spline.h"
#include "AdaptiveSpline.h"
#include "RAdaptiveSpline.h"

#include "Lorenz.h"
#include "ROde.h"

// Allows Splines to be returned from objects.  Note that this causes
// a copy, so copy constructors are required if there are any pointers
// that would get cleaned up on copy.
RCPP_EXPOSED_CLASS(Spline)
RCPP_EXPOSED_CLASS(AdaptiveSpline)

RCPP_MODULE(tree) {
  Rcpp::class_<Spline>("Spline")
    .constructor()
    .method("init", &Spline::init)
    .method("eval", &Spline::r_eval)
    .method("xy",   &Spline::r_get_xy)
    .property("size", &Spline::size)
    ;

  Rcpp::class_<AdaptiveSpline>("AdaptiveSpline")
    .constructor()
    .derives<Spline>("Spline")
    .method("construct_spline", &AdaptiveSpline::construct_spline)
    ;

  Rcpp::class_<RAdaptiveSpline>("RAdaptiveSpline")
    .derives<AdaptiveSpline>("AdaptiveSpline")
    .constructor<SEXP,SEXP,double,double>()
    .method("target", &RAdaptiveSpline::target)
    ;

  Rcpp::class_<Lorenz>("Lorenz")
    .constructor<double,double,double>()
    .method("derivs", &Lorenz::r_derivs)
    .property("size", &Lorenz::size)
    // ODE solving
    .method("set_state",  &Lorenz::ode_set_state)
    .method("get_state",  &Lorenz::ode_get_state)
    .method("get_time",   &Lorenz::ode_get_time)
    .method("step",       &Lorenz::ode_step)
    .method("step_fixed", &Lorenz::ode_step_fixed)
    .method("advance",    &Lorenz::ode_advance)
    .method("run",        &Lorenz::ode_r_run)
    ;

  Rcpp::class_<ROde>("ROde")
    .constructor<SEXP,SEXP,SEXP>()
    .method("derivs",     &ROde::r_derivs)
    .property("size",     &ROde::size)
    // ODE solving
    .method("set_state",  &ROde::ode_set_state)
    .method("get_state",  &ROde::ode_get_state)
    .method("get_time",   &ROde::ode_get_time)
    .method("step",       &ROde::ode_step)
    .method("step_fixed", &ROde::ode_step_fixed)
    .method("advance",    &ROde::ode_advance)
    .method("run",        &ROde::ode_r_run)
    ;
}
