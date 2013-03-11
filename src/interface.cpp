#include <Rcpp.h>

#include "Spline.h"
#include "AdaptiveSpline.h"
#include "RAdaptiveSpline.h"

#include "lorenz.h"
#include "ode_r.h"

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

  Rcpp::class_<ode::Lorenz>("Lorenz")
    .constructor<double,double,double>()
    .method("derivs",     &ode::Lorenz::r_derivs)
    .property("size",     &ode::Lorenz::size)
    // ODE solving
    .method("set_state",  &ode::Lorenz::ode_set_state)
    .method("get_state",  &ode::Lorenz::ode_get_state)
    .method("get_time",   &ode::Lorenz::ode_get_time)
    .method("step",       &ode::Lorenz::ode_step)
    .method("step_fixed", &ode::Lorenz::ode_step_fixed)
    .method("advance",    &ode::Lorenz::ode_advance)
    .method("run",        &ode::Lorenz::ode_r_run)
    ;

  Rcpp::class_<ode::OdeR>("OdeR")
    .constructor<SEXP,SEXP,SEXP>()
    .method("derivs",     &ode::OdeR::r_derivs)
    .property("size",     &ode::OdeR::size)
    // ODE solving
    .method("set_state",  &ode::OdeR::ode_set_state)
    .method("get_state",  &ode::OdeR::ode_get_state)
    .method("get_time",   &ode::OdeR::ode_get_time)
    .method("step",       &ode::OdeR::ode_step)
    .method("step_fixed", &ode::OdeR::ode_step_fixed)
    .method("advance",    &ode::OdeR::ode_advance)
    .method("run",        &ode::OdeR::ode_r_run)
    ;
}
