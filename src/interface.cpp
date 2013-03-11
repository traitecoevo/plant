#include <Rcpp.h>

#include "Spline.h"
#include "AdaptiveSpline.h"
#include "RAdaptiveSpline.h"

#include "Evolve.h"

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

  Rcpp::class_<Evolve>("Evolve")
    .constructor()
    .method("set_state",  &Evolve::set_state)
    .method("get_state",  &Evolve::get_state)
    .method("get_time",   &Evolve::get_time)
    .method("step",       &Evolve::step)
    .method("step_fixed", &Evolve::step_fixed)
    .method("advance",    &Evolve::advance)
    // R interface only...
    .method("run",        &Evolve::r_run)
    .method("derivs",     &Evolve::r_derivs)
    ;
}
