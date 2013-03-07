#include <Rcpp.h>

#include "Spline.h"
#include "AdaptiveSpline.h"
#include "RAdaptiveSpline.h"

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
}
