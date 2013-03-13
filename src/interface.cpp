#include <Rcpp.h>

#include "spline.h"
#include "adaptive_spline.h"
#include "adaptive_spline_r.h"

#include "lorenz.h"
#include "ode_r.h"

#include "strategy.h"
#include "parameters.h"
#include "plant.h"

// Allows Splines to be returned from objects.  Note that this causes
// a copy, so copy constructors are required if there are any pointers
// that would get cleaned up on copy.
RCPP_EXPOSED_CLASS(spline::Spline)
RCPP_EXPOSED_CLASS(spline::AdaptiveSpline)
RCPP_EXPOSED_CLASS(model::Strategy)

RCPP_MODULE(tree) {
  Rcpp::class_<spline::Spline>("Spline")
    .constructor()
    .method("init",   &spline::Spline::init)
    .method("eval",   &spline::Spline::r_eval)
    .method("xy",     &spline::Spline::r_get_xy)
    .property("size", &spline::Spline::size)
    ;

  Rcpp::class_<spline::AdaptiveSpline>("AdaptiveSpline")
    .constructor()
    .derives<spline::Spline>("Spline")
    .method("construct_spline", &spline::AdaptiveSpline::construct_spline)
    ;

  Rcpp::class_<spline::AdaptiveSplineR>("AdaptiveSplineR")
    .derives<spline::AdaptiveSpline>("AdaptiveSpline")
    .constructor<SEXP,SEXP,double,double>()
    .method("target", &spline::AdaptiveSplineR::target)
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

  Rcpp::class_<model::Strategy>("Strategy")
    .constructor()
    .method("get_params", &model::Strategy::get_params)
    .method("set_params", &model::Strategy::set_params)
    ;

  Rcpp::class_<model::Parameters>("Parameters")
    .constructor()
    .property("size",         &model::Parameters::size)
    .method("get_params",     &model::Parameters::get_params)
    .method("get_strategy",   &model::Parameters::get_strategy)
    .method("get_strategies", &model::Parameters::get_strategies)
    .method("set_params",     &model::Parameters::set_params)
    .method("add_strategy",   &model::Parameters::add_strategy)
    .method("set_strategy",   &model::Parameters::set_strategy)
    ;

  Rcpp::class_<model::Plant>("Plant")
    .constructor<model::Strategy>()
    ;
}
