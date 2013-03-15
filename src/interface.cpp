#include <Rcpp.h>

#include "spline.h"
#include "adaptive_spline.h"
#include "adaptive_spline_r.h"

#include "lorenz.h"
#include "ode_r.h"

#include "lookup.h"

#include "strategy.h"
#include "parameters.h"
#include "plant.h"

#include "functor.h"
#include "find_root.h"
#include "integrator.h"

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

  Rcpp::class_<ode::test::Lorenz>("Lorenz")
    .constructor<double,double,double>()
    .method("derivs",     &ode::test::Lorenz::r_derivs)
    .property("size",     &ode::test::Lorenz::size)
    // ODE solving
    .method("set_state",  &ode::test::Lorenz::ode_set_state)
    .method("get_state",  &ode::test::Lorenz::ode_get_state)
    .method("get_time",   &ode::test::Lorenz::ode_get_time)
    .method("step",       &ode::test::Lorenz::ode_step)
    .method("step_fixed", &ode::test::Lorenz::ode_step_fixed)
    .method("advance",    &ode::test::Lorenz::ode_advance)
    .method("run",        &ode::test::Lorenz::ode_r_run)
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

  Rcpp::class_<util::Lookup>("Lookup")
    // No constructor, because class contains pure virtual method.
    .method("get_parameters", &util::Lookup::get_parameters)
    .method("set_parameters", &util::Lookup::set_parameters)
    ;

  Rcpp::class_<model::Strategy>("Strategy")
    .derives<util::Lookup>("Lookup")
    .constructor()
    ;

  Rcpp::class_<model::Parameters>("Parameters")
    .constructor()
    .derives<util::Lookup>("Lookup")
    .property("size",         &model::Parameters::size)
    .method("get_strategy",   &model::Parameters::get_strategy)
    .method("get_strategies", &model::Parameters::get_strategies)
    .method("add_strategy",   &model::Parameters::add_strategy)
    .method("set_strategy",   &model::Parameters::set_strategy)
    ;

  // Pending: get_values, set_values, get_rates (for ODE)
  //          derivs (for R ODE)
  Rcpp::class_<model::Plant>("Plant")
    .constructor<model::Strategy>()
    .method("set_mass_leaf",        &model::Plant::set_mass_leaf)
    // Leaf distribution (external interface)
    .method("leaf_area_above",      &model::Plant::leaf_area_above)
    .method("q",                    &model::Plant::q)
    .method("Q",                    &model::Plant::Q)
    .method("Qp",                   &model::Plant::Qp)
    // R specific access
    .property("parameters",         &model::Plant::r_get_parameters)
    .property("vars_size",          &model::Plant::r_get_vars_size)
    .property("vars_phys",          &model::Plant::r_get_vars_phys)
    .method("assimilation_leaf",    &model::Plant::assimilation_leaf)
    .method("compute_assimilation", &model::Plant::r_compute_assimilation)
    .method("compute_assimilation_x", &model::Plant::r_compute_assimilation_x)
    ;

  // Misc functions
  Rcpp::function("test_functor",    &util::test::test_functor);
  Rcpp::function("test_find_root",  &util::test::test_find_root);
  Rcpp::function("test_find_value", &util::test::test_find_value);
  Rcpp::function("test_integrator", &util::test::test_integrator);
}
