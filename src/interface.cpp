#include <Rcpp.h>

#include "spline.h"
#include "adaptive_spline.h"
#include "adaptive_spline_r.h"

#include "ode_target.h"
#include "lorenz.h"
#include "ode_r.h"

#include "lookup.h"

#include "strategy.h"
#include "parameters.h"
#include "plant.h"
#include "cohort_discrete.h"
#include "cohort.h"
#include "patch.h"

#include "functor.h"
#include "find_root.h"
#include "integrator.h"

RCPP_MODULE(tree) {
  Rcpp::class_<spline::Spline>("Spline")
    .constructor()
    .method("init",   &spline::Spline::init)
    .method("eval",   &spline::Spline::r_eval)
    .property("xy",   &spline::Spline::r_get_xy)
    .property("size", &spline::Spline::size)
    ;

  Rcpp::class_<spline::AdaptiveSpline>("AdaptiveSpline")
    .constructor()
    .derives<spline::Spline>("Spline")
    .method("construct_spline", &spline::AdaptiveSpline::construct_spline)
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

  Rcpp::class_<ode::OdeTarget>("OdeTarget")
    .method("derivs",         &ode::OdeTarget::r_derivs)
    .property("ode_size",     &ode::OdeTarget::ode_size)
    .method("ode_values_set", &ode::OdeTarget::r_ode_values_set)
    .property("ode_values",   &ode::OdeTarget::r_ode_values)
    .property("ode_rates",    &ode::OdeTarget::r_ode_rates)
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
    .constructor<Rcpp::List>()
    ;

  Rcpp::class_<model::Parameters>("Parameters")
    .constructor()
    .derives<util::Lookup>("Lookup")
    .property("size",         &model::Parameters::size)
    .method("get_strategy",   &model::Parameters::r_get_strategy)
    .method("get_strategies", &model::Parameters::r_get_strategies)
    .method("add_strategy",   &model::Parameters::add_strategy)
    ;

  Rcpp::class_<model::Plant>("Plant")
    .constructor<model::Strategy>()
    .derives<ode::OdeTarget>("OdeTarget")
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
    .method("compute_vars_phys",    &model::Plant::r_compute_vars_phys)
    .method("germination_probability", &model::Plant::r_germination_probability)
    .method("offspring",            &model::Plant::offspring)
    .method("died",                 &model::Plant::r_died)
    .property("name",               &model::Plant::r_name)
    ;

  Rcpp::class_<model::CohortDiscrete>("CohortDiscrete")
    .derives<model::Plant>("Plant")
    .constructor<model::Strategy>()
    .constructor<model::Strategy, int>()
    .property("n_individuals", 
	      &model::CohortDiscrete::r_n_individuals,
	      &model::CohortDiscrete::r_set_n_individuals)
    ;

  Rcpp::class_<model::Cohort>("Cohort")
    .constructor<model::Strategy>()
    .derives<ode::OdeTarget>("OdeTarget")
    ;

  Rcpp::class_<model::PatchBase>("PatchBase")
    .derives<ode::OdeTarget>("OdeTarget")
    .property("size",             &model::PatchBase::r_size)
    .property("height_max",       &model::PatchBase::r_height_max)
    .method("canopy_openness",    &model::PatchBase::r_canopy_openness)
    .method("compute_light_environment",
	    &model::PatchBase::r_compute_light_environment)
    .property("light_environment",
	      &model::PatchBase::r_light_environment)
    .method("compute_vars_phys",  &model::PatchBase::r_compute_vars_phys)
    .property("age",              &model::PatchBase::r_age)
    .method("germination",        &model::PatchBase::r_germination)
    .method("get_plants",         &model::PatchBase::r_get_plants)
    .method("add_seeds",          &model::PatchBase::r_add_seeds)
    .method("get_mass_leaf",      &model::PatchBase::r_get_mass_leaf)
    .method("set_mass_leaf",      &model::PatchBase::r_set_mass_leaf)
    .method("clear",              &model::PatchBase::r_clear)
    .method("step",               &model::PatchBase::r_step)
    .method("step_deterministic", &model::PatchBase::step_deterministic)
    .method("step_stochastic",    &model::PatchBase::r_step_stochastic)
    ;

  Rcpp::class_< model::Patch<model::Plant> >("Patch")
    .constructor<model::Parameters>()
    .derives<model::PatchBase>("PatchBase")
    ;

  // Misc functions
  Rcpp::function("test_functor",    &util::test::test_functor);
  Rcpp::function("test_find_root",  &util::test::test_find_root);
  Rcpp::function("test_find_value", &util::test::test_find_value);
  Rcpp::function("test_integrator", &util::test::test_integrator);
  Rcpp::function("test_adaptive_spline", 
		 &spline::test::test_adaptive_spline);
  Rcpp::function("test_plant",      &model::test::test_plant);
}
