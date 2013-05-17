#include <Rcpp.h>

#include "spline.h"
#include "adaptive_spline.h"
#include "adaptive_spline_r.h"
#include "multi_spline.h"

#include "ode_target.h"
#include "lorenz.h"
#include "ode_r.h"

#include "lookup.h"

#include "strategy.h"
#include "parameters.h"
#include "control.h"

#include "disturbance.h"
#include "environment.h"

#include "plant.h"
#include "cohort_discrete.h"
#include "plant_spline.h"
#include "plant_approx.h"

#include "cohort_top.h"

#include "patch.h"

#include "metacommunity.h"

#include "functor.h"
#include "find_root.h"
#include "integrator.h"
#include "gradient.h"

RCPP_MODULE(tree) {
  Rcpp::class_<spline::Spline>("Spline")
    .constructor()
    .method("init",   &spline::Spline::init)
    .method("eval",   &spline::Spline::r_eval)
    .property("x",    &spline::Spline::r_get_x)
    .property("y",    &spline::Spline::r_get_y)
    .property("xy",   &spline::Spline::r_get_xy)
    .property("size", &spline::Spline::size)
    .property("min",  &spline::Spline::min)
    .property("max",  &spline::Spline::max)
    ;

  Rcpp::class_<spline::AdaptiveSpline>("AdaptiveSpline")
    .method("construct_spline", &spline::AdaptiveSpline::construct_spline)
    .method("set_control",      &spline::AdaptiveSpline::set_control)
    .method("eval_target",      &spline::AdaptiveSpline::eval_target)
    ;

  Rcpp::class_<spline::MultiSpline>("MultiSpline")
    .constructor<int>()
    .method("init",      &spline::MultiSpline::r_init)
    .method("init_self", &spline::MultiSpline::init_self)
    .method("add_point", &spline::MultiSpline::r_add_point)
    .method("reset",     &spline::MultiSpline::reset)
    .method("eval",      &spline::MultiSpline::r_eval)
    .method("eval_1",    &spline::MultiSpline::r_eval_1)
    .method("eval_r",    &spline::MultiSpline::r_eval_r)
    .property("x",       &spline::MultiSpline::r_get_x)
    .property("y",       &spline::MultiSpline::r_get_y)
    .property("size",    &spline::MultiSpline::size)
    .property("dim",     &spline::MultiSpline::dim)
    ;

  Rcpp::class_<ode::test::Lorenz>("Lorenz")
    .constructor<double,double,double>()
    .method("derivs",     &ode::test::Lorenz::r_derivs)
    .property("size",     &ode::test::Lorenz::size)
    // ODE solving
    .method("set_state",  &ode::test::Lorenz::ode_set_state)
    .property("state",    &ode::test::Lorenz::ode_get_state)
    .property("time",     &ode::test::Lorenz::ode_get_time)
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
    .property("state",    &ode::OdeR::ode_get_state)
    .property("time",     &ode::OdeR::ode_get_time)
    .method("step",       &ode::OdeR::ode_step)
    .method("step_fixed", &ode::OdeR::ode_step_fixed)
    .method("advance",    &ode::OdeR::ode_advance)
    .method("run",        &ode::OdeR::ode_r_run)
    ;

  Rcpp::class_<util::Lookup>("Lookup")
    // No constructor, because class contains pure virtual method.
    .property("parameters",   &util::Lookup::get_parameters)
    .method("set_parameters", &util::Lookup::set_parameters)
    ;

  Rcpp::class_<model::Strategy>("Strategy")
    .derives<util::Lookup>("Lookup")
    .constructor()
    .constructor<Rcpp::List>()
    .field("control", &model::Strategy::control)
    ;

  Rcpp::class_<model::Control>("Control")
    .derives<util::Lookup>("Lookup")
    .constructor()
    .constructor<Rcpp::List>()
    ;

  Rcpp::class_<model::Parameters>("Parameters")
    .constructor()
    .derives<util::Lookup>("Lookup")
    .property("size",         &model::Parameters::size)
    .method("[[",             &model::Parameters::r_at)
    .property("strategies",   &model::Parameters::r_get_strategies)
    .method("add_strategy",   &model::Parameters::add_strategy)
    .field("control",         &model::Parameters::control)
    ;

  Rcpp::class_<model::Disturbance>("Disturbance")
    .constructor()
    .constructor<double>()
    .method("survival_probability",
	    &model::Disturbance::survival_probability)
    ;

  Rcpp::class_<model::Plant>("Plant")
    .constructor<model::Strategy>()
    .derives<ode::OdeTarget>("OdeTarget")
    .method("set_mass_leaf",        &model::Plant::set_mass_leaf)
    .property("mass_leaf",          &model::Plant::mass_leaf)
    .property("height",             &model::Plant::get_height)
    .method("leaf_area_above",      &model::Plant::leaf_area_above)
    .method("compute_vars_phys",    &model::Plant::r_compute_vars_phys)
    .method("germination_probability", &model::Plant::r_germination_probability)
    .method("offspring",            &model::Plant::offspring)
    .method("died",                 &model::Plant::r_died)
    .property("survival_probability",  &model::Plant::survival_probability)
    .property("control",            &model::Plant::control)
    // R specific access
    .property("strategy",           &model::Plant::r_get_strategy)
    .property("vars_size",          &model::Plant::r_get_vars_size)
    .property("vars_phys",          &model::Plant::r_get_vars_phys)
    ;

  Rcpp::class_<model::PlantSpline>("PlantSpline")
    .constructor<model::Strategy,double,int>()
    .property("mass_leaf_max",   &model::PlantSpline::mass_leaf_max)
    .method("compute_vars_phys", &model::PlantSpline::r_compute_vars_phys)
    .method("ode_rates",         &model::PlantSpline::r_ode_rates)
    .property("plants",          &model::PlantSpline::r_get_plants)
    .property("plants_approx",   &model::PlantSpline::r_get_plants_approx)
    ;

  Rcpp::class_<model::CohortDiscrete>("CohortDiscrete")
    .derives<model::Plant>("Plant")
    .constructor<model::Strategy>()
    .constructor<model::Strategy, int>()
    .property("n_individuals", 
	      &model::CohortDiscrete::r_n_individuals,
	      &model::CohortDiscrete::r_set_n_individuals)
    ;

  Rcpp::class_<model::PlantApprox>("PlantApprox")
    .derives<model::Plant>("Plant")
    .constructor<model::Strategy, model::PlantSpline>()
    ;

  Rcpp::class_<model::CohortTop>("CohortTop")
    .constructor<model::Strategy>()
    .derives<model::Plant>("Plant")
    .method("compute_vars_phys_surv",
	    &model::CohortTop::r_compute_vars_phys_surv)
    .method("compute_initial_conditions",
	    &model::CohortTop::r_compute_initial_conditions)
    .method("growth_rate_gradient",
	    &model::CohortTop::r_growth_rate_gradient)
    .method("growth_rate_given_mass",
	    &model::CohortTop::r_growth_rate_given_mass)
    ;

  Rcpp::class_<model::SpeciesBase>("SpeciesBase")
    .derives<ode::OdeTarget>("OdeTarget")
    .property("size",          &model::SpeciesBase::size)
    .property("height_max",    &model::SpeciesBase::height_max)
    .method("leaf_area_above", &model::SpeciesBase::leaf_area_above)
    .method("compute_vars_phys", &model::SpeciesBase::r_compute_vars_phys)
    .method("add_seeds",       &model::SpeciesBase::add_seeds)
    .method("germination_probability", 
            &model::SpeciesBase::r_germination_probability)
    .method("clear",           &model::SpeciesBase::clear)
    .property("mass_leaf",     &model::SpeciesBase::r_get_mass_leaf)
    .method("set_mass_leaf",   &model::SpeciesBase::r_set_mass_leaf)
    .property("plants",        &model::SpeciesBase::r_get_plants)
    .property("n_individuals", &model::SpeciesBase::r_n_individuals)
    ;

  Rcpp::class_< model::Species<model::Plant> >("Species")
    .constructor<model::Strategy>()
    .derives<model::SpeciesBase>("SpeciesBase")
    .method("[[", &model::Species<model::Plant>::r_at)
    ;

  Rcpp::class_< model::Species<model::CohortDiscrete> >("SpeciesC")
    .constructor<model::Strategy>()
    .derives<model::SpeciesBase>("SpeciesBase")
    .method("[[", &model::Species<model::CohortDiscrete>::r_at)
    ;

  Rcpp::class_< model::Species<model::CohortTop> >("SpeciesCT")
    .constructor<model::Strategy>()
    .derives<model::SpeciesBase>("SpeciesBase")
    .method("[[", &model::Species<model::CohortTop>::r_at)
    ;


  Rcpp::class_<model::Environment>("Environment")
    .constructor<model::Disturbance, model::Control>()
    .method("canopy_openness", &model::Environment::canopy_openness)
    .method("patch_survival",  &model::Environment::patch_survival)
    .property("light_environment",
	      &model::Environment::get_light_environment,
	      &model::Environment::set_light_environment)
    .property("disturbance_regime",
	      &model::Environment::get_disturbance_regime)
    .property("age",
	      &model::Environment::get_age,
	      &model::Environment::set_age)
    ;

  Rcpp::class_<model::PatchBase>("PatchBase")
    .derives<ode::OdeTarget>("OdeTarget")
    .property("size",             &model::PatchBase::size)
    .property("height_max",       &model::PatchBase::r_height_max)
    .method("canopy_openness",    &model::PatchBase::r_canopy_openness)
    .method("compute_light_environment",
	    &model::PatchBase::r_compute_light_environment)
    .property("light_environment",
	      &model::PatchBase::r_light_environment)
    .method("compute_vars_phys",  &model::PatchBase::r_compute_vars_phys)
    .property("age",              &model::PatchBase::r_age)
    .method("germination",        &model::PatchBase::r_germination)
    .property("species",          &model::PatchBase::r_get_species)
    .method("add_seeds",          &model::PatchBase::r_add_seeds)
    .method("add_seedlings",      &model::PatchBase::r_add_seedlings)
    .property("mass_leaf",        &model::PatchBase::r_get_mass_leaf)
    .method("set_mass_leaf",      &model::PatchBase::r_set_mass_leaf)
    .property("n_individuals",    &model::PatchBase::r_n_individuals)
    .method("clear",              &model::PatchBase::clear)
    .method("step",               &model::PatchBase::r_step)
    .method("step_deterministic", &model::PatchBase::step_deterministic)
    .method("step_stochastic",    &model::PatchBase::r_step_stochastic)
    ;

  Rcpp::class_< model::Patch<model::Plant> >("Patch")
    .constructor<model::Parameters>()
    .derives<model::PatchBase>("PatchBase")
    .method("[[", &model::Patch<model::Plant>::r_at)
    ;

  Rcpp::class_< model::Patch<model::CohortDiscrete> >("PatchC")
    .constructor<model::Parameters>()
    .derives<model::PatchBase>("PatchBase")
    .method("[[", &model::Patch<model::CohortDiscrete>::r_at)
    ;

  Rcpp::class_<model::MetacommunityBase>("MetacommunityBase")
    .derives<ode::OdeTarget>("OdeTarget")
    .property("size",    &model::MetacommunityBase::size)
    .property("age",     &model::MetacommunityBase::r_age)
    .method("step",      &model::MetacommunityBase::r_step)
    .method("step_deterministic",
      	    &model::MetacommunityBase::step_deterministic)
    .method("step_stochastic",
      	    &model::MetacommunityBase::r_step_stochastic)
    // TODO: births & deaths?
    .method("clear",      &model::MetacommunityBase::clear)
    .method("add_seedlings", &model::MetacommunityBase::r_add_seedlings)
    .method("disperse", &model::MetacommunityBase::r_disperse)
    .property("n_individuals", &model::MetacommunityBase::r_n_individuals)
    .property("patches", &model::MetacommunityBase::r_get_patches)
    .property("mass_leaf",   &model::MetacommunityBase::r_get_mass_leaf)
    .method("set_mass_leaf", &model::MetacommunityBase::r_set_mass_leaf)
    ;

  Rcpp::class_< model::Metacommunity<model::Plant> >("Metacommunity")
    .constructor<model::Parameters>()
    .derives<model::MetacommunityBase>("MetacommunityBase")
    .method("[[", &model::Metacommunity<model::Plant>::r_at)
    ;

  Rcpp::class_< model::Metacommunity<model::CohortDiscrete> >("MetacommunityC")
    .constructor<model::Parameters>()
    .derives<model::MetacommunityBase>("MetacommunityBase")
    .method("[[", &model::Metacommunity<model::CohortDiscrete>::r_at)
    ;

  // Testing functions
  Rcpp::function("test_functor",    &util::test::test_functor);
  Rcpp::function("test_find_root",  &util::test::test_find_root);
  Rcpp::function("test_find_value", &util::test::test_find_value);
  Rcpp::function("test_integrator", &util::test::test_integrator);
  Rcpp::function("test_adaptive_spline", 
		 &spline::test::test_adaptive_spline);
  Rcpp::function("test_plant",      &model::test::test_plant);
  Rcpp::function("test_sum_double", &util::test::test_sum_double);
  Rcpp::function("test_sum_int",    &util::test::test_sum_int);
  Rcpp::function("test_to_rcpp_matrix",
		 &util::test::test_to_rcpp_matrix);
  Rcpp::function("test_from_rcpp_matrix",
		 &util::test::test_from_rcpp_matrix);
  Rcpp::function("test_gradient",   &util::test::test_gradient);
  Rcpp::function("test_gradient_richardson",
		 &util::test::test_gradient_richardson);

  // Useful functions
  Rcpp::function("set_sane_gsl_error_handling",
		 &util::set_sane_gsl_error_handling);
}
