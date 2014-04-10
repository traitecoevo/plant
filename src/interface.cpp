#include <Rcpp.h>

#include "interpolator.h"
#include "adaptive_interpolator.h"
#include "fake_light_environment.h"

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

#include "cohort_top.h"
#include "cohort_schedule.h"

#include "patch.h"

#include "metacommunity.h"
#include "ebt.h"
#include "ebt_mutant_runner.h"

#include "functor.h"
#include "find_root.h"
#include "quadrature.h"
#include "integration.h"
#include "gradient.h"
#include "state.h"

#ifdef __clang__
#pragma clang diagnostic push
// These I have no control over because they're Rcpp issues.
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wmissing-prototypes"
#endif
RCPP_MODULE(tree) {
#ifdef __clang__
#pragma clang diagnostic pop
#endif
  Rcpp::class_<interpolator::Interpolator>("Interpolator")
    .constructor()
    .constructor<bool,bool>()
    .property("type",     &interpolator::Interpolator::type)
    .method("init",       &interpolator::Interpolator::init)
    .method("eval",       &interpolator::Interpolator::r_eval)
    .method("deriv",      &interpolator::Interpolator::r_deriv)
    .property("x",        &interpolator::Interpolator::get_x)
    .property("y",        &interpolator::Interpolator::get_y)
    .property("xy",       &interpolator::Interpolator::r_get_xy)
    .property("size",     &interpolator::Interpolator::size)
    .property("min",      &interpolator::Interpolator::min)
    .property("max",      &interpolator::Interpolator::max)
    ;

  Rcpp::class_<interpolator::FakeLightEnvironment>("FakeLightEnvironment")
    .constructor<std::vector<double>, Rcpp::List>()
    .method("get",  &interpolator::FakeLightEnvironment::operator())
    ;

  Rcpp::class_<ode::test::Lorenz>("Lorenz")
    .constructor<double,double,double>()
    .method("derivs",        &ode::test::Lorenz::r_derivs)
    .property("size",        &ode::test::Lorenz::size)
    // ODE solving
    .method("set_state",     &ode::test::Lorenz::set_ode_state)
    .property("state",       &ode::test::Lorenz::ode_state)
    .property("time",        &ode::test::Lorenz::get_time)
    .property("times",       &ode::test::Lorenz::get_times)
    .method("step",          &ode::test::Lorenz::step)
    .method("step_fixed",    &ode::test::Lorenz::step_fixed)
    .method("step_to",       &ode::test::Lorenz::step_to)
    .method("advance",       &ode::test::Lorenz::advance)
    .method("advance_fixed", &ode::test::Lorenz::advance_fixed)
    .method("run",           &ode::test::Lorenz::r_run)
    ;

  Rcpp::class_<ode::OdeTarget>("OdeTarget")
    .method("derivs",         &ode::OdeTarget::r_derivs)
    .property("ode_size",     &ode::OdeTarget::ode_size)
    .property("ode_values",   &ode::OdeTarget::r_ode_values)
    .method("set_ode_values", &ode::OdeTarget::r_set_ode_values)
    .property("ode_rates",    &ode::OdeTarget::r_ode_rates)
    ;

  Rcpp::class_<ode::OdeR>("OdeR")
    .constructor<SEXP,SEXP,SEXP>()
    .constructor<SEXP,SEXP,SEXP,ode::OdeControl>()
    .property("size",        &ode::OdeR::size)
    // ODE solving
    .method("set_state",     &ode::OdeR::set_ode_state)
    .property("state",       &ode::OdeR::ode_state)
    .property("time",        &ode::OdeR::get_time)
    .property("times",       &ode::OdeR::get_times)
    .method("step",          &ode::OdeR::step)
    .method("step_fixed",    &ode::OdeR::step_fixed)
    .method("step_to",       &ode::OdeR::step_to)
    .method("advance",       &ode::OdeR::advance)
    .method("advance_fixed", &ode::OdeR::advance_fixed)
    .method("reset",         &ode::OdeR::reset)
    .method("derivs",        &ode::OdeR::r_derivs)
    .method("run",           &ode::OdeR::r_run)
    .property("control",     &ode::OdeR::r_control)
    .property("pars",        &ode::OdeR::r_pars)
    ;

  Rcpp::class_<util::Lookup>("Lookup")
    // NOTE: Parameters not set through .property because they can set
    // just a subset through this, which is poor semantics for '<-'.
    .property("parameters",   &util::Lookup::get_parameters)
    .method("set_parameters", &util::Lookup::set_parameters)
    .method("has_key",        &util::Lookup::has_key)
    ;

  Rcpp::class_<model::Strategy>("Strategy")
    .derives<util::Lookup>("Lookup")
    .constructor()
    .constructor<Rcpp::List>()
    .property("control", &model::Strategy::r_control,
	      &model::Strategy::set_control)
    .property("assimilation_fn",
	      &model::Strategy::r_assimilation_fn,
	      &model::Strategy::r_set_assimilation_fn)
    .property("integrator", &model::Strategy::r_integrator)
    .method("copy", &model::Strategy::r_copy)
    ;

  Rcpp::class_<model::Control>("Control")
    .derives<util::Lookup>("Lookup")
    .constructor()
    .constructor<Rcpp::List>()
    .property("ode_control", &model::Control::r_ode_control)
    ;

  Rcpp::class_<ode::OdeControl>("OdeControl")
    .derives<util::Lookup>("Lookup")
    .constructor()
    ;

  Rcpp::class_<model::Parameters>("Parameters")
    .constructor()
    .derives<util::Lookup>("Lookup")
    .method("clear",          &model::Parameters::clear)
    .property("size",         &model::Parameters::size)
    .property("n_residents",  &model::Parameters::n_residents)
    .property("n_mutants",    &model::Parameters::n_mutants)
    .method("[[",             &model::Parameters::r_at)
    .property("strategies",   &model::Parameters::r_get_strategies)
    .method("add_strategy",   &model::Parameters::add_strategy)
    .method("add_strategy_mutant",
	    &model::Parameters::add_strategy_mutant)
    .field("control",         &model::Parameters::control)
    .field("disturbance",     &model::Parameters::disturbance)
    .method("set_control_parameters",
	    &model::Parameters::set_control_parameters)
    .method("copy",           &model::Parameters::r_copy)
    .property("seed_rain",    &model::Parameters::r_seed_rain,
	      &model::Parameters::r_set_seed_rain)
    .property("is_resident",    &model::Parameters::r_is_resident,
	      &model::Parameters::r_set_is_resident)
    ;

  Rcpp::class_<model::Disturbance>("Disturbance")
    .constructor<double>()
    .method("density",         &model::Disturbance::density)
    .method("pr_survival",     &model::Disturbance::pr_survival)
    .method("pr_survival_conditional",
	    &model::Disturbance::pr_survival_conditional)
    .method("cdf",             &model::Disturbance::cdf)
    .property("mean_interval", &model::Disturbance::r_mean_interval)
    ;

  Rcpp::class_<model::Plant>("Plant")
    .constructor<model::Strategy>()
    .derives<ode::OdeTarget>("OdeTarget")
    .property("height",
	      &model::Plant::height,
	      &model::Plant::set_height)
    .property("heartwood_area",
          &model::Plant::heartwood_area,
          &model::Plant::set_area_heartwood)
    .property("heartwood_mass",
          &model::Plant::mass_heartwood,
          &model::Plant::set_mass_heartwood)

    .property("fecundity",          &model::Plant::fecundity)
    .property("leaf_area",          &model::Plant::leaf_area)
    .method("leaf_area_above",      &model::Plant::leaf_area_above)
    .method("compute_vars_phys",    &model::Plant::compute_vars_phys)
    .method("germination_probability", &model::Plant::germination_probability)
    .method("offspring",            &model::Plant::offspring)
    .method("died",                 &model::Plant::r_died)
    .property("survival_probability",  &model::Plant::survival_probability)
    .method("assimilation_given_height",
	    &model::Plant::assimilation_given_height)
    // R specific access
    .property("strategy",           &model::Plant::r_get_strategy)
    .property("vars_size",          &model::Plant::r_get_vars_size)
    .property("vars_phys",          &model::Plant::r_get_vars_phys)
    .property("vars_growth_decomp", &model::Plant::r_get_vars_growth_decomp)
    // State
    .property("state_size",         &model::Plant::state_size)
    .property("state",
	      &util::r_get_state<model::Plant>,
	      &util::r_set_state<model::Plant>)
    .method("copy",                 &model::Plant::r_copy)
    ;

  Rcpp::class_<model::CohortDiscrete>("CohortDiscrete")
    .derives<model::Plant>("Plant")
    .constructor<model::Strategy>()
    .constructor<model::Strategy, int>()
    .property("n_individuals",
	      &model::CohortDiscrete::get_n_individuals,
	      &model::CohortDiscrete::set_n_individuals)
    ;

  Rcpp::class_<model::CohortTop>("CohortTop")
    .constructor<model::Strategy>()
    .derives<model::Plant>("Plant")
    .method("compute_vars_phys",
	    &model::CohortTop::compute_vars_phys)
    .method("compute_initial_conditions",
	    &model::CohortTop::compute_initial_conditions)
    .method("growth_rate_gradient",
	    &model::CohortTop::r_growth_rate_gradient)
    .method("growth_rate_given_height",
	    &model::CohortTop::r_growth_rate_given_height)
    ;

  Rcpp::class_<model::CohortSchedule>("CohortSchedule")
    .constructor<size_t>()
    .property("size",       &model::CohortSchedule::size)
    .property("n_species",  &model::CohortSchedule::get_n_species)
    .method("expand",       &model::CohortSchedule::expand)
    .method("clear_times",  &model::CohortSchedule::r_clear_times)
    .method("set_times",    &model::CohortSchedule::r_set_times)
    .method("times",        &model::CohortSchedule::r_times)
    .method("reset",        &model::CohortSchedule::reset)
    .method("pop",          &model::CohortSchedule::pop)
    .property("next_event", &model::CohortSchedule::next_event)
    .property("remaining",  &model::CohortSchedule::remaining)
    .property("fixed_times",&model::CohortSchedule::fixed_times)
    .property("max_time",   &model::CohortSchedule::r_max_time,
	      &model::CohortSchedule::r_set_max_time)
    .property("ode_times",  &model::CohortSchedule::r_ode_times,
	      &model::CohortSchedule::r_set_ode_times)
    .method("clear_ode_times",
	    &model::CohortSchedule::r_clear_ode_times)
    .property("state",
	      &model::CohortSchedule::r_get_state,
	      &model::CohortSchedule::r_set_state)
    .property("all_times",
              &model::CohortSchedule::r_all_times,
              &model::CohortSchedule::r_set_all_times)
    .method("copy",         &model::CohortSchedule::r_copy)
    ;
  Rcpp::class_<model::CohortSchedule::Event>("CohortScheduleEvent")
    .constructor<double,int>()
    .property("species_index",
	      &model::CohortSchedule::Event::r_species_index)
    .field_readonly("times",
		    &model::CohortSchedule::Event::times)
    .property("time_introduction",
	      &model::CohortSchedule::Event::time_introduction)
    .property("time_end",
	      &model::CohortSchedule::Event::time_end)
    ;

  Rcpp::class_<model::SpeciesBase>("SpeciesBase")
    .derives<ode::OdeTarget>("OdeTarget")
    .property("size",          &model::SpeciesBase::size)
    .property("height_max",    &model::SpeciesBase::height_max)
    .method("leaf_area_above", &model::SpeciesBase::leaf_area_above)
    .property("seeds",         &model::SpeciesBase::seeds)
    .method("compute_vars_phys", &model::SpeciesBase::compute_vars_phys)
    .method("add_seeds",       &model::SpeciesBase::add_seeds)
    .method("germination_probability",
            &model::SpeciesBase::germination_probability)
    .method("clear",           &model::SpeciesBase::clear)
    .property("height",        &model::SpeciesBase::r_height,
	      &model::SpeciesBase::r_set_height)
    .property("plants",        &model::SpeciesBase::r_get_plants)
    .property("n_individuals", &model::SpeciesBase::r_n_individuals)
    .method("compute_assimilation_fn",
	    &model::SpeciesBase::compute_assimilation_fn)
    .method("rescale_assimilation_fn",
	    &model::SpeciesBase::rescale_assimilation_fn)
    .property("assimilation_fn",
	      &model::SpeciesBase::r_assimilation_fn)
    .property("strategy",      &model::SpeciesBase::r_strategy)
    .property("leaf_area",     &model::SpeciesBase::r_leaf_area)
    .property("vars_size",     &model::SpeciesBase::r_get_vars_size)
    .property("vars_phys",     &model::SpeciesBase::r_get_vars_phys)
    .property("state",
	      &model::SpeciesBase::r_get_state,
	      &model::SpeciesBase::r_set_state)
    .method("force_state",
	    &model::SpeciesBase::r_force_state)
    ;

  Rcpp::class_< model::Species<model::Plant> >("Species")
    .constructor<model::Strategy>()
    .derives<model::SpeciesBase>("SpeciesBase")
    .method("[[", &model::Species<model::Plant>::r_at)
    .property("seed", &model::Species<model::Plant>::r_seed)
    ;

  Rcpp::class_< model::Species<model::CohortDiscrete> >("SpeciesC")
    .constructor<model::Strategy>()
    .derives<model::SpeciesBase>("SpeciesBase")
    .method("[[", &model::Species<model::CohortDiscrete>::r_at)
    .property("seed", &model::Species<model::CohortDiscrete>::r_seed)
    ;

  Rcpp::class_< model::Species<model::CohortTop> >("SpeciesCT")
    .constructor<model::Strategy>()
    .derives<model::SpeciesBase>("SpeciesBase")
    .method("[[", &model::Species<model::CohortTop>::r_at)
    .property("seed", &model::Species<model::CohortTop>::r_seed)
    .method("leaf_area_error",
	    &model::Species<model::CohortTop>::leaf_area_error)
    ;

  Rcpp::class_<model::Environment>("Environment")
    .constructor<model::Parameters>()
    .method("canopy_openness",   &model::Environment::canopy_openness)
    .property("patch_survival",  &model::Environment::patch_survival)
    .method("patch_survival_conditional",
	    &model::Environment::patch_survival_conditional)
    .property("seed_rain_rate",  &model::Environment::seed_rain_rate)
    .property("disturbance_regime",
	      &model::Environment::get_disturbance_regime)
    .property("time",
	      &model::Environment::get_time,
	      &model::Environment::set_time)
    .property("light_environment",
	      &model::Environment::get_light_environment,
	      &model::Environment::set_light_environment)
    .property("seed_rain",       &model::Environment::r_get_seed_rain)
    .method("set_seed_rain_index",
	    &model::Environment::r_set_seed_rain_index)
    .property("state",
	      &model::Environment::r_get_state,
	      &model::Environment::r_set_state)
    ;

  Rcpp::class_<model::PatchBase>("PatchBase")
    .derives<ode::OdeTarget>("OdeTarget")
    .method("births",             &model::PatchBase::births)
    .method("deaths",             &model::PatchBase::deaths)
    .property("size",             &model::PatchBase::size)
    .property("time",             &model::PatchBase::get_time,
	      &model::PatchBase::set_time)
    .property("disturbance_regime",
	      &model::PatchBase::get_disturbance_regime)
    .property("height_max",       &model::PatchBase::r_height_max)
    .method("canopy_openness",    &model::PatchBase::canopy_openness)
    .method("leaf_area_above",    &model::PatchBase::leaf_area_above)
    .method("compute_light_environment",
	    &model::PatchBase::r_compute_light_environment)
    .property("environment",      &model::PatchBase::r_environment)
    .method("compute_vars_phys",  &model::PatchBase::r_compute_vars_phys)
    .method("germination",        &model::PatchBase::r_germination)
    .property("species",          &model::PatchBase::r_get_species)
    .method("add_seeds",          &model::PatchBase::r_add_seeds)
    .method("add_seedling",       &model::PatchBase::r_add_seedling)
    .method("add_seedlings",      &model::PatchBase::r_add_seedlings)
    .property("height",           &model::PatchBase::r_height,
	      &model::PatchBase::r_set_height)
    .property("n_individuals",    &model::PatchBase::r_n_individuals)
    .method("reset",              &model::PatchBase::reset)
    .property("parameters",
	      &model::PatchBase::r_parameters)
    .property("state",
	      &model::PatchBase::r_get_state,
	      &model::PatchBase::r_set_state)
    .method("force_state",
	    &model::PatchBase::r_force_state)
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

  Rcpp::class_< model::Patch<model::CohortTop> >("PatchCohortTop")
    .constructor<model::Parameters>()
    .derives<model::PatchBase>("PatchBase")
    .method("[[", &model::Patch<model::CohortTop>::r_at)
    ;

  Rcpp::class_<model::MetacommunityBase>("MetacommunityBase")
    .derives<ode::OdeTarget>("OdeTarget")
    .property("size",    &model::MetacommunityBase::size)
    .property("time",    &model::MetacommunityBase::get_time)
    .method("step",      &model::MetacommunityBase::r_step)
    .method("step_deterministic",
      	    &model::MetacommunityBase::step_deterministic)
    .method("step_stochastic",
      	    &model::MetacommunityBase::r_step_stochastic)
    .method("reset",      &model::MetacommunityBase::reset)
    .method("add_seedlings", &model::MetacommunityBase::r_add_seedlings)
    .method("disperse", &model::MetacommunityBase::r_disperse)
    .property("n_individuals", &model::MetacommunityBase::r_n_individuals)
    .property("patches", &model::MetacommunityBase::r_get_patches)
    .property("height",   &model::MetacommunityBase::r_height,
	      &model::MetacommunityBase::r_set_height)
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

  Rcpp::class_<model::EBT>("EBT")
    .derives<ode::OdeTarget>("OdeTarget")
    .constructor<model::Parameters>()
    .method("fitness",           &model::EBT::r_fitness)
    .property("fitnesses",       &model::EBT::fitnesses)
    .method("fitness_error",     &model::EBT::r_fitness_error)
    .method("fitness_cohort",    &model::EBT::r_fitness_cohort)
    .method("leaf_area_error",   &model::EBT::r_leaf_area_error)
    .property("patch",           &model::EBT::r_patch)
    .property("cohort_schedule", &model::EBT::r_cohort_schedule,
	      &model::EBT::r_set_cohort_schedule)
    .property("time",            &model::EBT::get_time)
    .method("reset",             &model::EBT::reset)
    .property("complete",        &model::EBT::complete)
    .method("run",               &model::EBT::run)
    .method("run_next",          &model::EBT::run_next)
    .property("ode_times",       &model::EBT::r_ode_times)
    .property("parameters",      &model::EBT::r_parameters)
    .method("times",             &model::EBT::r_times)
    .method("set_times",         &model::EBT::r_set_times)
    .property("state",           &model::EBT::r_get_state,
	      &model::EBT::r_set_state)
    ;

  Rcpp::class_<model::EBTMutantRunner>("EBTMutantRunner")
    .derives<model::EBT>("EBT")
    .constructor<model::Parameters, interpolator::FakeLightEnvironment>()
    ;

  Rcpp::class_<util::RFunctionWrapper>("RFunctionWrapper")
    .constructor<SEXP, SEXP>()
    .method("target", &util::RFunctionWrapper::target)
    ;

  Rcpp::class_<integration::QK>("QK")
    .constructor<size_t>()
    .method("integrate",         &integration::QK::r_integrate)
    .property("last_area",       &integration::QK::get_last_area)
    .property("last_error",      &integration::QK::get_last_error)
    .property("last_area_abs",   &integration::QK::get_last_area_abs)
    .property("last_area_asc",   &integration::QK::get_last_area_asc)
    .method("integrate_vector_x",&integration::QK::integrate_vector_x)
    .method("integrate_vector",  &integration::QK::integrate_vector)
    ;

  Rcpp::class_<integration::QAG>("QAG")
    .constructor<size_t, size_t, double, double>()
    .constructor<size_t>()
    .method("integrate",         &integration::QAG::r_integrate)
    .method("integrate_with_intervals",
	    &integration::QAG::r_integrate_with_intervals)
    .property("last_area",       &integration::QAG::get_last_area)
    .property("last_error",      &integration::QAG::get_last_error)
    .property("last_iterations", &integration::QAG::get_last_iterations)
    .property("last_intervals",  &integration::QAG::get_last_intervals)
    .property("is_adaptive",     &integration::QAG::is_adaptive)
    ;

  // Template helper functions
  Rcpp::function("species",       &model::species);
  Rcpp::function("patch",         &model::patch);
  Rcpp::function("metacommunity", &model::metacommunity);

  // Testing functions
  Rcpp::function("test_functor",    &util::test::test_functor);
  Rcpp::function("test_find_root",  &util::test::test_find_root);
  Rcpp::function("test_find_value", &util::test::test_find_value);
  Rcpp::function("test_adaptive_interpolator",
		 &interpolator::test::test_adaptive_interpolator);
  Rcpp::function("test_plant",      &model::test::test_plant);
  Rcpp::function("compute_assimilation_fn",
		 &model::test::compute_assimilation_fn);
  Rcpp::function("test_sum_double", &util::test::test_sum_double);
  Rcpp::function("test_sum_int",    &util::test::test_sum_int);
  Rcpp::function("test_to_rcpp_integer_matrix",
		 &util::test::test_to_rcpp_integer_matrix);
  Rcpp::function("test_from_rcpp_integer_matrix",
		 &util::test::test_from_rcpp_integer_matrix);
  Rcpp::function("test_to_rcpp_numeric_matrix",
		 &util::test::test_to_rcpp_numeric_matrix);
  Rcpp::function("test_from_rcpp_numeric_matrix",
		 &util::test::test_from_rcpp_numeric_matrix);
  Rcpp::function("test_gradient",   &util::test::test_gradient);
  Rcpp::function("test_gradient_richardson",
		 &util::test::test_gradient_richardson);

  // Useful functions
  Rcpp::function("set_sane_gsl_error_handling",
		 &util::set_sane_gsl_error_handling);
  Rcpp::function("trapezium",
		 &util::trapezium< std::vector<double>, std::vector<double> >);
  Rcpp::function("trapezium_vector",
		 &util::trapezium_vector< std::vector<double>,
		 std::vector<double> >);
  Rcpp::function("local_error_integration",
		 &util::local_error_integration);
}
