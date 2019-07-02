// -*-c++-*-
#ifndef _PLANT_H_
#define _PLANT_H_

#include <plant/util.h>

#include <plant/qk.h>
#include <plant/qag.h>
#include <plant/interpolator.h>
#include <plant/adaptive_interpolator.h>

#include <plant/ode_control.h>
#include <plant/ode_step.h>
#include <plant/ode_solver.h>
#include <plant/ode_runner.h>

#include <plant/disturbance.h>
#include <plant/environment.h>

#include <plant/control.h>
#include <plant/ff16_strategy.h>
#include <plant/parameters.h>
#include <plant/cohort_schedule.h>

// new physiology
// #include <plant/ff16r_strategy.h>

// Getting more serious down here.
#include <plant/plant.h>
// #include <plant/plant_plus.h>

#include <plant/cohort.h>
#include <plant/species.h>
#include <plant/patch.h>
#include <plant/scm.h>

// Stochastic model
#include <plant/stochastic_species.h>
#include <plant/stochastic_patch.h>
#include <plant/stochastic_patch_runner.h>

#include <plant/plant_runner.h>

// Purely for testing
#include <plant/lorenz.h>

// Include this early on.  It can be either after classes have been
// declared (but before Rcpp has been loaded) or first.  This file will
// attempt to provide declarations for the classes and namespaces that
// you use, but this might be fragile.
#include <plant/RcppR6_pre.hpp>

// Anything after this point is OK to include Rcpp.h.  This is
// probably where the meat of the included material goes if your
// classes directly use Rcpp types.  Otherwise you can just declare
// them earlier up.

#include <Rcpp.h>
#include <plant/ode_r.h>

// This line can safely be the last line in the file, but may go any
// point after RcppR6_pre.hpp is included.
#include <plant/RcppR6_post.hpp>
#include <plant/util_post_rcpp.h>
#include <plant/get_state.h>
#include <plant/plant_tools.h>

#endif
