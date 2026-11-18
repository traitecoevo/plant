// -*-c++-*-
#ifndef _PLANT_H_
#define _PLANT_H_

// First, before any header that defines arithmetic: see fp_contract.h.
#include <plant/numerics/fp_contract.h>

#include <plant/util.h>

#include <plant/numerics/qk.h>
#include <plant/numerics/qag.h>


#include <odelia/ode_control.hpp>
#include <odelia/ode_step.hpp>
#include <odelia/ode_solver.hpp>

#include <plant/models/environment.h>
#include <plant/models/resource_spline.h>

#include <plant/control.h>
#include <plant/models/model.h>
#include <plant/parameters.h>
#include <plant/solvers/schedule.h>
#include <plant/solvers/events.h>

// Disturbance regimes for meta-population
#include <plant/disturbances/disturbance_regime.h>
#include <plant/disturbances/no_disturbance.h>
#include <plant/disturbances/weibull_disturbance.h>

// Specific models
#include <plant/models/ff16_model.h>
#include <plant/models/tf24_model.h>
#include <plant/models/tf24f_model.h>
#include <plant/models/tf24t_model.h>
#include <plant/models/tf24w_model.h>
#include <plant/models/k93_model.h>
#include <plant/models/es20_model.h>
#include <plant/models/leaf_model.h>

// Getting more serious down here.
#include <plant/community/individual.h>
#include <plant/community/internals.h>

#include <plant/community/node.h>
#include <plant/community/species.h>
#include <plant/community/patch.h>
#include <plant/solvers/deterministic_patch_runner.h>

// Stochastic model
#include <plant/community/stochastic_species.h>
#include <plant/community/stochastic_patch.h>
#include <plant/solvers/stochastic_patch_runner.h>

#include <plant/solvers/individual_runner.h>

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
// #include <plant/ode_solver/ode_r.h>

// This line can safely be the last line in the file, but may go any
// point after RcppR6_pre.hpp is included.
#include <plant/RcppR6_post.hpp>
#include <plant/util_post_rcpp.h>

// gperftools profiler
// Uncomment next line if you want to use the profiler.
// #include "gperftools/profiler.h"

#endif
