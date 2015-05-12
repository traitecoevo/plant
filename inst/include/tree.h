// -*-c++-*-
#ifndef _TREE_H_
#define _TREE_H_

#include <tree/util.h>

#include <tree/qk.h>
#include <tree/qag.h>
#include <tree/interpolator.h>
#include <tree/adaptive_interpolator.h>

#include <tree/ode_control.h>
#include <tree/ode_step.h>
#include <tree/ode_solver.h>
#include <tree/ode_runner.h>

#include <tree/disturbance.h>
#include <tree/environment.h>

#include <tree/control.h>
#include <tree/ffw16_strategy.h>
#include <tree/parameters.h>
#include <tree/cohort_schedule.h>

// Getting more serious down here.
#include <tree/ffw16_plant_plus.h>
#include <tree/plant.h>

#include <tree/cohort.h>
#include <tree/species.h>
#include <tree/patch.h>
#include <tree/ebt.h>

#include <tree/plant_runner.h>

// Purely for testing
#include <tree/lorenz.h>

namespace tree {
// Use this to manually switch between the minimal plant type and the
// full one.  Eventually we'll get this stuff exposed nicely.
typedef Parameters<FFW16_Strategy> FFW16_Parameters;
typedef Plant<FFW16_Strategy>      FFW16_Plant;
typedef Cohort<FFW16_Plant>        FFW16_Cohort;
typedef Species<FFW16_Plant>       FFW16_Species;
typedef Patch<FFW16_Plant>         FFW16_Patch;
typedef EBT<FFW16_Plant>           FFW16_EBT;
}

// Include this early on.  It can be either after classes have been
// declared (but before Rcpp has been loaded) or first.  This file will
// attempt to provide declarations for the classes and namespaces that
// you use, but this might be fragile.
#include <tree/RcppR6_pre.hpp>

// Anything after this point is OK to include Rcpp.h.  This is
// probably where the meat of the included material goes if your
// classes directly use Rcpp types.  Otherwise you can just declare
// them earlier up.

#include <Rcpp.h>
#include <tree/ode_r.h>

// This line can safely be the last line in the file, but may go any
// point after RcppR6_pre.hpp is included.
#include <tree/RcppR6_post.hpp>
#include <tree/util_post_rcpp.h>
#include <tree/tree_utils.h>
#include <tree/tree_tools.h>

#endif
