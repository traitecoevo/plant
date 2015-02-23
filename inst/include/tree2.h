// -*-c++-*-
#ifndef _TREE2_H_
#define _TREE2_H_

#include <tree2/util.h>

#include <tree2/qk.h>
#include <tree2/qag.h>
#include <tree2/interpolator.h>
#include <tree2/adaptive_interpolator.h>

#include <tree2/ode_control.h>
#include <tree2/ode_step.h>
#include <tree2/ode_solver.h>
#include <tree2/ode_runner.h>

#include <tree2/disturbance.h>
#include <tree2/environment.h>

#include <tree2/control.h>
#include <tree2/strategy.h>
#include <tree2/parameters.h>
#include <tree2/cohort_schedule.h>

// Getting more serious down here.
#include <tree2/plant.h>
#include <tree2/plant_minimal.h>

#include <tree2/cohort.h>
#include <tree2/species.h>
#include <tree2/patch.h>
#include <tree2/ebt.h>

#include <tree2/plant_runner.h>

// Purely for testing
#include <tree2/lorenz.h>

namespace tree2 {
// Use this to manually switch between the minimal plant type and the
// full one.  Eventually we'll get this stuff exposed nicely.
//   typedef tree2::Plant PlantFF;
typedef tree2::PlantMinimal<tree2::Strategy> PlantFF;
}

// Include this early on.  It can be either after classes have been
// declared (but before Rcpp has been loaded) or first.  This file will
// attempt to provide declarations for the classes and namespaces that
// you use, but this might be fragile.
#include <tree2/RcppR6_pre.hpp>

// Anything after this point is OK to include Rcpp.h.  This is
// probably where the meat of the included material goes if your
// classes directly use Rcpp types.  Otherwise you can just declare
// them earlier up.

#include <Rcpp.h>
#include <tree2/ode_r.h>

// This line can safely be the last line in the file, but may go any
// point after RcppR6_pre.hpp is included.
#include <tree2/RcppR6_post.hpp>
#include <tree2/util_post_rcpp.h>
#include <tree2/tree_utils.h>
#include <tree2/tree_tools.h>

#endif
