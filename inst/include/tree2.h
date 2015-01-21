// -*-c++-*-
#ifndef _TREE2_H_
#define _TREE2_H_

// RcppR6 isn't going to do a good job of forward declaring these, so
// I'll need to.

namespace tree2 {
template <typename T> class Cohort;
template <typename T> class Species;
template <typename T> class Patch;
template <typename T> class EBT;
}
namespace ode {
template <typename T> class Runner;
}

#include <tree2/util.h>

// Include this early on.  It can be either after classes have been
// declared (but before Rcpp has been loaded) or first.  This file will
// attempt to provide declarations for the classes and namespaces that
// you use, but this might be fragile.
#include <tree2/RcppR6_pre.hpp>

// Hmm - I think we can actually include Rcpp here by this point
// anyway.  No big loss.

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
#include <tree2/cohort.h>
#include <tree2/species.h>
#include <tree2/patch.h>
#include <tree2/ebt.h>

// Anything after this point is OK to include Rcpp.h.  This is
// probably where the meat of the included material goes if your
// classes directly use Rcpp types.  Otherwise you can just declare
// them earlier up.

#include <Rcpp.h>
#include <tree2/lorenz.h>
#include <tree2/ode_r.h>

// This line can safely be the last line in the file, but may go any
// point after RcppR6_pre.hpp is included.
#include <tree2/RcppR6_post.hpp>
#include <tree2/util_rcpp.h>

#endif
