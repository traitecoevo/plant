##' @useDynLib tree
##' @import Rcpp
##' @import methods

##' @export CohortDiscrete
##' @export CohortSchedule
##' @export CohortScheduleEvent
##' @export CohortTop
##' @export Control
##' @export Disturbance
##' @export EBT
##' @export EBTMutantRunner
##' @export Environment
##' @export FakeLightEnvironment
##' @export Interpolator
##' @export Lookup
##' @export Lorenz
##' @export Metacommunity
##' @export MetacommunityBase
##' @export MetacommunityC
##' @export OdeR
##' @export OdeTarget
##' @export Parameters
##' @export Patch
##' @export PatchBase
##' @export PatchC
##' @export PatchCohortTop
##' @export Plant
##' @export QAG
##' @export QK
##' @export RFunctionWrapper
##' @export Species
##' @export SpeciesBase
##' @export SpeciesC
##' @export SpeciesCT
##' @export Strategy
##' @export compute_assimilation_fn
##' @export local_error_integration
##' @export make.reference
##' @export set_sane_gsl_error_handling
##' @export test_adaptive_interpolator
##' @export test_find_root
##' @export test_find_value
##' @export test_from_rcpp_integer_matrix
##' @export test_from_rcpp_numeric_matrix
##' @export test_functor
##' @export test_gradient
##' @export test_gradient_richardson
##' @export test_plant
##' @export test_sum_double
##' @export test_sum_int
##' @export test_to_rcpp_integer_matrix
##' @export test_to_rcpp_numeric_matrix
##' @export trapezium
##' @export trapezium_vector

.onLoad <- function(libname, pkgname){
  ## loadRcppModules() # this breaks devtools
  ## But the line below causes CHECK to fail because it can't find
  ## set_sane_gsl_error_handling.
  loadModule("tree", TRUE)
  set_sane_gsl_error_handling()
}
