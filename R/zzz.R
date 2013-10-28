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
##' @export Environment
##' @export Lookup
##' @export Lorenz
##' @export Metacommunity
##' @export MetacommunityBase
##' @export MetacommunityC
##' @export MultiSpline
##' @export OdeR
##' @export OdeTarget
##' @export Parameters
##' @export Patch
##' @export PatchBase
##' @export PatchC
##' @export PatchCohortTop
##' @export Plant
##' @export PlantApprox
##' @export PlantSpline
##' @export QAG
##' @export QK
##' @export RFunctionWrapper
##' @export Species
##' @export SpeciesBase
##' @export SpeciesC
##' @export SpeciesCT
##' @export Spline
##' @export Strategy
##' @export compute_assimilation_spline
##' @export make.reference
##' @export metacommunity
##' @export patch
##' @export set_sane_gsl_error_handling
##' @export species
##' @export test_adaptive_spline
##' @export test_find_root
##' @export test_find_value
##' @export test_from_rcpp_matrix
##' @export test_functor
##' @export test_gradient
##' @export test_gradient_richardson
##' @export test_plant
##' @export test_sum_double
##' @export test_sum_int
##' @export test_to_rcpp_matrix
##' @export trapezium

.onLoad <- function(libname, pkgname){
  ## loadRcppModules() # this breaks devtools
  ## But the line below causes CHECK to fail because it can't find
  ## set_sane_gsl_error_handling.
  loadModule("tree", TRUE)
  set_sane_gsl_error_handling()
}
