## Load the Rcpp module.  Automatically added to a load hook, so no
## need to add via .onLoad().
## These lines get processed by Roxygen
##' @useDynLib tree
##' @import Rcpp
##' @import methods
##' @export AdaptiveSpline
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
##' @export SeedRain
##' @export Species
##' @export SpeciesBase
##' @export SpeciesC
##' @export SpeciesCT
##' @export Spline
##' @export Strategy
##' @export make.reference
##' @export metacommunity
##' @export patch
##' @export seed_rain
##' @export set_sane_gsl_error_handling
##' @export species
##' @export test_adaptive_spline
##' @export test_find_root
##' @export test_find_value
##' @export test_from_rcpp_matrix
##' @export test_functor
##' @export test_gradient
##' @export test_gradient_richardson
##' @export test_integrator
##' @export test_plant
##' @export test_sum_double
##' @export test_sum_int
##' @export test_to_rcpp_matrix
##' @export trapezium

## and this will then load the Rcpp module.
loadModule("tree", TRUE)
