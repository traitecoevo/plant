library(methods)
library(Rcpp)
library(testthat)

options(warn=1)

if ( packageVersion("Rcpp") < package_version("0.10.0") )
  error("Needs Rcpp at least version 0.10.0")

if ( !exists("Spline") ) {
  dyn.load(file.path("../src", "tree.so"))

  tree_module <- Module("tree", "tree")

  ## Cause GSL errors to generate R errors.
  tree_module$set_sane_gsl_error_handling()
  
  Spline          <- tree_module$Spline
  AdaptiveSpline  <- tree_module$AdaptiveSpline
  AdaptiveSplineR <- tree_module$AdaptiveSplineR
  MultiSpline     <- tree_module$MultiSpline

  Lorenz <- tree_module$Lorenz
  OdeR   <- tree_module$OdeR

  Strategy    <- tree_module$Strategy
  Control     <- tree_module$Control
  Parameters  <- tree_module$Parameters
  PlantSpline <- tree_module$PlantSpline # poorly named...

  Plant          <- tree_module$Plant
  CohortDiscrete <- tree_module$CohortDiscrete
  PlantApprox    <- tree_module$PlantApprox

  ## Slightly different to the other Individuals...
  CohortTop      <- tree_module$CohortTop

  Species <- tree_module$Species
  SpeciesC <- tree_module$SpeciesC
  SpeciesCT <- tree_module$SpeciesCT

  Patch  <- tree_module$Patch
  PatchC <- tree_module$PatchC

  Metacommunity <- tree_module$Metacommunity
  MetacommunityC <- tree_module$MetacommunityC
}

source("../R/falster.R", chdir=TRUE)

## New expect_that helper functions; test that a number is in a range,
## or that a range contains a number.
is_within_interval <- function(lower, upper) {
  if ( missing(upper) && length(lower) == 2 ) {
    upper <- lower[[2]]
    lower <- lower[[1]]
  }
  if ( lower >= upper )
    stop("lower must be smaller than upper")
  err <- paste0("is not within range [", lower, ", ", upper, "]")
  function(actual) {
    expectation(actual > lower && actual < upper, err)
  }
}

contains <- function(value) {
  function(range) {
    if ( length(range) != 2 )
      stop("Expected a vector of length 2")
    if ( range[[1]] >= range[[2]] )
      stop("Expected that range[1] is smaler than range[2]")
    expectation(value > range[[1]] && value < range[[2]],
                paste("does not contain", value))
  }
}

is_greater_than <- function(value) {
  function(actual)
    expectation(actual > value, paste("is not greater than", value))
}

is_less_than <- function(value) {
  function(actual)
    expectation(actual < value, paste("is not less than", value))
}


## This is here to stop a really weird bug (not our fault(?)) around
## finalisers.
gc()
