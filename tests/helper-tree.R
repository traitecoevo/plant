library(methods)
library(Rcpp)
library(testthat)

options(warn=1)

if ( packageVersion("Rcpp") < package_version("0.10.0") )
  error("Needs Rcpp at least version 0.10.0")

if ( !exists("Spline") ) {
  ## TODO: Sort out how we run this.  Package may help.
  dyn.load(file.path("../src", "tree.so"))
  .Call("set_sane_gsl_error_handling", PACKAGE="tree")
  
  tree_module <- Module("tree", "tree")
  
  Spline          <- tree_module$Spline
  AdaptiveSpline  <- tree_module$AdaptiveSpline
  AdaptiveSplineR <- tree_module$AdaptiveSplineR
  MultiSpline     <- tree_module$MultiSpline

  Lorenz <- tree_module$Lorenz
  OdeR   <- tree_module$OdeR

  Strategy   <- tree_module$Strategy
  Parameters <- tree_module$Parameters

  Plant          <- tree_module$Plant
  CohortDiscrete <- tree_module$CohortDiscrete
  Cohort         <- tree_module$Cohort
  PlantSpline    <- tree_module$PlantSpline

  Patch  <- tree_module$Patch
  PatchC <- tree_module$PatchC
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
