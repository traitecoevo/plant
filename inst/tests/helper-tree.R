library(tree)
library(testthat)
options(warn=1)

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

## This makes a pretend light environment over the plant height,
## slightly concave up, whatever.
test.environment <- function(height, n=101, light.env=NULL,
                             n.strategies=1) {
  hh <- seq(0, height, length=n)
  if (is.null(light.env))
    light.env <- function(x)
      exp(x/(height*2)) - 1 + (1 - (exp(.5) - 1))/2
  ee <- light.env(hh)
  env <- new(Spline)
  env$init(hh, ee)

  parameters <- new(Parameters)
  for (i in seq_len(n.strategies))
    parameters$add_strategy(new(Strategy))

  ret <- new(Environment, parameters)
  ret$light_environment <- env
  attr(ret, "light.env") <- light.env
  ret
}

## This is here to stop a weird problem with class finalisation.  A
## fairly harmless error message at this point, but hopefully not
## something that will stick around.  Seems related to
## https://stat.ethz.ch/pipermail/r-devel/2011-December/062917.html
gc()
