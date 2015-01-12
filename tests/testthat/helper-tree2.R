## New expect_that helper functions; test that a number is in a range,
## or that a range contains a number.
is_within_interval <- function(lower, upper) {
  if ( missing(upper) && length(lower) == 2 ) {
    upper <- lower[[2]]
    lower <- lower[[1]]
  }
  if (lower >= upper) {
    stop("lower must be smaller than upper")
  }
  err <- paste0("is not within range [", lower, ", ", upper, "]")
  function(actual) {
    expectation(actual > lower && actual < upper, err)
  }
}

contains <- function(value) {
  function(range) {
    if (length(range) != 2L) {
      stop("Expected a vector of length 2")
    }
    if (range[[1]] >= range[[2]]) {
      stop("Expected that range[1] is smaler than range[2]")
    }
    expectation(value > range[[1]] && value < range[[2]],
                paste("does not contain", value))
  }
}

is_greater_than <- function(value) {
  function(actual) {
    expectation(actual > value, paste("is not greater than", value))
  }
}

is_less_than <- function(value) {
  function(actual) {
    expectation(actual < value, paste("is not less than", value))
  }
}

is_at_most <- function(value) {
  function(actual) {
    expectation(actual <= value, paste("is greater than", value))
  }
}

is_at_least <- function(value) {
  function(actual) {
    expectation(actual >= value, paste("is less than", value))
  }
}

## This makes a pretend light environment over the plant height,
## slightly concave up, whatever.
test_environment <- function(height, n=101, light_env=NULL,
                             n_strategies=1, seed_rain=0) {
  if (length(seed_rain) == 1) {
    seed_rain <- rep(seed_rain, length.out=n_strategies)
  }
  hh <- seq(0, height, length=n)
  if (is.null(light_env)) {
    light_env <- function(x) {
      exp(x/(height*2)) - 1 + (1 - (exp(.5) - 1))/2
    }
  }
  ee <- light_env(hh)
  env <- Interpolator()
  env$init(hh, ee)

  parameters <- Parameters()
  parameters$strategies <- rep(list(Strategy()), n_strategies)
  parameters$seed_rain <- seed_rain

  ret <- Environment(parameters)
  ret$light_environment <- env
  attr(ret, "light_env") <- light_env
  ret
}
