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
