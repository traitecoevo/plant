
get_list_of_strategy_types <- function() {
  list(FF16=FF16_Strategy,
       FFdev=FFdev_Strategy)
}

get_list_of_hyperpar_functions <- function() {
  list(FF16=FF16_hyperpar,
       FFdev=FFdev_hyperpar)
}

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
    testthat::expectation(actual > lower && actual < upper, err)
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
    testthat::expectation(value > range[[1]] && value < range[[2]],
                          paste("does not contain", value))
  }
}

is_greater_than <- function(value) {
  function(actual) {
    testthat::expectation(actual > value, paste("is not greater than", value))
  }
}

is_less_than <- function(value) {
  function(actual) {
    testthat::expectation(actual < value, paste("is not less than", value))
  }
}

is_at_most <- function(value) {
  function(actual) {
    testthat::expectation(actual <= value, paste("is greater than", value))
  }
}

is_at_least <- function(value) {
  function(actual) {
    testthat::expectation(actual >= value, paste("is less than", value))
  }
}

## This makes a pretend light environment over the plant height,
## slightly concave up, whatever.
test_environment <- function(height, n=101, light_env=NULL,
                             n_strategies=1, seed_rain=0) {
  if (length(seed_rain) == 1) {
    seed_rain <- rep(seed_rain, length.out=n_strategies)
  }
  hh <- seq(0, height, length.out=n)
  if (is.null(light_env)) {
    light_env <- function(x) {
      exp(x/(height*2)) - 1 + (1 - (exp(.5) - 1))/2
    }
  }
  ee <- light_env(hh)
  env <- Interpolator()
  env$init(hh, ee)

  parameters <- FF16_Parameters()
  parameters$strategies <- rep(list(FF16_Strategy()), n_strategies)
  parameters$seed_rain <- seed_rain
  parameters$is_resident <- rep(TRUE, n_strategies)

  ret <- make_environment(parameters)
  ret$light_environment <- env
  attr(ret, "light_env") <- light_env
  ret
}

test_ode_make_system <- function(obj) {
  make_derivs <- function(obj) {
    if (is.null(obj$set_ode_state)) {
      function(y, t) {
        obj$ode_state <- y
        obj$ode_rates
      }
    } else {
      function(y, t) {
        obj$set_ode_state(y, t)
        obj$ode_rates
      }
    }
  }
  ## Hmm, this is causing all sorts of trouble...
  make_state <- function(obj) {
    function() {
      obj$ode_state
    }
  }
  time <- if (is.null(obj$time)) 0.0 else obj$time
  sys <- OdeR(make_derivs(obj), make_state(obj), time)
}

test_ode_make_solver <- function(sys) {
  OdeRunner(class(sys)[[1]])(sys)
}

skip_if_no_plant_ml_python <- function() {
  skip_if_not_installed("plant.ml")
  ok <- !inherits(try(plant.ml:::ensure_python_loaded(), silent=TRUE),
                  "try-error")
  if (ok) {
    return()
  }
  skip("python packages missing")
}
