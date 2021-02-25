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

# ! Important the whitespace in the following function is used by the strategy scaffolder
get_list_of_strategy_types <- function() {
  list(
    FF16=FF16_Strategy,
    FF16r=FF16r_Strategy,
    K93=K93_Strategy
    )
}

get_list_of_environment_types <- function() {
  list(
    FF16="FF16_Env",
    FF16r="FF16_Env",
    K93="K93_Env"
    )
}

# ! Important the whitespace in the following function is used by the strategy scaffolder
get_list_of_hyperpar_functions <- function() {
  list(
    FF16=FF16_hyperpar,
    FF16r=FF16r_hyperpar,
    K93=K93_hyperpar
    )
}

test_environment <- function(type, ...) {
  switch(type,
    FF16=FF16_test_environment(...),
    FF16r=FF16_test_environment(...),
    K93=K93_test_environment(...),
    stop("Unknown type ", type))
}

fixed_environment<- function(type, ...) {
  switch(type,
    FF16=FF16_fixed_environment(...),
    FF16r=FF16_fixed_environment(...),
    K93=K93_fixed_environment(...),
    stop("Unknown type ", type))
}
