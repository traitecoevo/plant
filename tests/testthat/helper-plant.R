
source(
  file.path(
    rprojroot::find_testthat_root_file(),
    "FF16_reference", "make_reference_plant.R")
)

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

# NOTE: TF24f (#525) is intentionally excluded from these generic cross-strategy
# lists. It is an experimental variant with a 6th ODE state (tracked root-collar
# psi) and relies on birth-initialisation via Node::compute_initial_conditions, so
# it does not satisfy the generic invariants these lists drive (5-state layout;
# growth of a standalone Individual stepped a few times, which bypasses that
# initialisation). TF24f is covered via the SCM in test-strategy-tf24f.R.

# ! Important the whitespace in the following function is used by the strategy scaffolder
get_list_of_strategy_types <- function() {
  list(
    FF16=FF16_Strategy,
    TF24=TF24_Strategy,
    K93=K93_Strategy
    )
}

get_list_of_environment_types <- function() {
  list(
    FF16="FF16_Env",
    TF24="TF24_Env",
    K93="K93_Env"
    )
}

# ! Important the whitespace in the following function is used by the strategy scaffolder
get_list_of_hyperpar_functions <- function() {
  list(
    FF16=FF16_hyperpar,
    TF24=TF24_hyperpar,
    K93=K93_hyperpar
    )
}

