## Support functions for growing an individual to a given size or time, which
## run() on an Individual calls.

## Grows `individual` to each of `sizes` (increasing) of the state `size_name`.
## Returns time (NA for a size not reached), state (the ODE state at those
## times), individual (the plant at each), trajectory (every step taken) and
## env. `time_max` stops an unreachable size from running forever; `warn` warns
## when sizes are dropped, and `filter` drops them from the result.
##
## One run reaches every size in turn: each size's time is located on the
## interpolant of the step that crosses it, and the run lands there and goes
## on from there (advance_to_size(), on odelia's Solver::advance_until). So a
## size's time, state and individual describe the same instant, and the
## trajectory holds the step landed on each size.
individual_grow_to_sizes <- function(individual, sizes, size_name, env,
                                     time_max=Inf, warn=TRUE, filter=FALSE) {
  if (length(sizes) == 0L || is.unsorted(sizes)) {
    stop("sizes must be non-empty and sorted")
  }
  ode_names <- individual$ode_names
  j <- match(size_name, ode_names)
  if (is.na(j)) {
    stop(sprintf("'%s' is not one of the ODE states: %s", size_name,
                 paste(ode_names, collapse = ", ")))
  }
  if (individual$ode_state[[j]] > sizes[[1]]) {
    stop("Individual already bigger than smallest target size")
  }
  strategy_name <- individual$strategy_name
  runner <- OdeRunner(strategy_name)(
    IndividualRunner(strategy_name, environment_type(strategy_name))(individual, env))
  ## odelia refuses an infinite time_max; the largest double is no limit.
  until <- if (is.finite(time_max)) time_max else .Machine$double.xmax

  n <- length(sizes)
  time <- rep(NA_real_, n)
  state <- matrix(NA_real_, n, length(ode_names), dimnames = list(NULL, ode_names))
  at_size <- vector("list", n)
  start <- c(runner$time, runner$state)

  for (i in seq_len(n)) {
    if (i > 1L && sizes[[i]] == sizes[[i - 1L]]) {
      t <- time[[i - 1L]]
    } else if (runner$state[[j]] >= sizes[[i]]) {
      ## Reached already: where a run starts on the size, or lands a hair past
      ## one barely above the last.
      t <- runner$time
    } else {
      t <- tryCatch(runner$advance_to_size(j - 1L, sizes[[i]], until),
                    error = function(e) {
                      integration_failed(runner, e, n - i + 1L, warn)
                      NULL
                    })
      if (is.null(t)) {
        break
      }
      if (is.na(t)) {
        if (warn) {
          warning(sprintf("Time reached time_max, %d larger sizes dropped",
                          n - i + 1L), immediate.=TRUE)
        }
        break
      }
    }
    time[[i]] <- t
    state[i, ] <- runner$state
    at_size[i] <- list(runner$object$individual)
  }

  ret <- list(time=time,
              state=state,
              individual=at_size,
              trajectory=rbind(start, runner$steps_taken(), deparse.level = 0),
              env=env)
  colnames(ret$trajectory) <- c("time", ode_names)

  if (filter) {
    i <- !vlapply(ret$individual, is.null)
    if (!all(i)) {
      ret$time  <- ret$time[i]
      ret$state <- ret$state[i, , drop=FALSE]
      ret$individual <- ret$individual[i]
    }
  }
  ret
}

vlapply <- function(X, FUN, ...) {
  vapply(X, FUN, logical(1), ...)
}

vnapply <- function(X, FUN, ...) {
  vapply(X, FUN, numeric(1), ...)
}

## Grows `individual` and reports it at each of `times`: time, state (the ODE
## state at each) and individual (the plant at each), and env.
individual_grow_to_times <- function(individual, times, env) {
  if (any(times < 0.0)) {
    stop("Times must be positive")
  }
  n <- length(times)
  if (n == 0L) {
    stop("At least one time must be given")
  }

  y0 <- individual$ode_state
  t0 <- 0.0
  i <- 1L
  t_next <- times[[i]]
  strategy_name <- individual$strategy_name

  ir1 <- IndividualRunner(strategy_name, environment_type(strategy_name))(individual, env)
  ir2 <- IndividualRunner(strategy_name, environment_type(strategy_name))(individual, env)

  runner <- OdeRunner(strategy_name)(ir1)
  runner_detail <- OdeRunner(strategy_name)(ir2)

  ## `runner` steps freely; for each requested time a step crosses,
  ## `runner_detail` re-integrates from the start of that step to land on it
  ## exactly.
  state <- matrix(NA, n, length(y0))
  colnames(state) <- individual$ode_names

  individual <- vector("list", n)
  while (i <= n) {
    runner$step()
    t1 <- runner$time
    y1 <- runner$state
    while (t_next < t1 && i <= n) {
      runner_detail$set_state(y0, t0)
      ## odelia's Solver advances adaptively over a vector of times whose first
      ## element must be the current time; this steps from t0 to t_next.
      runner_detail$advance_adaptive(c(t0, t_next))
      state[i, ] <- runner_detail$state
      individual[[i]] <- runner_detail$object$individual
      i <- i + 1L
      t_next <- times[i] # allows out-of-bounds extraction
    }
    t0 <- t1
    y0 <- y1
  }

  list(time=times, state=state, individual=individual, env=env)
}

## internal funciton to grab the internal state of the ode runner:
get_individual_internals_fun <- function (individual) {
  get(paste0(individual$strategy_name, '_oderunner_individual_internals')) #!
}

## Warns that the run to `sizes` stopped on error `e`, with `dropped` sizes not
## reached.
integration_failed <- function(runner, e, dropped, warn) {
  if (warn) {
    msg <- c("Stopping early as integration failed with error: ",
             conditionMessage(e),
             sprintf("  %d larger sizes dropped", dropped))
    warning(paste(msg, collapse="\n"), immediate.=TRUE)
  }
  ## ⚠️ BOTH HALVES, BECAUSE A PARKED RATE IS WHAT A DEAD PLANT LOOKS LIKE.
  ## Every strategy's mortality_dt returns exactly zero once the accumulated
  ## hazard reaches establishment_failure_hazard, so the rate alone cannot
  ## separate "it died" from "it failed for some other reason" -- it reads
  ## zero for both. The hazard is the fact and the rate is the proxy, so a
  ## plant is taken to have died if either says so.
  rates <- runner$object$individual$ode_rates
  states <- runner$object$individual$ode_state
  if (rates[[2]] < 1e-10 && states[[2]] < 700) {
    warning("Integration may have failed for reasons other than mortality",
            immediate.=TRUE)
  }
}

#' Compute the whole plant light compensation point for a single
#' plant: the light level at which whole-plant carbon gain is zero.
#' @title Whole plant light compensation point
#' @param p An \code{Individual}, with strategy, height, etc set.
#' @param ... Additional arguments that are ignored
#' @export
#' @author Rich FitzJohn
resource_compensation_point <- function(p, ...) {
  p$resource_compensation_point()
}

#' `optimise_individual_rate_at_size_by_trait` solves for the maximum of
#' some rate (e.g. growth rate) at a specified height within
#' the interval of the bounds of a given trait
#' @param type The type of model to use (e.g. "FF16"). Defaults to "FF16"
#' @param bounds A vector giving the lower and upper bounds of the trait
#' @param log_scale Should the trait be optimised on a log scale? Defaults to TRUE
#' @param tol The tolerance for the optimisation
#' @param size  The size of the individual to optimise the rate at
#' @param size_name The name of the size variable specified by \code{size}
#' @param rate The name of the rate to optimise. Defaults to \code{size_name}
#' @param params The parameters of the model. Defaults to plant's testing set,
#'   \code{\link{default_parameters}(type, quiet = TRUE)}, so that the
#'   trade-offs apply to the trait being optimised.
#' @param env The environment of the model
#' @param hyperpars The trade-off function; defaults to the one \code{params}
#'   carries (\code{param_hyperpar(params)})
#' @param set_state_directly If TRUE, set the state directly to the size, otherwise grows the plant to that size. Defaults to FALSE
#' @export
#' @rdname optimise_individual_rate_at_size_by_trait
#' @author Isaac Towers, Daniel Falster and Andrew O'Reilly-Nugent

optimise_individual_rate_at_size_by_trait <- function(
    type = "FF16",
    bounds, log_scale = TRUE, tol = 1e-3,
    size = 1, size_name = "height",
    rate = size_name,
    params = default_parameters(type, quiet = TRUE),
    env = Environment(type),
    hyperpars = param_hyperpar(params),
    set_state_directly = FALSE) {
  
   # can't handle situations yet where bounds are outside of positive growth, not working for K93
#  bounds <- check_bounds(bounds)

  traits <- rownames(bounds)

  if (log_scale) {
    bounds[bounds[, 1] == -Inf, 1] <- 0
    bounds <- log(bounds)
    ff <- exp
  } else {
    ff <- I
  }

  ## Define function to optimise
  f <- function(x) {
    # create a strategy object
    s <- generate_strategy(params, ff(trait_matrix(x, rownames(bounds))), hyperpar = hyperpars, birth_rate = 1)[[1]]

    # Create an individual object
    types <- extract_RcppR6_template_types(params, "Parameters")
    indv <- do.call("Individual", types)(s) # equiavlent to calling Individual<TF24,TF24_Env> or FF16_individual(s)
    
    if(set_state_directly & size_name == "height"){
      # set inidividual at specified size
      indv$set_state(size_name, size)
      # compute rates given environment
      indv$compute_rates(env)
      #filter ode rate based on size name
      res <- indv$ode_rates[indv$ode_names == size_name]
      
    } else{
      res <- individual_grow_to_sizes(indv,
                                       sizes = size, size_name = size_name, env,
                                       time_max = 100, warn = TRUE, filter = TRUE) 
      #check if there was positive growth (indicated by the presence of postive numbers of rows in res$stae)
      if(nrow(res$state) == 0){
        res = NA
      } else{
      res <- res$individual[[1]]$ode_rates[res$individual[[1]]$ode_names == size_name]
      }

    # turn NA (non-positive growth) to 0, allows optimiser to get finite value but at minimum and thus avoided
    if (is.na(res)) {
      res <- 0
    }
    return(res)
    }
  }
  

  #solve for the trait value which maximise size growth
  ret <- solve_max_worker(bounds, f, tol = 1e-6, outcome = paste0(size_name, "_growth_rate"))

  #exponentiate the optimum trait value
  if (log_scale) {
    ret <- exp(ret)
  }
  
  #if the optimum growth rate is 0, means no postive growth, so set set optimum trait value to NA
  if (attr(ret, paste0(size_name, "_growth_rate")) == 0){
    ret[1] = NA
  }

  
  return(ret)
}

solve_max_worker <- function(bounds, f, tol = 1e-3, outcome) {
  if (length(rownames(bounds)) == 1L) {
    if (!all(is.finite(bounds))) {
      stop("Starting value did not have finite fitness; finite bounds required")
    }
    ## The suppressWarnings here is for warnings like:
    ##
    ## Warning message:
    ## In optimise(f, interval = bounds, maximum = TRUE, tol = tol) :
    ##   NA/Inf replaced by maximum positive value
    ##
    ## which is probably the desired behaviour here.
    out <- suppressWarnings(stats::optimise(f, interval = bounds, maximum = TRUE, tol = tol))
    # browser()
    ret <- out$maximum
    attr(ret, outcome) <- out$objective
  } else {
    ## This is not very well tested, and the tolerance is not useful:
    out <- stats::optim(rowMeans(bounds), f,
      method = "L-BFGS-B",
      lower = bounds[, "lower"], upper = bounds[, "upper"],
      control = list(fnscale = -1, factr = 1e10)
    )

    ret <- out$value
    attr(ret, outcome) <- out$par
  }
  return(ret)
}

