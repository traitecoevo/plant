## Support functions for growing an individual to a given size or time, which
## run() on an Individual calls. The integration itself is C++
## (grow_to_sizes() and grow_to_times() in solvers/individual_runner.h); these
## check the arguments and shape the result.

## Grows `individual` to each of `sizes` (increasing) of the state `size_name`.
## Returns time (NA for a size not reached), state (the ODE state at those
## times), individual (the plant at each), trajectory (every step taken) and
## env. `time_max` stops an unreachable size from running forever; `warn` warns
## when sizes are dropped, and `filter` drops them from the result.
##
## One run reaches every size in turn: each size's time is located on the
## interpolant of the step that crosses it, and the run lands there and goes
## on from there. So a size's time, state and individual describe the same
## instant, and the trajectory holds the step landed on each size.
individual_grow_to_sizes <- function(individual, sizes, size_name, env,
                                     time_max=Inf, warn=TRUE, filter=FALSE) {
  if (length(sizes) == 0L || anyNA(sizes) || is.unsorted(sizes)) {
    stop("sizes must be non-empty and sorted")
  }
  if (!all(is.finite(sizes))) {
    stop("sizes must be finite")
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
  if (!is.numeric(time_max) || length(time_max) != 1L || is.na(time_max) ||
      time_max < 0) {
    stop("time_max must be a single non-negative number")
  }
  ## odelia refuses an infinite time_max; the largest double is no limit.
  until <- if (time_max == Inf) .Machine$double.xmax else time_max

  res <- individual$grow_to_sizes(env, sizes, j - 1L, until)
  dropped <- sum(is.na(res$time))
  if (!is.null(res$error)) {
    integration_failed(res$last, res$error, dropped, warn)
  } else if (res$timed_out && warn) {
    warning(sprintf("Time reached time_max, %d larger sizes dropped", dropped),
            immediate.=TRUE)
  }

  colnames(res$state) <- ode_names
  colnames(res$trajectory) <- c("time", ode_names)
  ret <- list(time=res$time,
              state=res$state,
              individual=res$individual,
              trajectory=res$trajectory,
              env=env)

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
## state at each) and individual (the plant at each), and env. The run lands on
## each time in turn.
individual_grow_to_times <- function(individual, times, env) {
  if (length(times) == 0L) {
    stop("At least one time must be given")
  }
  if (anyNA(times) || any(times < 0.0)) {
    stop("Times must be positive")
  }
  if (is.unsorted(times)) {
    stop("times must be sorted")
  }
  res <- individual$grow_to_times(env, times)
  colnames(res$state) <- individual$ode_names
  list(time=times, state=res$state, individual=res$individual, env=env)
}

## Warns that the run to `sizes` stopped on error `message`, with `dropped`
## sizes not reached; `individual` is the plant where it stopped.
integration_failed <- function(individual, message, dropped, warn) {
  if (warn) {
    msg <- c("Stopping early as integration failed with error: ",
             message,
             sprintf("  %d larger sizes dropped", dropped))
    warning(paste(msg, collapse="\n"), immediate.=TRUE)
  }
  ## ⚠️ BOTH HALVES, BECAUSE A PARKED RATE IS WHAT A DEAD PLANT LOOKS LIKE.
  ## Every strategy's mortality_dt returns exactly zero once the accumulated
  ## hazard reaches establishment_failure_hazard, so the rate alone cannot
  ## separate "it died" from "it failed for some other reason" -- it reads
  ## zero for both. The hazard is the fact and the rate is the proxy, so a
  ## plant is taken to have died if either says so.
  rates <- individual$ode_rates
  states <- individual$ode_state
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

