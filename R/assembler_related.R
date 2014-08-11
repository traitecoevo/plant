## Functions *related* to assembly, but not actually doing it.

##' Compute max growth rate of a given set of values for a trait.
##' This is the log of per-capita seed production (i.e., fitness).
##'
##' Only works in one dimension
##' @title Compute Max Growth Rate
##' @param trait Name of the trait (e.g., \code{"lma"})
##' @param values Values to compute maximum growth rate for
##' @param p Parameters object to use.  Importantly, the
##' \code{strategy_default} element gets used here.
##' @param schedule \code{CohortSchedule} to use, or \code{NULL} to
##' generate a hopefully reasonable schedule.
##' @author Rich FitzJohn
##' @export
max_growth_rate <- function(trait, values, p, schedule=NULL) {
  log(landscape_empty(trait, values, p, schedule))
}

##' Compute the carrying capacity (equilibrium per-capita seed
##' production) for a set of values of a trait.  Each is considered in
##' isolation.
##'
##' @title Carrying Capacity
##' @param trait Name of the trait (e.g., \code{"lma"})
##' @param values Values to compute maximum growth rate for
##' @param p Parameters object to use.  Importantly, the
##' \code{strategy_default} element gets used here.
##' @param seed_rain Initial seed rain (optional)
##' @param parallel Use multiple processors?
##' @author Rich FitzJohn
##' @export
carrying_capacity <- function(trait, values, p, seed_rain=1,
                              parallel=FALSE) {
  f <- function(x) {
    carrying_capacity1(trait, x, p, seed_rain)
  }
  if (parallel) {
    unlist(parallel::mclapply(values, f, mc.preschedule=FALSE))
  } else {
    unlist(lapply(values, f))
  }
}

##' Compute region of positive fitness.  This will be the values where
##' fitness is approximately zero.
##'
##' @title Compute Region of Positive Fitnes
##' @param trait Name of the trait (e.g., \code{"lma"})
##' @param bounds 2D vector specifing range within which to search
##' @param p Parameters object to use.  Importantly, the
##' \code{strategy_default} element gets used here.
##' @param value Initial value - must have positive fitness itself!
##' If not given, then the value from the default strategy within
##' \code{p} is used.
##' @param log_scale Is the parameter naturally on a log scale?  If
##' so, this will greatly speed things up.
##' @param dx Amount to step the trait.  If \code{log_scale} is
##' \code{TRUE}, this is on a log scale.
##' @export
##' @author Rich FitzJohn
viable_fitness <- function(trait, p, bounds=NULL, value=NULL,
                           log_scale=TRUE, dx=1) {
  if(length(trait) > 1)
    stop("Doesn't yet support multiple traits")
  if (is.null(value)) {
    value <- p$strategy_default$parameters[[trait]]
  }
  if (is.null(bounds)) {
    if (log_scale) {
      bounds <- c(1e-5, 1e3)
    } else {
      bounds <- c(-Inf, Inf)
    }
  }
  if(value > bounds[2] || value < bounds[1])
    stop("Value does not lie within bounds")
  if (log_scale) {
    f <- function(x) {
      max_growth_rate(trait, exp(x), p)
    }
    out <- exp(positive(f, log(value), dx, lower=log(bounds[1]), upper=log(bounds[2])))
  } else {
    f <- function(x) {
      max_growth_rate(trait, x, p)
    }
    out <- positive(f, value, dx, lower=bounds[1], upper=bounds[2])
  }
  bounds <- rbind(x=out)
  colnames(bounds) <- c("lower", "upper")
  rownames(bounds) <- trait
  bounds
}

carrying_capacity1 <- function(trait, value, p, seed_rain=1) {
  p <- p$copy()
  p$clear()
  ## Use
  p <- expand_parameters(trait, value, p)
  p$seed_rain <- seed_rain
  if (is.null(schedule)) {
    schedule <- default_cohort_schedule(p)
  }
  res <- equilibrium_seed_rain(p)
  mean(res$seed_rain)
}

positive_bracket <- function(f, x, dx, lower=-Inf, upper=Inf) {
  factor <- 2
  fx <- f(x)
  if (fx < 0) {
    stop("Don't yet support doing this with no positive values")
  }

  L <- U <- x
  fL <- fU <- fx
  dL <- dU <- dx

  bracket <- function(x, dx, bound) {
    cleanup <- function(x, x_next, fx, fx_next) {
      ## This is *approximately* what uniroot will do anyway.  We'll
      ## get a usable negative number.
      if (fx_next == -Inf) {
        fx_next <- .Machine$double.xmin/2
      }
      if (dx < 0) {
        x <- c(x_next, x)
        fx <- c(fx_next, fx)
      } else {
        x <- c(x, x_next)
        fx <- c(fx, fx_next)
      }
      list(x=x, fx=fx)
    }
    hit_bounds <- FALSE
    repeat {
      x_next <- x + dx
      if ((dx < 0 && x_next < bound) || (dx > 0 && x_next > bound)) {
        x_next <- bound
        hit_bounds <- TRUE
      }
      fx_next <- f(x_next)
      if (fx_next < 0 || hit_bounds) {
        return(cleanup(x, x_next, fx, fx_next))
      } else {
        x <- x_next
        fx <- fx_next
        dx <- dx * factor
      }
    }
  }

  list(lower=bracket(x, -dx, lower), upper=bracket(x, dx, upper))
}

positive <- function(f, x, dx, lower=-Inf, upper=Inf, eps=1e-3) {
  b <- positive_bracket(f, x, dx, lower, upper)

  # Find lower root. If no root exists within that range, take
  # lower bound
  if(prod(b$lower$fx[1:2]) < 0){
        lower <- uniroot(f, b$lower$x,
                   f.lower=b$lower$fx[[1]], f.upper=b$lower$fx[[2]],
                   tol=eps)$root
  } else {
    lower <- b$lower$x[[1]]
  }
  # Find upper root. If no root exists within that range, take
  # upper bound
  if(prod(b$upper$fx[1:2]) < 0){
    upper <- uniroot(f, b$upper$x,
                   f.lower=b$upper$fx[[1]], f.upper=b$upper$fx[[2]],
                   tol=eps)$root
  } else {
    upper <- b$upper$x[[2]]
  }

  c(lower, upper)
}

#' Solve for equilbrium community with given values of specified trait
#'
#' Find seed rain and cohort schedule for equilbrium community with
#' given traits. This is point at which each resident seed returns
#' on average a single seed.
#' @param trait name of trait
#' @param values vector of trait values for community
#' @param p Parameters object to use.  Importantly, the
#' \code{strategy_default} element gets used here.
#' @param seed_rain=NULL vector of initial seed rains for community
#' @param verbose=TRUE prints details to screen
#' @author Daniel Falster
#' @export
get_equilibrium_community <- function(trait, values, p, seed_rain = NULL, verbose = TRUE) {
    p <- p$copy()
    p$clear()
    if (verbose)
        p$set_control_parameters(equilibrium_verbose())
    p <- expand_parameters(trait, values, p, FALSE)
    if (!is.null(seed_rain)) {
        if (length(values) != length(seed_rain))
            stop("incorrect vector lengths")
        p$seed_rain <- seed_rain
    }
    res <- equilibrium_seed_rain(p)
    p$seed_rain <- mean(res$seed_rain)
    list(p = p, schedule = res$schedule)
}

#' Find evolutionary attractor for single species and trait.
#'
#' Find evolutionary attractor for single species and trait.
#' This is point at which selection gradient equals zero.
#' Currently solved using \code{uniroot}
#' @param trait name of trait
#' @param interval  a vector containing the end-points of the
#' interval to be searched for the root.
#' @param p Parameters object to use.  Importantly, the
#' \code{strategy_default} element gets used here.
#' @param tol=1e-4 the desired accuracy (convergence tolerance).
#' @param ... set verbpse=TRUE for verbose output, seed_rain=value
#' gives starting values when solving for demographic equilibrium
#' @author Daniel Falster
#' @export
#' @return Same as \code{uniroot} function: A list with at least
#' four components:
#' ‘root’ and ‘f.root’ give the location of the root and the
#' value of the function evaluated at that point. ‘iter’ and
#' ‘estim.prec’ give the number of iterations
find_singularity_1D <- function(trait, interval, p, tol = 1e-04, ...) {

    f <- function(x) selection_gradient(trait, x, p, ...)
    root <- uniroot(f, interval, tol = tol)
}

#' Returns selection gradient in single species communities
#' at a range of trait values.
#'
#' Returns selection gradient in single species communities
#' at a range of trait values. This is derivative of fitness
#' with respect to trait value. The algorithm first solves
#' for the demographic attractor, using \code{seed_rain} as
#' a starting value.
#' @param trait name of trait.
#' @param values vector of trait values.
#' @param p Parameters object to use.  Importantly, the
#' \code{strategy_default} element gets used here.
#' @param dx=1e-4 Interval over which derivative is calculated
#' @param log_scale=TRUE Determines whether derivative is taken
#' with respect to raw or log-transformed x values. The latter
#' is useful when x is log-normally distributed.
#' @param seed_rain=NULL (optional) Starting values for seed rain
#' when solving for demographic euqilibrium
#' @param verbose=TRUE (optional) Print values to screen
#' @author Daniel Falster
#' @export
#' @seealso selection_gradient1
#' @return a vector of values with same length as \code{values}.
selection_gradient <- function(trait, values, p, dx = 1e-04, log_scale = TRUE, seed_rain = NULL,
    verbose = TRUE) {
    N <- length(values)
    ret <- rep(NA, N)
    for (i in seq_len(N)) ret[i] <- selection_gradient1(trait, values[i], p, dx,
        log_scale, seed_rain, verbose)
    ret
}

#' Returns selection gradient in a single-species community
#' with given trait value
#'
#' Returns selection gradient in a single-species community
#' with given trait value. This is derivative of fitness
#' with respect to trait value. The algorithm first solves
#' for the demographic attractor, using \code{seed_rain} as
#' a starting value.
#' @param trait name of trait
#' @param value value of trait where derivative is calculated
#' @param p Parameters object to use.  Importantly, the
#' \code{strategy_default} element gets used here.
#' @param dx=1e-4 Interval over which derivative is calculated
#' @param log_scale=TRUE Determines whether derivative is taken
#' with respect to raw or log-transformed x values. The latter
#' is useful when x is log-normally distributed.
#' @param seed_rain=NULL (optional) Starting values for seed rain
#' when solving for demographic euqilibrium
#' @param verbose=TRUE (optional) Print values to screen
#' @author Daniel Falster
#' @export
#' @seealso selection_gradient
#' @return a single value. The equilibirum seed rain is included as
#' an attribute.
selection_gradient1 <- function(trait, value, p, dx = 1e-04, log_scale = TRUE, seed_rain = NULL,
    verbose = TRUE) {

    if (verbose)
        message(sprintf("Calculating selection gradient for %s = %f", trait, value))
    res <- get_equilibrium_community(trait, value, p, seed_rain = seed_rain, verbose = verbose)

    f <- function(x) landscape(trait, x, res$p, res$schedule)
    ret <- gradient_fd(f, value, dx, log_scale)
    attr(ret, "seed_rain") <- res$p$seed_rain

    if (verbose)
        message(sprintf("selection gradient = %f", ret))
    ret
}

##' Find point of maximum fitness in empty landscape
##' within specified range.
##'
##' @title Find point of maximum fitness within some range.
##' @param trait Name of the trait (e.g., \code{"lma"})
##' @param bounds 2D vector specifing range within which to search
##' @param p Parameters object to use.  Importantly, the
##' \code{strategy_default} element gets used here.
##' @param log_scale Is the parameter naturally on a log scale?  If
##' so, this will greatly speed things up.
##' @export
##' @author Daniel Falster, Rich FitzJohn
max_fitness <- function(trait, p, bounds=NULL,
                           log_scale=TRUE) {
  if(length(trait) > 1)
    stop("Doesn't yet support multiple traits")
  if (log_scale) {
    if (is.null(bounds))
      bounds <- c(1E-5, 1E3)
    f <- function(x) max_growth_rate(trait, exp(x), p)
  } else {
    if (is.null(bounds))
      bounds <- c(-Inf, Inf)
    f <- function(x) max_growth_rate(trait, x, p)
  }
  out <- optimise(f, interval = bounds, maximum = TRUE, tol = 1E-3)
  if (log_scale)
    out$maximum <- exp(out$maximum)
  out
}
