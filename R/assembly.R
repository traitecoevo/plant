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
##' @param value Initial value - must have positive fitness itself!
##' @param p Parameters object to use.  Importantly, the
##' \code{strategy_default} element gets used here.
##' @param log_scale Is the parameter naturally on a log scale?  If
##' so, this will greatly speed things up.
##' @param dx Amount to step the trait.  If \code{log_scale} is
##' \code{TRUE}, this is on a log scale.
##' @export
##' @author Rich FitzJohn
viable_fitness <- function(trait, value, p, log_scale=TRUE, dx=1) {
  if (log_scale) {
    f <- function(x) {
      max_growth_rate(trait, exp(x), p)
    }
    exp(positive(f, log(value), dx))
  } else {
    f <- function(x) {
      max_growth_rate(trait, x, p)
    }
    positive(f, value, dx)
  }
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
    repeat {
      x_next <- x + dx
      if (dx < 0 && x_next < bound) {
        stop("Too low")
      } else if (dx > 0 && x_next > bound) {
        stop("Too high")
      }
      fx_next <- f(x_next)
      if (fx_next < 0) {
        if (dx < 0) {
          x <- c(x_next, x)
          fx <- c(fx_next, fx)
        } else {
          x <- c(x, x_next)
          fx <- c(fx, fx_next)
        }
        break
      } else {
        x <- x_next
        fx <- fx_next
        dx <- dx * factor
      }
    }
    list(x=x, fx=fx)
  }

  list(lower=bracket(x, -dx, lower), upper=bracket(x, dx, upper))
}

positive <- function(f, x, dx, lower=-Inf, upper=Inf, eps=1e-3) {
  b <- positive_bracket(f, x, dx, lower, upper)
  lower <- uniroot(f, b$lower$x,
                   f.lower=b$lower$fx[[1]], f.upper=b$lower$fx[[2]],
                   tol=eps)$root
  upper <- uniroot(f, b$upper$x,
                   f.lower=b$upper$fx[[1]], f.upper=b$upper$fx[[2]],
                   tol=eps)$root
  c(lower, upper)
}
