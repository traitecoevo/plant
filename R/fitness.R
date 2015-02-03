##' Construct a fitness landscape.
##'
##' @title Fitness Landscape
##' @param trait_matrix A matrix of traits corresponding to mutants to
##' introduce into the light environment constructed by the residents
##' in \code{p}.
##' @param p Parameters object.  Needs to contain residents with their
##' incoming seed rain.
##' @return Vector with the output seed rain.  Mutants have an
##' arbitrary seed rain of one, so this is the rate of seed
##' production per capita.
##' @author Rich FitzJohn
##' @export
fitness_landscape <- function(trait_matrix, p, raw_seed_rain=FALSE) {
  n_residents <- length(p$strategies)
  if (n_residents == 0L &&
      length(p$cohort_schedule_ode_times) == 0L) {
    p$cohort_schedule_ode_times <-
      c(p$cohort_schedule_times_default, p$cohort_schedule_max_time)
  }

  p_with_mutants <- expand_parameters(trait_matrix, p)

  ebt <- run_ebt(p_with_mutants,
                 use_ode_times=length(p$cohort_schedule_ode_times) > 0)
  seed_rain <- ebt$seed_rains
  if (n_residents > 0L) {
    seed_rain <- seed_rain[-seq_len(n_residents)]
  }
  if (raw_seed_rain) {
    seed_rain
  } else {
    log(seed_rain)
  }
}

##' Compute max growth rate of a given set of values for a trait.
##' This is the log of per-capita seed production (i.e., fitness).
##'
##' Only works in one dimension
##' @title Compute Max Growth Rate
##' @param trait_matrix A matrix of traits to compute maximum growth
##' rate at (see \code{\link{trait_matrix}})
##' @param p Parameters object to use.  Importantly, the
##' \code{strategy_default} element gets used here.
##' @author Rich FitzJohn
##' @export
max_growth_rate <- function(trait_matrix, p) {
  fitness_landscape(trait_matrix, remove_residents(p))
}

##' Compute the carrying capacity (equilibrium per-capita seed
##' production) for a set of values of a trait.  Each is considered in
##' isolation.
##'
##' @title Carrying Capacity
##' @param trait_matrix A matrix of traits
##' @param p Parameters object to use.  Importantly, the
##' \code{strategy_default} element gets used here.
##' @param seed_rain Initial seed rain (optional)
##' @param parallel Use multiple processors?
##' @author Rich FitzJohn
##' @export
carrying_capacity <- function(trait_matrix, p, seed_rain=1,
                              parallel=FALSE) {
  f <- function(x) {
    p <- expand_parameters(x, remove_residents(p))
    p$seed_rain <- seed_rain
    equilibrium_seed_rain(p)$seed_rain
  }
  unlist(loop(matrix_to_list(trait_matrix), f, parallel=parallel))
}

##' Find point of maximum fitness in empty fitness landscape within a
##' specified range.
##'
##' @title Find point of maximum fitness within some range.
##' @param trait Name of the trait (e.g., \code{"lma"})
##' @param bounds Two element vector specifing range within which to
##' search
##' @param p Parameters object to use.  Importantly, the
##' \code{strategy_default} element gets used here.
##' @param log_scale Is the parameter naturally on a log scale?  If
##' so, this will greatly speed things up.
##' @export
##' @author Daniel Falster, Rich FitzJohn
max_fitness <- function(bounds, p, log_scale=TRUE, tol=1e-3) {
  bounds <- check_bounds(bounds)
  traits <- rownames(bounds)

  if (length(traits) > 1L) {
    stop("Doesn't yet support multiple traits")
  }

  if (log_scale) {
    bounds <- log(bounds)
    f <- function(x) max_growth_rate(exp(trait_matrix(x, traits)), p)
  } else {
    f <- function(x) max_growth_rate(trait_matrix(x, traits), p)
  }

  if (length(traits) == 1L) {
    out <- optimise(f, interval=bounds, maximum=TRUE, tol=tol)
    ret <- if (log_scale) exp(out$maximum) else out$maximum
    fitness <- out$objective
  } else {
    ## This is not very well tested, and the tolerance is not useful:
    out <- optim(rowMeans(bounds), f, method="L-BFGS-B",
                 lower=bounds[,1], upper=bounds[,2],
                 control=list(fnscale=-1, factr=1e10))
    ret <- if (log_scale) exp(out$ret) else out$ret
    fitness <- out$value
  }

  attr(ret, "fitness") <- fitness
  ret
}

##' Compute region of positive fitness.  This will be the values where
##' fitness is approximately zero.
##'
##' @title Compute Region of Positive Fitnes
##' @param bounds Matrix of bounds; two columns corresponding to the
##' lower and upper limits, each row corresponds to a trait (the name
##' will be used).
##' @param p Parameters object to use.  Importantly, the
##' \code{strategy_default} element gets used here.
##' @param x Initial value - If not given, then the value from the
##' default strategy within \code{p} is used.
##' @param log_scale Is the parameter naturally on a log scale?  If
##' so, this will greatly speed things up.  Can be a vector of length
##' `\code{nrow(bounds)}`
##' @param dx Amount to step the trait.  If \code{log_scale} is
##' \code{TRUE}, this is on a log scale.
##' @export
##' @author Rich FitzJohn
viable_fitness <- function(bounds, p, x=NULL, log_scale=TRUE, dx=1) {
  bounds <- check_bounds(bounds)
  traits <- rownames(bounds)
  if (is.null(x)) {
    x <- unlist(p$strategy_default[traits])
  }
  x <- check_point(x, bounds)
  w <- max_growth_rate(x, p)

  if (w < 0) {
    message("Starting value had negative fitness, looking for max")
    x <- max_fitness(bounds, p, bounds, log_scale)
    w <- attr(x, "fitness")
    message(sprintf("\t...found max fitness at %2.5f (w=%2.5f)", x, w))
  }

  if (log_scale) {
    bounds <- log(bounds)
    x <- log(x)
    f <- function(x) {
      max_growth_rate(exp(trait_matrix(x, traits)), p)
    }
  } else {
    f <- function(x) {
      max_growth_rate(trait_matrix(x, traits), p)
    }
  }

  if (length(traits) == 1L) {
    out <- positive_1d(f, x, dx, lower=bounds[,1], upper=bounds[,2])
    ret <- rbind(out, deparse.level=0)
  } else {
    ## TODO: It seems a crying shame to throws all this information
    ## away.  I'd like to create some sort of polygon around the
    ## positive points but need to get the bilinear interpolation
    ## done.  Until then, this is lost.
    out <- positive_2d(f, x, lower=bounds[,1], upper=bounds[,2])
    ret <- t(apply(out$points[out$labels,], 2, range))
  }
  if (log_scale) {
    ret <- exp(ret)
  }
  colnames(ret) <- c("lower", "upper")
  rownames(ret) <- traits
  ret
}
