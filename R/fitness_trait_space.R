## These functions involve "trait space" as a concept, and as such
## care about hyperparameters; ways of converting things like lma into
## a set of parameters.  This is handled by strategy_list, which all
## of these functions use.

##' Find point of maximum fitness in empty fitness landscape within a
##' specified range.
##'
##' @title Find point of maximum fitness within some range.
##' @param bounds Two element vector specifing range within which to
##' search
##' @param p Parameters object to use.  Importantly, the
##' \code{strategy_default} element gets used here.
##' @param log_scale Is the parameter naturally on a log scale?  If
##' so, this will greatly speed things up.
##' @param tol Tolerance used in the optimisation
##' @importFrom stats optimise optim 
##' @export
##' @author Daniel Falster, Rich FitzJohn
max_fitness <- function(bounds, p, log_scale=TRUE, tol=1e-3) {
  bounds <- check_bounds(bounds)
  traits <- rownames(bounds)

  if (log_scale) {
    bounds[bounds[,1] == -Inf, 1] <- 0
    bounds <- log(bounds)
    f <- function(x) max_growth_rate(exp(trait_matrix(x, traits)), p)
  } else {
    f <- function(x) max_growth_rate(trait_matrix(x, traits), p)
  }

  if (length(traits) == 1L) {
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
    out <- suppressWarnings(optimise(f, interval=bounds, maximum=TRUE, tol=tol))
    ret <- out$maximum
    fitness <- out$objective
  } else {
    ## This is not very well tested, and the tolerance is not useful:
    out <- optim(rowMeans(bounds), f, method="L-BFGS-B",
                 lower=bounds[, "lower"], upper=bounds[, "upper"],
                 control=list(fnscale=-1, factr=1e10))
    ret <- out$par
    fitness <- out$value
  }

  if (log_scale) {
    ret <- exp(ret)
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
    plant_log_viable("Starting value had negative fitness, looking for max")
    x <- max_fitness(bounds, p, log_scale)
    w <- attr(x, "fitness")
    plant_log_viable(sprintf("\t...found max fitness at %s (w=%2.5f)",
                             paste(formatC(x), collapse=", "), w))
    if (w < 0) {
      return(NULL)
    }
  }

  if (log_scale) {
    bounds[bounds[,1] == -Inf, 1] <- 0
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
