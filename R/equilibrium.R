##' Determine the equilibrium seed rain for a set of strategies.
##'
##' This does not use any clever root finding algorithm to try and
##' speed up convergence, but simply iterates the solver.  Something
##' like the secant method would speed this up greatly and help in
##' divergent cases.
##'
##' @title Determine Equilibrium Seed Rain
##' @param p Parameters object
##' @param schedule Schedule to run EBT with, or \code{NULL} for a
##' default schedule (should generally be a reasonable choice as it's
##' only the starting seed for the usual schedule building with
##' \code{\link{build_schedule}}).
##' @author Rich FitzJohn
##' @export
equilibrium_seed_rain <- function(p, schedule=NULL) {
  if (p$size == 0) {
    stop("Need at least one species")
  }
  p <- p$copy() # don't modify what we're given
  schedule <- schedule_or_null(schedule, p)
  ## Might revert back to default periodically.
  schedule_default <- schedule$copy()

  control <- p$control$parameters
  eps <- control$equilibrium_eps
  progress <- control$equilibrium_progress
  large_seed_rain_change <- control$equilibrium_large_seed_rain_change
  history <- list()
  for (i in seq_len(control$equilibrium_nsteps)) {
    schedule <- build_schedule(p, schedule)
    res <- list(seed_rain = attr(schedule, "seed_rain", exact=TRUE),
                schedule  = schedule$copy())
    seed_rain <- res[["seed_rain"]]
    change <- seed_rain[,"out"] - seed_rain[,"in"]

    if (progress) {
      history <- c(history, list(res))
    }

    p$seed_rain <- seed_rain[,"out"]
    if (any(abs(change) > large_seed_rain_change)) {
      schedule <- schedule_default$copy()
    }

    if (control$equilibrium_verbose) {
      message(sprintf("*** %d: {%s} -> {%s} (delta = {%s})", i,
                      paste(prettyNum(seed_rain[,"in"]), collapse=","),
                      paste(prettyNum(seed_rain[,"out"]), collapse=","),
                      paste(prettyNum(change), collapse=",")))
    }

    if (eps > 0 && all(abs(change) < eps)) {
      if (control$equilibrium_verbose) {
        message(sprintf("Reached target accuracy (delta %2.5e < %2.5e eps)",
                        max(abs(change)), eps))
      }
      break
    }
  }
  if (progress) {
    attr(res, "progress") <- history
  }
  attr(res, "schedule") <- schedule
  res
}

##' Control parameters for \code{\link{equilibrium_seed_rain}} that
##' make progress noisier.  This is just a convenience function.
##'
##' @title Noisy Parameters for Equilibrium Finding
##' @export
##' @examples
##' p <- new(Parameters)
##' p$set_control_parameters(equilibrium_verbose())
equilibrium_verbose <- function() {
  list(schedule_verbose=TRUE,
       equilibrium_verbose=TRUE,
       equilibrium_progress=TRUE)
}
