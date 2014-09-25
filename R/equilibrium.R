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
##' @export
##' @rdname equilibrium_verbose
equilibrium_quiet <- function() {
  list(schedule_verbose=FALSE,
       equilibrium_verbose=FALSE,
       equilibrium_progress=FALSE)
}


## TODO: This is not yet used in the equlibrium search above, but it
## should be.
make_equilibrium_runner <- function(p, schedule_default=NULL,
                                    schedule_initial=NULL) {
  drop_schedule <- function(x) {
    attr(x, "schedule") <- NULL
    x
  }
  pretty_collapse <- function(x, collapse=", ") {
    paste(prettyNum(x), collapse=collapse)
  }
  p <- p$copy()
  schedule_default <- schedule_or_null(schedule_default, p)

  control <- p$control$parameters
  large_seed_rain_change <- control$equilibrium_large_seed_rain_change
  progress <- control$equilibrium_progress

  if (is.null(schedule_initial)) {
    last_schedule <- schedule_default$copy()
  } else {
    last_schedule <- schedule_initial$copy()
  }
  i <- 1L
  last <- NULL
  history <- NULL

  function(seed_rain_in,
           force_new_schedule=FALSE, force_old_schedule=FALSE) {
    p$seed_rain <- seed_rain_in
    if (force_old_schedule) {
      schedule <- last_schedule
    } else {
      if (force_new_schedule ||
          any(abs(seed_rain_in - p$seed_rain) > large_seed_rain_change)) {
        schedule <- schedule_default$copy()
      } else {
        schedule <- last_schedule$copy()
      }
      schedule <- build_schedule(p, schedule)
    }
    last_schedule <<- schedule

    ans <- attr(schedule, "seed_rain", exact=TRUE)
    last <<- list(seed_rain=drop_schedule(ans),
                  schedule=schedule$copy())
    if (progress) {
      history <<- c(history, list(last))
    }

    if (control$equilibrium_verbose) {
      message(sprintf("*** %d: {%s} -> {%s} (delta = {%s})", i,
                      pretty_collapse(ans[,"in"]),
                      pretty_collapse(ans[,"out"]),
                      pretty_collapse(ans[,"out"] - ans[,"in"])))
    }
    i <<- i + 1L

    attr(ans, "schedule") <- schedule$copy()
    ans
  }
}

equilibrium_runner_cleanup <- function(runner) {
  e <- environment(runner)
  res <- e$last
  attr(res, "progress") <- e$history
  res
}

equilibrium_seed_rain_iterate <- function(p, schedule_default=NULL,
                                          schedule_initial=NULL) {
  runner <- make_equilibrium_runner(p, schedule_default, schedule_initial)

  control <- p$control$parameters
  eps <- control$equilibrium_eps
  seed_rain <- p$seed_rain

  for (i in seq_len(control$equilibrium_nsteps)) {
    ans <- runner(seed_rain)
    seed_rain <- ans[,"out"]
    change <- ans[,"out"] - ans[,"in"]
    if (eps > 0 && all(abs(change) < eps)) {
      if (control$equilibrium_verbose) {
        message(sprintf("Reached target accuracy (delta %2.5e < %2.5e eps)",
                        max(abs(change)), eps))
      }
      break
    }
  }

  equilibrium_runner_cleanup(runner)
}

equilibrium_seed_rain_solve <- function(p, keep, schedule_default=NULL,
                                        schedule_initial=NULL,
                                        method="nleqslv") {
  keep <- recycle_simple(keep, p$size)

  runner <- make_equilibrium_runner(p, schedule_default, schedule_initial)
  target <- equilibrium_seed_rain_solve_target(runner, keep)

  tol <- p$control$parameters$equilibrium_eps
  maxit <- p$control$parameters$equilibrium_nsteps

  msg <- NULL
  if (method == "nleqslv") {
    control <- list(xtol=tol, ftol=tol, maxit=maxit)
    sol <- nleqslv::nleqslv(p$seed_rain, target, global="none",
                            control=control)
    if (sol$termcd > 2 || sol$termcd < 0) {
      msg <- sprintf("Equilibrium search has likely failed: code=%d, msg: %s",
                     sol$termcd, sol$message)
    }
  } else if (method == "dfsane") {
    control <- list(tol=tol, maxit=maxit)
    sol <- BB::dfsane(p$seed_rain, target, control=control, quiet=TRUE)
    if (sol$convergence != 0) {
      msg <- sprintf("Equilibrium search has likely failed: code=%d, msg: %s",
                     sol$convergence, sol$message)
    }
  } else {
    stop("Unknown method ", method)
  }
  if (!is.null(msg)) {
    warning(msg, immediate.=TRUE)
  }

  res <- equilibrium_runner_cleanup(runner)
  attr(res, "sol") <- sol
  res
}

equilibrium_seed_rain_solve_target <- function(runner, keep) {
  eps <- 1e-10
  force(runner)
  force(keep)
  function(x, ...) {
    ## Avoid negative seed rains:
    x[x < eps & keep] <- eps
    x[x < 0.0] <- 0.0
    xout <- unname(runner(x)[,"out"])

    xout[keep] <- log(xout[keep] / x[keep])
    i <- !keep & x > 0
    xout[i] <- log(xout[i] / x[i]) * x[i]
    ## !keep & x <= 0 will now  be zero
    if (any(xout[!keep & x <= 0] != 0.0)) {
      warning("This is not expected", immediate.=TRUE)
    }
    xout
  }

}
