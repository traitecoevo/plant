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

equilibrium_runner_cleanup <- function(runner, converged=TRUE) {
  e <- environment(runner)
  res <- e$last
  attr(res, "progress") <- e$history
  attr(res, "converged") <- converged
  res
}

equilibrium_seed_rain2 <- function(p, schedule_default=NULL,
                                   schedule_initial=NULL) {
  solver <- get_equilibrium_solver(p)
  if (p$control$parameters$equilibrium_verbose) {
    message("Solving seed rain using ", solver)
  }
  if (solver == "iteration") {
    equilibrium_seed_rain_iteration(p, schedule_default, schedule_initial)
  } else if (solver == "nleqslv") {
    equilibrium_seed_rain_solve_robust(p, schedule_default, schedule_initial,
                                       solver)
  } else if (solver == "runsteady") {
    equilibrium_seed_rain_runsteady(p, schedule_default, schedule_initial)
  } else if (solver == "hybrid") {
    equilibrium_seed_rain_hybrid(p, schedule_default, schedule_initial)
  } else {
    stop("Unknown solver ", solver)
  }
}

equilibrium_seed_rain_hybrid <- function(p, schedule_default=NULL,
                                         schedule_initial=NULL) {
  p <- p$copy()
  control <- p$control$parameters
  ## Things to help with the solver.
  solver <- "nleqslv"
  try_keep <- TRUE
  logN <- TRUE

  ans_it <- equilibrium_seed_rain_iteration(p, schedule_default,
                                            schedule_initial)
  p$seed_rain <- ans_it$seed_rain[,"out"]
  if (!isTRUE(attr(ans_it, "converged"))) {
    message("Iteration approach did not succeed: this may end badly")
  }

  ans_sol <- try(equilibrium_seed_rain_solve(p, schedule_default,
                                             ans_it$schedule,
                                             solver, try_keep, logN))

  if (inherits(ans_sol, "try-error") || !isTRUE(attr(ans_sol, "converged"))) {
    ret <- ans_it
  } else {
    ret <- ans_sol
    if (any(ans_sol$seed_rain[,"out"] == 0.0)) {
      message("Checking species driven to extinction")
      ## Add these species back at extremely low density and make sure
      ## that this looks like a legit extinction.
      y_in <- ans_sol$seed_rain[,"out"]
      i <- y_in <= 0.0
      y_in[i] <- control$equilibrium_extinct_seed_rain
      p2 <- p$copy()
      p2$seed_rain <- y_in
      y_out <- run_ebt(p2, ans_sol$schedule)$seed_rains
      if (any(y_out[i] > y_in[i])) {
        message("Solver drove viable species extinct: rejecting")
        ret <- ans_it
      }
    }
  }
  ret
}

equilibrium_seed_rain_iteration <- function(p, schedule_default=NULL,
                                            schedule_initial=NULL) {
  runner <- make_equilibrium_runner(p, schedule_default, schedule_initial)

  control <- p$control$parameters
  eps <- control$equilibrium_eps
  seed_rain <- p$seed_rain

  for (i in seq_len(control$equilibrium_nsteps)) {
    ans <- runner(seed_rain)
    seed_rain <- ans[,"out"]
    achange <- ans[,"out"] - ans[,"in"]
    rchange <- 1 - ans[,"out"] / ans[,"in"]
    converged <- eps > 0 && all(abs(achange) < eps | abs(rchange) < eps)
    if (converged) {
      if (control$equilibrium_verbose) {
        fmt <- "Reached target accuracy (delta %2.5e, %2.5e < %2.5e eps)"
        message(sprintf(fmt, max(abs(achange)), max(abs(rchange)), eps))
      }
      break
    }
  }

  equilibrium_runner_cleanup(runner, converged)
}

equilibrium_seed_rain_solve <- function(p, schedule_default=NULL,
                                        schedule_initial=NULL,
                                        solver="nleqslv",
                                        try_keep=TRUE,
                                        logN=FALSE) {
  seed_rain <- p$seed_rain
  runner <- make_equilibrium_runner(p, schedule_default, schedule_initial)
  if (try_keep) {
    ans <- runner(seed_rain)
    keep <- unname(ans[,"out"] >= ans[,"in"])
    message("Keeping species: ", paste(which(keep), collapse=", "))
  } else {
    keep <- rep(FALSE, p$size)
  }
  target <- equilibrium_seed_rain_solve_target(runner, keep, logN)

  tol <- p$control$parameters$equilibrium_eps
  ## NOTE: Hard coded minimum of 100 steps here.
  maxit <- max(100,
               p$control$parameters$equilibrium_nsteps)
  x0 <- if (logN) log(p$seed_rain) else p$seed_rain
  sol <- nlsolve(x0, target, tol=tol, maxit=maxit, solver=solver)
  res <- equilibrium_runner_cleanup(runner)
  attr(res, "sol") <- sol
  res
}

## Simple version of the solve function that does not try any clever
## business with 'keep'.  As such it's only good when close to the
## true root.
equilibrium_seed_rain_solve_simple <- function(p, schedule_default=NULL,
                                               schedule_initial=NULL,
                                               solver="nleqslv") {
  browser()
  runner <- make_equilibrium_runner(p, schedule_default, schedule_initial)
  target <- equilibrium_seed_rain_solve_target_simple(runner)

  tol <- p$control$parameters$equilibrium_eps
  maxit <- p$control$parameters$equilibrium_nsteps
  sol <- nlsolve(p$seed_rain, target, tol=tol, maxit=maxit, solver=solver)
  tmp <- nleqslv::nleqslv(p$seed_rain, target)

  res <- equilibrium_runner_cleanup(runner)
  attr(res, "sol") <- sol
  res
}

## Simple version of the solve function that does not try any clever
## business with 'keep'.  As such it's only good when close to the
## true root.
equilibrium_seed_rain_solve_log <- function(p, schedule_default=NULL,
                                            schedule_initial=NULL,
                                            solver="nleqslv") {
  runner <- make_equilibrium_runner(p, schedule_default, schedule_initial)
  target <- equilibrium_seed_rain_solve_target_log(runner)

  tol <- p$control$parameters$equilibrium_eps
  maxit <- p$control$parameters$equilibrium_nsteps
  sol <- nlsolve(log(p$seed_rain), target, tol=tol, maxit=maxit, solver=solver)
  res <- equilibrium_runner_cleanup(runner)
  attr(res, "sol") <- sol
  res
}



equilibrium_seed_rain_solve_robust <- function(p,
                                               schedule_default=NULL,
                                               schedule_initial=NULL,
                                               solver="nleqslv") {
  verbose <- p$control$parameters$equilibrium_verbose
  message_verbose <- function(...) {
    if (verbose) {
      message(...)
    }
  }

  fit <- try(equilibrium_seed_rain_solve(p, schedule_default,
                                         schedule_initial,
                                         solver))
  if (failed(fit)) {
    message_verbose("Falling back on iteration")
    fit <- equilibrium_seed_rain_iteration(p, schedule_default,
                                           schedule_initial)

    message_verbose("Trying again with the solver")
    p2 <- p$copy()
    p2$seed_rain <- unname(fit$seed_rain[,"out"])
    fit2 <- try(equilibrium_seed_rain_solve(p2, schedule_default,
                                            fit$schedule,
                                            solver))

    if (failed(fit2)) {
      message_verbose("Solver failed again: using iteration version")
    } else {
      fit <- fit2
      message_verbose("Solver worked on second attempt")
    }
  }
  fit
}

## This stuff is way uglier than it really should be.
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
           force_new_schedule=FALSE, force_old_schedule=FALSE,
           t=NULL) {
    seed_rain_last <- p$seed_rain
    p$seed_rain <- seed_rain_in
    if (force_old_schedule) {
      schedule <- last_schedule
      stop("This does not work!") # run
    } else {
      if (force_new_schedule ||
          any(abs(seed_rain_in - seed_rain_last) > large_seed_rain_change)) {
        schedule <- schedule_default$copy()
      } else {
        schedule <- last_schedule$copy()
      }
      schedule <- build_schedule(p, schedule)
    }
    last_schedule <<- schedule

    ans <- attr(schedule, "seed_rain", exact=TRUE)
    last <<- list(seed_rain=drop_schedule(ans),
                  schedule=schedule$copy(),
                  t=t)

    if (progress) {
      history <<- c(history, list(last[c("seed_rain", "t")]))
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

equilibrium_seed_rain_solve_target <- function(runner, keep, logN) {
  eps <- 1e-10
  force(runner)
  force(keep)
  force(logN)
  function(x, ...) {
    if (logN) {
      x <- exp(x)
    }
    ## Avoid negative seed rains:
    x[x < eps & keep] <- eps
    x[x < 0.0] <- 0.0
    xout <- unname(runner(x)[,"out"])

    xout[keep] <- xout[keep] / x[keep] - 1.0
    i <- !keep & x > 0
    xout[i] <- xout[i] - x[i]

    ## !keep & x <= 0 will now  be zero
    if (any(xout[!keep & x <= 0] != 0.0)) {
      warning("This is not expected", immediate.=TRUE)
    }
    xout
  }
}

## Simplified version of the ODE solver target:
equilibrium_seed_rain_solve_target_simple <- function(runner) {
  eps <- 1e-10
  force(runner)
  function(x, ...) {
    pos <- x > eps
    if (!any(pos)) {
      message("All species extinct?")
    }
    x[!pos] <- 0.0
    res <- runner(x)
    ret <- rep(0, length(x))
    ret[pos] <- unname(res[pos,"out"] - res[pos,"in"])
    ret
  }
}

## Simplified version of the ODE solver target:
equilibrium_seed_rain_solve_target_log <- function(runner) {
  eps <- 1e-10
  force(runner)
  function(lx, ...) {
    x <- exp(lx)
    pos <- x > eps
    if (!any(pos)) {
      message("All species extinct?")
    }
    x[!pos] <- 0.0
    res <- runner(x)
    ret <- rep(0, length(x))
    ret[pos] <- unname(res[pos,"out"] - res[pos,"in"])
    ret
  }
}

equilibrium_seed_rain_runsteady <- function(p, schedule_default=NULL,
                                            schedule_initial=NULL) {
  control <- p$control$parameters
  eps <- control$equilibrium_eps
  ode_tol <- control$equilibrium_runsteady_tol
  logN <- TRUE # for now, always work with d (log N) / dt

  runner <- make_equilibrium_runner(p, schedule_default, schedule_initial)
  f <- make_target_runsteady(runner, logN)
  x0 <- if (logN) log(p$seed_rain) else p$seed_rain
  ans <- rootSolve::runsteady(x0, func=f, parms=NULL, mf=22,
                              rtol=ode_tol, atol=ode_tol, stol=eps)
  equilibrium_runner_cleanup(runner, attr(ans, "steady"))
}

## We want
##
##   dN / dt = (exp(m) - 1) N
##
## Where m is fitness.  Because exp(m) = out / in
##
##   dN / dt = (out / in - 1) N
##   dN / dt = out - in
##
## which is raw seed rain.
make_target_runsteady <- function(f, logN=FALSE) {
  force(f)
  force(logN)
  eps <- 1e-10
  function(t, x, ...) {
    if (logN) {
      x <- exp(x)
    }
    pos <- x > eps
    if (!any(pos)) {
      message("All species extinct?")
    }
    x[!pos] <- 0.0
    res <- f(x, t=t)
    ret <- rep(0, length(x))
    if (logN) {
      ret[pos] <- unname(res[pos,"out"] / res[pos,"in"]) - 1.0
    } else {
      ret[pos] <- unname(res[pos,"out"] - res[pos,"in"])
    }
    list(ret)
  }
}

equilibrium_solvers <- function() {
  c("iteration", "nleqslv", "runsteady", "hybrid")
}

equilibrium_solver_name <- function(i) {
  equilibrium_solvers()[[i]]
}
equilibrium_solver_code <- function(name) {
  i <- match(name, equilibrium_solvers())
  if (is.na(i)) {
    stop("Solver not found")
  }
  i
}

##' @export
set_equilibrium_solver <- function(name, p) {
  code <- equilibrium_solver_code(name)
  p$set_control_parameters(list(equilibrium_solver=code))
}
##' @export
get_equilibrium_solver <- function(p) {
  equilibrium_solver_name(p$control$parameters$equilibrium_solver)
}

##' Check low-abundance strategies for viability.
##'
##' @title Check low-abundance strategies for viability
##' @param p A Parameters object
##' @export
check_inviable <- function(p) {
  ## eps_test: *Relative* value to use for determining what
  ## "low abundance" means.  Species that have a seed rain of less than
  ## `eps_test * max(p$seed_rain)` will be tested.  By default
  ##' this is 1 100th of the maximum seed rain.
  eps <- p$control$parameters$equilibrium_extinct_seed_rain
  eps_test <- p$control$parameters$equilibrium_inviable_test
  seed_rain <- p$seed_rain
  ## TODO: This should take and pass in the cohort schedule here.
  eq <- make_equilibrium_runner(p)
  res <- eq(seed_rain)

  test <- which(res[,"out"] < res[,"in"] &
                seed_rain < max(seed_rain) * eps_test)
  test <- test[order(seed_rain[test])]

  drop <- logical(length(seed_rain))
  ret <- res

  for (i in test) {
    message("Testing species ", i)
    x <- ret[,"out"]
    x[i] <- eps
    res <- eq(x)
    if (diff(res[i,]) < 0) {
      message("\t...removing")
      drop[i] <- TRUE
      ret[i,"out"] <- 0.0
      ret <- res
    }
  }

  attr(ret, "drop") <- drop
  ret
}
