##' Run system to seed rain equilibrium
##'
##' @title Run system to seed rain equilibrium
##' @param p A \code{\link{Parameters}} object
##' @return A Parameters object, with seed rain and cohort schedule
##' elements set.
##' @export
##' @author Rich FitzJohn
equilibrium_seed_rain <- function(p) {
  solver <- p$control$equilibrium_solver_name
  if (p$control$equilibrium_verbose) {
    message("eq> Solving seed rain using ", solver)
  }
  switch(solver,
         iteration=equilibrium_seed_rain_iteration(p),
         nleqslv=equilibrium_seed_rain_solve_robust(p, solver),
         dfsane=equilibrium_seed_rain_solve_robust(p, solver),
         hybrid=equilibrium_seed_rain_hybrid(p),
         stop("Unknown solver ", solver))
}

## This is the simplest solver: it simply iterates the outgoing seed
## rain as incoming seed rain.  No attempt at projection is made.
equilibrium_seed_rain_iteration <- function(p) {
  check <- function(x_in, x_out, eps, verbose) {
    achange <- x_out - x_in
    rchange <- 1 - x_out / x_in
    ## eps > 0 && # <- this was a precondition - seems odd.
    converged <- all(abs(achange) < eps | abs(rchange) < eps)
    if (verbose && converged) {
      fmt <- "eq> Reached target accuracy (delta %2.5e, %2.5e < %2.5e eps)"
      message(sprintf(fmt, max(abs(achange)), max(abs(rchange)), eps))
    }
    converged
  }

  eps <- p$control$equilibrium_eps
  verbose <- p$control$equilibrium_verbose
  seed_rain <- p$seed_rain
  runner <- make_equilibrium_runner(p)

  for (i in seq_len(p$control$equilibrium_nsteps)) {
    seed_rain_out <- runner(seed_rain)
    converged <- check(seed_rain, seed_rain_out, eps, verbose)
    seed_rain <- seed_rain_out
    if (converged) {
      break
    }
  }

  equilibrium_runner_cleanup(runner, converged)
}

equilibrium_seed_rain_solve <- function(p, solver="nleqslv") {
  try_keep <- p$control$equilibrium_solver_try_keep
  logN <- p$control$equilibrium_solver_logN
  min_seed_rain <- 1e-10 # should also be in the controls?
  message_verbose <- make_message_verbose(p$control$equilibrium_verbose)

  message_verbose("eq> Solving seed rain using %s", solver)
  seed_rain <- p$seed_rain
  runner <- make_equilibrium_runner(p)

  ## First, we exclude species that have seed rains below some minimum
  ## level.
  to_drop <- seed_rain < min_seed_rain
  if (any(to_drop)) {
    i_keep <- which(!to_drop)
    message_verbose("eq> Species %s extinct: excluding from search",
                    paste(which(to_drop), collapse=" & "))
    seed_rain_full <- seed_rain
    seed_rain <- seed_rain_full[i_keep]

    runner_full <- runner
    runner <- function(x) {
      x_full <- rep(0, length(seed_rain_full))
      x_full[i_keep] <- x
      runner_full(x_full)[i_keep]
    }
  }

  ## Then see if any species should be retained:
  if (try_keep) {
    ans <- runner(seed_rain)
    keep <- unname(ans >= seed_rain)
    message_verbose("eq> Keeping species: %s",
                    paste(which(!to_drop)[keep], collapse=", "))
  } else {
    keep <- rep(FALSE, length(p$strategies))
  }

  ## TODO: This is annoying, but turns out to be a problem for getting
  ## the solution working nicely.
  max_seed_rain <- pmax(seed_rain * 100, 10000)
  target <- equilibrium_seed_rain_solve_target(runner, keep, logN,
                                               min_seed_rain, max_seed_rain,
                                               p$control$equilibrium_verbose)
  x0 <- if (logN) log(seed_rain) else seed_rain

  tol <- p$control$equilibrium_eps
  ## NOTE: Hard coded minimum of 100 steps here.
  maxit <- max(100,
               p$control$equilibrium_nsteps)
  sol <- nlsolve(x0, target, tol=tol, maxit=maxit, solver=solver)
  res <- equilibrium_runner_cleanup(runner, attr(sol, "converged"))
  attr(res, "sol") <- sol
  res
}

equilibrium_seed_rain_solve_robust <- function(p, solver="nleqslv") {
  message_verbose <- make_message_verbose(p$control$equilibrium_verbose)

  fit <- try(equilibrium_seed_rain_solve(p, solver))
  if (failed(fit)) {
    message_verbose("eq> Falling back on iteration")
    fit_it <- equilibrium_seed_rain_iteration(p)
    ## Here, we might want to pick a set of values that look good from
    ## the previous set.

    message_verbose("eq> Trying again with the solver")
    p$seed_rain <- unname(fit_it$seed_rain)
    fit2 <- try(equilibrium_seed_rain_solve(p, solver))

    if (failed(fit2)) {
      message_verbose("eq> Solver failed again: using iteration version")
      fit <- fit_it
    } else {
      message_verbose("eq> Solver worked on second attempt")
      fit <- fit2
    }
  }
  fit
}

## The idea is to use rounds of iteration to try and push the
## system into the basin of attraction of the stable equilibrium.  The
## final approach is slow so use a root-finding approach there.
## However, if we are *not* in the basin of attraction the root finder
## will happily select zero seed rains for species that are not
## zeros.  So after running a round with the solver, check any species
## that were zeroed to make sure they're really dead.
equilibrium_seed_rain_hybrid <- function(p) {
  attempts <- p$control$equilibrium_nattempts
  message_verbose <- make_message_verbose(p$control$equilibrium_verbose)

  ## Then expand this out so that we can try alternating solvers
  solver <- rep(c("nleqslv", "dfsane"), length.out=attempts)

  for (i in seq_len(attempts)) {
    ans_it <- equilibrium_seed_rain_iteration(p)
    p$seed_rain <- ans_it$seed_rain
    converged_it <- isTRUE(attr(ans_it, "converged"))
    message_verbose("eq> Iteration %d %s",
                    i, if (converged_it) "converged" else "did not converge")

    ans_sol <- try(equilibrium_seed_rain_solve(p, solver[[i]]))
    converged_sol <- isTRUE(attr(ans_sol, "converged"))
    message_verbose("eq> Solve %d %s",
                    i, if (converged_sol) "converged" else "did not converge")

    if (converged_sol) {
      if (any(ans_sol$seed_rain == 0.0)) {
        message_verbose("eq> Checking species driven to extinction")
        ## Add these species back at extremely low density and make sure
        ## that this looks like a legit extinction.
        y_in <- ans_sol$seed_rain
        i <- y_in <= 0.0
        y_in[i] <- p$control$equilibrium_extinct_seed_rain

        p_check <- p
        p_check$seed_rain <- y_in
        y_out <- run_ebt(p_check)$seed_rains
        if (any(y_out[i] > y_in[i])) {
          message_verbose("eq> Solver drove viable species extinct: rejecting")
          next
        }
      }
      message_verbose("eq> Accepting solution via solver")
      return(ans_sol)
    }
  }

  message_verbose("eq> Repeated rounds failed to find optimum")
  ans_it
}

## Support code:
make_equilibrium_runner <- function(p) {
  pretty_num_collapse <- function(x, collapse=", ") {
    paste0("{", paste(prettyNum(x), collapse=collapse), "}")
  }

  p <- validate(p)

  large_seed_rain_change <- p$control$equilibrium_large_seed_rain_change

  i <- 1L
  last_seed_rain <- p$seed_rain
  default_schedule_times <- rep(list(p$cohort_schedule_times_default),
                                length(p$seed_rain))
  last_schedule_times <- p$cohort_schedule_times
  history <- NULL

  function(seed_rain_in) {
    if (any(abs(seed_rain_in - last_seed_rain) > large_seed_rain_change)) {
      p$cohort_schedule_times <- default_schedule_times
    }

    p$seed_rain <- seed_rain_in

    p_new <- build_schedule(p)
    seed_rain_out <- attr(p_new, "seed_rain_out", exact=TRUE)

    ## These all write up to the containing environment:
    p <<- p_new
    last_schedule_times <<- p_new$cohort_schedule_times
    last_seed_rain      <<- seed_rain_in
    history <<- c(history, list(c("in"=seed_rain_in, out=seed_rain_out)))

    if (p$control$equilibrium_verbose) {
      message(sprintf("eq> %d: %s -> %s (delta = %s)", i,
                      pretty_num_collapse(seed_rain_in),
                      pretty_num_collapse(seed_rain_out),
                      pretty_num_collapse(seed_rain_out - seed_rain_in)))
    }
    i <<- i + 1L

    attr(seed_rain_out, "schedule_times") <- last_schedule_times
    seed_rain_out
  }
}

equilibrium_runner_cleanup <- function(runner, converged=TRUE) {
  ## This is super gross.
  e <- environment(runner)
  if (is.function(e$runner_full)) {
    runner <- e$runner_full
    e <- environment(runner)
  }

  p <- e$p
  p$seed_rain <- as.numeric(e$last_seed_rain)
  p$cohort_schedule_times <- e$last_schedule_times
  attr(p, "progress") <- rbind_list(e$history)
  attr(p, "converged") <- converged
  p
}

## Another layer of runner for the solver code:
equilibrium_seed_rain_solve_target <- function(runner, keep, logN,
                                               min_seed_rain, max_seed_rain,
                                               verbose) {
  force(runner)
  force(keep)
  force(logN)
  force(min_seed_rain)
  force(max_seed_rain)
  message_verbose <- make_message_verbose(verbose)
  function(x, ...) {
    if (logN) {
      x <- exp(x)
    }
    ## Avoid negative seed rains:
    x[x < min_seed_rain & keep] <- min_seed_rain
    x[x < min_seed_rain & !keep] <- 0.0
    if (!any(x > 0)) {
      message_verbose("eq> All species extinct?")
    }
    too_high <- x > max_seed_rain
    if (any(too_high)) {
      message_verbose("eq> Truncating seed rain of species %s",
                      paste(which(too_high), collapse=", "))
      if (length(max_seed_rain) == 1L) {
        max_seed_rain <- rep(max_seed_rain, length.out=length(x))
      }
      x[too_high] <- max_seed_rain[too_high]
    }
    xout <- unname(runner(x))

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

##' Check low-abundance strategies for viability.
##'
##' @title Check low-abundance strategies for viability
##' @param p A Parameters object
##' @param schedule_default Schedule to use by default
##' @param schedule_initial Schedule to start from
##' @export
check_inviable <- function(p) {
  ## eps_test: *Relative* value to use for determining what
  ## "low abundance" means.  Species that have a seed rain of less than
  ## `eps_test * max(p$seed_rain)` will be tested.  By default
  ##  this is 1 100th of the maximum seed rain.
  ## TODO: don't do anything if we don't have at least 2 species?
  eps <- p$control$equilibrium_extinct_seed_rain
  ## TODO: This was p$control$equilibrium_inviable_test, but I think
  ## that birth seed rain actually makes more sense?  It's fractional
  ## though so who knows.
  eps_test <- 1e-2
  seed_rain <- p$seed_rain
  ## NOTE: We don't actually run to equilibrium here; this is just
  ## because it's a useful way of doing incoming -> outgoing seed
  ## rain.
  runner <- make_equilibrium_runner(p)
  seed_rain_out <- runner(seed_rain)

  test <- which(seed_rain_out < seed_rain &
                seed_rain < max(seed_rain_out) * eps_test)
  test <- test[order(seed_rain_out[test])]

  drop <- logical(length(seed_rain_out))

  for (i in test) {
    message("check_inviable> Testing species ", i)
    x <- seed_rain_out
    x[i] <- eps
    res <- runner(x)
    if (res[[i]] < eps) {
      message("check_inviable>\t...removing")
      drop[[i]] <- TRUE
      res[[i]] <- 0.0
      seed_rain_out <- res
    }
  }

  ## It's possible that things slip through and get driven extinct by
  ## the time that they reach here.
  drop <- drop | seed_rain_out < eps

  attr(seed_rain_out, "drop") <- drop
  seed_rain_out
}
