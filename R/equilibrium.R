##' Run system to offspring arrival equilibrium
##'
##' @title Run system to offspring arrival equilibrium
##' @param p A \code{Parameters} object
##' @return A Parameters object, with offspring arrival and cohort schedule
##' elements set.
##' @export
##' @author Rich FitzJohn
equilibrium_offspring_arriving <- function(p) {
  solver <- p$control$equilibrium_solver_name
  plant_log_info(sprintf("Solving offspring arrival using %s", solver),
                 routine="equilibrium", stage="start", solver=solver)
  switch(solver,
         iteration=equilibrium_offspring_arriving_iteration(p),
         nleqslv=equilibrium_offspring_arriving_solve_robust(p, solver),
         dfsane=equilibrium_offspring_arriving_solve_robust(p, solver),
         hybrid=equilibrium_offspring_arriving_hybrid(p),
         stop("Unknown solver ", solver))
}

## This is the simplest solver: it simply iterates the outgoing offspring
## produced as incoming offspring arrival.  No attempt at projection is made.
equilibrium_offspring_arriving_iteration <- function(p) {
  check <- function(x_in, x_out, eps, verbose) {
    achange <- x_out - x_in
    rchange <- 1 - x_out / x_in
    ## eps > 0 && # <- this was a precondition - seems odd.
    converged <- all(abs(achange) < eps | abs(rchange) < eps)
    if (converged) {
      achange <- max(abs(achange))
      rchange <- max(abs(rchange))
      fmt <- "Reached target accuracy (delta %2.5e, %2.5e < %2.5e eps)"
      plant_log_eq(sprintf(fmt, achange, rchange, eps),
                   stage="converged", achange=achange, rchange=rchange)
    }
    converged
  }

  eps <- p$control$equilibrium_eps
  verbose <- p$control$equilibrium_verbose
  offspring_arriving <- p$offspring_arriving
  runner <- make_equilibrium_runner(p)

  for (i in seq_len(p$control$equilibrium_nsteps)) {
    offspring_produced <- runner(offspring_arriving)
    converged <- check(offspring_arriving, offspring_produced, eps, verbose)
    offspring_arriving <- offspring_produced
    if (converged) {
      break
    }
  }

  equilibrium_runner_cleanup(runner, converged)
}

equilibrium_offspring_arriving_solve <- function(p, solver="nleqslv") {
  try_keep <- p$control$equilibrium_solver_try_keep
  logN <- p$control$equilibrium_solver_logN
  min_offspring_arriving <- 1e-10 # TODO: should also be in the controls?

  plant_log_eq(paste("Solving offspring arrival using", solver),
               stage="start", solver=solver)

  offspring_arriving <- p$offspring_arriving
  runner <- make_equilibrium_runner(p)

  ## First, we exclude species that have offspring arrivals below some minimum
  ## level.
  to_drop <- offspring_arriving < min_offspring_arriving
  if (any(to_drop)) {
    i_keep <- which(!to_drop)
    msg <- sprintf("Species %s extinct: excluding from search",
                   paste(which(to_drop), collapse=" & "))
    plant_log_eq(msg, stage="drop species", drop=which(to_drop))
    offspring_arriving_full <- offspring_arriving
    offspring_arriving <- offspring_arriving_full[i_keep]

    runner_full <- runner
    runner <- function(x) {
      x_full <- rep(0, length(offspring_arriving_full))
      x_full[i_keep] <- x
      runner_full(x_full)[i_keep]
    }
  }

  ## Then see if any species should be retained:
  if (try_keep) {
    ans <- runner(offspring_arriving)
    keep <- unname(ans >= offspring_arriving)

    msg <- sprintf("Keeping species %s",
                   paste(which(!to_drop)[keep], collapse=", "))
    plant_log_eq(msg, stage="keep species", keep=which(!to_drop)[keep])
  } else {
    keep <- rep(FALSE, length(p$strategies))
  }

  ## TODO: This is annoying, but turns out to be a problem for getting
  ## the solution working nicely.
  max_offspring_arriving <- pmax(offspring_arriving * 100, 10000)
  target <- equilibrium_offspring_arriving_solve_target(runner, keep, logN,
                                               min_offspring_arriving, max_offspring_arriving,
                                               p$control$equilibrium_verbose)
  x0 <- if (logN) log(offspring_arriving) else offspring_arriving

  tol <- p$control$equilibrium_eps
  ## NOTE: Hard coded minimum of 100 steps here.
  maxit <- max(100,
               p$control$equilibrium_nsteps)
  sol <- nlsolve(x0, target, tol=tol, maxit=maxit, solver=solver)
  res <- equilibrium_runner_cleanup(runner, attr(sol, "converged"))
  attr(res, "sol") <- sol
  res
}

## NOTE: I don't know that this is used?  Perhaps it is?
equilibrium_offspring_arriving_solve_robust <- function(p, solver="nleqslv") {
  fit <- try(equilibrium_offspring_arriving_solve(p, solver))
  if (failed(fit)) {
    plant_log_eq("Falling back on iteration")
    fit_it <- equilibrium_offspring_arriving_iteration(p)
    ## Here, we might want to pick a set of values that look good from
    ## the previous set.

    plant_log_eq("Trying again with the solver")
    p$offspring_arriving <- unname(fit_it$offspring_arriving)
    fit2 <- try(equilibrium_offspring_arriving_solve(p, solver))

    if (failed(fit2)) {
      plant_log_eq("Solver failed again: using iteration version")
      fit <- fit_it
    } else {
      plant_log_eq("Solver worked on second attempt")
      fit <- fit2
    }
  }
  fit
}

## The idea is to use rounds of iteration to try and push the
## system into the basin of attraction of the stable equilibrium.  The
## final approach is slow so use a root-finding approach there.
## However, if we are *not* in the basin of attraction the root finder
## will happily select zero offspring arrivals for species that are not
## zeros.  So after running a round with the solver, check any species
## that were zeroed to make sure they're really dead.
equilibrium_offspring_arriving_hybrid <- function(p) {
  attempts <- p$control$equilibrium_nattempts

  ## Then expand this out so that we can try alternating solvers
  solver <- rep(c("nleqslv", "dfsane"), length.out=attempts)

  for (i in seq_len(attempts)) {
    ans_it <- equilibrium_offspring_arriving_iteration(p)
    p$offspring_arriving <- ans_it$offspring_arriving
    converged_it <- isTRUE(attr(ans_it, "converged"))
    msg <- sprintf("Iteration %d %s",
                   i, if (converged_it) "converged" else "did not converge")
    plant_log_eq(msg, step="iteration", converged=converged_it, iteration=i)

    ans_sol <- try(equilibrium_offspring_arriving_solve(p, solver[[i]]))
    converged_sol <- isTRUE(attr(ans_sol, "converged"))
    msg <- sprintf("Solve %d %s",
                    i, if (converged_sol) "converged" else "did not converge")
    plant_log_eq(msg, step="solve", converged=converged_sol, iteration=i)

    if (converged_sol) {
      if (any(ans_sol$offspring_arriving == 0.0)) {
        plant_log_eq("Checking species driven to extinction")
        ## Add these species back at extremely low density and make sure
        ## that this looks like a legit extinction.
        y_in <- ans_sol$offspring_arriving
        i <- y_in <= 0.0
        y_in[i] <- p$control$equilibrium_extinct_offspring_arriving

        p_check <- p
        p_check$offspring_arriving <- y_in
        y_out <- run_scm(p_check)$all_offspring_produced
        if (any(y_out[i] > y_in[i])) {
          plant_log_eq("Solver drove viable species extinct: rejecting")
          next
        }
      }
      plant_log_eq("Accepting solution via solver")
      return(ans_sol)
    }
  }

  ## This one should be a warning?
  plant_log_eq("Repeated rounds failed to find optimum")
  ans_it
}

## Support code:
make_equilibrium_runner <- function(p) {
  pretty_num_collapse <- function(x, collapse=", ") {
    paste0("{", paste(prettyNum(x), collapse=collapse), "}")
  }

  p <- validate(p)

  large_offspring_arriving_change <- p$control$equilibrium_large_offspring_arriving_change

  i <- 1L
  last_offspring_arriving <- p$offspring_arriving
  default_schedule_times <- rep(list(p$cohort_schedule_times_default),
                                length(p$offspring_arriving))
  last_schedule_times <- p$cohort_schedule_times
  history <- NULL

  function(offspring_arriving) {
    if (any(abs(offspring_arriving - last_offspring_arriving) > large_offspring_arriving_change)) {
      p$cohort_schedule_times <- default_schedule_times
    }

    p$offspring_arriving <- last_offspring_arriving

    p_new <- build_schedule(p)
    offspring_produced <- attr(p_new, "offspring_produced", exact=TRUE)

    ## These all write up to the containing environment:
    p <<- p_new
    last_schedule_times <<- p_new$cohort_schedule_times
    offspring_arriving      <<- offspring_arriving
    history <<- c(history, list(c("in"=offspring_arriving, out=offspring_produced)))

    msg <- sprintf("eq> %d: %s -> %s (delta = %s)", i,
                   pretty_num_collapse(offspring_arriving),
                   pretty_num_collapse(offspring_produced),
                   pretty_num_collapse(offspring_produced - offspring_arriving))
    plant_log_eq(msg,
                 stage="runner",
                 iteration=i,
                 offspring_arriving=offspring_arriving,
                 offspring_produced=offspring_produced)
    i <<- i + 1L

    attr(offspring_produced, "schedule_times") <- last_schedule_times
    offspring_produced
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
  p$offspring_arriving <- as.numeric(e$last_offspring_arriving)
  p$cohort_schedule_times <- e$last_schedule_times
  attr(p, "progress") <- rbind_list(e$history)
  attr(p, "converged") <- converged
  p
}

## Another layer of runner for the solver code:
equilibrium_offspring_arriving_solve_target <- function(runner, keep, logN,
                                               min_offspring_arriving, max_offspring_arriving,
                                               verbose) {
  force(runner)
  force(keep)
  force(logN)
  force(min_offspring_arriving)
  force(max_offspring_arriving)
  function(x, ...) {
    if (logN) {
      x <- exp(x)
    }
    ## TODO: most of the plant_log_eq things here should be DEBUG not INFO?
    ## Avoid negative offspring arrivals:
    x[x < min_offspring_arriving & keep] <- min_offspring_arriving
    x[x < min_offspring_arriving & !keep] <- 0.0
    if (!any(x > 0)) {
      plant_log_eq("All species extinct?")
    }
    too_high <- x > max_offspring_arriving
    if (any(too_high)) {
      plant_log_eq(sprintf("Truncating offspring arrival of species %s",
                           paste(which(too_high), collapse=", ")))
      if (length(max_offspring_arriving) == 1L) {
        max_offspring_arriving <- rep(max_offspring_arriving, length.out=length(x))
      }
      x[too_high] <- max_offspring_arriving[too_high]
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
##' @export
check_inviable <- function(p) {
  ## eps_test: *Relative* value to use for determining what
  ## "low abundance" means.  Species that have a offspring arrival of less than
  ## `eps_test * max(p$offspring_arriving)` will be tested.  By default
  ##  this is 1 100th of the maximum offspring arrival.
  ## TODO: don't do anything if we don't have at least 2 species?
  eps <- p$control$equilibrium_extinct_offspring_arriving
  ## TODO: This was p$control$equilibrium_inviable_test, but I think
  ## that birth offspring arrival actually makes more sense?  It's fractional
  ## though so who knows.
  eps_test <- 1e-2
  offspring_arriving <- p$offspring_arriving
  ## NOTE: We don't actually run to equilibrium here; this is just
  ## because it's a useful way of doing incoming -> outgoing offspring
  ## rain.
  runner <- make_equilibrium_runner(p)
  offspring_produced <- runner(offspring_produced)

  test <- which(offspring_produced < offspring_arriving &
                offspring_arriving < max(offspring_produced) * eps_test)
  test <- test[order(offspring_produced[test])]

  drop <- logical(length(offspring_produced))

  for (i in test) {
    plant_log_inviable(paste("Testing species", i),
                       stage="testing", species=i)
    x <- offspring_produced
    x[i] <- eps
    res <- runner(x)
    if (res[[i]] < eps) {
      plant_log_inviable(paste("Removing species", i),
                         stage="removing", species=i)
      drop[[i]] <- TRUE
      res[[i]] <- 0.0
      offspring_produced <- res
    }
  }

  ## It's possible that things slip through and get driven extinct by
  ## the time that they reach here.
  drop <- drop | offspring_produced < eps

  attr(offspring_produced, "drop") <- drop
  offspring_produced
}
