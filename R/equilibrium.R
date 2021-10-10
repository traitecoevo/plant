##' Run system to offspring arrival equilibrium
##'
##' @title Run system to offspring arrival equilibrium
##' @param p A \code{Parameters} object
##' @param ctrl Control object
##' @return A Parameters object, with offspring arrival and cohort schedule
##' elements set.
##' @export
##' @author Rich FitzJohn
equilibrium_birth_rate <- function(p, ctrl) {
  solver <- ctrl$equilibrium_solver_name
  plant_log_info(sprintf("Solving offspring arrival using %s", solver),
                 routine="equilibrium", stage="start", solver=solver)
  switch(solver,
         iteration=equilibrium_birth_rate_iteration(p),
         nleqslv=equilibrium_birth_rate_solve_robust(p, solver),
         dfsane=equilibrium_birth_rate_solve_robust(p, solver),
         hybrid=equilibrium_birth_rate_hybrid(p),
         stop("Unknown solver ", solver))
}

## This is the simplest solver: it simply iterates the outgoing offspring
## produced as incoming offspring arrival.  No attempt at projection is made.
equilibrium_birth_rate_iteration <- function(p, ctrl) {
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

  eps <- ctrl$equilibrium_eps
  verbose <- ctrl$equilibrium_verbose
  birth_rate <- p$birth_rate
  runner <- make_equilibrium_runner(p)

  for (i in seq_len(ctrl$equilibrium_nsteps)) {
    net_reproduction_ratios <- runner(birth_rate)
    converged <- check(birth_rate, net_reproduction_ratios, eps, verbose)
    birth_rate <- net_reproduction_ratios
    if (converged) {
      break
    }
  }

  equilibrium_runner_cleanup(runner, converged)
}

equilibrium_birth_rate_solve <- function(p, ctrl = scm_base_control(),
                                         solver="nleqslv") {
  try_keep <- ctrl$equilibrium_solver_try_keep
  logN <- ctrl$equilibrium_solver_logN
  min_offspring_arriving <- 1e-10 # TODO: should also be in the controls?

  plant_log_eq(paste("Solving offspring arrival using", solver),
               stage="start", solver=solver)

  birth_rate <- p$birth_rate
  runner <- make_equilibrium_runner(p)

  ## First, we exclude species that have offspring arrivals below some minimum
  ## level.
  to_drop <- birth_rate < min_offspring_arriving
  if (any(to_drop)) {
    i_keep <- which(!to_drop)
    msg <- sprintf("Species %s extinct: excluding from search",
                   paste(which(to_drop), collapse=" & "))
    plant_log_eq(msg, stage="drop species", drop=which(to_drop))
    offspring_arriving_full <- birth_rate
    birth_rate <- offspring_arriving_full[i_keep]

    runner_full <- runner
    runner <- function(x) {
      x_full <- rep(0, length(offspring_arriving_full))
      x_full[i_keep] <- x
      runner_full(x_full)[i_keep]
    }
  }

  ## Then see if any species should be retained:
  if (try_keep) {
    ans <- runner(birth_rate)
    keep <- unname(ans >= birth_rate)

    msg <- sprintf("Keeping species %s",
                   paste(which(!to_drop)[keep], collapse=", "))
    plant_log_eq(msg, stage="keep species", keep=which(!to_drop)[keep])
  } else {
    keep <- rep(FALSE, length(p$strategies))
  }

  ## TODO: This is annoying, but turns out to be a problem for getting
  ## the solution working nicely.
  max_offspring_arriving <- pmax(birth_rate * 100, 10000)
  target <- equilibrium_birth_rate_solve_target(runner, keep, logN,
                                               min_offspring_arriving, max_offspring_arriving,
                                               ctrl$equilibrium_verbose)
  x0 <- if (logN) log(birth_rate) else birth_rate

  tol <- ctrl$equilibrium_eps
  ## NOTE: Hard coded minimum of 100 steps here.
  maxit <- max(100,
               p$control$equilibrium_nsteps)
  sol <- nlsolve(x0, target, tol=tol, maxit=maxit, solver=solver)
  res <- equilibrium_runner_cleanup(runner, attr(sol, "converged"))
  attr(res, "sol") <- sol
  res
}

## NOTE: I don't know that this is used?  Perhaps it is?
equilibrium_birth_rate_solve_robust <- function(p, solver="nleqslv") {
  fit <- try(equilibrium_birth_rate_solve(p, solver))
  if (failed(fit)) {
    plant_log_eq("Falling back on iteration")
    fit_it <- equilibrium_birth_rate_iteration(p)
    ## Here, we might want to pick a set of values that look good from
    ## the previous set.

    plant_log_eq("Trying again with the solver")
    p$birth_rate <- unname(fit_it$birth_rate)
    fit2 <- try(equilibrium_birth_rate_solve(p, solver))

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
equilibrium_birth_rate_hybrid <- function(p, ctrl) {
  attempts <- ctrl$equilibrium_nattempts

  ## Then expand this out so that we can try alternating solvers
  solver <- rep(c("nleqslv", "dfsane"), length.out=attempts)

  for (i in seq_len(attempts)) {
    ans_it <- equilibrium_birth_rate_iteration(p)
    p$birth_rate <- ans_it$birth_rate
    converged_it <- isTRUE(attr(ans_it, "converged"))
    msg <- sprintf("Iteration %d %s",
                   i, if (converged_it) "converged" else "did not converge")
    plant_log_eq(msg, step="iteration", converged=converged_it, iteration=i)

    ans_sol <- try(equilibrium_birth_rate_solve(p, solver[[i]]))
    converged_sol <- isTRUE(attr(ans_sol, "converged"))
    msg <- sprintf("Solve %d %s",
                    i, if (converged_sol) "converged" else "did not converge")
    plant_log_eq(msg, step="solve", converged=converged_sol, iteration=i)

    if (converged_sol) {
      if (any(ans_sol$birth_rate == 0.0)) {
        plant_log_eq("Checking species driven to extinction")
        ## Add these species back at extremely low density and make sure
        ## that this looks like a legit extinction.
        y_in <- ans_sol$birth_rate
        i <- y_in <= 0.0
        y_in[i] <- ctrl$equilibrium_extinct_offspring_arriving

        p_check <- p
        p_check$birth_rate <- y_in
        y_out <- run_scm(p_check)$net_reproduction_ratios
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
make_equilibrium_runner <- function(p, ctrl) {
  pretty_num_collapse <- function(x, collapse=", ") {
    paste0("{", paste(prettyNum(x), collapse=collapse), "}")
  }

  p <- validate(p)

  large_offspring_arriving_change <- ctrl$equilibrium_large_offspring_arriving_change

  i <- 1L
  last_offspring_arriving <- p$birth_rate
  default_schedule_times <- rep(list(p$cohort_schedule_times_default),
                                length(p$birth_rate))
  last_schedule_times <- p$cohort_schedule_times
  history <- NULL

  function(birth_rate) {
    if (any(abs(birth_rate - last_offspring_arriving) > large_offspring_arriving_change)) {
      p$cohort_schedule_times <- default_schedule_times
    }

    p$birth_rate <- last_offspring_arriving

    p_new <- build_schedule(p, ctrl)
    net_reproduction_ratios <- attr(p_new, "net_reproduction_ratios", exact=TRUE)

    ## These all write up to the containing environment:
    p <<- p_new
    last_schedule_times <<- p_new$cohort_schedule_times
    birth_rate      <<- birth_rate
    history <<- c(history, list(c("in"=birth_rate, "out"=net_reproduction_ratios)))

    msg <- sprintf("eq> %d: %s -> %s (delta = %s)", i,
                   pretty_num_collapse(birth_rate),
                   pretty_num_collapse(net_reproduction_ratios),
                   pretty_num_collapse(net_reproduction_ratios - birth_rate))
    plant_log_eq(msg,
                 stage="runner",
                 iteration=i,
                 birth_rate=birth_rate,
                 net_reproduction_ratios=net_reproduction_ratios)
    i <<- i + 1L

    attr(net_reproduction_ratios, "schedule_times") <- last_schedule_times
    net_reproduction_ratios
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
  p$birth_rate <- as.numeric(e$last_offspring_arriving)
  p$cohort_schedule_times <- e$last_schedule_times
  attr(p, "progress") <- rbind_list(e$history)
  attr(p, "converged") <- converged
  p
}

## Another layer of runner for the solver code:
equilibrium_birth_rate_solve_target <- function(runner, keep, logN,
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
##' @param ctrl Control object
##' @export
check_inviable <- function(p, ctrl) {
  ## eps_test: *Relative* value to use for determining what
  ## "low abundance" means.  Species that have a offspring arrival of less than
  ## `eps_test * max(p$birth_rate)` will be tested.  By default
  ##  this is 1 100th of the maximum offspring arrival.
  ## TODO: don't do anything if we don't have at least 2 species?
  eps <- ctrl$equilibrium_extinct_offspring_arriving
  ## TODO: This was ctrl$equilibrium_inviable_test, but I think
  ## that birth offspring arrival actually makes more sense?  It's fractional
  ## though so who knows.
  eps_test <- 1e-2
  birth_rate <- p$birth_rate
  ## NOTE: We don't actually run to equilibrium here; this is just
  ## because it's a useful way of doing incoming -> outgoing offspring
  ## rain.
  runner <- make_equilibrium_runner(p)
  net_reproduction_ratios <- runner(net_reproduction_ratios)

  test <- which(net_reproduction_ratios < birth_rate &
                birth_rate < max(net_reproduction_ratios) * eps_test)
  test <- test[order(net_reproduction_ratios[test])]

  drop <- logical(length(net_reproduction_ratios))

  for (i in test) {
    plant_log_inviable(paste("Testing species", i),
                       stage="testing", species=i)
    x <- net_reproduction_ratios
    x[i] <- eps
    res <- runner(x)
    if (res[[i]] < eps) {
      plant_log_inviable(paste("Removing species", i),
                         stage="removing", species=i)
      drop[[i]] <- TRUE
      res[[i]] <- 0.0
      net_reproduction_ratios <- res
    }
  }

  ## It's possible that things slip through and get driven extinct by
  ## the time that they reach here.
  drop <- drop | net_reproduction_ratios < eps

  attr(net_reproduction_ratios, "drop") <- drop
  net_reproduction_ratios
}
