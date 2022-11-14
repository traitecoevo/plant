##' Run system to offspring arrival equilibrium
##'
##' @title Run system to offspring arrival equilibrium
##' @param p A \code{Parameters} object
##' @param ctrl Control object
##' @return A Parameters object, with offspring arrival and node schedule
##' elements set.
##' @export
##' @author Rich FitzJohn
equilibrium_birth_rate <- function(p, ctrl) {
  solver <- ctrl$equilibrium_solver_name
  plant_log_info(sprintf("Solving offspring arrival using %s", solver),
                 routine = "equilibrium", stage = "start", solver = solver)
  switch(solver,
         iteration = equilibrium_birth_rate_iteration(p, ctrl = ctrl),
         nleqslv = equilibrium_birth_rate_solve_robust(p, ctrl = ctrl, solver = solver),
         dfsane = equilibrium_birth_rate_solve_robust(p, ctrl = ctrl, solver = solver),
         hybrid = equilibrium_birth_rate_hybrid(p, ctrl = ctrl),
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
  
  birth_rate <- purrr::map_dbl(p$strategies, ~ purrr::pluck(., birth_rate_y))

  runner <- make_equilibrium_runner(p, ctrl = ctrl)
  
  for (i in seq_len(ctrl$equilibrium_nsteps)) {
    offspring_production <- runner(birth_rate)
    converged <- check(birth_rate, offspring_production, eps, verbose)
    birth_rate <- offspring_production
    if (converged) {
      break
    }
  }

  # TODO: revisit 'gross' behaviour in cleanup utility
  equilibrium_runner_cleanup(runner, converged)
}

equilibrium_birth_rate_solve <- function(p, ctrl = scm_base_control(),
                                         solver="nleqslv") {
  try_keep <- ctrl$equilibrium_solver_try_keep
  logN <- ctrl$equilibrium_solver_logN
  min_offspring_arriving <- 1e-10 # TODO: should also be in the controls?

  plant_log_eq(paste("Solving offspring arrival using", solver),
               stage="start", solver=solver)

  birth_rate <- sapply(p$strategies, function(s) s$birth_rate_y, simplify = TRUE)
  runner <- make_equilibrium_runner(p,ctrl =ctrl)

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
    
    f <- function(s, br){
      s$birth_rate_y <- br
      return(s)
    }
    
    #this may or may not work (original function below)     
    #p$birth_rate <- unname(fit_it$birth_rate)

    p$strategies <- mapply(f, p$strategies, fit_it$birth_rate, SIMPLIFY = FALSE)
    
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
    eq_solution_iteration <- equilibrium_birth_rate_iteration(p, ctrl = ctrl)
    
    converged_it <- isTRUE(attr(eq_solution_iteration, "converged"))
    msg <- sprintf("Iteration %d %s",
                   i, if (converged_it) "converged" else "did not converge")
    plant_log_eq(msg, step="iteration", converged = converged_it, iteration=i)

    eq_solution <- try(
      equilibrium_birth_rate_solve(
        eq_solution_iteration, 
        ctrl = ctrl, 
        solver = solver[[i]]
        )
      )

    converged_sol <- isTRUE(attr(eq_solution, "converged"))

    msg <- sprintf("Solve %d %s",
                    i, if (converged_sol) "converged" else "did not converge")
    plant_log_eq(msg, step="solve", converged=converged_sol, iteration=i)

    if (converged_sol) {

      # check species with zero eq. birth rate are truly unviable.
      extinct = purrr::map_lgl(eq_solution$strategies, function(s) s$birth_rate_y == 0.0)

      if (any(extinct)) {
        plant_log_eq("Checking species driven to extinction")
        
        ## Add extinct species back at extremely low density and make sure
        ## that this looks like a legit extinction.
        p_check <- eq_solution
        p_check$strategies[extinct]$birth_rate_y <- ctrl$equilibrium_extinct_birth_rate

        res <- run_scm(p_check)

        # `next` breaks the loop iterating over solutions and does not return `eq_solution`
        if (any(res$offspring_production[extinct] > ctrl$equilibrium_extinct_birth_rate)) {
          plant_log_eq("Solver drove viable species extinct: rejecting")
          next
        }
      }
      plant_log_eq("Accepting solution via solver")
      return(eq_solution)
    }
  }

  ## This one should be a warning?
  plant_log_eq("Repeated rounds failed to find optimum; returning solution from equilibrium_birth_rate_iteration")
  return(eq_solution_iteration)
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
