## Support functions for growing a plant to a given size, etc.  Used
## in dfalster/TraitGrowthTrajectories

##' Grow a plant up to particular sizes.
##'
##' @title Grow plant to given size
##' @param plant A \code{Plant} object.
##' @param sizes A vector of sizes to grow the plant to
##' @param size_name The name of the size variable within
##' \code{Plant$vars_phys} (e.g., height).
##' @param env An \code{Environment} object.
##' @param time_max Time to run the ODE out for -- only exists to
##' prevent an infinite loop (say, on an unreachable size).
##' @return A list with elements \code{time} (the time that a given
##' size was reached), \code{state} (the \emph{ode state} at these
##' times, as a matrix) and \code{plant} a list of plants grown to the
##' appropriate size.  Note that if only a single size is given,
##' a list of length 1 is returned.
##' @export
grow_plant_to_size <- function(plant, sizes, size_name, env, time_max=Inf) {
  obj <- grow_plant_bracket(plant, sizes, size_name, env, time_max)
  res <- lapply(seq_along(sizes), function(i)
                grow_plant_bisect(obj$runner, sizes[[i]], size_name,
                                  obj$t0[[i]], obj$t1[[i]], obj$y0[i,]))
  state <- t(sapply(res, "[[", "state"))
  colnames(state) <- colnames(obj$state)
  list(time=sapply(res, "[[", "time"),
       state=state,
       plant=lapply(res, "[[", "plant"))
}

grow_plant_bracket <- function(plant, sizes, size_name, env,
                               time_max=Inf) {
  if (length(sizes) == 0L || is.unsorted(sizes)) {
    stop("sizes must be non-empty and sorted")
  }
  if (plant$internals[[size_name]] > sizes[[1]]) {
    stop("Plant already bigger than smallest target size")
  }

  runner <- OdeRunner("PlantRunner")(PlantRunner(plant, env))
  i <- 1L
  n <- length(sizes)
  j <- integer(n)
  state <- list(list(time=runner$time, state=runner$state))

  while (i <= n) {
    runner$step()
    state <- c(state, list(list(time=runner$time, state=runner$state)))
    while (i <= n && oderunner_plant_size(runner)[[size_name]] > sizes[[i]]) {
      j[[i]] <- length(state) - 1L
      i <- i + 1L
    }
    if (runner$time > time_max) {
      stop("Time exceeded time_max")
    }
  }

  t <- sapply(state, "[[", "time")
  m <- t(sapply(state, "[[", "state"))
  colnames(m) <- runner$object$plant$ode_names
  list(t0=t[j], t1=t[j + 1L],
       y0=m[j,,drop=FALSE], y1=m[j + 1L,,drop=FALSE],
       time=t, state=m, index=j, runner=runner)
}

grow_plant_bisect <- function(runner, size, size_name, t0, t1, y0) {
  f <- function(t1) {
    runner$set_state(y0, t0)
    runner$step_to(t1)
    oderunner_plant_size(runner)[[size_name]] - size
  }
  ## NOTE: if we get access to the previously computed distances here
  ## we can save a little time.  But I don't think that's going to
  ## work well.
  root <- uniroot(f, lower=t0, upper=t1)
  list(time=root$root, state=runner$state, plant=runner$object$plant)
}
