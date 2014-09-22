
## Need to replicate the little runner function:
schedule_or_null <- tree:::schedule_or_null
make_equilibrium_runner <- function(p, schedule=NULL) {
  p <- p$copy()
  schedule_default <- schedule_or_null(schedule, p)

  control <- p$control$parameters
  large_seed_rain_change <- control$equilibrium_large_seed_rain_change
  last_schedule <- schedule_default$copy()

  function(seed_rain_in,
           force_new_schedule=FALSE, force_old_schedule=FALSE) {
    last_seed_rain_in <- p$seed_rain
    message(sprintf("seed_rain: %s, last: %s",
                    paste(prettyNum(seed_rain_in, digits=3), collapse=", "),
                    paste(prettyNum(last_seed_rain_in, digits=3), collapse=", ")))
    p$seed_rain <- seed_rain_in
    if (force_old_schedule) {
      schedule <- last_schedule
    } else {
      if (force_new_schedule ||
          any(abs(seed_rain_in - last_seed_rain_in) > large_seed_rain_change)) {
        schedule <- schedule_default$copy()
      } else {
        schedule <- last_schedule$copy()
      }
      schedule <- build_schedule(p, schedule)
    }
    last_schedule <<- schedule
    
    seed_rain_out <- attr(schedule, "seed_rain", exact=TRUE)
    attr(seed_rain_out, "schedule") <- schedule$copy()
    seed_rain_out
  }
}

drop_schedule <- function(x) {
  attr(x, "schedule") <- NULL
  x
}

equilibrium <- function(p, schedule=NULL) {
  p <- p$copy()
  target <- make_equilibrium_runner(p, schedule)

  control <- p$control$parameters
  eps <- control$equilibrium_eps
  progress <- control$equilibrium_progress
  history <- list()
  seed_rain <- p$seed_rain

  for (i in seq_len(control$equilibrium_nsteps)) {
    ans <- target(seed_rain)
    seed_rain <- ans[,"out"]
    change <- ans[,"out"] - ans[,"in"]
    res <- list(seed_rain=drop_schedule(ans),
                schedule=attr(ans, "schedule"))
    if (progress) {
      history <- c(history, list(res))
    }
    if (control$equilibrium_verbose) {
      message(sprintf("*** %d: {%s} -> {%s} (delta = {%s})", i,
                      paste(prettyNum(ans[,"in"]), collapse=","),
                      paste(prettyNum(ans[,"out"]), collapse=","),
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
  res
}

## This is a wrapper for nleqslv() to deal with negative seed rains.
## I'm going to treat everything with zero seed rain as not being
## absent in the population and therefore having zero fitness.  This
## induces a discontinuitiy though, so we might pay for this later.
##
## Another approach would be to encourage the function away from the
## trivial equilibrium by multiplying (seed_out - seed_in) by
## something that approaches 1/x for x -> Inf and 0 for x -> 0, and is
## monotonic.  With those conditions, I think we'll stay away from the
## trivial root iff it is unstable.
make_target <- function(f) {
  force(f)
  function(x, ...) {
    neg <- x < 0
    if (all(neg)) {
      rep(0.0, length(x))
    } else {
      x[neg] <- 0.0
      res <- f(x, ...)
      ret <- rep(0, length(x))
      pos <- !neg
      ret[pos] <- log(res[pos,"out"] / res[pos,"in"])
      ret
    }
  }
}

make_target_ode <- function(f) {
  force(f)
  g <- make_target(f)
  function(x) {
    g(x) * pmax(x, 0.0)
  }
}

make_pars <- function(pars, time_disturbance) {
  p <- ebt_base_parameters()
  p$set_control_parameters(list(equilibrium_nsteps=20))
  p$strategy_default <- new(Strategy, pars)
  p$disturbance <- new(Disturbance, time_disturbance)
  p
}

## This *probably* could be a useful member function for the
## Parameters object.  It gets used in the assembler_community object
## (see `to_parameters`)
add_strategy <- function(p, pars) {
  s <- p$strategy_default$copy()
  s$set_parameters(pars)
  p$add_strategy(s)
}
