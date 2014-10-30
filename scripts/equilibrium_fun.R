make_equilibrium_runner <- tree:::make_equilibrium_runner
equilibrium <- tree:::equilibrium_seed_rain_iteration

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
##
## Also need to get the true zero value done better here; should
## use same flag that indicates that 
make_target <- function(f) {
  force(f)
  function(x, ...) {
    pos <- x > 0
    if (!any(pos)) {
      rep(0.0, length(x))
    } else {
      x[!pos] <- 0.0
      res <- f(x, ...)
      ret <- rep(0, length(x))
      ret[pos] <- log(res[pos,"out"] / res[pos,"in"])
      ret
    }
  }
}

make_target_ode <- function(f) {
  g <- make_target(f)
  function(x) {
    g(x) * pmax(x, 0.0)
  }
}

make_target2 <- function(f, allow_zero) {
  g <- make_target(f)
  force(allow_zero)
  function(x) {
    if (length(x) != length(allow_zero)) {
      stop("unexpected length")
    }
    x[x == 0 & allow_zero] <- 1e-10
    res <- g(x)
    res[allow_zero] <- res[allow_zero] * pmax(x[allow_zero], 0.0)
    res
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
