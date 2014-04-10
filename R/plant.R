# Things to do with plants

##' Grow a plant up to a certain size from a seed.
##'
##' This is a little utility and is likely to fail.
##'
##' To stop at a height of 10, the \code{distance.from.target}
##' function should be \code{function(plant) plant$height - 10}, so
##' that the return value is negative when the plant is smaller than
##' the target size.
##' @title Grow Plant To Given Size
##' @param plant A \code{Plant} object.
##' @param env An \code{Environment} object.
##' @param distance.from.target A function that takes \code{plant} as
##' its only argument and computes the distance from a target.  The
##' function must return negative values when the plant is smaller
##' than the target.  See Details.
##' @param max.time Time to run the ODE out for -- only exists to
##' prevent an infinite loop (say, on an unreachable size).
##' @return A plant grown to the appropriate size.  This is a
##' \emph{copy} of the input plant.
##' @author Rich FitzJohn
##' @export
grow.plant.to.size <- function(plant, env, distance.from.target,
                               max.time=1000) {
  # derivatives function (finds plant and env through local scope)
  derivs <- function(t, y, pars) {
    plant$set_ode_values(t, y)
    plant$compute_vars_phys(env)
    plant$ode_rates
  }

  # Work off a copy of the plant so that the input plant is not
  # modified:
  plant <- plant$copy()

  # Set up the ode solver.
  solver <- new(OdeR, derivs, new.env(), NULL)
  solver$set_state(plant$ode_values, 0)

  # Bracket the solution:
  # Starting size must be smaller than the target size:
  if (distance.from.target(plant) > 0) {
    stop("Plant already bigger than target")
  }

  while (solver$time < max.time && distance.from.target(plant) < 0) {
    t0 <- solver$time
    y0 <- solver$state
    solver$step()
  }

  if (distance.from.target(plant) < 0) {
    stop("Could not find plant size big enough by max_time")
  }

  t1 <- solver$time

  # Isolate the root within this interval [t0, t1]:

  # Helper for uniroot: stepping locally and non-adaptively from t0 to
  # t, setting up the plant with ode variables y0.  These are going to
  # be the values at the end of bracketing the solution
  target <- function(t, solver, t0, y0) {
    solver$reset()
    solver$set_state(y0, t0)
    solver$step_to(t)
    distance.from.target(plant)
  }
  ans <- uniroot(target, c(t0, t1), solver, t0, y0)

  # And return the plant
  plant
}
