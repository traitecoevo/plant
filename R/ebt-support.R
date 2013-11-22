## EBT support functions.

##' Generate a suitable set of default cohort introduction times,
##' biased so that introductions are more closely packed at the
##' beginning of time, become increasingly spread out.
##'
##' The reason for the stepped distribution is to keep step sizes as
##' series of doublings.  Doing this limits the range of possible
##' introduction times from an infinite set of possible values to a
##' very limited subset of values (based on combinations of 1, 0.5,
##' 0.25, 0.125 etc).  The reason for doing this is to minimise the
##' number of unique introduction times across all species. The ODE
##' stepper needs to stop at each point where a cohort is introduced.
##' If each species was selecting a bunch of points that was
##' essentially unique (compared to those selected for all other
##' species), the number of unique cohort introductions times could
##' get very large, requiring more ODE steps.
##'
##' @title Generate Default Cohort Introduction Times
##' @param max.time Time to generate introduction times up to (the
##' last introduction time will be at least \code{max.time}).
##' @param multiplier The rate of increase of step size with time.
##' The greater the number the faster step size will increase.
##' @param min.step.size The smallest gap between introduction times
##' (must be greater than zero, and will be the first introduction
##' time).
##' @param max.step.size The largest gap between introduction times
##' (may be infinite).
##' @return Vector of introduction times.
##' @export
##' @author Rich FitzJohn, adapted from original C++ code by Daniel
##' S. Falster.
cohort.introduction.times <- function(max.time, multiplier=0.2,
                                      min.step.size=1e-5,
                                      max.step.size=2.0) {
  if (min.step.size <= 0)
    stop("The minimum step size must be greater than zero")
  times <- numeric(0)
  dt <- time <- 0
  while (time <= max.time) {
    times <- c(times, time)
    dt <- 2^floor(log2(time * multiplier))
    time <- time + max(min(dt, max.step.size), min.step.size)
  }
  times
}

##' Parameters for running the simulations more quickly (but less
##' accurately) than the defaults.  Used in a number of places.
##'
##' @title Fast Control Defaults
##' @return A `list` of values to be passed into
##' \code{Control$set_parameters}
##' @author Rich FitzJohn
##' @export
fast.control <- function() {
  ctrl <- list()
  ctrl$environment_light_rescale_usually <- TRUE
  ctrl$environment_light_tol <- 1e-4
  ctrl$plant_assimilation_rule <- 21
  ctrl$plant_assimilation_over_distribution <- FALSE
  ctrl$plant_assimilation_tol <- 1e-4
  ctrl$ode_tol_rel <- 1e-4
  ctrl$ode_tol_abs <- 1e-4
  ctrl$ode_step_size_max <- 5
  ctrl$cohort_gradient_direction <- -1
  ctrl$cohort_gradient_richardson <- FALSE
  ctrl
}

##' Generate a default schedule of cohort introduction times.  At
##' present this is just linear from 0..max.t.
##'
##' Note that this is basically useless for practical purposes, but
##' will merge with \code{cohort.introduction.times} at some point.
##'
##' @title Simple Linear Cohort Schedule
##' @param nt Number of cohorts to introduce
##' @param max.t Time that the simulation stops at
##' @return A \code{CohortSchedule} object
##' @author Rich FitzJohn
##' @export
default.schedule <- function(nt, max.t) {
  times <- seq(0, max.t, length=nt + 1)[-(nt + 1)]
  sched <- new(CohortSchedule, 1)
  sched$set_times(times, 1)
  sched$max_time <- max.t
  sched
}
