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
##' will merge with \code{\link{cohort.introduction.times}} at some
##' point.
##'
##' @title Simple Linear Cohort Schedule
##' @param nt Number of cohorts to introduce
##' @param max.t Time that the simulation stops at
##' @return A \code{CohortSchedule} object
##' @author Rich FitzJohn
##' @export
default.schedule <- function(nt, max.t) {
  schedule.from.times(seq(0, max.t, length.out=nt + 1))
}

##' Generate a CohortSchedule object from a series of times.  The last
##' time is ommited but will be set as the \code{$max_time}.
##'
##' Only works for a single species at the moment.
##'
##' @title Create CohortSchedule From Vector of Times
##' @param times A numeric vector of increasing times
##' @return A \code{\link{CohortSchedule}} object.
##' @author Rich FitzJohn
##' @export
schedule.from.times <- function(times) {
  if (any(diff(times) <= 0.0))
    stop("Times must be strictly increasing")
  sched <- new(CohortSchedule, 1)
  sched$set_times(times[-length(times)], 1)
  sched$max_time <- last(times)
  sched
}

##' Run the EBT model, given a Parameters and CohortSchedule
##'
##' This is mostly a simple wrapper around some of the EBT functions.
##' Not sure if this is how we will generally want to do this.
##' Consider this function liable to change.
##'
##' @title Run the EBT, Collecting Output
##' @param p A Parameters object
##' @param sched A CohortSchedule Object
##' @export
##' @author Rich FitzJohn
run.ebt.collect <- function(p, sched) {
  get.state <- function(ebt)
    list(time=ebt$time,
         species=ebt$state$patch$species,
         light.env=get.light.env(ebt))

  ebt <- new(EBT, p)
  ebt$cohort_schedule <- sched
  res <- list(get.state(ebt))

  while (ebt$cohort_schedule$remaining > 0) {
    ebt$run_next()
    st <- get.state(ebt)
    if (st$time > last(res)$time) {
      res <- c(res, list(st))
    } else {
      res[[length(res)]] <- st
    }
  }
  pad.list.to.array <- tree:::pad.list.to.array

  time <- sapply(res, "[[", "time")
  light.env <- lapply(res, "[[", "light.env")
  ## The aperm() here means that dimensions are
  ## [variable,time,cohort], so that taking species[[1]]["height",,]
  ## gives a matrix that has time down rows and cohorts across columns
  ## (so is therefore plottable with matplot)
  species <- lapply(res, "[[", "species")
  species <- lapply(seq_along(species[[1]]), function(i)
                    aperm(pad.list.to.array(lapply(species, "[[", i)),
                          c(1, 3, 2)))
  ## Drop the boundary condition; we do this mostly because it cannot
  ## be compared against the reference output, which does not contain
  ## this.  This does have the nice property of giving a non-square
  ## matrix, so the difference between time and cohort becomes a
  ## little more obvious.
  species <- lapply(species, function(m) m[,,-dim(m)[[3]]])

  list(time=time, species=species, light.env=light.env,
       fitness=ebt$fitnesses)
}

get.light.env <- function(ebt) {
  light.env <- ebt$patch$environment$light_environment$xy
  colnames(light.env) <- c("height", "canopy.openness")
  light.env
}

##' Run the EBT, returning the EBT object for interrogation
##'
##' This is the simplest way of using the EBT, probably.
##' @title Run EBT
##' @param p Parameters object
##' @param sched CohortSchedule object
##' @return A \code{EBT} object.
##' @author Rich FitzJohn
##' @export
run.ebt <- function(p, sched) {
  ebt <- new(EBT, p)
  ebt$cohort_schedule <- sched
  ebt$run()
  ebt
}
