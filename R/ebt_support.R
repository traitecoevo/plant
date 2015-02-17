##' @title Fast Control Defaults
##' @return A Control object with parameters set.
##' @author Rich FitzJohn
##' @export
##' @param base An optional \code{Control} object.  If omitted, the
##' defaults are used.
fast_control <- function(base=Control()) {
  base$environment_light_rescale_usually <- TRUE
  base$environment_light_tol <- 1e-4

  base$plant_assimilation_adaptive <- FALSE
  base$plant_assimilation_rule <- 21
  base$plant_assimilation_over_distribution <- FALSE
  base$plant_assimilation_tol <- 1e-4

  base$ode_tol_rel <- 1e-4
  base$ode_tol_abs <- 1e-4
  base$ode_step_size_max <- 5

  base$cohort_gradient_direction <- -1
  base$cohort_gradient_richardson <- FALSE

  base
}

##' Control parameters for \code{\link{equilibrium_seed_rain}} that
##' make progress noisier.  This is just a convenience function.
##'
##' @title Noisy Parameters for Equilibrium Finding
##' @export
##' @param base An optional \code{Control} object.  If omitted, the
##' defaults are used.
##' @examples
##' p <- new(Parameters)
##' p$set_control_parameters(equilibrium_verbose())
equilibrium_verbose <- function(base=Control()) {
  base$schedule_verbose <- TRUE
  base$equilibrium_verbose <- TRUE
  base
}
##' @export
##' @rdname equilibrium_verbose
equilibrium_quiet <- function(base=Control()) {
  base$schedule_verbose <- FALSE
  base$equilibrium_verbose <- FALSE
  base
}

##' Run the EBT, returning the EBT object for interrogation
##'
##' This is the simplest way of using the EBT, probably.
##' @title Run EBT
##' @param p Parameters object
##' @param use_ode_times Should ODE times be used?
##' @return A \code{EBT} object.
##' @author Rich FitzJohn
##' @export
run_ebt <- function(p, use_ode_times=FALSE) {
  ebt <- EBT(p)
  if (use_ode_times) {
    ebt$use_ode_times <- TRUE
  }
  ebt$run()
  ebt
}

##' Hopefully sensible set of parameters for use with the EBT.  Turns
##' accuracy down a bunch, makes it noisy.
##' @title Sensible, fast (ish) EBT parameters
##' @author Rich FitzJohn
##' @export
ebt_base_parameters <- function() {
  ctrl <- equilibrium_verbose(fast_control())
  ctrl$schedule_eps=0.005
  Parameters(patch_area=1.0, control=ctrl)
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
##' @author Rich FitzJohn
run_ebt_collect <- function(p) {
  ebt <- EBT(p)
  res <- list(ebt$state)

  while (!ebt$complete) {
    ebt$run_next()
    res <- c(res, list(ebt$state))
  }

  time <- sapply(res, "[[", "time")
  light_env <- lapply(res, "[[", "light_env")
  species <- lapply(res, "[[", "species")
  ## The aperm() here means that dimensions are
  ## [variable,time,cohort], so that taking species[[1]]["height",,]
  ## gives a matrix that has time down rows and cohorts across columns
  ## (so is therefore plottable with matplot)
  species <- lapply(seq_along(species[[1]]), function(i)
                    aperm(pad_list_to_array(lapply(species, "[[", i)),
                          c(1, 3, 2)))
  ## Drop the boundary condition; we do this mostly because it cannot
  ## be compared against the reference output, which does not contain
  ## this.  This does have the nice property of giving a non-square
  ## matrix, so the difference between time and cohort becomes a
  ## little more obvious.
  species <- lapply(species, function(m) m[,,-dim(m)[[3]]])

  list(time=time, species=species, light_env=light_env,
       seed_rain=ebt$seed_rains)
}

run_ebt_error <- function(p) {
  ebt <- EBT(p)
  n_spp <- length(p$strategies)

  lai_error <- rep(list(NULL), n_spp)
  while (!ebt$complete) {
    added <- ebt$run_next()
    for (idx in added) {
      lai_error[[idx]] <-
        c(lai_error[[idx]], list(ebt$area_leaf_error(idx)))
    }
  }

  lai_error <- lapply(lai_error, function(x) rbind_list(pad_matrix(x)))
  seed_rain_error <- ebt$seed_rain_error
  f <- function(m) {
    suppressWarnings(apply(m, 2, max, na.rm=TRUE))
  }
  total <- lapply(seq_len(n_spp), function(idx)
                  f(rbind(lai_error[[idx]], seed_rain_error[[idx]])))

  list(seed_rain=ebt$seed_rains,
       err=list(lai=lai_error, seed_rain=seed_rain_error, total=total),
       ode_times=ebt$ode_times)
}
