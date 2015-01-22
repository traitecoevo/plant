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
  base$schedule_verbose=TRUE
  base$equilibrium_verbose=TRUE
  base$equilibrium_progress=TRUE
  base
}
##' @export
##' @rdname equilibrium_verbose
equilibrium_quiet <- function(base=Control()) {
  base$schedule_verbose=FALSE
  base$equilibrium_verbose=FALSE
  base$equilibrium_progress=FALSE
  base
}

##' Run the EBT, returning the EBT object for interrogation
##'
##' This is the simplest way of using the EBT, probably.
##' @title Run EBT
##' @param p Parameters object
##' @return A \code{EBT} object.
##' @author Rich FitzJohn
##' @export
run_ebt <- function(p, sched) {
  ebt <- new(EBT, p)
  ebt$run()
  ebt$update_ode_times()
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
