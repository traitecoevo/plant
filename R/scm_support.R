##' Sets reasonable defaults for fast numerical calculations
##' @title Fast Control Defaults
##' @return A Control object with parameters set.
##' @author Rich FitzJohn
##' @export
##' @param base An optional \code{Control} object.  If omitted, the
##' defaults are used.
fast_control <- function(base=Control()) {
  base$function_integration_rule <- 21

  base$ode_tol_rel <- 1e-4
  base$ode_tol_abs <- 1e-4
  base$ode_step_size_max <- 5

  base$node_gradient_direction <- -1
  base$node_gradient_richardson <- FALSE

  base
}

##' Hopefully sensible set of parameters for use with the SCM.  Turns
##' accuracy down a bunch, makes it noisy, sets up the
##' hyperparameterisation that we most often use.
##' @title Sensible, fast (ish) SCM control settings
##' @author Rich FitzJohn
##' @export
scm_base_control <- function() {
  ctrl <- fast_control()
  ctrl$schedule_eps <- 0.005
  return(ctrl)
}


##' Basic default settings for a given strategy, environment only really
##' used for templating initially and will be overloaded later by passing
##' an environment to the SCM API (suggesting perhaps the template could be
##' removed).
##' @title Basic default parameters for a given strategy
##' @author Rich FitzJohn
##' @param type Any strategy name as a string, e.g.: \code{"FF16"}.
##' @param env And environment object
##' @export
scm_base_parameters <- function(type = NA, env = environment_type(type)) {
  
   Parameters(type, env)(patch_area=1.0)
}


##' Run the SCM, returning the SCM object for interrogation
##'
##' This is the simplest way of using the SCM, probably.
##' @title Run SCM
##' @param p Parameters object
##' @param env Environment object (defaults to FF16_Environment)
##' @param ctrl Control object
##' @param use_ode_times Should ODE times be used?
##' @return A \code{SCM} object.
##' @author Rich FitzJohn
##' @export
run_scm <- function(p, env = NULL, 
                    ctrl = scm_base_control(), use_ode_times=FALSE, collect = FALSE) {

  types <- extract_RcppR6_template_types(p, "Parameters")
  
  if(is.null(env))
    env <- Environment(types[[1]])

  scm <- do.call('SCM', types)(p, env, ctrl)
  if (use_ode_times) {
    scm$use_ode_times <- TRUE
  }
  if(collect) {
    scm$collect <- TRUE
  }
  scm$run()
  scm
}


##' Run the SCM model, given a Parameters and NodeSchedule
##'
##' This is mostly a simple wrapper around some of the SCM functions.
##' Not sure if this is how we will generally want to do this.
##' Consider this function liable to change.
##'
##' @title Run the SCM, Collecting Output
##' @param p A \code{Parameters} object
##' @param env Environment object (defaults to FF16_Environment)
##' @param ctrl Control object
##' competition_effect)
##' @author Rich FitzJohn
##' @export
run_scm_collect <- function(p, env = NULL, 
                            ctrl = scm_base_control()) {
  
  scm <- run_scm(p, env, ctrl, collect = TRUE)

  results <- lapply(scm$history, "[[", "state") |> tidy_patch()

  results[["offspring_production"]] <- scm$offspring_production
  results[["net_reproduction_ratios"]] <- scm$net_reproduction_ratios
  
  results[["p"]] <- p

  results
}

run_scm_error <- function(p, env = Environment(parameters = p),
                          ctrl = scm_base_control()) {
  types <- extract_RcppR6_template_types(p, "Parameters")
  scm <- do.call('SCM', types)(p, env, ctrl)
  n_spp <- length(p$strategies)

  lai_error <- rep(list(NULL), n_spp)
  while (!scm$complete) {
    added <- scm$run_next()
    for (idx in added) {
      lai_error[[idx]] <-
        c(lai_error[[idx]], list(scm$competition_effect_error(idx)))
    }
  }

  rbind_list <- function(x) do.call("rbind", as.list(x))

  lai_error <- lapply(lai_error, function(x) rbind_list(pad_matrix(x)))
  average_fecundity_error <- scm$average_fecundity_error
  f <- function(m) {
    suppressWarnings(apply(m, 2, max, na.rm=TRUE))
  }
  total <- lapply(seq_len(n_spp), function(idx)
                  f(rbind(lai_error[[idx]], average_fecundity_error[[idx]])))
  list(offspring_production=scm$offspring_production,
       err=list(lai=lai_error, offspring_production=average_fecundity_error, total=total),
       ode_times=scm$ode_times)
}

