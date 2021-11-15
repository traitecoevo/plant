##' Sets reasonable defaults for fast numerical calculations
##' @title Fast Control Defaults
##' @return A Control object with parameters set.
##' @author Rich FitzJohn
##' @export
##' @param base An optional \code{Control} object.  If omitted, the
##' defaults are used.
fast_control <- function(base=Control()) {
  base$assimilator_adaptive_integration <- FALSE
  base$assimilator_integration_rule <- 21
  base$assimilator_integration_tol <- 1e-4

  base$ode_tol_rel <- 1e-4
  base$ode_tol_abs <- 1e-4
  base$ode_step_size_max <- 5

  base$cohort_gradient_direction <- -1
  base$cohort_gradient_richardson <- FALSE

  base
}

##' Control parameters for \code{\link{equilibrium_birth_rate}} that
##' make progress noisier.  This is just a convenience function.
##'
##' @title Noisy Parameters for Equilibrium Finding
##' @export
##' @param base An optional \code{Control} object.  If omitted, the
##' defaults are used.
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

##' Hopefully sensible set of parameters for use with the SCM.  Turns
##' accuracy down a bunch, makes it noisy, sets up the
##' hyperparameterisation that we most often use.
##' @title Sensible, fast (ish) SCM control settings
##' @author Rich FitzJohn
##' @export
scm_base_control <- function() {
  ctrl <- equilibrium_verbose(fast_control())
  ctrl$schedule_eps <- 0.005
  ctrl$equilibrium_eps <- 1e-3
  return(ctrl)
}


##' Basic default settings for a given strategy, environment only really
##' used for templating initially and will be overloaded later by passing
##' an environment to the SCM API (suggesting perhaps the template could be
##' removed).
##' @title Basic default parameters for a given strategy
##' @author Rich FitzJohn
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
##' @param state A optional State object matching the strategies in \code{p}
##' @param use_ode_times Should ODE times be used?
##' @return A \code{SCM} object.
##' @author Rich FitzJohn
##' @export
run_scm <- function(p, env = make_environment(parameters = p),
                       ctrl = scm_base_control(), 
                       state=NULL, 
                       use_ode_times=FALSE) {
  scm <- make_scm(p, env, ctrl, state)
  if (use_ode_times) {
    scm$use_ode_times <- TRUE
  }
  scm$run()
  scm
}


##' Run the SCM model, given a Parameters and CohortSchedule
##'
##' This is mostly a simple wrapper around some of the SCM functions.
##' Not sure if this is how we will generally want to do this.
##' Consider this function liable to change.
##'
##' @title Run the SCM, Collecting Output
##' @param p A \code{Parameters} object
##' @param env Environment object (defaults to FF16_Environment)
##' @param ctrl Control object
##' @param state An optional State object matching the strategies in \code{p}
##' @param include_competition_effect Include total leaf area (will change; see
##' issue #138)
##' @author Rich FitzJohn
##' @export
run_scm_collect <- function(p, env = make_environment(parameters = p), 
                            ctrl = scm_base_control(),
                            state = NULL,
                            include_competition_effect=FALSE) {
  collect_default <- function(scm) {
    scm$state
  }
  collect_competition_effect <- function(scm) {
    ret <- scm$state
    competition_effect <- numeric(length(ret$species))
    for (i in seq_along(ret$species)) {
      ## ret$species[[i]] <- rbind(
      ##   ret$species[[i]]
      ##   competition_effect=c(scm$patch$species[[i]]$competition_effects, 0.0))
      competition_effect[i] <- scm$patch$species[[i]]$compute_competition(0.0)
    }
    ret$competition_effect <- competition_effect
    ret
  }

  collect <- if (include_competition_effect) collect_competition_effect else collect_default
  
  scm <- make_scm(p, env, ctrl, state)
  res <- list(collect(scm))

  while (!scm$complete) {
    scm$run_next()
    res <- c(res, list(collect(scm)))
  }

  time <- sapply(res, "[[", "time")
  env <- lapply(res, "[[", "env")
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

  patch_density <- scm$patch$density(time)

  ret <- list(time=time, species=species,
              env=env,
              net_reproduction_ratios=scm$net_reproduction_ratios,
              patch_density=patch_density,
              p=p)

  if (include_competition_effect) {
    ret$competition_effect <- do.call("rbind", lapply(res, "[[", "competition_effect"))
  }

  ret
}

##' Functions for reconstructing a Patch from an SCM
##' @title Reconstruct a patch
##' @param state State object created by \code{scm_state}
##' @param p Parameters object
##' @param env Environment object (defaults to FF16_Environment)
##' @param ctrl Control object
##' @export
make_patch <- function(state, p, env = make_environment(parameters = p),
                       ctrl = scm_base_control()) {
  types <- extract_RcppR6_template_types(p, "Parameters")
  n <- viapply(state$species, ncol)
  patch <- do.call('Patch', types)(p, env, ctrl)
  patch$set_state(state$time, unlist(state$species), n, state$env)
  patch
}

##' Functions for reconstructing a Patch from an SCM
##' @title Reconstruct a patch
##' @param p Parameters object
##' @param state An optional State object matching the strategies in \code{p}
##' @export
make_scm <- function(p, env, ctrl, state=NULL) {
  types <- extract_RcppR6_template_types(p, "Parameters")
  scm <- do.call('SCM', types)(p, env, ctrl)
  
  if(!is.null(state)) {
    n_str <- length(p$strategies)
    n_spp = length(state$species)
    
    if(n_spp != n_str)
      stop("State object has more species than strategies defined in Parameters")
  
    # need to append cohort times to enable integration of net fecundity
    if(state$time != 0)
      warning("Solver must start from 0, resetting initial state time")

    times <- scm$cohort_schedule$all_times

    initial_cohorts <- sapply(times, function(t) sum(t == 0), simplify = F)
    new_cohorts <- mapply(function(s, i) max(0, ncol(s) - i), 
                          state$species, initial_cohorts, SIMPLIFY = F)
    
    new_times <- mapply(function(i, t) c(rep(0, i), t), 
                        new_cohorts, times, SIMPLIFY = F)

    # this introduces one more ind. than necessary, but if we
    # overwrite the oldest cohorts first then we can just start at t1
    scm$set_cohort_schedule_times(new_times)
    scm$run_next()

    # next step starts from first non-zero time
    start_time <- sapply(times, function(t) min(t[t>0]))
        
    # add as many new cohorts as required to fit `state` object
    scm$set_state(min(start_time), unlist(state$species), 
                  n = mapply(`+`, new_cohorts, initial_cohorts))
  }  
  
  return(scm)  
}
  
##' @rdname make_patch
##' @param i Index to extract from \code{x}
##' @param x Result of running \code{\link{run_scm_collect}}
##' @export
scm_state <- function(i, x) {
  f_sp <- function(el) {
    el <- el[, i, ]
    el[, !is.na(el[1, ]), drop=FALSE]
  }
  list(time=x$time[[i]], species=lapply(x$species, f_sp),
       env=x$env[[i]])
}

##' @export
##' @rdname make_patch
scm_patch <- function(i, x) {
  make_patch(scm_state(i, x), x$p)
}

run_scm_error <- function(p, env = make_environment(parameters = p),
                          ctrl = scm_base_control(), state = NULL) {
  scm <- make_scm(p, env, ctrl, state)
  n_spp <- length(p$strategies)

  lai_error <- rep(list(NULL), n_spp)
  while (!scm$complete) {
    added <- scm$run_next()
    for (idx in added) {
      lai_error[[idx]] <-
        c(lai_error[[idx]], list(scm$competition_effect_error(idx)))
    }
  }

  lai_error <- lapply(lai_error, function(x) rbind_list(pad_matrix(x)))
  net_reproduction_ratio_errors <- scm$net_reproduction_ratio_errors
  f <- function(m) {
    suppressWarnings(apply(m, 2, max, na.rm=TRUE))
  }
  total <- lapply(seq_len(n_spp), function(idx)
                  f(rbind(lai_error[[idx]], net_reproduction_ratio_errors[[idx]])))

  # schedule is needed to update parameters, not sure why ode_times is carried through
  # saving min height here is lazy, but I'm not sure where else build_schedule can get it
  list(net_reproduction_ratios=scm$net_reproduction_ratios,
       err=list(lai=lai_error, net_reproduction_ratios=net_reproduction_ratio_errors, total=total),
       schedule=scm$cohort_schedule$all_times,
       min_heights = sapply(scm$state$species, function(s) min(s["height", ])),
       ode_times=scm$ode_times)
}

##' Helper function for creating parameter objects suitable for an
##' assembly.
##' @title Helper function for creating parameter objects
##' @param ... Named set of parameters
##' @param pars A list of parameters
##' @param base_parameters_fn Function for creating base parameter set (default scm_base_parameters)
##' @param base_control_fn Function for creating base Control object (default scm_base_control)
##' @param make_hyperpar_fn Function for creating hyperparameterisation (default make_FF16_hyperpar)
##' @export
assembly_parameters <- function(..., pars=NULL, type = NA,
                                base_parameters_fn = scm_base_parameters,
                                base_control_fn = scm_base_control,
                                make_hyperpar_fn = make_FF16_hyperpar) {

  p <- base_parameters_fn(type)
  ctrl <- base_control_fn()

  ## These are nice to have:
  ctrl$equilibrium_solver_name <- "hybrid"
  ctrl$equilibrium_nsteps <- 60

  if (is.null(pars)) {
    pars <- list(...)
  } else if (length(list(...)) > 0L) {
    stop("Do not provide both ... and pars")
  }

  if (length(pars) > 0L) {
    assert_named_if_not_empty(pars)

    excl <- c("strategy_default", "hyperpar")
    pos <- setdiff(c(names(formals(make_hyperpar_fn)),
                     names(p),
                     names(p$strategy_default)),
                   excl)
    unk <- setdiff(names(pars), pos)
    if (length(unk) > 0L) {
      stop("Unknown parameters: ", paste(unk, collapse=", "))
    }

    nms_hyper <- intersect(names(pars), names(formals(make_hyperpar_fn)))
    p                  <- modify_list(p,                  pars)
    p$strategy_default <- modify_list(p$strategy_default, pars)
  }
  p
}

scm_to_internals <- function(obj, use_environment=TRUE) {
  dat <- lapply(seq_along(obj$time), function(i)
    patch_to_internals(scm_patch(i, obj), use_environment))
  f <- function(i) {
    aperm(pad_list_to_array(lapply(dat, function(x) t(x[[i]]))), c(1, 3, 2))
  }
  lapply(seq_along(dat[[1]]), f)
}

patch_to_internals <- function(x, use_environment=TRUE) {
  env <- if (use_environment) x$environment else NULL
  lapply(x$species, species_to_internals, env)
}


species_to_internals <- function(sp, environment=NULL) {
  # Aggregate and extract plants
  sp_p <- lapply(sp$cohorts, function(x) x$plant )
  new_names <- c(sp_p[[1]]$ode_names, paste0(sp_p[[1]]$ode_names, '_dt'), sp_p[[1]]$aux_names)
  ints <- do.call("rbind", lapply(sp_p, function(x) c(x$internals$states, x$internals$rates, x$internals$auxs)))
  colnames(ints) <- new_names
  cbind(ints,
        log_density=sp$log_densities,
        offspring_produced_survival_weighted=sp$net_reproduction_ratio_by_cohort)
}

##' Create a function that allows integrating aggregate properties of
##' the SCM system.
##'
##' The workflow here is to run an SCM to create an SCM by running
##' \code{run_scm}, or a set of data from \code{run_scm_collect} and
##' then reconstitute all the intermediate bits of data so that an any
##' variable that \code{PlantPlus} tracks can be integrated out.
##' Because the pre-processing step is reasonably slow, this function
##' returns a function that takes a variable name and integrates it.
##'
##' @title Integrate SCM variables
##' @param obj An object from \code{run_scm} or \code{run_scm_collect}
##' @export
make_scm_integrate <- function(obj) {
  ## TODO: This needs to be made to work with the output of
  ## run_scm_collect, which means that it could possibly work with the
  ## output over time, which would be cool; that might help with some
  ## of the stuff from the "emergent" vignette.
  if (inherits(obj, "SCM")) {
    internals <- patch_to_internals(obj$patch)
    n <- length(internals)
    sched <- obj$cohort_schedule
    a <- lapply(seq_len(n), sched$times)
    pa <- lapply(a, obj$patch$density)
  } else {
    internals <- patch_to_internals(scm_patch(length(obj$time), obj))
    n <- length(internals)
    a <- obj$p$cohort_schedule_times
    if(obj$p$patch_type == 'meta-population')
      d = Weibull_Disturbance_Regime(obj$p$max_patch_lifetime)
    else
      d = No_Disturbance()
    pa = lapply(a, d$density)
  }

  if (n == 0L) {
    stop("This just isn't going to work out (plant internals are empty)")
  }

  pos <- colnames(internals[[1]])

  f1 <- function(i, name, error=FALSE) {
    x <- a[[i]]
    y <- pa[[i]] * internals[[i]][, name]
    total <- trapezium(x, y)
    if (error) {
      local_error_integration(x, y, total)
    } else {
      total
    }
  }
  function(name, error=FALSE) {
    if (!(name %in% pos)) {
      stop("Unknown variable: ", name)
    }
    vnapply(seq_len(n), f1, name, error)
  }
}
