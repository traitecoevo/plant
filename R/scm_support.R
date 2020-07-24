##' Sets reasonable defaults for fast numerical calculations
##' @title Fast Control Defaults
##' @return A Control object with parameters set.
##' @author Rich FitzJohn
##' @export
##' @param base An optional \code{Control} object.  If omitted, the
##' defaults are used.
fast_control <- function(base=Control()) {
  base$environment_rescale_usually <- TRUE
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

##' Run the SCM, returning the SCM object for interrogation
##'
##' This is the simplest way of using the SCM, probably.
##' @title Run SCM
##' @param p Parameters object
##' @param use_ode_times Should ODE times be used?
##' @return A \code{SCM} object.
##' @author Rich FitzJohn
##' @export
run_scm <- function(p, use_ode_times=FALSE) {
  types <- extract_RcppR6_template_types(p, "Parameters")
  scm <- do.call('SCM', types)(p)
  if (use_ode_times) {
    scm$use_ode_times <- TRUE
  }
  scm$run()
  scm
}

##' Hopefully sensible set of parameters for use with the SCM.  Turns
##' accuracy down a bunch, makes it noisy, sets up the
##' hyperparameterisation that we most often use.
##' @title Sensible, fast (ish) SCM parameters
##' @author Rich FitzJohn
##' @param type Name of model (defaults to FF16 but any strategy name is valid).
##' @export
scm_base_parameters <- function(type="FF16", env=sprintf("%s_Env", type)) {
  ctrl <- equilibrium_verbose(fast_control())
  ctrl$schedule_eps <- 0.005
  ctrl$equilibrium_eps <- 1e-3
  Parameters(type, env)(patch_area=1.0, control=ctrl, hyperpar=hyperpar(type))
}

##' Run the SCM model, given a Parameters and CohortSchedule
##'
##' This is mostly a simple wrapper around some of the SCM functions.
##' Not sure if this is how we will generally want to do this.
##' Consider this function liable to change.
##'
##' @title Run the SCM, Collecting Output
##' @param p A \code{Parameters} object
##' @param include_competition_effect Include total leaf area (will change; see
##' issue #138)
##' @author Rich FitzJohn
##' @export
run_scm_collect <- function(p, include_competition_effect=FALSE) {
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
  types <- extract_RcppR6_template_types(p, "Parameters")

  make_environment("FF16", p)
  scm <- do.call('SCM', types)(p)
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

  patch_density <- scm$patch$environment$disturbance_regime$density(time)

  ret <- list(time=time, species=species,
              env=env,
              seed_rain=scm$seed_rains,
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
##' @export
make_patch <- function(state, p) {
  types <- extract_RcppR6_template_types(p, "Parameters")
  n <- viapply(state$species, ncol)
  patch <- do.call('Patch', types)(p)
  patch$set_state(state$time, unlist(state$species), n, state$env)
  patch
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

run_scm_error <- function(p) {
  types <- extract_RcppR6_template_types(p, "Parameters")
  scm <- do.call('SCM', types)(p)
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
  seed_rain_error <- scm$seed_rain_error
  f <- function(m) {
    suppressWarnings(apply(m, 2, max, na.rm=TRUE))
  }
  total <- lapply(seq_len(n_spp), function(idx)
                  f(rbind(lai_error[[idx]], seed_rain_error[[idx]])))

  list(seed_rain=scm$seed_rains,
       err=list(lai=lai_error, seed_rain=seed_rain_error, total=total),
       ode_times=scm$ode_times)
}

##' Set a suitable hyperparameter function for chosen physiological model
##' @title Hyperparameters for FF16 physiological model
##' @param type Any strategy name as a string, e.g.: \code{"FF16"}.
##' @rdname Hyperparameter_functions
##' @export
# if you update this function (even syntactic changes) update the function update_smc_support in the scaffolder
make_hyperpar <- function(type) {
  switch(type,
         FF16r=make_FF16_hyperpar,
         K93=make_FF16_hyperpar,
         FF16=make_FF16_hyperpar,
         stop("Unknown type ", type))
}

##' @rdname Hyperparameter_functions
##' @export
# if you update this function (even syntactic changes) update the function update_smc_support in the scaffolder
hyperpar <- function(type) {
  switch(type,
         FF16r=FF16_hyperpar,
         K93=K93_hyperpar,
         FF16=FF16_hyperpar,
         stop("Unknown type ", type))
}

##' Helper function for creating parameter objects suitable for an
##' assembly.
##' @title Helper function for creating parameter objects
##' @param ... Named set of parameters
##' @param pars A list of parameters
##' @param base_parameters_fn Function for creating base parameter set (default scm_base_parameters)
##' @param make_hyperpar_fn Function for creating hyperparameterisation (default make_FF16_hyperpar)
##' @export
assembly_parameters <- function(..., pars=NULL, base_parameters_fn = scm_base_parameters,
                                  make_hyperpar_fn = make_FF16_hyperpar) {

  p <- base_parameters_fn()

  ## These are nice to have:
  p$control$equilibrium_solver_name <- "hybrid"
  p$control$equilibrium_nsteps <- 60

  if (is.null(pars)) {
    pars <- list(...)
  } else if (length(list(...)) > 0L) {
    stop("Do not provide both ... and pars")
  }

  if (length(pars) > 0L) {
    assert_named_if_not_empty(pars)

    excl <- c("control", "strategy_default", "hyperpar")
    pos <- setdiff(c(names(formals(make_hyperpar_fn)),
                     names(p),
                     names(p$control),
                     names(p$strategy_default)),
                   excl)
    unk <- setdiff(names(pars), pos)
    if (length(unk) > 0L) {
      stop("Unknown parameters: ", paste(unk, collapse=", "))
    }

    nms_hyper <- intersect(names(pars), names(formals(make_hyperpar_fn)))
    p$hyperpar <- do.call("make_hyperpar_fn", pars[nms_hyper])
    p                  <- modify_list(p,                  pars)
    p$control          <- modify_list(p$control,          pars)
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
        seeds_survival_weighted=sp$seeds)
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
    pa <- lapply(a, obj$patch$environment$disturbance_regime$density)
  } else {
    internals <- patch_to_internals(scm_patch(length(obj$time), obj))
    n <- length(internals)
    a <- obj$p$cohort_schedule_times
    pa <- lapply(a, Disturbance(obj$p$disturbance_mean_interval)$density)
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

make_environment<- function(type, ...) {
  switch(type,
    FF16=FF16_make_environment(...),
    FF16r=FF16r_make_environment(...),
    stop("Unknown type ", type))
}
