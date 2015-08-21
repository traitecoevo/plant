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
  type <- extract_RcppR6_template_type(p, "Parameters")
  ebt <- EBT(type)(p)
  if (use_ode_times) {
    ebt$use_ode_times <- TRUE
  }
  ebt$run()
  ebt
}

##' Hopefully sensible set of parameters for use with the EBT.  Turns
##' accuracy down a bunch, makes it noisy, sets up the
##' hyperparameterisation that we most often use.
##' @title Sensible, fast (ish) EBT parameters
##' @author Rich FitzJohn
##' @param type Name of model (defaults to FFW16 but FFdev also valid)
##' @export
ebt_base_parameters <- function(type="FFW16") {
  ctrl <- equilibrium_verbose(fast_control())
  ctrl$schedule_eps <- 0.005
  ctrl$equilibrium_eps <- 1e-3
  Parameters(type)(patch_area=1.0, control=ctrl, hyperpar=hyperpar(type))
}

##' Run the EBT model, given a Parameters and CohortSchedule
##'
##' This is mostly a simple wrapper around some of the EBT functions.
##' Not sure if this is how we will generally want to do this.
##' Consider this function liable to change.
##'
##' @title Run the EBT, Collecting Output
##' @param p A \code{Parameters} object
##' @param include_area_leaf Include total leaf area (will change; see
##' issue #138)
##' @author Rich FitzJohn
##' @export
run_ebt_collect <- function(p, include_area_leaf=FALSE) {
  collect_default <- function(ebt) {
    ebt$state
  }
  collect_area_leaf <- function(ebt) {
    ret <- ebt$state
    area_leaf <- numeric(length(ret$species))
    for (i in seq_along(ret$species)) {
      ## ret$species[[i]] <- rbind(
      ##   ret$species[[i]]
      ##   area_leaf=c(ebt$patch$species[[i]]$area_leafs, 0.0))
      area_leaf[i] <- ebt$patch$species[[i]]$area_leaf_above(0.0)
    }
    ret$area_leaf <- area_leaf
    ret
  }
  collect <- if (include_area_leaf) collect_area_leaf else collect_default
  type <- extract_RcppR6_template_type(p, "Parameters")

  ebt <- EBT(type)(p)
  res <- list(collect(ebt))

  while (!ebt$complete) {
    ebt$run_next()
    res <- c(res, list(collect(ebt)))
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

  patch_density <- ebt$patch$environment$disturbance_regime$density(time)

  ret <- list(time=time, species=species,
              light_env=light_env,
              seed_rain=ebt$seed_rains,
              patch_density=patch_density,
              p=p)

  if (include_area_leaf) {
    ret$area_leaf <- do.call("rbind", lapply(res, "[[", "area_leaf"))
  }

  ret
}

##' Functions for reconstructing a Patch from an EBT
##' @title Reconstruct a patch
##' @param state State object created by \code{ebt_state}
##' @param p Parameters object
##' @export
make_patch <- function(state, p) {
  type <- extract_RcppR6_template_type(p, "Parameters")
  n <- viapply(state$species, ncol)
  patch <- Patch(type)(p)
  patch$set_state(state$time, unlist(state$species), n, state$light_env)
  patch
}

##' @rdname make_patch
##' @param i Index to extract from \code{x}
##' @param x Result of running \code{\link{run_ebt_collect}}
##' @export
ebt_state <- function(i, x) {
  f_sp <- function(el) {
    el <- el[, i, ]
    el[, !is.na(el[1, ]), drop=FALSE]
  }
  list(time=x$time[[i]], species=lapply(x$species, f_sp),
       light_env=x$light_env[[i]])
}

##' @export
##' @rdname make_patch
ebt_patch <- function(i, x) {
  make_patch(ebt_state(i, x), x$p)
}

run_ebt_error <- function(p) {
  type <- extract_RcppR6_template_type(p, "Parameters")
  ebt <- EBT(type)(p)
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

##' Hyperparameters for plant
##' @title Hyperparameters for plant
##' @param B_kl_slope Slope of lma / leaf turnover log-log relationship
##' @param lma_0 LMA value...
##' @param B_kl_base ...
##' @param rho_0 ...
##' @param B_dI_base ...
##' @param B_dI_slope ...
##' @param B_ks_base ...
##' @param B_ks_slope ...
##' @param narea_0 ...
##' @param c_ext ...
##' @export
##' @rdname FFW16_hyperpar
make_FFW16_hyperpar <- function(B_kl_slope=1.71,
                               lma_0=0.1978791,
                               B_kl_base=0.4565855,
                               rho_0=608.0,
                               B_dI_base=0.01,
                               B_dI_slope=0.0,
                               B_ks_base=0.2,
                               B_ks_slope=0.0,
                               narea_0=1.87e-3,
                               c_ext=0.5,
                               latitude=0) {
  force(B_kl_slope)
  force(lma_0)
  force(B_kl_base)
  force(rho_0)
  force(B_dI_base)
  force(B_dI_slope)
  force(B_ks_slope)
  force(narea_0)
  force(c_ext)
  force(latitude)

  function(m, s, filter=TRUE) {
    with_default <- function(name, default_value=s[[name]]) {
      rep_len(if (name %in% colnames(m)) m[, name] else default_value,
              nrow(m))
    }
    lma       <- with_default("lma")
    rho       <- with_default("rho")
    omega     <- with_default("omega")
    narea     <- with_default("narea", narea_0)

    ## lma / leaf turnover relationship:
    k_l   <- B_kl_base * (lma / lma_0) ^ (-B_kl_slope)

    ## rho / mortality relationship:
    d_I  <- B_dI_base * (rho / rho_0) ^ (-B_dI_slope)

    ## rho / wood turnover relationship:
    k_s  <- B_ks_base *  (rho / rho_0) ^ (-B_ks_slope)

    ## rho / sapwood respiration relationship:

    ## Respiration rates are per unit mass, so this next line has the
    ## effect of holding constant the respiration rate per unit volume.
    ## So respiration rates per unit mass vary with rho, respiration
    ## rates per unit volume don't.
    r_s <- 4012.0 / rho
    r_b <- 2.0 * r_s # bark respiration follows from sapwood

    ## omega / accessory cost relationship
    a_f3 <- 3.0 * omega

    ## narea / photosynthesis / respiration
    ## Photosynthesis per mass leaf N [mol CO2 / kgN / yr]

    assimilian_rectangular_hyperbolae <- function(I, Amax, theta, quantum) {
      x <- quantum * I + Amax
      (x - (x^2 - 4 * theta * quantum * I * Amax)^0.5)/(2 * theta)
    }

    approximate_annual_assimilation <- function(pars, E = seq(0, 1, by = 0.02)) {

      # Only integrate over half year, as solar path is symmetrical
      D <- seq(0, 365/2, length.out = 10000)
      I <- PAR_given_solar_angle(solar_angle(D, latitude = abs(pars$latitude)))

      AA <- NA * E

      for (i in seq_len(length(E))) {
        AA[i] <- 2*trapezium(D, assimilian_rectangular_hyperbolae(
                                pars$c_ext * I * E[i],
                                pars$Amax, pars$theta, pars$QY)
                          )
      }

      data <- data.frame(E = E, AA = AA)
      fit <- nls(AA ~ p1 * E/(p2 + E), data, start = list(p1 = 100, p2 = 0.2))
      coef(fit)
    }

    # This needed in case narea has length zero, in which case trapezium fails
    a_p1  <- a_p2  <- 0 * narea
    if(length(narea) > 0 ) {
      NUE <- 5120.738 * 24 * 3600/1e+06
      # (mol CO2/ kg N /day )
      pars <- list(latitude = latitude, Amax = narea * NUE, theta = 0.5, QY = 0.04, c_ext = c_ext)
      y <- approximate_annual_assimilation(pars)
      a_p1  <- y[["p1"]]
      a_p2  <- y[["p2"]]
    }

    ## Respiration per mass leaf N [mol CO2 / kgN / yr]
    ## = (6.66e-4 * (365*24*60*60))
    ## Obtained from global average of ratio of dark respiration rate to
    ## leaf nitrogen content using the GLOPNET dataset
    c_RN  <-  21000
    ## Respiration rates are per unit mass, so convert to mass-based
    ## rate by dividing with lma
    ## So respiration rates per unit mass vary with lma, while
    ## respiration rates per unit area don't.
    r_l  <- c_RN * narea / lma

    extra <- cbind(k_l,                   # lma
                   d_I, k_s, r_s, r_b, # rho
                   a_f3,                 # omega
                   a_p1, a_p2,            # narea   
                   r_l)                  # lma, narea

    overlap <- intersect(colnames(m), colnames(extra))
    if (length(overlap) > 0L) {
      stop("Attempt to overwrite generated parameters: ",
           paste(overlap, collapse=", "))
    }

    ## Filter extra so that any column where all numbers are with eps
    ## of the default strategy are not replaced:
    if (filter) {
      if (nrow(extra) == 0L) {
        extra <- NULL
      } else {
        pos <- diff(apply(extra, 2, range)) == 0
        if (any(pos)) {
          eps <- sqrt(.Machine$double.eps)
          x1 <- extra[1, pos]
          x2 <- unlist(s[names(x1)])
          drop <- abs(x1 - x2) < eps & abs(1 - x1/x2) < eps
          if (any(drop)) {
            keep <- setdiff(colnames(extra), names(drop)[drop])
            extra <- extra[, keep, drop=FALSE]
          }
        }
      }
    }

    cbind(m, extra)
  }
}
##' @rdname FFW16_hyperpar
##' @export
##' @param m A trait matrix
##' @param s A default strategy
##' @param filter Logical, indicating if generated parameters that are
##' the same as the default should be removed.
FFW16_hyperpar <- make_FFW16_hyperpar()

##' @rdname FFW16_hyperpar
##' @export
FFdev_hyperpar <- make_FFW16_hyperpar()

##' @rdname FFW16_hyperpar
##' @param type Either \code{"FFW16"} or \code{"FFdev"}.
##' @export
make_hyperpar <- function(type) {
  switch(type,
         FFW16=make_FFW16_hyperpar,
         FFdev=make_FFW16_hyperpar,
         stop("Unknown type ", type))
}

##' @rdname FFW16_hyperpar
##' @export
hyperpar <- function(type) {
  switch(type,
         FFW16=FFW16_hyperpar,
         FFdev=FFW16_hyperpar,
         stop("Unknown type ", type))
}

##' Helper function for creating parameter objects suitable for an
##' assembly.
##' @title Helper function for creating parameter objects
##' @param ... Named set of parameters
##' @param pars A list of parameters
##' @export
assembly_parameters <- function(..., pars=NULL) {
  p <- plant::ebt_base_parameters()

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
    pos <- setdiff(c(names(formals(make_FFW16_hyperpar)),
                     names(p),
                     names(p$control),
                     names(p$strategy_default)),
                   excl)
    unk <- setdiff(names(pars), pos)
    if (length(unk) > 0L) {
      stop("Unknown parameters: ", paste(unk, collapse=", "))
    }

    nms_hyper <- intersect(names(pars), names(formals(make_FFW16_hyperpar)))
    p$hyperpar <- do.call("make_FFW16_hyperpar", pars[nms_hyper])
    p                  <- modify_list(p,                  pars)
    p$control          <- modify_list(p$control,          pars)
    p$strategy_default <- modify_list(p$strategy_default, pars)
  }
  p
}

ebt_to_internals <- function(obj, use_environment=TRUE) {
  dat <- lapply(seq_along(obj$time), function(i)
    patch_to_internals(ebt_patch(i, obj), use_environment))
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
  sp_pp <- lapply(sp$cohorts, function(x)
    plant_to_plant_plus(x$plant, environment))
  ints <- do.call("rbind", lapply(sp_pp, function(x) unlist(x$internals)))
  cbind(ints,
        log_density=sp$log_densities,
        seeds_survival_weighted=sp$seeds)
}

##' Create a function that allows integrating aggregate properties of
##' the EBT system.
##'
##' The workflow here is to run an EBT to create an EBT by running
##' \code{run_ebt}, or a set of data from \code{run_ebt_collect} and
##' then reconstitute all the intermediate bits of data so that an any
##' variable that \code{PlantPlus} tracks can be integrated out.
##' Because the pre-processing step is reasonably slow, this function
##' returns a function that takes a variable name and integrates it.
##'
##' @title Integrate EBT variables
##' @param ebt An object from \code{run_ebt} or \code{run_ebt_collect}
##' @export
make_ebt_integrate <- function(obj) {
  ## TODO: This needs to be made to work with the output of
  ## run_ebt_collect, which means that it could possibly work with the
  ## output over time, which would be cool; that might help with some
  ## of the stuff from the "emergent" vignette.
  if (inherits(obj, "EBT")) {
    internals <- patch_to_internals(obj$patch)
    n <- length(internals)
    sched <- obj$cohort_schedule
    a <- lapply(seq_len(n), sched$times)
    pa <- lapply(a, obj$patch$environment$disturbance_regime$density)
  } else {
    internals <- patch_to_internals(ebt_patch(length(obj$time), obj))
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
