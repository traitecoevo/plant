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
  ebt <- FFW16_EBT(p)
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
##' @export
ebt_base_parameters <- function() {
  ctrl <- equilibrium_verbose(fast_control())
  ctrl$schedule_eps <- 0.005
  ctrl$equilibrium_eps <- 1e-3
  FFW16_Parameters(patch_area=1.0, control=ctrl,
                   hyperpar=FFW16_hyperpar)
}

##' Run the EBT model, given a Parameters and CohortSchedule
##'
##' This is mostly a simple wrapper around some of the EBT functions.
##' Not sure if this is how we will generally want to do this.
##' Consider this function liable to change.
##'
##' @title Run the EBT, Collecting Output
##' @param p A \code{\link{FFW16_Parameters}} object
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

  ebt <- FFW16_EBT(p)
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

  ret <- list(time=time, species=species, light_env=light_env,
              seed_rain=ebt$seed_rains)

  if (include_area_leaf) {
    ret$area_leaf <- do.call("rbind", lapply(res, "[[", "area_leaf"))
  }

  ret
}

run_ebt_error <- function(p) {
  ebt <- FFW16_EBT(p)
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

##' Hyperparameters for tree
##' @title Hyperparameters for tree
##' @param B4 Slope of lma / leaf turnover log-log relationship
##' @param lma_0 LMA value...
##' @param k_l_0 ...
##' @param rho_0 ...
##' @param d0_0 ...
##' @param d1 ...
##' @param k_s_0 ...
##' @param B5 ...
##' @param narea_0 ...
##' @export
##' @rdname FFW16_hyperpar
make_FFW16_hyperpar <- function(B4=1.71,
                               lma_0=0.1978791,
                               k_l_0=0.4565855,
                               rho_0=608.0,
                               d0_0=0.01,
                               d1=0.0,
                               k_s_0=0.2,
                               B5=0.0,
                               narea_0=1.87e-3) {
  force(B4)
  force(lma_0)
  force(k_l_0)
  force(rho_0)
  force(d0_0)
  force(d1)
  force(B5)
  function(m, s, filter=TRUE) {
    with_default <- function(name, default_value=s[[name]]) {
      rep_len(if (name %in% colnames(m)) m[, name] else default_value,
              nrow(m))
    }
    lma       <- with_default("lma")
    rho       <- with_default("rho")
    mass_seed <- with_default("mass_seed")
    narea     <- with_default("narea", narea_0)

    ## lma / leaf turnover relationship:
    k_l   <- k_l_0 * (lma / lma_0) ^ (-B4)

    ## rho / mortality relationship:
    c_d0  <- d0_0 * (rho / rho_0) ^ (-d1)

    ## rho / wood turnover relationship:
    k_s  <- k_s_0 *  (rho / rho_0) ^ (-B5)

    ## rho / sapwood respiration relationship:

    ## Respiration rates are per unit mass, so this next line has the
    ## effect of holding constant the respiration rate per unit volume.
    ## So respiration rates per unit mass vary with rho, respiration
    ## rates per unit volume don't.
    c_Rs <- 4012.0 / rho
    c_Rb <- 2.0 * c_Rs # bark respiration follows from sapwood

    ## mass_seed / accessory cost relationship
    c_acc <- 3.0 * mass_seed

    ## narea / photosynthesis / respiration
    ## Photosynthesis per mass leaf N [mol CO2 / kgN / yr]
    ## TODO: more transparency needed for this calculation. Current
    ## value comes from offline model
    ## TODO: this is where we improve photosynthesis model
    ## TODO: Move control of this into the arguments
    c_PN  <- 80406.42
    c_p1  <- c_PN * narea

    ## Respiration per mass leaf N [mol CO2 / kgN / yr]
    ## = (6.66e-4 * (365*24*60*60))
    ## Obatined from global average of ratio of dark respiration rate to
    ## leaf nitrogen content using the GLOPNET dataset
    c_RN  <-  21000
    ## Respiration rates are per unit mass, so convert to mass-based
    ## rate by dividing with lma
    ## So respiration rates per unit mass vary with lma, while
    ## respiration rates per unit area don't.
    c_Rl  <- c_RN * narea / lma

    extra <- cbind(k_l,                   # lma
                   c_d0, k_s, c_Rs, c_Rb, # rho
                   c_acc,                 # mass_seed
                   c_p1, c_Rl)            # narea

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
