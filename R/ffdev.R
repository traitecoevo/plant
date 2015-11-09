## We can probably actually do better than this with an S3 method on
## the actual strategy?  That would need to be organised by the
## templating though and that's stretched to the limit.

##' Create a FFdev Plant or Cohort
##' @title Create a FFdev Plant or Cohort
##' @param s A \code{\link{FFdev_Strategy}} object
##' @export
##' @rdname FFdev
##' @examples
##' pl <- FFdev_Plant()
##' pl$height
FFdev_Plant <- function(s=FFdev_Strategy()) {
  Plant("FFdev")(s)
}

##' @export
##' @rdname FFdev
FFdev_Cohort <- function(s=FFdev_Strategy()) {
  Cohort("FFdev")(s)
}

##' @export
##' @rdname FFdev
FFdev_Species <- function(s=FFdev_Strategy()) {
  Species("FFdev")(s)
}

##' @export
##' @rdname FFdev
##' @param ... Arguments!
FFdev_Parameters <- function(...) {
  Parameters("FFdev")(...)
}

##' @export
##' @rdname FFdev
##' @param p A \code{Parameters<FFdev>} object
FFdev_Patch <- function(p) {
  Patch("FFdev")(p)
}

##' @export
##' @rdname FFdev
FFdev_SCM <- function(p) {
  SCM("FFdev")(p)
}

##' @export
##' @rdname FFdev
FFdev_StochasticSpecies <- function(s=FFdev_Strategy()) {
  StochasticSpecies("FFdev")(s)
}

##' @export
##' @rdname FFdev
FFdev_StochasticPatch <- function(p) {
  StochasticPatch("FFdev")(p)
}

##' @export
##' @rdname FFdev
FFdev_StochasticPatchRunner <- function(p) {
  StochasticPatchRunner("FFdev")(p)
}

##' @export
##' @rdname FFdev
FFdev_PlantPlus <- function(s=FFdev_Strategy()) {
  PlantPlus("FFdev")(s)
}


##' Hyperparameters for FFdev physiological model
##' @title Hyperparameters for FFdev physiological model
##' @param lma_0 Central (mean) value for leaf mass per area
##' @param B_kl1 Rate of leaf turnover at phi_0
##' @param B_kl2 Scaling slope for phi in leaf turnover
##' @param rho_0 Central (mean) value for wood density
##' @param B_dI1 Rate of instantaneous mortality at rho_0
##' @param B_dI2 Scaling slope for wood density in intrinsic mortality
##' @param B_ks1 Rate of sapwood turnover at rho_0
##' @param B_ks2 Scaling slope for rho in sapwood turnover
##' @param B_rs1 CO_2 respiration per unit sapwood volume
##' @param B_rb1 CO_2 respiration per unit sapwood volume
##' @param B_f1 Cost of seed accessories per unit seed mass
##' @param narea_0 central (mean) value for nitrogen per leaf area
##' @param B_lf1 Potential CO_2 photosynthesis per unit narea
##' @param B_lf2 Curvature of leaf photosynthetic light response curve
##' @param B_lf3 Quantum yield of leaf photosynthetic light response curve
##' @param B_lf4 CO_2 respiration per unit leaf nitrogen
##' @param k_I light extinction coefficient
##' @param latitude degrees from equator (0-90), used in solar model
##' @export
##' @rdname FFdev_hyperpar
make_FFdev_hyperpar <- function(
                                lma_0=0.1978791,
                                B_kl1=0.4565855,
                                B_kl2=1.71,
                                rho_0=608.0,
                                B_dI1=0.01,
                                B_dI2=0.0,
                                B_ks1=0.2,
                                B_ks2=0.0,
                                B_rs1=4012.0,
                                B_rb1=2.0*4012.0,
                                B_f1 =3.0,
                                narea_0=1.87e-3,
                                B_lf1=5120.738 * 24 * 3600 / 1e+06,
                                B_lf2=0.5,
                                B_lf3=0.04,
                                B_lf4=21000,
                                k_I=0.5,
                                latitude=0) {
  assert_scalar <- function(x, name=deparse(substitute(x))) {
    if (length(x) != 1L) {
      stop(sprintf("%s must be a scalar", name), call. = FALSE)
    }
  }
  assert_scalar(lma_0)
  assert_scalar(B_kl1)
  assert_scalar(B_kl2)
  assert_scalar(rho_0)
  assert_scalar(B_dI1)
  assert_scalar(B_dI2)
  assert_scalar(B_ks1)
  assert_scalar(B_ks2)
  assert_scalar(B_rs1)
  assert_scalar(B_rb1)
  assert_scalar(B_f1)
  assert_scalar(narea_0)
  assert_scalar(B_lf1)
  assert_scalar(B_lf2)
  assert_scalar(B_lf3)
  assert_scalar(B_lf4)
  assert_scalar(k_I)
  assert_scalar(latitude)

  ## TODO: k_I should actually be in default parameter set, so perhaps don't pass into function?

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
    k_l   <- B_kl1 * (lma / lma_0) ^ (-B_kl2)

    ## rho / mortality relationship:
    d_I  <- B_dI1 * (rho / rho_0) ^ (-B_dI2)

    ## rho / wood turnover relationship:
    k_s  <- B_ks1 *  (rho / rho_0) ^ (-B_ks2)

    ## rho / sapwood respiration relationship:

    ## Respiration rates are per unit mass, so this next line has the
    ## effect of holding constant the respiration rate per unit volume.
    ## So respiration rates per unit mass vary with rho, respiration
    ## rates per unit volume don't.
    r_s <- B_rs1 / rho
    # bark respiration follows from sapwood
    r_b <- B_rb1 / rho

    ## omega / accessory cost relationship
    a_f3 <- B_f1 * omega

    ## narea / photosynthesis / respiration
    ## Photosynthesis per mass leaf N [mol CO2 / kgN / yr]
    assimilation_rectangular_hyperbolae <- function(I, Amax, theta, QY) {
      x <- QY * I + Amax
      (x - sqrt(x^2 - 4 * theta * QY * I * Amax)) / (2 * theta)
    }

    approximate_annual_assimilation <- function(narea, latitude) {
      E <- seq(0, 1, by=0.02)
      ## Only integrate over half year, as solar path is symmetrical
      D <- seq(0, 365/2, length.out = 10000)
      I <- PAR_given_solar_angle(solar_angle(D, latitude = abs(latitude)))

      Amax <- narea * B_lf1
      theta <- B_lf2
      QY <- B_lf3

      AA <- NA * E

      for (i in seq_len(length(E))) {
        AA[i] <- 2 * trapezium(D, assimilation_rectangular_hyperbolae(
                                    k_I * I * E[i], Amax, theta, QY))
      }
      if(all(diff(AA) < 1E-8)) {
        # line fitting will fail if all have are zero, or potentially same value
        ret <- c(last(AA), 0)
        names(ret) <- c("p1","p2")
      } else {
        fit <- nls(AA ~ p1 * E/(p2 + E), data.frame(E = E, AA = AA), start = list(p1 = 100, p2 = 0.2))
        ret <- coef(fit)
      }
      ret
    }

    # This needed in case narea has length zero, in which case trapezium fails
    a_p1 <- a_p2 <- 0 * narea
    ## TODO: Remove the 0.5 hardcoded default for k_I here, and deal
    ## with this more nicely.
    if (length(narea) > 0 || k_I != 0.5) {
      i <- match(narea, unique(narea))
      y <- vapply(unique(narea), approximate_annual_assimilation,
                  numeric(2), latitude)
      a_p1  <- y["p1", i]
      a_p2  <- y["p2", i]
    }

    ## Respiration rates are per unit mass, so convert to mass-based
    ## rate by dividing with lma
    ## So respiration rates per unit mass vary with lma, while
    ## respiration rates per unit area don't.
    r_l  <- B_lf4 * narea / lma

    extra <- cbind(k_l,                # lma
                   d_I, k_s, r_s, r_b, # rho
                   a_f3,               # omega
                   a_p1, a_p2,         # narea
                   r_l)                # lma, narea

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

    if (!is.null(extra)) {
      m <- cbind(m, extra)
    }
    m
  }
}

##' @title Hyperparameters for FFdev physiological model
##' @rdname FFdev_hyperpar
##' @export
FFdev_hyperpar <- make_FFdev_hyperpar()
