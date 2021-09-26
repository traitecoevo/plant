# Built from  R/ff16.R on Mon Jul 19 11:01:04 2021 using the scaffolder, from the strategy:  FF16
##' Create a FF16w Plant or Cohort
##' @title Create a FF16w Plant or Cohort
##' @param s A \code{\link{FF16w_Strategy}} object
##' @export
##' @rdname FF16w
##' @examples
##' pl <- FF16w_Individual()
##' pl$height
FF16w_Individual <- function(s=FF16w_Strategy()) {
  Individual("FF16w", "FF16_Env")(s)
}

##' @export
##' @rdname FF16w
FF16w_Cohort <- function(s=FF16w_Strategy()) {
  Cohort("FF16w", "FF16_Env")(s)
}

##' @export
##' @rdname FF16w
FF16w_Species <- function(s=FF16w_Strategy()) {
  Species("FF16w", "FF16_Env")(s)
}

##' @export
##' @rdname FF16w
FF16w_Parameters <- function() {
  Parameters("FF16w","FF16_Env")()
}

##' @export
##' @rdname FF16w
##' @param p A \code{Parameters<FF16w,FF16_Env>} object
FF16w_Patch <- function(p) {
  Patch("FF16w", "FF16_Env")(p)
}

##' @export
##' @rdname FF16w
FF16w_SCM <- function(p) {
  SCM("FF16w", "FF16_Env")(p)
}

##' @export
##' @rdname FF16w
FF16w_StochasticSpecies <- function(s=FF16w_Strategy()) {
  StochasticSpecies("FF16w", "FF16_Env")(s)
}

##' @export
##' @rdname FF16w
FF16w_StochasticPatch <- function(p) {
  StochasticPatch("FF16w", "FF16_Env")(p)
}

##' @export
##' @rdname FF16w
FF16w_StochasticPatchRunner <- function(p) {
  StochasticPatchRunner("FF16w", "FF16_Env")(p)
}


## Helper:
##' @export
##' @rdname FF16_Environment
##' @param infil_rate rate of water entering the first layer
##' @param n_layers the number of layers
##' @param init starting conditions
FF16w_make_environment <- function(canopy_light_tol = 1e-4,
                                   canopy_light_nbase = 17,
                                   canopy_light_max_depth = 16,
                                   canopy_rescale_usually = TRUE,
                                   soil_number_of_depths = 1,
                                   soil_initial_state = 0.0,
                                   soil_infiltration_rate = 0.0) {
  
  if (soil_number_of_depths > 0 &&
      soil_number_of_depths != length(soil_initial_state))
    stop("Not enough starting points for all layers")
  
  FF16_Environment(
    canopy_light_tol,
    canopy_light_nbase,
    canopy_light_max_depth,
    canopy_rescale_usually,
    soil_number_of_depths,
    soil_initial_state,
    soil_infiltration_rate
  )
}
  

##' Construct a fixed environment for FF16w strategy
##'
##' @param e Value of environment (deafult  = 1.0)
##' @param ctrl Control object
##' @param height_max = 150.0 maximum possible height in environment
##' @rdname FF16_Environment
##'
##' @export
FF16w_fixed_environment <- function(e=1.0, height_max = 150.0) {
  env <- FF16w_make_environment()
  env$set_fixed_environment(e, height_max)
  env
}


##' This makes a pretend light environment over the plant height,
##' slightly concave up, whatever.
##' @title Create a test environment for FF16w startegy
##' @param height top height of environment object
##' @param n number of points
##' @param light_env function for light environment in test object
##' @param n_strategies number of strategies for test environment
##' @param birth_rate birth_rate for test environment
##' @export
##' @rdname FF16w_test_environment
##' @examples
##' environment <- FF16w_test_environment(10)
FF16w_test_environment <- function(height,
                                   n = 101,
                                   light_env = NULL,
                                   n_strategies = 1) {
  
  hh <- seq(0, height, length.out = n)
  if (is.null(light_env)) {
    light_env <- function(x) {
      exp(x / (height * 2)) - 1 + (1 - (exp(.5) - 1)) / 2
    }
  }
  ee <- light_env(hh)
  interpolator <- Interpolator()
  interpolator$init(hh, ee)
  
  ret <- FF16w_make_environment()
  ret$canopy$canopy_interpolator <- interpolator
  attr(ret, "light_env") <- light_env
  ret
}

##' Hyperparameters for FF16w physiological model
##' @title Hyperparameters for FF16w physiological model
##' @param lma_0 Central (mean) value for leaf mass per area [kg /m2]
##' @param B_kl1 Rate of leaf turnover at lma_0 [/yr]
##' @param B_kl2 Scaling slope for phi in leaf turnover [dimensionless]
##' @param rho_0 Central (mean) value for wood density [kg /m3]
##' @param B_dI1 Rate of instantaneous mortality at rho_0 [/yr]
##' @param B_dI2 Scaling slope for wood density in intrinsic mortality [dimensionless]
##' @param B_ks1 Rate of sapwood turnover at rho_0 [/yr]
##' @param B_ks2 Scaling slope for rho in sapwood turnover [dimensionless]
##' @param B_rs1 CO_2 respiration per unit sapwood volume [mol / yr / m3 ]
##' @param B_rb1 CO_2 respiration per unit sapwood volume [mol / yr / m3 ]
##' @param B_f1 Cost of seed accessories per unit seed mass [dimensionless]
##' @param narea nitrogen per leaf area [kg / m2]
##' @param narea_0 central (mean) value for nitrogen per leaf area [kg / m2]
##' @param B_lf1 Potential CO_2 photosynthesis at average leaf nitrogen [mol / d / m2]
##' @param B_lf2 Curvature of leaf photosynthetic light response curve [dimensionless]
##' @param B_lf3 Quantum yield of leaf photosynthetic light response curve [dimensionless]
##' @param B_lf4 CO_2 respiration per unit leaf nitrogen [mol / yr / kg]
##' @param B_lf5 Scaling exponent for leaf nitrogen in maximum leaf photosynthesis [dimensionless]
##' @param k_I light extinction coefficient [dimensionless]
##' @param latitude degrees from equator (0-90), used in solar model [deg]
##' @importFrom stats coef nls
##' @export
##' @rdname FF16w_hyperpar
make_FF16w_hyperpar <- function(
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
                                narea=1.87e-3,
                                narea_0=1.87e-3,
                                B_lf1=5120.738 * 1.87e-3 * 24 * 3600 / 1e+06,
                                B_lf2=0.5,
                                B_lf3=0.04,
                                B_lf4=21000,
                                B_lf5=1,
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
  assert_scalar(narea)
  assert_scalar(narea_0)
  assert_scalar(B_lf1)
  assert_scalar(B_lf2)
  assert_scalar(B_lf3)
  assert_scalar(B_lf4)
  assert_scalar(B_lf5)
  assert_scalar(k_I)
  assert_scalar(latitude)

  function(m, s, filter=TRUE) {
    with_default <- function(name, default_value=s[[name]]) {
      rep_len(if (name %in% colnames(m)) m[, name] else default_value,
              nrow(m))
    }
    lma       <- with_default("lma")
    rho       <- with_default("rho")
    omega     <- with_default("omega")
    narea     <- with_default("narea", narea)

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

    ## Narea, photosynthesis, respiration

    assimilation_rectangular_hyperbolae <- function(I, Amax, theta, QY) {
      x <- QY * I + Amax
      (x - sqrt(x^2 - 4 * theta * QY * I * Amax)) / (2 * theta)
    }

    ## Photosynthesis  [mol CO2 / m2 / yr]
    approximate_annual_assimilation <- function(narea, latitude) {
      E <- seq(0, 1, by=0.02)
      ## Only integrate over half year, as solar path is symmetrical
      D <- seq(0, 365/2, length.out = 10000)
      I <- PAR_given_solar_angle(solar_angle(D, latitude = abs(latitude)))

      Amax <- B_lf1 * (narea/narea_0) ^  B_lf5
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

##' Hyperparameter function for FF16w physiological model
##' @title Hyperparameter function for FF16w physiological model
##' @param m A matrix of trait values, as returned by \code{trait_matrix}
##' @param s A strategy object
##' @param filter A flag indicating whether to filter columns. If TRUE, any numbers
##' that are within eps of the default strategy are not replaced.
##' @export
FF16w_hyperpar <- make_FF16w_hyperpar()
