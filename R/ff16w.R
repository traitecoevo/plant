# Built from  R/ff16.R on Mon Jul 19 11:01:04 2021 using the scaffolder, from the strategy:  FF16
##' Create a FF16w Plant or Node
##' @title Create a FF16w Plant or Node
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
FF16w_Node <- function(s=FF16w_Strategy()) {
  Node("FF16w", "FF16_Env")(s)
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
##'
##' @param soil_number_of_depths the number of soil layers
##' @param canopy_light_tol 
##' @param canopy_light_nbase 
##' @param canopy_light_max_depth 
##' @param canopy_rescale_usually 
##' @param soil_initial_state initial soil moisture (m3 m^-3)
##' @param vpd Vapour pressure deficit (kPa). Can be a constant value or a list (x = time, y = VPD)
##' @param co2 CO2 partial pressure (Pa). Can be a constant value or a list (x = time, y = co2)
##' @param rainfall Rainfall (m). Can be a constant value or a list (x = time, y = rainfall)
FF16w_make_environment <- function(canopy_light_tol = 1e-4, 
                                   canopy_light_nbase = 17,
                                   canopy_light_max_depth = 16, 
                                   canopy_rescale_usually = TRUE,
                                   soil_number_of_depths = 1,
                                   soil_initial_state = 0.0,
                                   rainfall = 1,
                                   vpd = 1,
                                   co2 = 40) {
  
  if(soil_number_of_depths < 1)
    stop("FF16w Environment must have at least one soil layer")
  
  e <- FF16_Environment(canopy_rescale_usually, 
                        soil_number_of_depths)
  
  e$canopy <- Canopy(canopy_light_tol, 
                     canopy_light_nbase, 
                     canopy_light_max_depth)

  # there might be a better way to skip this if using defaults
  if(sum(soil_initial_state) > 0.0) {
    if(soil_number_of_depths != length(soil_initial_state))
      stop("Not enough starting points for all layers")
    
    e$set_soil_water_state(soil_initial_state)
  }
    
  drivers <- ExtrinsicDrivers()
  if (is.list(rainfall)) {
    drivers$set_variable("rainfall", rainfall$x, rainfall$y)
  } else if (is.numeric(rainfall)) {
    drivers$set_constant("rainfall", rainfall)
    drivers$set_extrapolate("rainfall", FALSE)
  } else {
    stop("Invalid type in birth_rate - need either a list with x, y control points or a numeric")
  }
  
  #TO DO: This probably wasn't the nicest way to add VPD, would recommend we think about how to do this better? - Isaac

  if (is.list(vpd)) {
    drivers$set_variable("vpd", vpd$x, vpd$y)
  } else if (is.numeric(vpd)) {
    drivers$set_constant("vpd", vpd)
    drivers$set_extrapolate("vpd", FALSE)
  } else {
    stop("Invalid type in birth_rate - need either a list with x, y control points or a numeric")
  }
  
  if (is.list(co2)) {
    drivers$set_variable("co2", co2$x, co2$y)
  } else if (is.numeric(co2)) {
    drivers$set_constant("co2", co2)
    drivers$set_extrapolate("co2", FALSE)
  } else {
    stop("Invalid type in birth_rate - need either a list with x, y control points or a numeric")
  }
  
  e$extrinsic_drivers <- drivers
  
  return(e)
}
  

##' Construct a fixed environment for FF16w strategy
##'
##' @param e Value of environment (deafult  = 1.0)
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
##'
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
##' @param B_lf1 Beta coefficient for empirical relationship between narea ~ lma [g/m2] (Dong et al. 2022)
##' @param B_lf2 Beta coefficient for empirical relationship between narea ~ vcmax [umol / m2 / s] (Dong et al. 2022)
##' @param B_lf4 CO_2 respiration per unit leaf nitrogen [mol / yr / kg]
##' @param k_I light extinction coefficient [dimensionless]
##' @param a_lf1 intercept for empirical relationship between narea and vcmax, lma (Dong et al. 2022)
##' @param B_Hv1 p50 at K_s_0 [-MPa]
##' @param B_Hv2 Scaling slope for K_s in p50 [dimensionless]
##' @param latitude degrees from equator (0-90), used in solar model [deg]
##' @param K_s_0 Central (mean) value for maximum sapwood conductivity [kg /m2 / s / MPA]

##'
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
                                B_hks1 = 0.2,
                                B_hks2 = 0.0,
                                B_ks2=0.0,
                                B_rs1=4012.0,
                                B_rb1=2.0*4012.0,
                                B_f1 =3.0,
                                a_lf1=0.535, 
                                B_lf1=0.009, 
                                B_lf2=0.004,
                                B_lf3=0.0008,
                                B_lf4=21000,
                                B_lf5= 40000,
                                k_I=0.5,
                                latitude=0,
                                B_Hv1 = 1.731347,
                                B_Hv2 = -0.7246377,
                                K_s_0 = 2) {
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
  assert_scalar(B_hks1)
  assert_scalar(B_hks2)
  assert_scalar(B_ks2)
  assert_scalar(B_rs1)
  assert_scalar(B_rb1)
  assert_scalar(B_f1)
  assert_scalar(a_lf1)
  assert_scalar(B_lf1)
  assert_scalar(B_lf2)
  assert_scalar(B_lf3)
  assert_scalar(B_lf4)
  assert_scalar(B_lf5)
  assert_scalar(B_Hv1)
  assert_scalar(B_Hv2)
  assert_scalar(K_s_0)
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
    K_s     <- with_default("K_s")
    vcmax     <- with_default("vcmax")
    c     <- with_default("c")
    jmax     <- with_default("jmax")
    

    ## lma / leaf turnover relationship:
    k_l   <- B_kl1 * (lma / lma_0) ^ (-B_kl2)

    ## rho / mortality relationship:
    d_I  <- B_dI1 * (rho / rho_0) ^ (-B_dI2)

    ## rho / wood turnover relationship:
    k_s  <- B_ks1 *  (rho / rho_0) ^ (-B_ks2)

    ## rho / moisture-wood turnover relationship:
    hk_s  <- B_hks1 *  (rho / rho_0) ^ (-B_hks2)
    
    ## hard coded model parameters for now
    ## p_50 sapwood specific conductivity turnover:
    p_50 <- B_Hv1*(K_s/K_s_0)^(B_Hv2)

    ## sensitivity parameter hydraulic vulnerability curve, water potential at 37% conductivity remaining (return unitless):
    b <- p_50/((-log(1-50/100))^(1/c))

    ## water potential at critical xylem failure (95%) (return -MPa):
    psi_crit <- b*(log(1/0.05))^(1/c)

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

    ## n_area from structural (lma) and metabolic (vcmax) N (Dong et al. 2022)

    narea_s <- (a_lf1 + B_lf1*lma*1000)/1000
    narea_p <- (B_lf2*vcmax + B_lf3*jmax)/1000
    
    ## Respiration rates are per unit mass, so convert to mass-based
    ## rate by dividing with lma
    ## So respiration rates per unit mass vary with lma, while
    ## respiration rates per unit area don't.
    r_ls  <- B_lf4 * narea_s / lma
    r_lp  <- B_lf5 * narea_p / lma
    
    r_l <- r_ls + r_lp
    
    extra <- cbind(k_l,                 # lma
                   d_I, k_s, r_s, r_b, hk_s,  # rho
                   a_f3,               # omega
                   r_l,                # lma, narea
                   p_50, b, psi_crit)  # K_s, c              

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
