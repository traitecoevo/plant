##' Create a FF16 Individual
##' @title Create a FF16 Individual
##' @param s A \code{\link{FF16_Strategy}} object
##' @export
##' @rdname FF16_individual
##' @examples
##' pl <- FF16_Individual()
##' pl$height
FF16_Individual <- function(s=FF16_Strategy()) {
  Individual("FF16", "FF16_Env")(s)
}

##' @title Setup an a model system with default or specified parameters
##'
##' @description Setup an a model system with default or specified parameters. 
##' This function enables you initialize a model system. Use the model name to start different models. 
##' @param ... Arguments to be passed to the model constructor. These include
##' 
##'   *`patch_area`: Area of idnividfual patch. Only relevant for stochastic model. Default is 1.0m2.
##'   *`max_patch_lifetime`: The maximum time in years we want to simulate
##'   *`strategies`: A list of stratgies to simulate. The default is an empty list.
##'   *`strategy_default`: Values for the default startegy. The default values are those specified in the C++ code for the model.
##'   *`node_schedule_times_default`: Default vector of times at which to introduce nodes. The default is chosen to have close spacing at the start of the simulation.
##'   *`node_schedule_times`: A list with each element containing the vector of times we want to introduce nodes for each strategy. The default is an empty list.
##'   *`ode_times`: A vector of patch ages we want the ode solver to stop at
##' @export
##' @rdname FF16_Parameters
##' @examples
##' p1 <- FF16_Parameters()
##' p2 <- FF16_Parameters(max_patch_lifetime = 10.0, patch_area = 1.0, strategies = list(FF16_Strategy()), strategy_default = FF16_Strategy(), node_schedule_times_default = node_schedule_times_default(10.0), node_schedule_times = list(node_schedule_times_default(10.0)), ode_times = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
FF16_Parameters <- function(...) {
  Parameters("FF16","FF16_Env")(...)
}

## Helper functions to create an XXX_Environment object. Useful for running individuals
##' @title create a model Environment object
##' @param light_availability_spline_tol Error tolerance of adpative spline method. Deafult is 1e-4.
##' @param light_availability_spline_nbase Parameter used in adaptive spline method. Default is 17.
##' @param light_availability_spline_max_depth Parameter used in adaptive spline method. Default is 16.
##' @inheritParams FF16_Environment
##' @export
##' @rdname FF16_make_environment
FF16_make_environment <- function(light_availability_spline_tol = 1e-4, 
                                  light_availability_spline_nbase = 17,
                                  light_availability_spline_max_depth = 16, 
                                  light_availability_spline_rescale_usually = TRUE) {
  
  e <- FF16_Environment(light_availability_spline_rescale_usually, 
                        soil_number_of_depths = 0)
  
  # Shading defaults have lower tolerance which are overwritten for speed
  e$light_availability <- ResourceSpline(light_availability_spline_tol, 
                     light_availability_spline_nbase, 
                     light_availability_spline_max_depth, 
                     light_availability_spline_rescale_usually)
  
  return(e)
}

##' Construct a fixed environment for a model
##'
##' @param e Value of environment (deafult  = 1.0)
##' @param height_max maximum possible height in environment
##' @param ... Additional parameters to be passed to \code{XXXX_make_environment}, where XXXX referes to model name.
##' @rdname FF16_fixed_environment
##'
##' @export
FF16_fixed_environment <- function(e=1.0, height_max = 150.0, ...) {
  env <- FF16_make_environment(...)
  env$set_fixed_environment(e, height_max)
  env
}


##' This makes a pretend light environment over the plant height,
##' slightly concave up, whatever.
##' @title Create a test environment for FF16 startegy
##' @param height top height of environment object
##' @param n number of points
##' @param light_env function for light environment in test object
##' @param n_strategies number of strategies for test environment
##' @rdname FF16_test_environment
##' @examples
##' environment <- plant:::FF16_test_environment(10)
FF16_test_environment <- function(height, n=101, light_env=NULL,
                                  n_strategies=1) {
  
  hh <- seq(0, height, length.out=n)
  if (is.null(light_env)) {
    light_env <- function(x) {
      exp(x/(height*2)) - 1 + (1 - (exp(.5) - 1))/2
    }
  }
  ee <- light_env(hh)
  interpolator <- Interpolator()
  interpolator$init(hh, ee)

  ret <- FF16_make_environment()
  ret$light_availability$spline <- interpolator
  attr(ret, "light_env") <- light_env
  ret
}

##' Generates a report on stand grown with FF16 strategy
##'
##' Builds a detailed report on stand grown with FF16 strategy, based on the template Rmd file provided.  The reports are
##' rendered as html files and saved in the specified output folder.
##'
##' @param results results of runnning \code{run_scm_collect}
##' @param output_file name of output file
##' @param overwrite logical value to determine whether to overwrite existing report
##' @param target_ages Patches ages at which to make plots
##' @param input_file report script (.Rmd) file to build study report
##' @param quiet An option to suppress printing during rendering from knitr, pandoc command line and others.
##'
##' @rdname FF16_generate_stand_report
##' @return html file of the rendered report located in the specified output folder.
##' @export
FF16_generate_stand_report <- function(results,
                                    output_file = "FF16_report.html",
                                    overwrite = FALSE,
                                    target_ages = NA,
                                    input_file = system.file("reports", "FF16_report.Rmd", package = "plant"),
                                    quiet = TRUE) {
  

  output_dir <- dirname(output_file)
  
  if (!file.exists(output_dir)) {
    dir.create(output_dir, FALSE, TRUE)
  }
  
  #output_file <- basename(output_file)

  if (overwrite | !file.exists(output_file)) {
    # knit and render. Note, call render directly
    # in preference to knit, then render, as leaflet widget
    # requires this to work
    result <-
      rmarkdown::render(
        input_file,
        output_dir = output_dir,
        output_file = output_file,
        quiet = quiet,
        params = list(
          results = results,
          target_ages = target_ages
        )
    )

    # remove temporary Rmd
    message(sprintf("Report for FF16 stand saved at %s", output_file))
  } else {
    message(sprintf("Report for FF16 stand already exists at %s", output_file))
  }
}

##' Hyperparameters for FF16 physiological model
##' @title Hyperparameters for FF16 physiological model
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
##' @rdname make_FF16_hyperpar
make_FF16_hyperpar <- function(
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

    ## Check for infitinte values - these cause issues
    if(any(is.infinite(extra))) {
      stop("Attempt to use infinite value in derived parameters: ",
           paste(colnames(extra)[is.infinite(extra)], collapse=", "))
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

##' Hyperparameter function for FF16 physiological model
##' @title Hyperparameter function for FF16 physiological model
##' @param m A matrix of trait values, as returned by \code{trait_matrix}
##' @param s A strategy object
##' @param filter A flag indicating whether to filter columns. If TRUE, any numbers
##' that are within eps of the default strategy are not replaced.
##' @rdname FF16_hyperpar
##' @export
FF16_hyperpar <- make_FF16_hyperpar()
