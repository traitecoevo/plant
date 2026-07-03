# Built from  R/TF24.R on Mon Feb 12 09:52:27 2024 using the scaffolder, from the strategy:  TF24
##' @title Create a TF24 Plant or Node
##' @param s A \code{\link{TF24_Strategy}} object
##' @export
##' @rdname TF24_Individual
##' @examples
##' pl <- TF24_Individual()
##' pl$height
TF24_Individual <- function(s=TF24_Strategy()) {
  Individual("TF24", "TF24_Env")(s)
}

##' @title Setup an a TF24 system with  default or specified parameters
##' @description Setup an a model system with default or specified parameters.
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
##' @rdname TF24_Parameters
##' @examples
##' p1 <- TF24_Parameters()
##' p2 <- TF24_Parameters(max_patch_lifetime = 10.0)
TF24_Parameters <- function(...) {
  Parameters("TF24","TF24_Env")(...)
}

##' Generates a report on stand grown with TF24 strategy
##'
##' Builds a detailed report on stand grown with TF24 strategy, based on the template Rmd file provided.  The reports are
##' rendered as html files and saved in the specified output folder.
##'
##' @inheritParams FF16_generate_stand_report
##' @rdname TF24_generate_stand_report
##' @return html file of the rendered report located in the specified output folder.
##' @export
TF24_generate_stand_report <- function(results,
                                    output_file = "TF24_report.html",
                                    overwrite = FALSE,
                                    target_ages = NA,
                                    input_file = system.file("reports", "TF24_report.Rmd", package = "plant"),
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
    message(sprintf("Report for TF24 stand saved at %s", output_file))
  } else {
    message(sprintf("Report for TF24 stand already exists at %s", output_file))
  }
}

##' Hyperparameters for TF24 physiological model
##' @title Hyperparameters for TF24 physiological model
##' @param lma_0 Central (mean) value for leaf mass per area [kg /m2]
##' @param B_kl1 Rate of leaf turnover at lma_0 [/yr]
##' @param B_kl2 Scaling slope for phi in leaf turnover [dimensionless]
##' @param rho_0 Central (mean) value for wood density [kg /m3]
##' @param B_dI1 Rate of instantaneous mortality at rho_0 [/yr]
##' @param B_dI2 Scaling slope for wood density in intrinsic mortality [dimensionless]
##' @param B_hks1 Intercept for the g1_TF24 ~ rho relationship at rho_0 [dimensionless]
##' @param B_hks2 Scaling slope for rho in the g1_TF24 relationship [dimensionless]
##' @param B_ks1 Rate of sapwood turnover at rho_0 [/yr]
##' @param B_ks2 Scaling slope for rho in sapwood turnover [dimensionless]
##' @param B_rs1 CO_2 respiration per unit sapwood volume [mol / yr / m3 ]
##' @param B_rb1 CO_2 respiration per unit sapwood volume [mol / yr / m3 ]
##' @param B_f1 Cost of seed accessories per unit seed mass [dimensionless]
##' @param B_lf1 Beta coefficient for empirical relationship between narea_ls ~ lma [g/m2] (Dong et al. 2022)
##' @param B_lf2 Beta coefficient for empirical relationship between narea_lp ~ vcmax [g/m2] (Dong et al. 2022)
##' @param B_lf3 Beta coefficient for empirical relationship between narea_lp ~ jmax [umol / m2 / s] (Dong et al. 2022)
##' @param B_lf4 CO_2 respiration per unit structural leaf nitrogen [mol / yr / kg]
##' @param B_lf5 CO_2 respiration per unit photosynthetic leaf nitrogen [mol / yr / kg]
##' @param a_lf1 intercept for empirical relationship between narea and vcmax, lma (Dong et al. 2022)
##' @param B_Hv1 p50 at K_s = 1 [-MPa]
##' @param B_Hv2 Scaling slope for K_s in p50 [dimensionless]
##' @param B_c1 Shape parameter c of the vulnerability curve at p_50 = 0 [dimensionless]
##' @param B_c2 Scaling slope for p_50 in the vulnerability-curve shape parameter c [dimensionless]
##' @param latitude degrees from equator (0-90), used in solar model [deg]
##' @export
##' @rdname make_TF24_hyperpar
make_TF24_hyperpar <- function(lma_0=0.1978791,
                                B_kl1=0.4565855,
                                B_kl2=1.71,
                                rho_0=608.0,
                                B_dI1=0.01,
                                B_dI2=0.0,
                                B_hks1=7.5,
                                B_hks2=0.0,
                                B_ks1=0.2,
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
                                latitude=0,
                                B_Hv1 = 0.4607063,
                                B_Hv2 = -0.2,
                                B_c1 = 2.04,
                                B_c2 = 0) {


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
  assert_scalar(B_hks1)
  assert_scalar(B_hks2)
  assert_scalar(B_ks1)
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
  assert_scalar(B_c1)
  assert_scalar(B_c2)
  assert_scalar(latitude)

  function(m, s, filter=TRUE) {
    with_default <- function(name, default_value=s$pars[[name]]) {
      rep_len(if (name %in% colnames(m)) m[, name] else default_value,
              nrow(m))
    }
    lma       <- with_default("lma")
    rho       <- with_default("rho")
    omega     <- with_default("omega")
    K_s     <- with_default("K_s")
    vcmax_25     <- with_default("vcmax_25")
    jmax_25     <- with_default("jmax_25")


    ## lma / leaf turnover relationship:
    k_l   <- B_kl1 * (lma / lma_0) ^ (-B_kl2)

    ## rho / mortality relationship:
    d_I  <- B_dI1 * (rho / rho_0) ^ (-B_dI2)

    ## Reuse the legacy hk_s parameterisation to derive g1_TF24:
    g1_TF24 <- B_hks1 * (rho / rho_0) ^ (-B_hks2)

    ## rho / wood turnover relationship:
    k_s  <- B_ks1 *  (rho / rho_0) ^ (-B_ks2)

    ## TODO: Convert the p50 ks function back to a mean centred function using K_s_0

    ## p_50 sapwood specific conductivity turnover:
    if (any(K_s <= 0, na.rm = TRUE)) {
      stop("K_s must be > 0 for p_50 derivation", call. = FALSE)
    }
    p_50 <- 10^(B_Hv1 + B_Hv2 * log10(K_s))

    ## p_50 shape parameter trade off
    c <- B_c1 * exp(-B_c2 * p_50)
    ## scale parameter b of the vulnerability curve exp(-(psi/b)^c): the water
    ## potential at 1/e (~37%) conductivity remaining, solved here from the 50%
    ## loss-of-conductivity point p_50 [-MPa]:
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

    narea_ls <- (a_lf1 + B_lf1*lma*1000)/1000
    narea_lp <- (B_lf2*vcmax_25 + B_lf3*jmax_25)/1000
    
    
    nmass_ls <- narea_ls / lma
    nmass_lp <- narea_lp / lma
    nmass_l <- nmass_ls + nmass_lp
    ## Respiration rates are per unit mass, so convert to mass-based
    ## rate by dividing with lma
    ## So respiration rates per unit mass vary with lma, while
    ## respiration rates per unit area don't.
    r_ls  <- B_lf4 * nmass_ls
    r_lp  <- B_lf5 * nmass_lp
    
    r_l <- r_ls + r_lp

    extra <- cbind(k_l,                        # lma
                   d_I, g1_TF24, k_s, r_s, r_b, # rho
                   a_f3,                        # omega
                   r_l, nmass_l,                # lma, narea
                   c, p_50, b, psi_crit)        # K_s

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
          x2 <- unlist(s$pars[names(x1)])
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

##' Hyperparameter function for TF24 physiological model
##' @title Hyperparameter function for TF24 physiological model
##' @param m A matrix of trait values, as returned by \code{trait_matrix}
##' @param s A strategy object
##' @param filter A flag indicating whether to filter columns. If TRUE, any numbers
##' that are within eps of the default strategy are not replaced.
##' @export
##' @rdname TF24_hyperpar
TF24_hyperpar <- make_TF24_hyperpar()

#' @export
#' @importFrom rlang .data
#' @rdname expand_state
TF24_expand_state <- function(results) {
  data <- split(results$species, results$species$species)

  for (i in seq_len(results$n_spp)) {
    s <- results$p$strategies[[i]]
    d <- data[[i]]
    # Derived size/mass columns are computed by the strategy's own C++
    # allometry functions (see src/strategy_expand.cpp) rather than being
    # re-derived here, so the formulas live in exactly one place.
    allom <- TF24_strategy_expand_allometry(s, d$height, d$area_heartwood,
                                            d$mass_heartwood)
    data[[i]] <- dplyr::bind_cols(d, allom)
  }

  results$species <- data %>% dplyr::bind_rows()

  results
}
