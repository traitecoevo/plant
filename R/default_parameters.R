## plant's testing trade-off functions, the Parameters slot they travel in, and
## the parameter table. The recommended, living trade-off functions and the
## curated species/site parameter sets live in the fieldtrial package; what is
## here is frozen so that plant's tests and goldens stay self-contained.

##' Plant's testing parameter set: today's defaults plus frozen trade-offs
##'
##' \code{default_parameters(model)} is \code{\link{scm_base_parameters}(model)}
##' with plant's \emph{testing} trade-off function attached as
##' \code{p$hyperpar}. It exists so that every test and golden in plant is
##' self-contained: the trade-offs it attaches are the ones plant's defaults
##' were tuned with, frozen as they stood when they moved out of the model
##' files, and they are \strong{not} a recommended parameterisation of any
##' species or site. Curated species and site parameter sets, and the
##' recommended, living trade-off functions, are in the \pkg{fieldtrial}
##' package. Changing the frozen copy here means re-blessing plant's goldens.
##'
##' A \code{Parameters} object carries its own trade-offs in \code{p$hyperpar}:
##' an R closure satisfying the contract \code{f(m, s, filter = TRUE)} (a trait
##' matrix in, the matrix with derived parameter columns appended out, with
##' \code{filter = TRUE} dropping columns that sit at the strategy's defaults;
##' optionally carrying \code{attr(f, "input_traits")}, the traits it reads
##' that are not fields of \code{pars}), or \code{NULL}.
##' \code{\link{generate_strategy}}, \code{\link{add_strategies}},
##' \code{\link{add_mutant}}, \code{\link{stand_gradient}} and
##' \code{\link{fixed_point_fitness}} read it through \code{param_hyperpar(p)}
##' as their default, so a set built by \pkg{fieldtrial}, \pkg{regnans} or by
##' hand works through one API; passing \code{hyperpar =} explicitly still
##' overrides, and \code{hyperpar = NULL} is the explicit identity. The slot is a
##' field of the C++ \code{Parameters}, held opaquely, so it survives
##' \code{validate()}, \code{\link{run_scm}} (including \code{scm$parameters})
##' and a saved object.
##'
##' The bare \code{scm_base_parameters(model)} attaches no trade-offs: a trait
##' named in a trait matrix is written straight to \code{pars}, and
##' \code{generate_strategy} says so once per session.
##'
##' @param model A model name: \code{"FF16"}, \code{"TF24"}, \code{"TF24f"},
##'   \code{"TF24t"}, \code{"TF24w"} or \code{"K93"}.
##' @param quiet If \code{FALSE} (the default), warn that this is a testing set.
##' @return For \code{default_parameters}, a \code{Parameters} object with
##'   \code{hyperpar} set.
##' @seealso \code{\link{scm_base_parameters}}, \code{\link{parameter_table}},
##'   \code{\link{generate_strategy}}
##' @examples
##' p <- default_parameters("FF16", quiet = TRUE)
##' p <- add_strategies(p, trait_matrix(0.0825, "lma"))
##' p$strategies[[1]]$pars$k_l   # derived from lma by the attached trade-off
##'
##' p0 <- scm_base_parameters("FF16")
##' is.null(p0$hyperpar)          # TRUE: nothing attached
##' @export
##' @rdname default_parameters
default_parameters <- function(model, quiet = FALSE) {
  p <- scm_base_parameters(model)
  p$hyperpar <- default_hyperpar(model)
  if (!quiet) {
    warning("default_parameters(\"", model, "\") is plant's testing parameter ",
            "set: the current C++ defaults with the frozen trade-off functions ",
            "plant's tests were built on. It is not a recommended ",
            "parameterisation of any species or site; see the fieldtrial ",
            "package for curated sets. Pass quiet = TRUE to silence this.",
            call. = FALSE)
  }
  p
}

##' @param ... For \code{default_hyperpar}, coefficient overrides passed to the
##'   frozen maker (\code{B_kl1}, \code{B_Hv1}, \code{narea}, ...; see
##'   \code{\link{frozen_FF16_hyperpar}} and \code{\link{frozen_TF24_hyperpar}}
##'   for the lists). For tests and diagnostics only: with a coefficient moved
##'   the result is no longer the testing set. K93 accepts none.
##' @return For \code{default_hyperpar}, the frozen trade-off closure for
##'   \code{model}. TF24's serves TF24, TF24f, TF24t and TF24w, which share
##'   \code{TF24_Pars}; K93 has no trade-offs and gets the identity function.
##'   The closure carries \code{attr(f, "traits")} (every trait-matrix column
##'   it reads), \code{attr(f, "input_traits")} (those that are not fields of
##'   \code{pars}) and \code{attr(f, "coefficients")} (the coefficients it was
##'   built with, overrides applied).
##' @export
##' @rdname default_parameters
default_hyperpar <- function(model, ...) {
  maker <- switch(model,
    FF16 = frozen_FF16_hyperpar,
    TF24 = ,
    TF24f = ,
    TF24t = ,
    TF24w = frozen_TF24_hyperpar,
    K93 = NULL,
    stop("Unknown type ", model))
  if (is.null(maker)) {
    if (...length() > 0L) {
      stop("K93 has no trade-off coefficients to set", call. = FALSE)
    }
    f <- identity_hyperpar
    coefficients <- list()
  } else {
    f <- maker(...)
    ## The coefficients the closure was built with: the maker's defaults, with
    ## the overrides applied. Recorded so a caller (logpile's fingerprint, a
    ## parameter set's provenance) can read them off the function rather than
    ## off its formals, which is also the shape fieldtrial's closures carry.
    coefficients <- lapply(formals(maker), eval, envir = baseenv())
    overrides <- list(...)
    coefficients[names(overrides)] <- overrides
  }
  ## Every column the function reads from the trait matrix, `pars` fields
  ## included, so parameter_table() can show the traits. `input_traits` (set in
  ## the maker) stays the narrower list generate_strategy() needs: the ones that
  ## are not `pars` fields.
  attr(f, "traits") <- switch(model,
    FF16 = c("lma", "rho", "omega", "narea"),
    K93 = character(0),
    c("lma", "rho", "omega", "K_s", "vcmax_25", "jmax_25"))
  attr(f, "coefficients") <- coefficients
  f
}

##' @param parameters A \code{Parameters} object.
##' @return For \code{param_hyperpar}, the trade-off function \code{parameters}
##'   carries in \code{$hyperpar}, or the identity when it carries none.
##' @export
##' @rdname default_parameters
param_hyperpar <- function(parameters) {
  h <- parameters$hyperpar
  if (is.null(h)) {
    return(unattached_hyperpar)
  }
  if (!is.function(h)) {
    stop("`p$hyperpar` must be a trade-off function f(m, s, filter = TRUE), ",
         "or NULL", call. = FALSE)
  }
  h
}

## The identity trade-off: traits are the parameters. `hyperpar = NULL` means
## this explicitly; param_hyperpar() returns the `unattached` copy when a
## Parameters carries none, which is what lets generate_strategy() say so.
identity_hyperpar <- function(m, s, filter = TRUE) m
attr(identity_hyperpar, "identity") <- TRUE
unattached_hyperpar <- identity_hyperpar
attr(unattached_hyperpar, "unattached") <- TRUE

is_identity_hyperpar <- function(f) {
  is.function(f) && isTRUE(attr(f, "identity", exact = TRUE))
}

## Top-level numeric scalar fields of a strategy beyond the six every model
## carries (`pars`, `control`, the birth-rate trio and `collect_all_auxiliary`):
## TF24t's thermal parameters and TF24f's acclimation gain. These are the
## fields a trait matrix can reach besides `pars` (see generate_strategy()).
strategy_field_names <- function(strategy) {
  infrastructure <- c("pars", "control", "collect_all_auxiliary",
                      "birth_rate_x", "birth_rate_y", "is_variable_birth_rate")
  nms <- setdiff(names(strategy), infrastructure)
  nms[vapply(nms, function(n) {
    v <- strategy[[n]]
    is.numeric(v) && length(v) == 1L
  }, logical(1))]
}

##' The frozen FF16 trade-offs plant tests against
##'
##' The trade-off function FF16's defaults were tuned with, frozen when the
##' recommended trade-offs moved to the \pkg{fieldtrial} package. Reached
##' through \code{\link{default_hyperpar}("FF16", ...)}; not for direct use.
##' It derives \code{k_l} from \code{lma}; \code{d_I}, \code{k_s}, \code{r_s}
##' and \code{r_b} from \code{rho}; \code{a_f3} from \code{omega}; \code{a_p1}
##' and \code{a_p2} from \code{narea} and \code{latitude} by integrating a
##' solar-driven PAR day and fitting a saturating curve; and \code{r_l} from
##' \code{narea} and \code{lma}.
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
##' @param narea nitrogen per leaf area [kg / m2]
##' @param narea_0 central (mean) value for nitrogen per leaf area [kg / m2]
##' @param B_lf1 Potential CO_2 photosynthesis at average leaf nitrogen [mol / d / m2]
##' @param B_lf2 Curvature of leaf photosynthetic light response curve [dimensionless]
##' @param B_lf3 Quantum yield of leaf photosynthetic light response curve [dimensionless]
##' @param B_lf4 CO_2 respiration per unit leaf nitrogen [mol / yr / kg]
##' @param B_lf5 Scaling exponent for leaf nitrogen in maximum leaf photosynthesis [dimensionless]
##' @param latitude degrees from equator (0-90), used in solar model [deg]
##' @importFrom stats coef nls
##' @return A trade-off closure \code{f(m, s, filter = TRUE)} carrying
##'   \code{attr(f, "input_traits") = "narea"}.
##' @importFrom stats coef nls
##' @keywords internal
##' @rdname frozen_FF16_hyperpar
frozen_FF16_hyperpar <- function(
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
  assert_scalar(latitude)

  f <- function(m, s, filter=TRUE) {
    with_default <- function(name, default_value=s$pars[[name]]) {
      rep_len(if (name %in% colnames(m)) m[, name] else default_value,
              nrow(m))
    }
    lma       <- with_default("lma")
    rho       <- with_default("rho")
    omega     <- with_default("omega")
    narea     <- with_default("narea", narea)

    ## Light extinction coefficient: source from the strategy so the
    ## assimilation integral below stays consistent with the canopy
    ## model, rather than carrying a separate hard-coded default here.
    k_I       <- s$pars$k_I

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
        ret <- c(dplyr::last(AA), 0)
        names(ret) <- c("p1","p2")
      } else {
        fit <- nls(AA ~ p1 * E/(p2 + E), data.frame(E = E, AA = AA), start = list(p1 = 100, p2 = 0.2))
        ret <- coef(fit)
      }
      ret
    }

    # Guard against narea having length zero (e.g. an empty trait matrix),
    # in which case trapezium() fails.
    a_p1 <- a_p2 <- 0 * narea
    if (length(narea) > 0) {
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

  ## The traits this hyperpar reads from the MATRIX that are not fields of the
  ## strategy's `pars` (#636). `generate_strategy()` refuses an unknown trait
  ## name, and it cannot know about these from `pars` alone: `narea` is consumed
  ## here to derive `a_p1`, `a_p2` and `r_l`, and is never itself assigned to the
  ## strategy -- so without this declaration a legitimate `narea` sweep would be
  ## rejected.
  ##
  ## ⚠️ ONLY THE ONES READ FROM `m`, i.e. those with a `with_default()` call above.
  ## The nineteen CONSTRUCTOR arguments of this function (`B_kl1`, `B_ks1`,
  ## `B_lf4`, ...) are not traits: they parameterise the trade-offs once, when the
  ## hyperpar is built. Naming one in a trait matrix has never done anything, and
  ## listing them here would restore exactly the silence #636 removes --
  ## `tests/testthat/test-ode-individual-runner.R` had been passing nine of them
  ## for years, to no effect.
  attr(f, "input_traits") <- "narea"
  f
}

##' The frozen TF24 trade-offs plant tests against
##'
##' The trade-off function TF24's defaults were tuned with, frozen when the
##' recommended trade-offs moved to the \pkg{fieldtrial} package. Reached
##' through \code{\link{default_hyperpar}("TF24", ...)}, and shared by TF24f,
##' TF24t and TF24w; not for direct use. It derives \code{k_l} from
##' \code{lma}; \code{d_I}, \code{TF24_cost_scale}, \code{k_s}, \code{r_s} and
##' \code{r_b} from \code{rho}; \code{a_f3} from \code{omega}; \code{stem_P50},
##' \code{stem_c}, \code{stem_b} and \code{psi_crit} from \code{K_s}; and
##' \code{nmass_l} and \code{r_l} from \code{lma}, \code{vcmax_25} and
##' \code{jmax_25} (Dong et al. 2022): thirteen parameters. TF24's C++ default
##' \code{stem_P50} is this function's value at the default \code{K_s}, so a
##' strategy built from traits and the bare default share one curve.
##'
##' @param lma_0 Central (mean) value for leaf mass per area [kg /m2]
##' @param B_kl1 Rate of leaf turnover at lma_0 [/yr]
##' @param B_kl2 Scaling slope for phi in leaf turnover [dimensionless]
##' @param rho_0 Central (mean) value for wood density [kg /m3]
##' @param B_dI1 Rate of instantaneous mortality at rho_0 [/yr]
##' @param B_dI2 Scaling slope for wood density in intrinsic mortality [dimensionless]
##' @param B_hks1 Intercept for the TF24_cost_scale ~ rho relationship at rho_0 [dimensionless]
##' @param B_hks2 Scaling slope for rho in the TF24_cost_scale relationship [dimensionless]
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
##' @param B_Hv1 p50 at terminal-segment K_s = 1 [-MPa]. The relation is keyed
##'   on a tip conductivity rather than a whole-stem one, so the intercept
##'   describes a different physical quantity than it did before the stem path
##'   integral existed; it was re-anchored to hold stem_P50 fixed at the default
##'   parameters. Whether the Ks-p50 relation should be keyed on tip conductivity
##'   at all is open -- twig segments are what people usually measure.
##' @param B_Hv2 Scaling slope for K_s in p50 [dimensionless]
##' @param B_c1 Shape parameter stem_c of the vulnerability curve at stem_P50 = 0 [dimensionless]
##' @param B_c2 Scaling slope for stem_P50 in the vulnerability-curve shape parameter stem_c [dimensionless]
##' @param latitude degrees from equator (0-90), used in solar model [deg]
##' @return A trade-off closure \code{f(m, s, filter = TRUE)}.
##' @keywords internal
##' @rdname frozen_TF24_hyperpar
frozen_TF24_hyperpar <- function(lma_0=0.1978791,
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
                                # K_s is a TERMINAL-SEGMENT conductivity,
                                # 2.9782x smaller than the whole-stem value this
                                # relation was calibrated against. Left at
                                # 0.4607063 the same relation would have moved
                                # stem_P50 from 2.889 to 3.593 MPa -- buying a height
                                # exponent and silently selling the vulnerability
                                # curve. Shifted by -0.2*log10(2.9782) so stem_P50 is
                                # unchanged at the default.
                                B_Hv1 = 0.36591565341924093,
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

    ## Reuse the legacy hk_s parameterisation to derive TF24_cost_scale (the TF24
    ## hydraulic cost scale, named `g1_TF24` before #634):
    TF24_cost_scale <- B_hks1 * (rho / rho_0) ^ (-B_hks2)

    ## rho / wood turnover relationship:
    k_s  <- B_ks1 *  (rho / rho_0) ^ (-B_ks2)

    ## TODO: Convert the p50 ks function back to a mean centred function using K_s_0

    ## stem_P50 sapwood specific conductivity turnover:
    if (any(K_s <= 0, na.rm = TRUE)) {
      stop("K_s must be > 0 for stem_P50 derivation", call. = FALSE)
    }
    stem_P50 <- 10^(B_Hv1 + B_Hv2 * log10(K_s))

    ## stem_P50 shape parameter trade off
    stem_c <- B_c1 * exp(-B_c2 * stem_P50)
    ## scale parameter stem_b of the vulnerability curve exp(-(psi/b)^c): the water
    ## potential at 1/e (~37%) conductivity remaining, solved here from the 50%
    ## loss-of-conductivity point stem_P50 [MPa]. NOTE phylloptim derives this and
    ## psi_crit from (stem_P50, stem_c) itself, by these same formulas -- the two
    ## below are reported for the record and are not inputs to the leaf (#634).
    stem_b <- stem_P50/((-log(1-50/100))^(1/stem_c))

    ## water potential at critical xylem failure (95%) [MPa]:
    psi_crit <- stem_b*(log(1/0.05))^(1/stem_c)

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

    extra <- cbind(k_l,                              # lma
                   d_I, TF24_cost_scale, k_s, r_s, r_b, # rho
                   a_f3,                              # omega
                   r_l, nmass_l,                      # lma, narea
                   stem_c, stem_P50, stem_b, psi_crit) # K_s

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

##' Every parameter of a model, grouped by how a user reaches it
##'
##' One table of everything a model can be told, so the options are visible in
##' one place (\code{#698}). Rows are grouped:
##' \describe{
##'   \item{\code{trait}}{what the trade-off function reads from a trait matrix:
##'     the ecological traits. A trait that is not itself a field of
##'     \code{pars} (FF16's \code{narea}) has no default here.}
##'   \item{\code{derived}}{what the trade-off function derives from the
##'     traits. Naming one of these in a trait matrix is refused, because the
##'     trade-off would overwrite it; set the trait it comes from, or build the
##'     strategies with \code{hyperpar = NULL}.}
##'   \item{\code{parameter}}{the remaining fields of \code{pars}: settable as
##'     traits, with no trade-off attached to them.}
##'   \item{\code{strategy}}{numeric fields a strategy carries at its top level,
##'     beside \code{pars} (TF24t's thermal parameters, TF24f's acclimation
##'     settings). A trait matrix reaches these too.}
##'   \item{\code{control}}{\code{\link{Control}}'s numerical settings, which
##'     govern how the model is solved rather than what is solved.}
##' }
##' Which rows are \code{trait} and \code{derived} depends on the trade-off
##' function: with none attached every field of \code{pars} is a
##' \code{parameter}. The coefficients of the trade-offs themselves (slopes and
##' intercepts) are not model parameters and are not listed; they are documented
##' with the trade-off functions (\pkg{fieldtrial}'s coefficient tables, or
##' \code{\link{frozen_TF24_hyperpar}} for plant's frozen copy).
##'
##' @param x A model name (\code{"FF16"}, \code{"TF24"}, ...), or a
##'   \code{Parameters} object, whose own trade-offs are then read.
##' @param hyperpar The trade-off function to classify by. For a model name it
##'   defaults to plant's frozen testing trade-offs
##'   (\code{\link{default_hyperpar}}); for a \code{Parameters} object to
##'   \code{param_hyperpar(x)}. Pass \code{NULL} explicitly for no trade-offs
##'   (every \code{pars} field becomes a \code{parameter}).
##' @return A tibble with one row per parameter: \code{name}, \code{group},
##'   \code{default} (the strategy's or control's default, as text; read the
##'   exact double from \code{Strategy()$pars} or \code{Control()}),
##'   \code{description}, \code{symbol} and \code{unit} where
##'   \code{inst/docs} documents the name, \code{NA} otherwise.
##' @seealso \code{\link{default_parameters}}, \code{\link{generate_strategy}},
##'   \code{\link{trait_matrix}}
##' @examples
##' tab <- parameter_table("TF24")
##' table(tab$group)
##' tab[tab$group == "derived", "name"]
##' @export
parameter_table <- function(x, hyperpar = NULL) {
  explicit_hyperpar <- !missing(hyperpar)
  if (inherits(x, "Parameters")) {
    p <- x
    model <- extract_RcppR6_template_types(p, "Parameters")[[1]][[1]]
    if (!explicit_hyperpar) {
      hyperpar <- param_hyperpar(p)
    }
  } else if (is.character(x) && length(x) == 1L) {
    model <- x
    p <- scm_base_parameters(model)
    if (!explicit_hyperpar) {
      hyperpar <- default_hyperpar(model)
    }
  } else {
    stop("`x` must be a model name or a Parameters object", call. = FALSE)
  }
  if (is.null(hyperpar)) {
    hyperpar <- identity_hyperpar
  }
  if (!is.function(hyperpar)) {
    stop("`hyperpar` must be a trade-off function f(m, s, filter = TRUE), ",
         "or NULL", call. = FALSE)
  }
  s <- p$strategy_default
  pars <- s$pars
  ctrl <- s$control

  derived <- colnames(hyperpar(matrix(numeric(0), 1L, 0L), s, filter = FALSE))
  if (is.null(derived)) {
    derived <- character(0)
  }
  traits <- union(attr(hyperpar, "traits", exact = TRUE),
                  attr(hyperpar, "input_traits", exact = TRUE))
  if (is.null(traits)) {
    traits <- character(0)
  }
  fields <- strategy_field_names(s)

  as_text <- function(v) {
    if (is.null(v) || length(v) == 0L) return(NA_character_)
    if (is.character(v)) return(paste(v, collapse = ", "))
    paste(as.character(v), collapse = ", ")
  }
  row <- function(names, group, values) {
    data.frame(name = names, group = rep_len(group, length(names)),
               default = vapply(values, as_text, character(1)),
               stringsAsFactors = FALSE)
  }

  par_names <- names(pars)
  trait_in_pars <- intersect(traits, par_names)
  trait_consumed <- setdiff(traits, par_names)
  derived_in_pars <- intersect(derived, par_names)
  rest <- setdiff(par_names, c(trait_in_pars, derived_in_pars))

  tab <- rbind(
    row(trait_in_pars, "trait", pars[trait_in_pars]),
    row(trait_consumed, "trait", rep(list(NULL), length(trait_consumed))),
    row(derived_in_pars, "derived", pars[derived_in_pars]),
    row(rest, "parameter", pars[rest]),
    row(fields, "strategy", lapply(fields, function(n) s[[n]])),
    row(names(ctrl), "control", unclass(ctrl)[names(ctrl)]))
  rownames(tab) <- NULL

  docs <- parameter_docs()
  at <- match(tab$name, docs$Code)
  tab$description <- docs$Description[at]
  tab$symbol <- docs$Symbol[at]
  tab$unit <- docs$Unit[at]
  attr(tab, "model") <- model
  tibble::as_tibble(tab)
}

## The parameter descriptions plant documents, by code name. FF16's table
## covers most of the shared economic and allometric names; a name it does not
## list gets NA.
parameter_docs <- function() {
  path <- system.file("docs", "FF16", "FF16-params-core.csv", package = "plant",
                      mustWork = TRUE)
  docs <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  docs <- docs[!is.na(docs$Code) & nzchar(docs$Code), , drop = FALSE]
  docs$Unit[!nzchar(docs$Unit)] <- NA_character_
  docs$Symbol[!nzchar(docs$Symbol)] <- NA_character_
  docs
}

## The solar model the frozen FF16 trade-off integrates its assimilation day
## with. Internal: the recommended driver code lives in fieldtrial, with the
## declination corrected (this one uses sin(delta) as delta, about 0.6 degrees
## low at the solstices); correcting it here would move a_p1 and a_p2 and so
## every FF16 golden, which is why the frozen copy keeps it as it was.

# Returns incoming PAR on surface normal to incoming light, for given solar
# angle (in radians) S_atmos = radiation outside atmosphere W / m2/s tau:
# transmisivity of shortest atmospheric path PAR_frac Fraction of light that is
# Photosynthetically Active return value in mol PAR / m2/ day
PAR_given_solar_angle <- function(solar_angle, S_atmos = 1360, tau = 0.8, PAR_frac = 0.5) {
  
  # Adjust intensity for relative path-length through atmosphere
  
  # Would adjust for altitude here, but this function disabled
  altitude_adjustment <- 1
  
  M <- 1/cos(pi/2 - solar_angle) * altitude_adjustment
  
  mol_per_watt <- 4.6/1e+06
  sec_per_day <- 3600 * 24
  # Watts /m2/s
  ret <- (S_atmos * PAR_frac) * tau^M * mol_per_watt * sec_per_day
  ret[solar_angle <= 0] <- 0
  
  ret
}

# Returns solar angle = angle between the horizontal and the line connecting to
# the sun for given latitude and time of day, in radians
solar_angle <- function(decimal_day_time, latitude) {
  
  day <- floor(decimal_day_time)
  time <- (decimal_day_time - day) * 24
  
  radians <- function(x) {
    x/180 * pi
  }
  lat_radians <- radians(latitude)
  
  # solar declination (radians) - angle b/w the earth-sun line and the equatorial
  # plane
  delta <- 0.39785 * sin(4.869 + day/365 * 2 * pi + 0.03345 * sin(6.224 + day/365 * 
    2 * pi))
  
  # Hour Angle (degrees) - angular distance that earth has rotated in a day
  HrAng <- radians(180 - time * 15)
  
  # Sin of solar angle & Solar Angle - angle between the horizontal and the line
  # connecting to the sun
  SinA <- sin(lat_radians) * sin(delta) + cos(lat_radians) * cos(delta) * cos(HrAng)
  asin(SinA)
}
