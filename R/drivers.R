##' Extrinsic drivers are interpolated with a cubic spline, which is a poor fit
##' for intermittent forcing such as daily rainfall. Between mostly-zero control
##' points a cubic overshoots in both directions: a realistic daily rainfall
##' series with a ~10\% wet-day fraction evaluates *negative* at roughly 45\% of
##' points, and overshoots pulse peaks several-fold. Negative rainfall is
##' floored at zero inside the TF24 soil balance, but the distortion of pulse
##' magnitudes remains, so it is worth inspecting a series before committing to
##' a long run.
##'
##' This runs that inspection: it evaluates the driver on a fine grid and
##' compares the result against the control points you supplied.
##'
##' @title Check how a driver's control points survive interpolation
##' @param environment An \code{Environment} object with the driver already set
##'   (via \code{environment$extrinsic_drivers_set_variable(driver_name, x, y)}).
##' @param driver_name Name of the driver to check, e.g. \code{"rainfall"}.
##' @param x,y The control points that were supplied for this driver. Used to
##'   compute the intended integral and value range to compare against.
##' @param n_eval Number of points on the fine evaluation grid.
##' @param warn If \code{TRUE} (default), raise a warning when the interpolated
##'   series goes negative or its integral departs from the supplied series by
##'   more than \code{tol_integral}.
##' @param tol_integral Relative tolerance on the integral mismatch before
##'   warning.
##'
##' @return Invisibly, a list with the evaluated \code{min}, \code{max},
##'   \code{frac_negative}, the \code{integral_supplied} and
##'   \code{integral_evaluated} (trapezoid), their \code{integral_rel_error},
##'   and \code{negative_area} — the amount spuriously removed by undershoot.
##'   Printed as a short summary as a side effect.
##'
##' @examples
##' e <- TF24_Environment()
##' t <- seq(0, 2, length.out = 730)
##' # intermittent daily rainfall: wet season only, most days dry
##' set.seed(1)
##' y <- ifelse((t %% 1) > 0.5 & runif(length(t)) < 0.3, 5, 0)
##' e$extrinsic_drivers_set_variable("rainfall", t, y)
##' check_driver_interpolation(e, "rainfall", t, y)
##' @export
check_driver_interpolation <- function(environment, driver_name, x, y,
                                       n_eval = 20001, warn = TRUE,
                                       tol_integral = 0.01) {
  if (length(x) != length(y)) {
    stop("`x` and `y` must be the same length.")
  }
  if (length(x) < 2L) {
    stop("Need at least two control points to check an interpolation.")
  }

  grid <- seq(min(x), max(x), length.out = n_eval)
  evaluated <- environment$extrinsic_drivers_evaluate_range(driver_name, grid)

  trapz <- function(xs, ys) sum(diff(xs) * (utils::head(ys, -1) + utils::tail(ys, -1)) / 2)

  integral_supplied  <- trapz(x, y)
  integral_evaluated <- trapz(grid, evaluated)
  # Area under the negative excursions: water the driver removes rather than adds.
  negative_area <- -trapz(grid, pmin(evaluated, 0))

  rel_error <- if (integral_supplied != 0) {
    (integral_evaluated - integral_supplied) / abs(integral_supplied)
  } else {
    NA_real_
  }

  out <- list(
    driver             = driver_name,
    min                = min(evaluated),
    max                = max(evaluated),
    supplied_min       = min(y),
    supplied_max       = max(y),
    frac_negative      = mean(evaluated < 0),
    negative_area      = negative_area,
    integral_supplied  = integral_supplied,
    integral_evaluated = integral_evaluated,
    integral_rel_error = rel_error
  )

  message(sprintf(
    paste0("driver '%s': evaluated range [%.4g, %.4g] vs supplied [%.4g, %.4g]\n",
           "  negative at %.1f%% of points (area %.4g)\n",
           "  integral %.6g vs supplied %.6g (%+.2f%%)"),
    driver_name, out$min, out$max, out$supplied_min, out$supplied_max,
    100 * out$frac_negative, out$negative_area,
    out$integral_evaluated, out$integral_supplied, 100 * rel_error))

  if (warn) {
    if (out$min < 0 && out$supplied_min >= 0) {
      warning(sprintf(
        paste0("Driver '%s' interpolates negative (min %.4g) although every ",
               "supplied value is non-negative: the cubic spline is ",
               "undershooting. Consider supplying smoother control points or ",
               "aggregating to a coarser interval."),
        driver_name, out$min), call. = FALSE)
    }
    if (!is.na(rel_error) && abs(rel_error) > tol_integral) {
      warning(sprintf(
        paste0("Driver '%s' interpolated integral departs from the supplied ",
               "control points by %+.1f%% (tolerance %.1f%%)."),
        driver_name, 100 * rel_error, 100 * tol_integral), call. = FALSE)
    }
  }

  invisible(out)
}
