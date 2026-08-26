## Soil layer geometry for TF24.
##
## TF24_Environment's geometry setters are C++ methods reached one at a time, and
## they have an ordering constraint between them: changing the layers throws the
## water state away, so the state must be set last. That constraint is real and is
## not going away, but it does not belong in every caller's head -- so it is solved
## once, here, in set_tf24_soil().

##' Soil layer widths grading from a thin surface layer
##'
##' Builds a vector of layer widths for
##' \code{\link{TF24_Environment}}$\code{set_soil_layer_widths()}: \code{n} layers
##' spanning \code{depth}, the first of width \code{top} and each subsequent one a
##' fixed factor thicker. The motivating case is soil evaporation, which needs a
##' very thin layer at the top of the profile (#626) — a 30 cm surface layer
##' averages away the wetting and drying that drives it.
##'
##' The ratio is solved so the widths sum to \code{depth} exactly; the last width
##' is then taken as the remainder rather than computed, so rounding cannot leave
##' the column a hair shallower or deeper than asked for.
##'
##' @section What a thin layer costs:
##' Two things, neither of which \code{plant} corrects for you.
##'
##' **The column gets stiff.** The soil water balance's Jacobian diagonal is
##' \code{K'(theta)/dz}, so it scales as \code{1/dz}: measured 1.8e4 yr^-1 at the
##' 0.3 m default and 2.6e5 at 2 cm, where the characteristic time is about two
##' minutes. It only bites when the layer is wet — just after rain, which is
##' exactly when an evaporation layer matters. Below about 2.6 mm the stable step
##' falls under \code{Control()}'s \code{ode_step_size_min} and the step controller
##' cannot satisfy stability at all.
##'
##' **Runoff is recalibrated.** The infiltration term keys off the surface layer
##' alone, through \code{(theta_0/theta_sat)^b_infil} with \code{b_infil = 8}. A
##' thin layer 0 fills far faster and so holds that switch shut for longer:
##' replacing the uniform profile with a 2 cm top layer moved annual runoff from
##' 13.4\% to 20.4\% under an event-driven rainfall series, with no parameter
##' changed. \code{a_infil} and \code{b_infil} were calibrated against a 30 cm
##' surface layer.
##'
##' @param depth Total column depth, m.
##' @param n Number of layers.
##' @param top Width of the surface layer, m. Must be less than
##'   \code{depth / n}, since a graded profile thickens with depth.
##' @return A vector of \code{n} layer widths, summing to \code{depth}.
##' @seealso \code{\link{set_tf24_soil}}
##' @export
##' @examples
##' soil_widths_graded(depth = 1.5, n = 5, top = 0.02)
##' sum(soil_widths_graded(depth = 1.5, n = 5, top = 0.02)) # 1.5
soil_widths_graded <- function(depth = 1.5, n = 5, top = 0.02) {
  if (length(n) != 1L || n < 1 || n != round(n)) {
    stop("n must be a single positive whole number")
  }
  if (length(depth) != 1L || !is.finite(depth) || depth <= 0) {
    stop("depth must be a single positive number")
  }
  if (length(top) != 1L || !is.finite(top) || top <= 0) {
    stop("top must be a single positive number")
  }
  n <- as.integer(n)
  if (n == 1L) {
    return(depth)
  }
  ## top >= depth/n would need a ratio <= 1, i.e. a profile thinning with depth.
  ## Refused rather than served, because it is the sign of a transposed argument
  ## far more often than a deliberate choice -- and an inverted profile puts the
  ## stiffest layer where nothing needs the resolution.
  if (top >= depth / n) {
    stop("top (", format(top), ") must be less than depth / n (",
         format(depth / n), "): a graded profile thickens with depth")
  }

  ## Widths are top * r^(0..n-1), so sum = top * (r^n - 1)/(r - 1) = depth. No
  ## closed form for r, so solve it; the sum is strictly increasing in r, and
  ## r = 1 gives n*top < depth, so the root is bracketed above 1.
  total <- function(r) top * (r^n - 1) / (r - 1)
  upper <- 2
  while (total(upper) < depth) {
    upper <- upper * 2
    if (upper > 1e6) {
      stop("could not find a grading ratio for depth = ", format(depth),
           ", n = ", n, ", top = ", format(top))
    }
  }
  r <- stats::uniroot(function(r) total(r) - depth,
                      c(1 + 1e-12, upper), tol = .Machine$double.eps^0.75)$root

  out <- top * r^(seq_len(n) - 1L)
  ## The last width absorbs the residual, so sum(out) == depth to the bit. Without
  ## this the column comes out a rounding error off the requested depth, and
  ## `depth` is read by rooting_depth_max comparisons downstream.
  out[[n]] <- depth - sum(out[-n])
  if (out[[n]] <= 0) {
    stop("grading gave a non-positive final layer; try a larger `top` or ",
         "smaller `n`")
  }
  out
}

##' Set a TF24 soil profile in one call
##'
##' Sets the layer geometry, then the per-layer parameters, then the water state,
##' in that order — which is the only order that works, and the reason this
##' function exists.
##'
##' Changing the number of layers necessarily discards the water state and any
##' per-layer parameters, because both are per layer and there is no defensible
##' way to re-map them onto a different column. So a caller doing this by hand has
##' to know that \code{set_soil_layer_widths()} comes before
##' \code{set_soil_parameters()}, which comes before
##' \code{set_soil_water_state()}. Getting it wrong is silent: the profile is
##' simply the previous one, and both are valid geometries.
##'
##' @param env A \code{\link{TF24_Environment}}.
##' @param widths Layer widths, m, top down — e.g. from
##'   \code{\link{soil_widths_graded}}. A single number sets one layer of that
##'   width. If \code{NULL} the geometry is left alone.
##' @param theta Initial volumetric soil moisture. A single number is used for
##'   every layer; a vector must have one entry per layer. If \code{NULL} the
##'   layers start at half saturation, matching the constructor.
##' @param soil_moist_sat,K_sat,a_psi,n_psi Per-layer soil parameters, each either
##'   \code{NULL} (use the scalar field for every layer) or a vector with one
##'   entry per layer. If all four are \code{NULL} the layered-parameter path is
##'   not engaged at all.
##' @return \code{env}, modified in place and returned invisibly.
##' @seealso \code{\link{soil_widths_graded}}
##' @export
##' @examples
##' env <- Environment("TF24")
##' set_tf24_soil(env, soil_widths_graded(1.5, 5, top = 0.02), theta = 0.2)
##' env$get_soil_layer_widths()
set_tf24_soil <- function(env, widths = NULL, theta = NULL,
                          soil_moist_sat = NULL, K_sat = NULL,
                          a_psi = NULL, n_psi = NULL) {
  if (!is.null(widths)) {
    env$set_soil_layer_widths(as.numeric(widths))
  }
  n <- env$get_soil_number_of_depths()

  pars <- list(soil_moist_sat = soil_moist_sat, K_sat = K_sat,
               a_psi = a_psi, n_psi = n_psi)
  if (any(!vapply(pars, is.null, logical(1)))) {
    ## Passing the current count means set_soil_parameters() takes its same-n
    ## branch and leaves the geometry -- and so the widths just set -- alone.
    env$set_soil_parameters(n, soil_moist_sat, K_sat, a_psi, n_psi)
  }

  if (is.null(theta)) {
    theta <- rep(env$soil_moist_sat * 0.5, n)
  } else if (length(theta) == 1L) {
    theta <- rep(theta, n)
  }
  if (length(theta) != n) {
    stop("theta must be length 1 or ", n, " (the number of layers), not ",
         length(theta))
  }
  env$set_soil_water_state(theta)

  invisible(env)
}
