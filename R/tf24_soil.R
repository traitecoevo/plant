## Soil layer geometry for TF24.
##
## TF24_Environment's geometry setters are C++ methods reached one at a time, and
## they have an ordering constraint between them: changing the layers throws the
## water state away, so the state must be set last. That constraint is real and is
## not going away, but it does not belong in every caller's head -- so it is solved
## once, here, in set_tf24_soil().

##' Set a TF24 soil profile in one call
##'
##' Sets the layer geometry, then the per-layer parameters, then the water state,
##' in that order -- which is the only order that works, and the reason this
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
##' @param widths Layer widths, m, top down -- e.g. a profile grading from a
##'   thin surface layer, as \pkg{fieldtrial}'s \code{soil_widths_graded()}
##'   builds. Must not get thinner with depth (equal is fine). A single number
##'   sets one layer of that width. If \code{NULL} the geometry is left alone.
##' @param theta Initial volumetric soil moisture. A single number is used for
##'   every layer; a vector must have one entry per layer. If \code{NULL} the
##'   layers start at half saturation, matching the constructor.
##' @param soil_moist_sat,K_sat,a_psi,n_psi Per-layer soil parameters, each either
##'   \code{NULL} (use the scalar field for every layer) or a vector with one
##'   entry per layer. If all four are \code{NULL} the layered-parameter path is
##'   not engaged at all. \code{a_psi} is the retention-curve scale in
##'   MPa.
##' @return \code{env}, modified in place and returned invisibly.
##' @export
##' @examples
##' env <- Environment("TF24")
##' set_tf24_soil(env, c(0.02, 0.06, 0.18, 0.44, 0.80), theta = 0.2)
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
