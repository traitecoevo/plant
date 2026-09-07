## Model-specific event constructors.
##
## The event machinery in events.R is deliberately taxa- and model-agnostic: it
## knows about resource pools, not about water. Names that are only true of a
## particular model belong here, as thin wrappers that read correctly at their
## own call sites. Nothing in the generic layer depends on this file.

##' A rainfall pulse for TF24.
##'
##' TF24's resources are its soil layers, so a rainfall pulse is a
##' \code{\link{resource_pulse}} of water into one of them. This is the
##' motivating case for the event mechanism (#522): TF24 otherwise takes
##' rainfall as a smooth cubic spline, which smears out exactly the variability
##' that matters — large pulses.
##'
##' The surface layer accepts what its free capacity allows and sheds the rest,
##' which is not a corner case: a realistic dryland event (~13 mm) already
##' exceeds a moderately wet layer's capacity. \code{scm$event_log} reports how
##' much of each pulse landed and how much was shed.
##'
##' @param time Event time(s), in years of patch age.
##' @param depth Depth of water delivered, in m — so 13 mm is \code{0.013}.
##' @param layer Soil layer to deliver it to. The default, 1, is the surface,
##'   i.e. rain. A deeper layer is irrigation, or a water table rising into the
##'   profile.
##' @return Event rows, to be combined with \code{\link{events}}.
##' @seealso \code{\link{events}}, \code{\link{resource_pulse}}
##' @export
##' @examples
##' p <- scm_base_parameters("TF24")
##' p <- add_strategies(p, trait_matrix(1, "lma"))
##' ev <- events(
##'   events_default(p),
##'   rainfall_pulse(time = c(1.5, 3.2), depth = c(0.013, 0.050))
##' )
##' @param env Optionally, the \code{Environment} the run will use. When given,
##'   \code{rainfall_pulse()} warns if that environment still carries a non-zero
##'   continuous \code{rainfall} driver — pulses \emph{add to} the continuous
##'   forcing rather than replacing it, so supplying both silently double-counts
##'   the water. Set the driver to zero to model rainfall as pulses alone.
##' @rdname rainfall_pulse
rainfall_pulse <- function(time, depth, layer = 1, env = NULL) {
  if (!is.null(env)) warn_if_continuous_rainfall(env)
  ## Builds the rows directly rather than calling resource_pulse(), so that a
  ## length mismatch is reported against `depth` -- the argument the caller
  ## actually typed -- instead of the generic `amount`.
  event_rows("resource_pulse", time = time, target = "environment",
             target_index = layer, params = list(depth = depth))
}

## Pulses are added on top of whatever the `rainfall` driver is doing; nothing
## disables it. Silently double-counting the water is the easy mistake here, so
## say so rather than leaving it to be discovered in a water budget.
warn_if_continuous_rainfall <- function(env, times = seq(0, 105.32, length.out = 64)) {
  drivers <- tryCatch(env$extrinsic_drivers_get_names(),
                      error = function(e) character(0))
  if (!"rainfall" %in% drivers) return(invisible(FALSE))
  vals <- tryCatch(env$extrinsic_drivers_evaluate_range("rainfall", times),
                   error = function(e) 0)
  if (any(abs(vals) > 0)) {
    warning("This environment still has a non-zero continuous `rainfall` ",
            "driver. Rainfall pulses are added on top of it, not instead of ",
            "it, so the run will receive both. Set the driver to zero with ",
            "env$extrinsic_drivers_set_constant(\"rainfall\", 0) to model ",
            "rainfall as pulses alone.", call. = FALSE)
    return(invisible(TRUE))
  }
  invisible(FALSE)
}
