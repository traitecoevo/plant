##' @title Create a TF24 Plant or Node
##' @param s A \code{\link{TF24_Model}} object
##' @export
##' @rdname TF24_Individual
##' @examples
##' pl <- TF24_Individual()
##' pl$height
TF24_Individual <- function(s=TF24_Model()) {
  Individual("TF24", "TF24_Env")(s)
}

##' @title Setup an a TF24 system with  default or specified parameters
##' @description Setup an a model system with default or specified parameters.
##' @inheritParams FF16_Parameters
##' @export
##' @rdname TF24_Parameters
##' @examples
##' p1 <- TF24_Parameters()
##' p2 <- TF24_Parameters(max_patch_lifetime = 10.0)
TF24_Parameters <- function(...) {
  Parameters("TF24","TF24_Env")(...)
}

##' Height at which the stem hydraulic path is anchored [m].
##'
##' TF24's resistance depends on `theta_tip` and `K_s` only through their product, so
##' there is exactly one free scalar and the height-dependent model can be made
##' to agree with the older height-linear one at exactly **one** height -- never
##' two. `TF24_H_ANCHOR` is that height. It is the pivot of a rotation: below it
##' every plant is more resistant than under the linear model, above it every
##' plant is less.
##'
##' Anchored at 1 m, so resistance is unchanged at 1 m stature and falls
##' progressively above it -- reaching 0.21x by 60 m. Only sub-metre plants pay
##' more than they did. It is deliberately a fixed constant rather than being
##' read from `pars$hmat`: sweeping `hmat` as a trait must not silently move the
##' hydraulics with it.
##'
##' This is provenance for a default, not a model parameter, which is why it does
##' not live in `TF24_Pars`. Exposing it would add a fourth knob to the grouping
##' `L_tip^beta / (K_s * theta_tip)`, which is already identifiable only as a whole --
##' and inviting independent calibration of its members is exactly what
##' invariance criterion I5 forbids.
##' @export
TF24_H_ANCHOR <- 1

#' @export
#' @importFrom rlang .data
#' @rdname expand_state
TF24_expand_state <- function(results) {
  expand_state_by_allometry(results, TF24_strategy_expand_allometry)
}
