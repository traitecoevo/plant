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
##' TF24's resistance depends on `theta` and `K_s` only through their ratio, so
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
##' `theta * L_tip^beta / K_s`, which is already identifiable only as a whole --
##' and inviting independent calibration of its members is exactly what
##' invariance criterion I5 forbids.
##' @export
TF24_H_ANCHOR <- 1

##' Convert a whole-stem sapwood-specific conductivity to a terminal-segment one.
##'
##' `K_s` is a property of the terminal segment, but conductivity is usually
##' *reported* for a stem section at some height -- and it was defined that way
##' in TF24 before the stem path integral existed. This converts such a value
##' into the terminal-segment `K_s` that reproduces it at `H_anchor`.
##'
##' The ratio is obtained by calling the same `effective_path_length` the model
##' itself uses, at `beta` and at zero, rather than being written as a literal:
##' that way the `eta_c` convention (the path runs to the leaf-area-weighted mean
##' leaf height, not the apex) and the thin-crown approximation cancel out of the
##' ratio instead of being smuggled into a hand-computed constant.
##'
##' At the shipped defaults the ratio is 2.9782, so a whole-stem 1 becomes
##' 0.3357738.
##'
##' @param K_s Whole-stem sapwood-specific conductivity.
##' @param D_c Conduit widening exponent.
##' @param theta_c Huber-profile exponent.
##' @param L_tip Terminal segment length [m].
##' @param H_anchor Height at which the whole-stem value was defined [m].
##' @param eta Canopy shape parameter, used for the crown-centre coordinate.
##' @return The terminal-segment `K_s`.
##' @export
TF24_K_s_from_whole_stem <- function(K_s, D_c = 0.2, theta_c = 0,
                                    L_tip = 0.02,
                                    H_anchor = TF24_H_ANCHOR, eta = 12) {
  eta_c <- 1 - 2 / (1 + eta) + 1 / (1 + 2 * eta)
  L_top <- H_anchor * eta_c
  beta <- 2 * D_c + theta_c
  # R_L = (theta/K_s) * L_eff, so matching R_L at H_anchor across the two models
  # means K_s_tip / K_s = L_eff(beta) / L_eff(0). Widening SHORTENS the effective
  # path, so K_s must come DOWN to compensate -- the tip of a widening stem is
  # narrower, hence less conductive, than the uniform stem it replaces.
  K_s * test_stem_effective_path_length(L_top, L_tip, beta) /
    test_stem_effective_path_length(L_top, L_tip, 0)
}

#' @export
#' @importFrom rlang .data
#' @rdname expand_state
TF24_expand_state <- function(results) {
  expand_state_by_allometry(results, TF24_strategy_expand_allometry)
}
