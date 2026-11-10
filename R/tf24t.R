##' @title Create a TF24t Plant or Node
##' @description TF24t is TF24 with phylloptim's lasting leaf thermal-damage and
##' acclimation layer switched on (#566). It adds four ODE states after TF24's:
##' \code{acclim_thermostab}, \code{impaired},
##' \code{damaged} and \code{mean_leaf_temp}. Its thermal parameters
##' sit at the top level of \code{\link{TF24t_Strategy}}, with every rate
##' constant per day. The Penman-Monteith energy balance is always on, so the
##' \code{leaf_temp} driver is read as air temperature.
##' @param s A \code{\link{TF24t_Strategy}} object
##' @export
##' @rdname TF24t_Individual
##' @examples
##' pl <- TF24t_Individual()
##' pl$height
TF24t_Individual <- function(s=TF24t_Strategy()) {
  Individual("TF24t", "TF24_Env")(s)
}

##' @title Set up a TF24t system with default or specified parameters
##' @description Set up a model system with default or specified parameters.
##' @inheritParams FF16_Parameters
##' @export
##' @rdname TF24t_Parameters
##' @examples
##' p1 <- TF24t_Parameters()
##' p2 <- TF24t_Parameters(max_patch_lifetime = 10.0)
TF24t_Parameters <- function(...) {
  Parameters("TF24t","TF24_Env")(...)
}

#' @export
#' @importFrom rlang .data
#' @rdname expand_state
TF24t_expand_state <- function(results) {
  expand_state_by_allometry(results, TF24t_strategy_expand_allometry)
}
