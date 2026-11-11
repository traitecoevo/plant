##' @title Create a TF24w Plant or Node
##' @param s A \code{\link{TF24w_Model}} object
##' @export
##' @rdname TF24w_Individual
##' @examples
##' pl <- TF24w_Individual()
##' pl$height
TF24w_Individual <- function(s=TF24w_Model()) {
  Individual("TF24w", "TF24_Env")(s)
}

##' @title Set up a TF24w system with default or specified parameters
##' @description Set up a model system with default or specified parameters.
##' @inheritParams FF16_Parameters
##' @export
##' @rdname TF24w_Parameters
##' @examples
##' p1 <- TF24w_Parameters()
##' p2 <- TF24w_Parameters(max_patch_lifetime = 10.0)
TF24w_Parameters <- function(...) {
  Parameters("TF24w","TF24_Env")(...)
}

#' @export
#' @importFrom rlang .data
#' @rdname expand_state
TF24w_expand_state <- function(results) {
  expand_state_by_allometry(results, TF24w_strategy_expand_allometry)
}
