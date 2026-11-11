##' @title Create a TF24f Plant or Node
##' @param s A \code{\link{TF24f_Model}} object
##' @export
##' @rdname TF24f_Individual
##' @examples
##' pl <- TF24f_Individual()
##' pl$height
TF24f_Individual <- function(s=TF24f_Model()) {
  Individual("TF24f", "TF24_Env")(s)
}

##' @title Set up a TF24f system with default or specified parameters
##' @description Set up a model system with default or specified parameters.
##' @inheritParams FF16_Parameters
##' @export
##' @rdname TF24f_Parameters
##' @examples
##' p1 <- TF24f_Parameters()
##' p2 <- TF24f_Parameters(max_patch_lifetime = 10.0)
TF24f_Parameters <- function(...) {
  Parameters("TF24f","TF24_Env")(...)
}

#' @export
#' @importFrom rlang .data
#' @rdname expand_state
TF24f_expand_state <- function(results) {
  expand_state_by_allometry(results, TF24f_strategy_expand_allometry)
}
