##' Create a K93 Individual or Node
##' @title Create a K93 Individual or Node
##' @param s A \code{\link{K93_Strategy}} object
##' @export
##' @rdname K93_Individual
##' @examples
##' pl <- K93_Individual()
##' pl$height
K93_Individual <- function(s=K93_Strategy()) {
  Individual("K93", "K93_Env")(s)
}

##' @export
##' @rdname FF16_Parameters
K93_Parameters <- function() {
  Parameters("K93","K93_Env")()
}


#' @export
#' @importFrom rlang .data
#' @rdname expand_state
K93_expand_state <- function(results) {

  # currently not doing anything
  
  results
}
