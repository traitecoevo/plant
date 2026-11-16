##' Create a K93 Individual or Node
##' @title Create a K93 Individual or Node
##' @param s A \code{\link{K93_Model}} object
##' @export
##' @rdname K93_Individual
##' @examples
##' pl <- K93_Individual()
##' pl$height
K93_Individual <- function(s=K93_Model()) {
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
  ## K93's size is stem diameter, so it has no heartwood or leaf columns: its
  ## expansion is the stem's diameter and basal area, in m and m^2.
  data <- split(results$species, results$species$species)
  for (i in seq_len(results$n_spp)) {
    d <- data[[i]]
    stem <- K93_strategy_expand_allometry(results$p$strategies[[i]], d$height)
    data[[i]] <- dplyr::bind_cols(d, stem)
  }
  results$species <- dplyr::bind_rows(data)
  results
}
