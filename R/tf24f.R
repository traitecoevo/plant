##' @title Create a TF24f Plant or Node
##' @param s A \code{\link{TF24f_Strategy}} object
##' @export
##' @rdname TF24f_Individual
##' @examples
##' pl <- TF24f_Individual()
##' pl$height
TF24f_Individual <- function(s=TF24f_Strategy()) {
  Individual("TF24f", "TF24_Env")(s)
}

##' @title Set up a TF24f system with default or specified parameters
##' @description Set up a model system with default or specified parameters.
##' @param ... Arguments to be passed to the model constructor. These include
##'
##'   *`patch_area`: Area of an individual patch. Only relevant for the stochastic model. Default is 1.0m2.
##'   *`max_patch_lifetime`: The maximum time in years we want to simulate
##'   *`strategies`: A list of strategies to simulate. The default is an empty list.
##'   *`strategy_default`: Values for the default strategy. The default values are those specified in the C++ code for the model.
##'   *`node_schedule_times_default`: Default vector of times at which to introduce nodes. The default is chosen to have close spacing at the start of the simulation.
##'   *`node_schedule_times`: A list with each element containing the vector of times we want to introduce nodes for each strategy. The default is an empty list.
##'   *`ode_times`: A vector of patch ages we want the ode solver to stop at
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
  data <- split(results$species, results$species$species)

  for (i in seq_len(results$n_spp)) {
    s <- results$p$strategies[[i]]
    d <- data[[i]]
    # Derived size/mass columns are computed by the strategy's own C++
    # allometry functions (see src/strategy_expand.cpp) rather than being
    # re-derived here, so the formulas live in exactly one place.
    allom <- TF24f_strategy_expand_allometry(s, d$height, d$area_heartwood,
                                            d$mass_heartwood)
    data[[i]] <- dplyr::bind_cols(d, allom)
  }

  results$species <- data %>% dplyr::bind_rows()

  results
}
