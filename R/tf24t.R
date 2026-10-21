##' @title Create a TF24t Plant or Node
##' @description TF24t is TF24 with phylloptim's lasting leaf thermal-damage and
##' acclimation layer switched on (#566). It adds four ODE states after TF24's:
##' \code{acclim_thermostab}, \code{damage_recoverable},
##' \code{damage_permanent} and \code{mean_leaf_temp}. Its thermal parameters
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
##' @rdname TF24t_Parameters
##' @examples
##' p1 <- TF24t_Parameters()
##' p2 <- TF24t_Parameters(max_patch_lifetime = 10.0)
TF24t_Parameters <- function(...) {
  Parameters("TF24t","TF24_Env")(...)
}

##' @export
##' @rdname make_TF24_hyperpar
make_TF24t_hyperpar <- function(...) {
  ## TF24t_Strategy inherits TF24_Strategy and reuses TF24_Pars, so its
  ## hyperparameterisation is TF24's. The thermal parameters are top-level
  ## strategy fields, which a trait matrix does not reach.
  make_TF24_hyperpar(...)
}

##' Hyperparameter function for TF24t physiological model
##' @title Hyperparameter function for TF24t physiological model
##' @param m A matrix of trait values, as returned by \code{trait_matrix}
##' @param s A strategy object
##' @param filter A flag indicating whether to filter columns. If TRUE, any numbers
##' that are within eps of the default strategy are not replaced.
##' @export
##' @rdname TF24t_hyperpar
TF24t_hyperpar <- make_TF24t_hyperpar()

#' @export
#' @importFrom rlang .data
#' @rdname expand_state
TF24t_expand_state <- function(results) {
  data <- split(results$species, results$species$species)

  for (i in seq_len(results$n_spp)) {
    s <- results$p$strategies[[i]]
    d <- data[[i]]
    # Derived size/mass columns are computed by the strategy's own C++
    # allometry functions (see src/strategy_expand.cpp) rather than being
    # re-derived here, so the formulas live in exactly one place.
    allom <- TF24t_strategy_expand_allometry(s, d$height, d$area_heartwood,
                                            d$mass_heartwood)
    data[[i]] <- dplyr::bind_cols(d, allom)
  }

  results$species <- data %>% dplyr::bind_rows()

  results
}
