# Built from R/tf24f.R for the strategy: TF24w (#681)
##' @title Create a TF24w Plant or Node
##' @param s A \code{\link{TF24w_Strategy}} object
##' @export
##' @rdname TF24w_Individual
##' @examples
##' pl <- TF24w_Individual()
##' pl$height
TF24w_Individual <- function(s=TF24w_Strategy()) {
  Individual("TF24w", "TF24_Env")(s)
}

##' @title Set up a TF24w system with default or specified parameters
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
##' @rdname TF24w_Parameters
##' @examples
##' p1 <- TF24w_Parameters()
##' p2 <- TF24w_Parameters(max_patch_lifetime = 10.0)
TF24w_Parameters <- function(...) {
  Parameters("TF24w","TF24_Env")(...)
}

##' Generates a report on stand grown with TF24w strategy
##'
##' Builds a detailed report on stand grown with TF24w strategy, based on the template Rmd file provided.  The reports are
##' rendered as html files and saved in the specified output folder.
##'
##' @inheritParams FF16_generate_stand_report
##' @rdname TF24w_generate_stand_report
##' @return html file of the rendered report located in the specified output folder.
##' @export
TF24w_generate_stand_report <- function(results,
                                    output_file = "TF24w_report.html",
                                    overwrite = FALSE,
                                    target_ages = NA,
                                    input_file = system.file("reports", "TF24w_report.Rmd", package = "plant"),
                                    quiet = TRUE) {
  

  output_dir <- dirname(output_file)
  
  if (!file.exists(output_dir)) {
    dir.create(output_dir, FALSE, TRUE)
  }
  
  #output_file <- basename(output_file)

  if (overwrite | !file.exists(output_file)) {
    # knit and render. Note, call render directly
    # in preference to knit, then render, as leaflet widget
    # requires this to work
    result <-
      rmarkdown::render(
        input_file,
        output_dir = output_dir,
        output_file = output_file,
        quiet = quiet,
        params = list(
          results = results,
          target_ages = target_ages
        )
    )

    # remove temporary Rmd
    message(sprintf("Report for TF24w stand saved at %s", output_file))
  } else {
    message(sprintf("Report for TF24w stand already exists at %s", output_file))
  }
}

##' @export
##' @rdname make_TF24_hyperpar
##' @param ... Arguments passed to \code{\link{make_TF24_hyperpar}}
make_TF24w_hyperpar <- function(...) {
  ## TF24w is TF24 with a private root-zone soil column per cohort (#681):
  ## TF24w_Strategy inherits TF24_Strategy and reuses TF24_Environment +
  ## TF24_Pars, so its hyperparameterisation is identical to the parent.
  ## Delegate to keep TF24 the single source of truth.
  make_TF24_hyperpar(...)
}

##' Hyperparameter function for TF24w physiological model
##' @title Hyperparameter function for TF24w physiological model
##' @param m A matrix of trait values, as returned by \code{trait_matrix}
##' @param s A strategy object
##' @param filter A flag indicating whether to filter columns. If TRUE, any numbers
##' that are within eps of the default strategy are not replaced.
##' @export
##' @rdname TF24w_hyperpar
TF24w_hyperpar <- make_TF24w_hyperpar()

#' @export
#' @importFrom rlang .data
#' @rdname expand_state
TF24w_expand_state <- function(results) {
  data <- split(results$species, results$species$species)

  for (i in seq_len(results$n_spp)) {
    s <- results$p$strategies[[i]]
    d <- data[[i]]
    # Derived size/mass columns are computed by the strategy's own C++
    # allometry functions (see src/strategy_expand.cpp) rather than being
    # re-derived here, so the formulas live in exactly one place.
    allom <- TF24w_strategy_expand_allometry(s, d$height, d$area_heartwood,
                                            d$mass_heartwood)
    data[[i]] <- dplyr::bind_cols(d, allom)
  }

  results$species <- data %>% dplyr::bind_rows()

  results
}
