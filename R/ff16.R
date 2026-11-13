##' Create a FF16 Individual
##' @title Create a FF16 Individual
##' @param s A \code{\link{FF16_Model}} object
##' @export
##' @rdname FF16_individual
##' @examples
##' pl <- FF16_Individual()
##' pl$height
FF16_Individual <- function(s=FF16_Model()) {
  Individual("FF16", "FF16_Env")(s)
}

##' @title Setup an a model system with default or specified parameters
##'
##' @description Setup an a model system with default or specified parameters. 
##' This function enables you initialize a model system. Use the model name to start different models. 
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
##' @rdname FF16_Parameters
##' @examples
##' p1 <- FF16_Parameters()
##' p2 <- FF16_Parameters(max_patch_lifetime = 10.0, patch_area = 1.0, strategies = list(FF16_Model()), strategy_default = FF16_Model(), node_schedule_times_default = node_schedule_times_default(10.0), node_schedule_times = list(node_schedule_times_default(10.0)), ode_times = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
FF16_Parameters <- function(...) {
  Parameters("FF16","FF16_Env")(...)
}

# Todo -- generalise report generation

##' Generates a report on stand grown with FF16 strategy
##'
##' Builds a detailed report on stand grown with FF16 strategy, based on the template Rmd file provided.  The reports are
##' rendered as html files and saved in the specified output folder.
##'
##' @param results results of running \code{tidy(run(p, ..., record = "steps"))}
##' @param output_file name of output file
##' @param overwrite logical value to determine whether to overwrite existing report
##' @param target_ages Patches ages at which to make plots
##' @param input_file report script (.Rmd) file to build study report
##' @param quiet An option to suppress printing during rendering from knitr, pandoc command line and others.
##'
##' @rdname FF16_generate_stand_report
##' @return html file of the rendered report located in the specified output folder.
##' @export
FF16_generate_stand_report <- function(results,
                                    output_file = "FF16_report.html",
                                    overwrite = FALSE,
                                    target_ages = NA,
                                    input_file = system.file("reports", "FF16_report.Rmd", package = "plant"),
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
    message(sprintf("Report for FF16 stand saved at %s", output_file))
  } else {
    message(sprintf("Report for FF16 stand already exists at %s", output_file))
  }
}

#' @export
#' @importFrom rlang .data
#' @rdname expand_state
FF16_expand_state <- function(results) {
  expand_state_by_allometry(results, FF16_strategy_expand_allometry)
}
