
##' Set a suitable hyperparameter function for chosen physiological model
##' @title Hyperparameters for physiological model
##' @param type Any strategy name as a string, e.g.: \code{"FF16"}.
##' @param parameters A parameters object
##' @rdname Hyperparameter_functions
##' @export
# if you update this function (even syntactic changes) update the function update_smc_support in the scaffolder
make_hyperpar <- function(type) {
  switch(type,
         FF16=make_FF16_hyperpar,
         TF24=make_TF24_hyperpar,
         TF24f=make_TF24f_hyperpar,
         K93=make_K93_hyperpar,
         stop("Unknown type ", type))
}

##' @rdname Hyperparameter_functions
##' @export
param_hyperpar <- function(parameters) {
  type <- attr(parameters$strategy_default, "class")
  switch(type,
         FF16_Strategy=FF16_hyperpar,
         TF24_Strategy=TF24_hyperpar,
         TF24f_Strategy=TF24f_hyperpar,
         K93_Strategy=K93_hyperpar,
         stop("Unknown type ", type))
}


##' @rdname Hyperparameter_functions
##' @export
# if you update this function (even syntactic changes) update the function update_smc_support in the scaffolder
hyperpar <- function(type) {
  switch(type,
         FF16=FF16_hyperpar,
         TF24=TF24_hyperpar,
         TF24f=TF24f_hyperpar,
         K93=K93_hyperpar,
         stop("Unknown type ", type))
}

##' @rdname Environment
##' @export
environment_type <- function(type) {
  switch(type,
         FF16=sprintf("FF16_Env"),
         TF24=sprintf("TF24_Env"),
         TF24f=sprintf("TF24_Env"),
         K93=sprintf("K93_Env"),
         stop("Unknown type ", type))
}

##' Scientific version of a physiological model
##'
##' The scientific version increments only when a model's equations or default
##' parameters change the simulation output for identical inputs. It is
##' independent of the package \code{Version} (which also moves for refactors,
##' performance and interface changes). Downstream tools such as \pkg{logpile}
##' use it to decide when archived simulations must be re-run: reruns happen
##' when the scientific version changes, not on every software release.
##'
##' The version is returned as a string. It is usually a single integer
##' (\code{"1"}), but may be compound for a model defined relative to another:
##' \code{TF24f} is a fast approximation of \code{TF24}, so its version is
##' \code{"<TF24 version>.<approximation revision>"} (e.g. \code{"2.1"}); the
##' major component auto-tracks \code{TF24}, so a \code{TF24} scientific change
##' also invalidates \code{TF24f}.
##'
##' The number is authored in C++ (the \code{scientific_version} /
##' \code{approximation_revision} constants on each strategy class, see
##' \code{inst/include/plant/models/*_strategy.h}) and read here through the
##' compiled \code{strategy_scientific_version()}; there is no duplicated copy
##' in R.
##'
##' @param type Any strategy name as a string, e.g.: \code{"FF16"}.
##' @return For \code{model_version}, a version string (e.g. \code{"1"} or
##'   \code{"2.1"}). For \code{model_id}, a string of the form \code{"FF16@v1"}
##'   (model name and scientific version).
##' @rdname model_version
##' @export
# if you add a new strategy, add its `scientific_version` constant to the model
# header and a dispatch arm to strategy_scientific_version() in src/strategy_version.cpp
model_version <- function(type) {
  strategy_scientific_version(type)
}

##' @rdname model_version
##' @export
model_id <- function(type) {
  sprintf("%s@v%s", type, model_version(type))
}

##' Creates an environment object of specified type
##' @param type Any environment name as a string, e.g.: \code{"FF16_Env"}.
##' @rdname Environment
##' @export
Environment <- function(type = NULL) {

  switch(type,
         FF16=FF16_Environment(),
         FF16_Env=FF16_Environment(),
         TF24=TF24_Environment(),
         TF24f=TF24_Environment(),
         TF24_Env=TF24_Environment(),
         K93=K93_Environment(),
         K93_Env=K93_Environment(),
         stop("Unknown type ", type))
}

#' Add additional state variables to the species component in output of a model.
#'
#' @param results from `tidy_patch`
#' @return similar format to input, but with additional columns for additional state variables
#' @export
#' @importFrom rlang .data
#' @rdname expand_state
expand_state <- function(results) {
  type <- extract_RcppR6_template_types(results$p, "Parameters")[[1]][1]

  switch(type,
    FF16 = FF16_expand_state(results),
    TF24 = TF24_expand_state(results),
    TF24f = TF24f_expand_state(results),
    K93 = K93_expand_state(results),
    stop("Unknown type ", type)
  )
}

node_schedule_default <- function(p) {
  cl <- class(p)[[1]]
  switch(cl,
         "Parameters<FF16,FF16_Env>"=`node_schedule_default__Parameters___FF16__FF16_Env`,
         "Parameters<TF24,TF24_Env>"=`node_schedule_default__Parameters___TF24__TF24_Env`,
         "Parameters<TF24f,TF24_Env>"=`node_schedule_default__Parameters___TF24f__TF24_Env`,
         "Parameters<K93,K93_Env>"=`node_schedule_default__Parameters___K93__K93_Env`,
         stop("Unknown type: ", cl))(p)
}

make_node_schedule <- function(p) {
  cl <- class(p)[[1]]
  switch(cl,
         "Parameters<FF16,FF16_Env>"=`make_node_schedule__Parameters___FF16__FF16_Env`,
         "Parameters<TF24,TF24_Env>"=`make_node_schedule__Parameters___TF24__TF24_Env`,
         "Parameters<TF24f,TF24_Env>"=`make_node_schedule__Parameters___TF24f__TF24_Env`,
         "Parameters<K93,K93_Env>"=`make_node_schedule__Parameters___K93__K93_Env`,
                  stop("Unknown type: ", cl))(p)
}
