##' @rdname Environment
##' @export
environment_type <- function(type) {
  switch(type,
         FF16=sprintf("FF16_Env"),
         TF24=sprintf("TF24_Env"),
         TF24f=sprintf("TF24_Env"),
         TF24t=sprintf("TF24_Env"),
         TF24w=sprintf("TF24_Env"),
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
##' \code{TF24f} is an approximation of \code{TF24}, so its version is
##' \code{"<TF24 version>.<approximation revision>"} (e.g. \code{"2.1"}); the
##' major component auto-tracks \code{TF24}, so a \code{TF24} scientific change
##' also invalidates \code{TF24f}. \code{TF24t} is compound the same way,
##' \code{"<TF24 version>.<thermal revision>"}, the revision moving only for a
##' change to its thermal-damage layer.
##'
##' The number is authored in C++ (the \code{scientific_version},
##' \code{approximation_revision} and \code{thermal_revision} constants on each
##' model class, see \code{inst/include/plant/models/*_model.h}) and read
##' here through the
##' compiled \code{strategy_scientific_version()}; there is no duplicated copy
##' in R.
##'
##' @param type Any strategy name as a string, e.g.: \code{"FF16"}.
##' @return For \code{model_version}, a version string (e.g. \code{"1"} or
##'   \code{"2.1"}). For \code{model_id}, a string of the form \code{"FF16@v1"}
##'   (model name and scientific version).
##' @rdname model_version
##' @export
# if you add a new model, add its `scientific_version` constant to the model
# header and a dispatch arm to strategy_scientific_version() in src/model_exports.cpp
model_version <- function(type) {
  strategy_scientific_version(type)
}

##' @rdname model_version
##' @export
model_id <- function(type) {
  sprintf("%s@v%s", type, model_version(type))
}

##' Creates an environment object of specified type
##'
##' An environment is a live C++ object, which \code{saveRDS()} cannot keep. A
##' \code{\link{scenario}} therefore keeps its environment as a recipe: plain
##' values for everything a run reads from it (for TF24, the soil column, its
##' parameters, the soil water a run starts from and the drivers; FF16 and K93
##' read nothing a run does not reset). Given a recipe, \code{Environment()}
##' builds the environment it describes, so \code{Environment(res$scenario$env)}
##' gives back the one a result ran in, after \code{readRDS()} as well.
##' @param type Any environment name as a string, e.g.: \code{"FF16_Env"}, or
##'   a recipe, as a scenario's \code{env} holds one.
##' @rdname Environment
##' @export
Environment <- function(type = NULL) {
  if (inherits(type, "plant_environment_recipe")) {
    return(environment_from_recipe(type))
  }

  switch(type,
         FF16=FF16_Environment(),
         FF16_Env=FF16_Environment(),
         TF24=TF24_Environment(),
         TF24f=TF24_Environment(),
         TF24t=TF24_Environment(),
         TF24w=TF24_Environment(),
         TF24_Env=TF24_Environment(),
         K93=K93_Environment(),
         K93_Env=K93_Environment(),
         stop("Unknown type ", type))
}

## The recipe of an environment: its type, and the values its class's recipe()
## method reads, which set_recipe() writes back. A class without the pair
## (FF16, K93) holds nothing a run does not reset, so its recipe is its type.
## An environment that has a recipe already is returned as it is.
environment_recipe <- function(env) {
  if (inherits(env, "plant_environment_recipe")) {
    return(env)
  }
  type <- if (inherits(env, "R6")) {
    switch(class(env)[[1]],
           FF16_Environment = "FF16_Env",
           K93_Environment = "K93_Env",
           TF24_Environment = "TF24_Env",
           NULL)
  }
  if (is.null(type)) {
    stop("`env` must be an Environment, or a recipe of one as a scenario's ",
         "`env` holds", call. = FALSE)
  }
  check_environment_alive(env)
  values <- if (is.function(env$recipe)) env$recipe() else list()
  structure(c(list(type = type), values), class = "plant_environment_recipe")
}

environment_from_recipe <- function(recipe) {
  env <- Environment(recipe$type)
  values <- unclass(recipe)[setdiff(names(recipe), "type")]
  if (is.function(env$set_recipe)) {
    env$set_recipe(values)
  } else if (length(values) > 0L) {
    stop("A ", recipe$type, " environment's recipe holds its type alone, not ",
         paste(names(values), collapse = ", "), call. = FALSE)
  }
  env
}

## The environment a run is handed. A scenario holds a recipe; one built by
## hand, or saved before scenarios held recipes, may hold the environment.
live_environment <- function(env) {
  if (inherits(env, "plant_environment_recipe")) {
    return(environment_from_recipe(env))
  }
  check_environment_alive(env)
  env
}

check_environment_alive <- function(env) {
  ## Every environment has a time, and reading it fails only on a pointer that
  ## readRDS() left empty.
  alive <- !inherits(env, "R6") ||
    tryCatch({env$time; TRUE}, error = function(e) FALSE)
  if (!alive) {
    stop("This environment was read back with readRDS(), which keeps no live ",
         "object, so it must be built again. A scenario keeps its environment ",
         "as a recipe, which saveRDS() keeps; a scenario holding a dead ",
         "environment like this one was made by an older plant, and is made ",
         "again with scenario()", call. = FALSE)
  }
}

##' @export
print.plant_environment_recipe <- function(x, ...) {
  cat(sprintf("<plant_environment_recipe> %s", x$type))
  if (x$type == "TF24_Env") {
    varying <- names(Filter(function(d) is.null(d$value), x$drivers))
    cat(sprintf(", %d soil layers to %g m, drivers %s",
                length(x$layer_widths), x$depth,
                if (length(varying) == 0L) "all constant"
                else paste("varying:", paste(varying, collapse = ", "))))
  }
  cat("\n")
  invisible(x)
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
    TF24t = TF24t_expand_state(results),
    TF24w = TF24w_expand_state(results),
    K93 = K93_expand_state(results),
    stop("Unknown type ", type)
  )
}

schedule_default <- function(p) {
  parameters_export(p, "schedule_default")(p)
}

make_schedule <- function(p) {
  parameters_export(p, "make_schedule")(p)
}

# The C++ helpers exported once per <Model, Environment> pair, from
# src/model_exports.cpp, are named for the pair; this finds the one for `p`.
parameters_export <- function(p, prefix) {
  cl <- class(p)[[1]]
  types <- regmatches(cl, regexec("^Parameters<([^,>]+),([^,>]+)>$", cl))[[1]]
  f <- if (length(types) == 3L) {
    get0(sprintf("%s__Parameters___%s__%s", prefix, types[[2]], types[[3]]),
         envir = asNamespace("plant"), mode = "function", inherits = FALSE)
  }
  if (is.null(f)) {
    stop("Unknown type: ", cl)
  }
  f
}

# A run's species rows, one piece per species in strategy order. `species` is
# stored as character, so splitting on it directly orders "1", "10", "11", "2",
# ..., and from the tenth species on a piece would meet another's strategy.
split_by_species <- function(results) {
  levels <- as.character(seq_len(results$n_spp))
  split(results$species,
        factor(as.character(results$species$species), levels = levels))
}

# Each model's expand_state: the derived size and mass columns, computed by the
# strategy's own C++ allometry (src/strategy_expand.cpp) so the formulas live in
# one place, bound onto each species' rows.
expand_state_by_allometry <- function(results, allometry) {
  data <- split_by_species(results)
  for (i in seq_len(results$n_spp)) {
    d <- data[[i]]
    allom <- allometry(results$p$strategies[[i]], d$height, d$area_heartwood,
                       d$mass_heartwood)
    data[[i]] <- dplyr::bind_cols(d, allom)
  }
  results$species <- dplyr::bind_rows(data)
  results
}
