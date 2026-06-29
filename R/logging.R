## Log fields included in output:
##
## * time    (from logger)
## * level   (from logger)
## * routine {equilibrium, viable, inviable, ...}
plant_log_info <- function(message, routine=NULL, ...) {
  if (!is.null(routine)) message <- paste0(crayon::yellow(routine), "> ", message)
  logger::log_info(message, namespace="plant", .topenv=parent.frame())
}

plant_log_debug <- function(message, routine=NULL, ...) {
  if (!is.null(routine)) message <- paste0(crayon::yellow(routine), "> ", message)
  logger::log_debug(message, namespace="plant", .topenv=parent.frame())
}

plant_log_eq <- function(...) {
  plant_log_info(..., routine="equilibrium")
}

plant_log_viable <- function(...) {
  plant_log_info(..., routine="viable")
}

plant_log_inviable <- function(...) {
  plant_log_info(..., routine="inviable")
}

make_plant_layout <- function(colour=TRUE) {
  col_time <- if (colour) crayon::silver else identity
  function(level, msg, namespace, .logcall, .topcall, .topenv, .timestamp=Sys.time()) {
    sprintf("[%s] %s",
            col_time(format(.timestamp, "%H:%M:%OS3")),
            paste(msg, collapse=" "))
  }
}

##' Activate logging
##'
##' By default plant prints little information about its progress.
##' This can be modified by enabling logging.  Log entries include a
##' timestamp and, where applicable, a routine label indicating which
##' part of the model is running.
##'
##' "Schedule" events (splitting) are sent to the DEBUG stream,
##' everything else is sent to INFO.
##' @title Activate logging
##' @param file_name File to save output (default = "console")
##' @param colour Use colour in console output?
##' @param threshold Minimum log level to emit: "DEBUG", "INFO", etc.
##' @export
plant_log_console <- function(file_name="console", colour=TRUE, threshold="INFO") {
  logger::log_threshold(threshold, namespace="plant")
  if (file_name == "console") {
    logger::log_appender(logger::appender_console, namespace="plant")
  } else {
    logger::log_appender(logger::appender_file(file_name), namespace="plant")
  }
  logger::log_layout(make_plant_layout(colour), namespace="plant")
}

.onLoad <- function(libname, pkgname) {
  logger::log_threshold(logger::OFF, namespace = "plant")

  # Wrap Leaf() to error on partial argument name matches.
  # R's partial matching silently accepts any unambiguous prefix (e.g. `vcma`
  # resolves to `vcmax_25`), turning typos into wrong-value bugs.  We derive
  # valid_args from formals() so the check stays correct automatically whenever
  # RcppR6.R is regenerated with new parameters.
  ns <- asNamespace(pkgname)
  original_Leaf <- get("Leaf", envir = ns, inherits = FALSE)
  valid_args <- names(formals(original_Leaf))

  leaf_wrapper <- function() {
    supplied <- names(as.list(sys.call())[-1])
    supplied <- supplied[nzchar(supplied)]
    bad <- setdiff(supplied, valid_args)
    if (length(bad) > 0) {
      stop(sprintf(
        "Unknown argument(s) to Leaf(): %s\nCheck for misspelled parameter names.",
        paste(bad, collapse = ", ")
      ), call. = FALSE)
    }
    env <- environment()
    do.call(original_Leaf, lapply(valid_args, function(nm) get(nm, envir = env)))
  }
  formals(leaf_wrapper) <- formals(original_Leaf)
  assign("Leaf", leaf_wrapper, envir = ns)
}
