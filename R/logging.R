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
  logger::log_threshold(logger::OFF, namespace="plant")
}
