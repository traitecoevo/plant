##' Activate logging (deprecated)
##'
##' plant writes nothing to a log, so this has no effect. It is kept so that
##' existing scripts still run, and will be removed in a later release.
##' @title Activate logging (deprecated)
##' @param file_name,colour,threshold Ignored.
##' @export
plant_log_console <- function(file_name = "console", colour = TRUE,
                              threshold = "INFO") {
  .Deprecated(msg = "plant_log_console() has no effect: plant writes nothing to a log.")
  invisible(NULL)
}
