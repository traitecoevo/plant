## Because we work mostly with things that are just pointers to data
## structures, serialisation is hard, and currently unsolved.  These
## are hacks to get things working for now.
serialise_parameters <- function(p) {
  if (p$size > 0) {
    stop("Only working for empty parameters right now")
  }
  ret <- list()
  ret$disturbance <- p$disturbance$mean_interval
  ret$control <- p$control$parameters
  ret$strategy_default <- p$strategy_default$parameters
  ret$parameters <- p$parameters
  ret$species <- list()
  class(ret) <- "parameters_serialised"
  ret
}

unserialise_parameters <- function(obj) {
  if (!inherits(obj, "parameters_serialised")) {
    stop("Expected a serialised parameters object")
  }
  if (length(obj$species) > 0) {
    stop("Can't serialise non-empty parameters")
  }
  p <- new(Parameters)
  p$set_parameters(obj$parameters)
  p$set_control_parameters(obj$control)
  p$disturbance <- new(Disturbance, obj$disturbance)
  p$strategy_default <- new(Strategy, obj$strategy)
  p
}
