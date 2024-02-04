##' Create a list of Strategies or Plants by varying a single trait.
##'
##' @title Create a list of Strategies
##' @param x Values for the trait.  This must be a \emph{matrix}, with
##'   column names corresponding to entries in \code{Strategy} and
##'   rows representing different values.
##' @param hyperpar Hyperparameter function to use. By default links
##' to standard function for this strategy type.
##' @param parameters \code{Parameters} object containing a
##'   default strategy to modify.  Any hyperparameterisation included
##'   will be applied.
##' @param birth_rate_list List object with birth rates for each species in
##' x. Birth rates can take the form of a scalar (constant) or a vector.
##' In either case birth rates are set as \code{strategy$birth_rate_y}, however
##' varying birth rates will also have \code{strategy$birth_rate_x} and
##  \code{`is_variable_bithrate = TRUE`}
##'
##' @export
strategy_list <- function(x, parameters, hyperpar=param_hyperpar(parameters), birth_rate_list) {
  if (!is.matrix(x)) {
    stop("Invalid type x -- expected a matrix")
  }

  strategy <- parameters$strategy_default
  x <- hyperpar(x, strategy)

  trait_names <- colnames(x)
  f <- function(xi, br) {
    strategy[trait_names] <- xi
    if (is.list(br)) {
      strategy$birth_rate_x <- br$x
      strategy$birth_rate_y <- br$y
      strategy$is_variable_birth_rate <- TRUE
    } else if (is.numeric(br)) {
      strategy$birth_rate_y <- br
      strategy$is_variable_birth_rate <- FALSE
    } else {
      stop("Invalid type in birth_rate_list - need either a list with x, y control points or a numeric")
    }
    strategy
  }
  
  # insert custom traits and birth values into default strategy template
  strategies <- mapply(f, matrix_to_list(x), birth_rate_list, SIMPLIFY = FALSE)
  return(strategies)
}

##' @export
##' @rdname strategy_list
strategy_default <- function(parameters, hyperpar=param_hyperpar(parameters)) {
  strategy(trait_matrix(1, "a")[, -1, drop=FALSE], parameters, hyperpar)
}

##' @export
##' @rdname strategy_list
strategy <- function(x, parameters, hyperpar=param_hyperpar(parameters), birth_rate_list) {
  if (nrow(x) != 1L) {
    stop("Expected a single type")
  }
  strategy_list(x, parameters, hyperpar, birth_rate_list)[[1]]
}

##' @rdname strategy_list
##' @export
individual_list <- function(x, parameters, hyperpar=param_hyperpar(parameters), birth_rate_list) {

  if (!inherits(parameters, "Parameters")) {
    stop("parameters must be a 'Parameters' object")
  }
  types <- extract_RcppR6_template_types(parameters, "Parameters")
  lapply(strategy_list(x, parameters, hyperpar, birth_rate_list), do.call('Individual', types))
}

##' Helper function to create trait matrices suitable for
##' \code{\link{strategy_list}}.
##'
##' @title Create trait matrix
##' @param x Values
##' @param trait_name Name of a single trait
##' @export
##' @author Rich FitzJohn
trait_matrix <- function(x, trait_name) {
  m <- matrix(x, ncol=length(trait_name))
  colnames(m) <- trait_name
  m
}

##' The functions expand_parameters and mutant_parameters convert trait values into parametr objects for the model. By default, expand_parameters adds an extra strategy to existing.

##' @title Setup parameters to run resindets or mutants
##' @param trait_matrix A matrix of traits corresponding to the
##' new types to introduce.
##' @param p A \code{Parameters} object.
##' @param hyperpar Hyperparameter function to use. By default links
##' to standard function for this strategy type.
##' @param birth_rate_list List object with birth rates for each species in
##' @param keep_existing_strategies Should existing residents be retained
##' x. Birth rates can take the form of a scalar (constant) or a vector.
##' In either case birth rates are set as \code{strategy$birth_rate_y}, however
##' varying birth rates will also have \code{strategy$birth_rate_x} and
##  \code{`is_variable_bithrate = TRUE`}
##' @author Rich FitzJohn
##' @export
##' @rdname expand_parameters
expand_parameters <- function(trait_matrix, p, hyperpar=param_hyperpar(p), birth_rate_list = 1, keep_existing_strategies = TRUE) {

  if(nrow(trait_matrix) != length(birth_rate_list)) {
    stop("Must provide exactly one birth rate input for each species")
  }
  extra <- strategy_list(trait_matrix, p, hyperpar, birth_rate_list)
  n_extra <- length(extra)

  ret <- p <- validate(p) # Ensure times are set up correctly.

  ## Determine node introduction times
  if (length(p$strategies) == 0L) {
    times_new <- p$node_schedule_times_default
  } else {
    ## if residnets are present, use all unique times of all residents
    times_new <- unique(sort(unlist(p$node_schedule_times)))
  }
  
  if (keep_existing_strategies) {
    ret$strategies <- c(p$strategies, extra)
    ret$node_schedule_times <- c(p$node_schedule_times, 
                                  rep(list(times_new), n_extra))
  } else {
    ret$strategies <- extra
    ret$node_schedule_times <- rep(list(times_new), n_extra)
  }

  ## Clear this if it's present:
  attr(ret, "offspring_production") <- NULL

  ret
}

##' @export
##' @rdname expand_parameters
mutant_parameters <- function(..., keep_existing_strategies = FALSE) {
  expand_parameters(..., keep_existing_strategies = keep_existing_strategies)
}



remove_residents <- function(p) {
  if (length(p$strategies) > 0L) {
    p$strategies <- list()
    p$birth_rate <- numeric(0)
    p$node_schedule_times <- list()
  }
  p
}
