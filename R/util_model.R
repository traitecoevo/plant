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
  mapply(f, matrix_to_list(x), birth_rate_list, SIMPLIFY = FALSE)
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

##' Expand Parameters to include mutants.  All mutants get the same
##' schedule, equal to all the unique times that any resident was
##' introduced (if mutants) or the default schedule (if residents).
##' This results in more work than is really needed, but should be
##' reasonable most of the time.
##'
##' @title Expand Parameters to include mutants
##' @param trait_matrix A matrix of traits corresponding to the
##' new types to introduce.
##' @param p A \code{Parameters} object.
##' @param hyperpar Hyperparameter function to use. By default links 
##' to standard function for this strategy type. 
##' @param mutant Create new types as \emph{mutants}?  These will have
##' no effect on other plants in the community (i.e. have zero
##' density).
##' @author Rich FitzJohn
##' @export
expand_parameters <- function(trait_matrix, p, hyperpar=param_hyperpar(p), mutant=TRUE, birth_rate_list) {
  if (length(mutant) != 1L) {
    stop("mutant must be scalar")
  }
  if(nrow(trait_matrix) != length(birth_rate_list)) {
    stop("Must provide exactly one birth rate input for each species")
  }
  extra <- strategy_list(trait_matrix, p, hyperpar, birth_rate_list)
  n_extra <- length(extra)
  
  ret <- p <- validate(p) # Ensure times are set up correctly.
  ret$strategies <- c(p$strategies, extra)
  ret$is_resident <- c(p$is_resident, rep(!mutant, n_extra))
  
  ## Introduce mutants at all unique times:
  if (length(p$strategies) == 0L || !mutant) {
    times_new <- p$cohort_schedule_times_default
  } else {
    times_new <- unique(sort(unlist(p$cohort_schedule_times)))
  }
  ret$cohort_schedule_times <- c(p$cohort_schedule_times,
                                 rep(list(times_new), n_extra))
  
  ## Clear this if it's present:
  attr(ret, "net_reproduction_ratios") <- NULL
  
  ret
}

remove_residents <- function(p) {
  if (length(p$strategies) > 0L) {
    p$strategies <- list()
    p$is_resident <- logical(0)
    p$birth_rate <- numeric(0)
    p$cohort_schedule_times <- list()
  }
  p
}
