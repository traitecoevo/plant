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


##' @export
##' @rdname strategy_list
strategy_default <- function(parameters, hyperpar=param_hyperpar(parameters)) {
  strategy(trait_matrix(1, "a")[, -1, drop=FALSE], parameters, hyperpar)
}

##' @export
##' @rdname strategy_list
strategy <- function(x, parameters, hyperpar=param_hyperpar(parameters)) {
  if (nrow(x) != 1L) {
    stop("Expected a single type")
  }
  strategy_list(x, parameters, hyperpar)[[1]]
}

##' @rdname strategy_list
##' @export
individual_list <- function(x, parameters, hyperpar=param_hyperpar(parameters)) {
  
  if (!inherits(parameters, "Parameters")) {
    stop("parameters must be a 'Parameters' object")
  }
  types <- extract_RcppR6_template_types(parameters, "Parameters")
  lapply(strategy_list(x, parameters, hyperpar), do.call('Individual', types))
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
##' @param strategy A model of species traits to use as a template
##' @param trait_matrix A matrix of traits corresponding to the
##' new types to introduce.
##' @param hyperpar Hyperparameter function to use. By default links 
##' to standard function for this strategy type. 
##' @param mutant Create new types as \emph{mutants}?  These will have
##' no effect on other plants in the community (i.e. have zero
##' density).
##' @author Rich FitzJohn
##' @export
expand_traits <- function(strategy = FF16_Strategy(), 
                          trait_matrix, 
                          hyperpar_fn = NULL) {
  
  if (!is.matrix(trait_matrix)) {
    stop("Invalid type of trait_matrix -- expected a matrix")
  }
  
  # some strategies don't use hyperparameterisation
  # in which case this passes through
  if (is.null(hyperpar)) {
    hyperpar_fn = hyperpar(strategy)
  }
  
  # override defaults with hyper-parameterisation
  species_matrix <- hyperpar_fn(trait_matrix, strategy)
  trait_names <- colnames(species_matrix)
    
  f <- function(xi) {
    new_traits <- strategy
    new_traits[trait_names] <- xi
    new_traits
  }
  
  species_list <- lapply(matrix_to_list(species_matrix), f)
  return(species_list)
}

set_introduction_parameters <- function(species_list, strategy, control, mutant = FALSE) {
  if (length(mutant) != 1L) {
    stop("mutant must be scalar")
  }
  type <- class(strategy)
  species_parameters <- SpeciesParameters(type)()

  n_spp <- length(species_list)
  species_parameters$species <- species_list
  species_parameters$is_resident <- rep(!mutant, n_spp)
  species_parameters$birth_rate <- rep(1.0, n_spp)

  species_parameters$setup_cohort_schedule(control)  
  species_parameters$validate()
  
  return(species_parameters)
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
