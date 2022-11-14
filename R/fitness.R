##' Construct a fitness landscape.
##'
##' @title Fitness Landscape
##' @param trait_matrix A matrix of traits corresponding to mutants to
##' introduce into the light environment constructed by the residents
##' in \code{p}.
##' @param p Parameters object.  Needs to contain residents with their
##' incoming seed rain.
##' @param hyperpar Hyperparameter function to use. By default links
##' to standard function for this strategy type.
##' @param log_fitness Logical; if \code{TRUE} report per capita
##' seed rain rather than fitness.
##' @return Vector with the output seed rain.  Mutants have an
##' arbitrary seed rain of one, so this is the rate of seed
##' production per capita.
##' @author Rich FitzJohn
##' @export
fitness_landscape <- function(trait_matrix, p, hyperpar=param_hyperpar(p), log_fitness=TRUE) {
  n_residents <- length(p$strategies)
  
  birth_rates <- rep(0, nrow(trait_matrix))
  
  p_with_mutants <- expand_parameters(trait_matrix, p, hyperpar,
                                      mutant = TRUE,
                                      birth_rate_list = birth_rates)
  
  scm <- run_scm(p_with_mutants,
                 use_ode_times=length(p$node_schedule_ode_times) > 0)
  
  net_reproduction_ratios <- scm$net_reproduction_ratios
  if (n_residents > 0L) {
    net_reproduction_ratios <- net_reproduction_ratios[-seq_len(n_residents)]
  }
  if (log_fitness) {
    log(net_reproduction_ratios)
  } else {
    net_reproduction_ratios
  }
}

##' Compute max growth rate of a given set of values for a trait.
##' This is the log of per-capita seed production (i.e., fitness).
##'
##' Only works in one dimension
##' @title Compute Max Growth Rate
##' @param trait_matrix A matrix of traits to compute maximum growth
##' rate at (see \code{\link{trait_matrix}})
##' @param p Parameters object to use.  Importantly, the
##' \code{strategy_default} element gets used here.
##' @author Rich FitzJohn
##' @export
fundamental_fitness <- function(trait_matrix, p) {
  fitness_landscape(trait_matrix, remove_residents(p))
}
