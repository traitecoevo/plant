##' Construct a fitness landscape.
##'
##' @title Fitness Landscape
##' @param trait_matrix A matrix of traits corresponding to mutants to
##' introduce into the light environment constructed by the residents
##' in \code{p}.
##' @param p Parameters object.  Needs to contain residents with their
##' arriving offspring.
##' @param raw_offspring_produced Logical; if \code{TRUE} report per capita
##' offspring production rather than fitness.
##' @return Vector with the output offspring produced.  Mutants have an
##' arbitrary offspring of one, so this is the rate of offspring
##' production per capita.
##' @author Rich FitzJohn
##' @export
fitness_landscape <- function(trait_matrix, p, hyperpar=param_hyperpar(p), raw_offspring_produced=FALSE) {
  n_residents <- length(p$strategies)
  p_with_mutants <- expand_parameters(trait_matrix, p, hyperpar)
  scm <- run_scm(p_with_mutants,
                 use_ode_times=length(p$cohort_schedule_ode_times) > 0)
  net_reproduction_ratios <- scm$net_reproduction_ratios
  if (n_residents > 0L) {
    net_reproduction_ratios <- net_reproduction_ratios[-seq_len(n_residents)]
  }
  if (raw_offspring_produced) {
    net_reproduction_ratios
  } else {
    log(net_reproduction_ratios)
  }
}

##' Compute max growth rate of a given set of values for a trait.
##' This is the log of per-capita offspring production (i.e., fitness).
##'
##' Only works in one dimension
##' @title Compute Max Growth Rate
##' @param trait_matrix A matrix of traits to compute maximum growth
##' rate at (see \code{\link{trait_matrix}})
##' @param p Parameters object to use.  Importantly, the
##' \code{strategy_default} element gets used here.
##' @author Rich FitzJohn
##' @export
max_growth_rate <- function(trait_matrix, p) {
  fitness_landscape(trait_matrix, remove_residents(p))
}

##' Compute the carrying capacity (equilibrium per-capita offspring
##' production) for a set of values of a trait.  Each is considered in
##' isolation.
##'
##' @title Carrying Capacity
##' @param trait_matrix A matrix of traits
##' @param p Parameters object to use.  Importantly, the
##' \code{strategy_default} element gets used here.
##' @param birth_rate Initial offspring arrival (optional)
##' @param parallel Use multiple processors?
##' @author Rich FitzJohn
##' @export
carrying_capacity <- function(trait_matrix, p, birth_rate=1,
                              parallel=FALSE) {
  f <- function(x) {
    ## NOTE: This is not very pretty because matrix_to_list, which we
    ## use to iterate over this, does not preserve array-ness or
    ## names.
    p <- expand_parameters(trait_matrix(x, traits),
                           remove_residents(p),
                           mutant=FALSE)
    p$birth_rate <- birth_rate
    equilibrium_offspring_arriving(p)$birth_rate
  }
  traits <- colnames(trait_matrix)
  unlist(loop(matrix_to_list(trait_matrix), f, parallel=parallel))
}
