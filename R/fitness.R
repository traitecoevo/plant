##' Construct a fitness landscape.
##'
##' @title Fitness Landscape
##' @param trait_matrix A matrix of traits corresponding to mutants to
##' introduce into the light environment constructed by the residents
##' in \code{p}.
##' @param p Parameters object.  Needs to contain residents with their
##' incoming seed rain.
##' @param raw_seed_rain Logical; if \code{TRUE} report per capita
##' seed rain rather than fitness.
##' @return Vector with the output seed rain.  Mutants have an
##' arbitrary seed rain of one, so this is the rate of seed
##' production per capita.
##' @author Rich FitzJohn
##' @export
fitness_landscape <- function(trait_matrix, p, hyperpar=param_hyperpar(p), raw_seed_rain=FALSE) {
  n_residents <- length(p$strategies)
  p_with_mutants <- expand_parameters(trait_matrix, p, hyperpar)
  scm <- run_scm(p_with_mutants,
                 use_ode_times=length(p$cohort_schedule_ode_times) > 0)
  seed_rain <- scm$seed_rains
  if (n_residents > 0L) {
    seed_rain <- seed_rain[-seq_len(n_residents)]
  }
  if (raw_seed_rain) {
    seed_rain
  } else {
    log(seed_rain)
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
max_growth_rate <- function(trait_matrix, p) {
  fitness_landscape(trait_matrix, remove_residents(p))
}

##' Compute the carrying capacity (equilibrium per-capita seed
##' production) for a set of values of a trait.  Each is considered in
##' isolation.
##'
##' @title Carrying Capacity
##' @param trait_matrix A matrix of traits
##' @param p Parameters object to use.  Importantly, the
##' \code{strategy_default} element gets used here.
##' @param seed_rain Initial seed rain (optional)
##' @param parallel Use multiple processors?
##' @author Rich FitzJohn
##' @export
carrying_capacity <- function(trait_matrix, p, seed_rain=1,
                              parallel=FALSE) {
  f <- function(x) {
    ## NOTE: This is not very pretty because matrix_to_list, which we
    ## use to iterate over this, does not preserve array-ness or
    ## names.
    p <- expand_parameters(trait_matrix(x, traits),
                           remove_residents(p),
                           mutant=FALSE)
    p$seed_rain <- seed_rain
    equilibrium_seed_rain(p)$seed_rain
  }
  traits <- colnames(trait_matrix)
  unlist(loop(matrix_to_list(trait_matrix), f, parallel=parallel))
}
