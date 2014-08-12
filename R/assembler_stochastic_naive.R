##' @export
assembler_stochastic_naive <- function(community0,
                                       n_mutants=1L,
                                       n_immigrants=1L,
                                       vcv=NULL, vcv_p=0.001,
                                       seed_rain_eps=1e-3,
                                       compute_viable_fitness=FALSE,
                                       filename=NULL) {
  ## This is probably not a good idea and is causing some problems
  ## because we are having problems with recomputing bounds (below).
  ## Better would be to tidy this up on entry.
  if (is.null(vcv)) {
    vcv <- mutational_vcv_proportion(community0$bounds, vcv_p)
  }
  if (compute_viable_fitness) {
    community0$set_viable_bounds()
  }
  ## This makes the immigrants come from the good (viable) bounds
  births_sys <- make_births_stochastic_naive(n_mutants, vcv,
                                             n_immigrants, community0$bounds)
  deaths_sys <- make_deaths_stochastic_naive(seed_rain_eps)
  assembler(community0, births_sys, deaths_sys, filename)
}

## Support functions:
make_births_stochastic_naive <- function(n_mutants, vcv, n_immigrants, bounds, check_positive=TRUE) {
  mutation    <- make_mutation_stochastic_naive(n_mutants, vcv)
  immigration <- make_immigration_stochastic_naive(n_immigrants, bounds)
  function(sys) {
    to_add <- rbind(mutation(sys), immigration())
    # check seed production of new mutants if possible
    seed_production <- sys$make_landscape()
    if (check_positive && !is.null(seed_production) &&
        nrow(to_add) > 0) {
      R_new <- seed_production(to_add)
      to_add  <- to_add[R_new > 1,,drop=FALSE]
    }
    to_add
  }
}

## This code all derives from the successional-diversity project.
## There is considerable overlap with code in Revolve, but at this
## point it's not totally clear how general the code should be.

## Mutation: Draw (on average) n_mutants from the population with a
## mutational variance of vcv on the log scale.
make_mutation_stochastic_naive <- function(n_mutants, vcv) {
  n_traits <- ncol(vcv)
  blank <- matrix(nrow=0, ncol=n_traits)
  function(sys) {
    if (sys$size() == 0) {
      return(blank)
    }
    n <- rpois(1, n_mutants)
    if (n == 0) {
      return(blank)
    }
    traits <- sys$traits(TRUE)
    weights <- sys$seed_rain()
    i <- sample(length(weights), n, replace=TRUE, prob=weights)
    unname(exp(log(traits[i,,drop=FALSE]) + mvtnorm::rmvnorm(n, sigma=vcv)))
  }
}

## Immigration: Same from the bounds, in log space.
make_immigration_stochastic_naive <- function(n_immigrants, bounds) {
  lower <- log(bounds[,1])
  range <- log(bounds[,2]) - lower
  n_traits <- nrow(bounds)
  function() {
    n <- rpois(1, n_immigrants)
    if (n == 0) {
      m <- matrix(nrow=0, ncol=n_traits)
    } else {
      u <- t(lhs::randomLHS(n, n_traits))
      m <- exp(t(lower + range * u))
    }
    colnames(m) <- rownames(bounds)
    m
  }
}

## Deaths
make_deaths_stochastic_naive <- function(eps) {
  function(sys) {
    sys$drop(sys$seed_rain() < eps)
  }
}

mutational_vcv_proportion <- function(bounds, p=0.001) {
  p * diag(nrow(bounds)) * as.numeric(diff(t(log(bounds))))
}
