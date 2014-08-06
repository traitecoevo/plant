## Support functions:
make_births_stochastic_naive <- function(n_mutants, vcv, n_immigrants, bounds) {
  mutation    <- make_mutation_stochastic_naive(n_mutants, vcv)
  immigration <- make_immigration_stochastic_naive(n_immigrants, bounds)
  function(sys, must_grow=FALSE) {
    repeat {
      new_traits <- rbind(mutation(sys), immigration())
      if (nrow(new_traits) > 0 || !must_grow) {
        break
      }
    }
    new_traits
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
