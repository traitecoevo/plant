## Functions to support naive assembly.

## This code all derives from the successional-diversity project.
## There is considerable overlap with code in Revolve, but at this
## point it's not totally clear how general the code should be.

## Generator function for the "mutation" step.  There are two things
## here; mutation and immigrants.
##' @export
make_births <- function(n_mutants, vcv, n_immigrants, bounds,
                        seed_rain0) {
  mutation    <- make_mutation(n_mutants, vcv)
  immigration <- make_immigration(n_immigrants, bounds)
  function(sys, must_grow=FALSE) {
    repeat {
      new_traits <- rbind(mutation(sys), immigration())
      if (nrow(new_traits) > 0 || !must_grow) {
        break
      }
    }
    if (nrow(new_traits) > 0) {
      sys$traits <- rbind(sys$traits, new_traits)
      sys$seed_rain <- c(sys$seed_rain,
                         rep(seed_rain0, nrow(new_traits)))
    }
    sys
  }
}

## Mutation: Draw (on average) n_mutants from the population with a
## mutational variance of vcv on the log scale.
##' @export
make_mutation <- function(n_mutants, vcv) {
  n_traits <- ncol(vcv)
  blank <- matrix(nrow=0, ncol=n_traits)
  function(sys) {
    traits <- sys[["traits"]]
    weights <- sys[["seed_rain"]]
    assertthat::assert_that(is.matrix(traits)             &&
                            ncol(traits) == n_traits      &&
                            nrow(traits) == length(weights))
    if (length(weights) == 0) {
      return(blank)
    }
    n <- rpois(1, n_mutants)
    if (n == 0) {
      return(blank)
    }
    i <- sample(length(weights), n, replace=TRUE, prob=weights)
    exp(log(traits[i,,drop=FALSE]) + rmvnorm(n, sigma=vcv))
  }
}

## Immigration: Same from the bounds, in log space.
##' @export
make_immigration <- function(n_immigrants, bounds) {
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

## Deaths:
##' @export
make_deaths <- function(eps) {
  function(sys) {
    drop <- sys[["seed_rain"]] < eps
    if (any(drop)) {
      keep <- !drop
      sys$traits    <- sys$traits[keep,,drop=FALSE]
      sys$seed_rain <- sys$seed_rain[keep]
    }
    sys
  }
}
