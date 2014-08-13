##' @export
assembler_stochastic_naive <- function(community0, vcv,
                                       n_mutants=1L, n_immigrants=1L,
                                       seed_rain_eps=1e-3,
                                       compute_viable_fitness=TRUE,
                                       jump_to_attractor=FALSE,
                                       filename=NULL) {
  births_sys <- make_births_stochastic_naive(n_mutants, vcv, n_immigrants)
  deaths_sys <- make_deaths_stochastic_naive(seed_rain_eps)
  assembler(community0, births_sys, deaths_sys, filename,
            compute_viable_fitness, jump_to_attractor)
}

## Support functions:
make_births_stochastic_naive <- function(n_mutants, vcv, n_immigrants,
                                         check_positive=TRUE) {
  mutation    <- make_mutation_stochastic_naive(n_mutants, vcv)
  immigration <- make_immigration_stochastic_naive(n_immigrants)
  function(sys) {
    to_add <- rbind(mutation(sys), immigration(sys))
    if (check_positive && nrow(to_add) > 0) {
      seed_production <- sys$make_landscape()
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
make_immigration_stochastic_naive <- function(n_immigrants) {
  function(sys) {
    bounds <- sys$bounds
    n_traits <- length(sys$trait_names)
    if (is.null(bounds)) {
      ret <- matrix(nrow=0, ncol=n_traits)
    } else {
      lower <- log(bounds[,1])
      range <- log(bounds[,2]) - lower
      n <- rpois(1, n_immigrants)
      if (n == 0) {
        ret <- matrix(nrow=0, ncol=n_traits)
      } else {
        u <- t(lhs::randomLHS(n, n_traits))
        ret <- exp(t(lower + range * u))
      }
    }
    colnames(ret) <- sys$trait_names
    ret
  }
}

## Deaths
make_deaths_stochastic_naive <- function(eps) {
  function(sys) {
    sys$drop(sys$seed_rain() < eps)
  }
}

##' @export
mutational_vcv_proportion <- function(x, p=0.001) {
  if (inherits(x, "community")) {
    bounds <- x$bounds
  } else {
    bounds <- x
  }
  p * diag(nrow(bounds)) * as.numeric(diff(t(log(bounds))))
}
