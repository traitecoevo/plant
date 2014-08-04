## Generator function for the "mutation" step.  There are two things
## here; mutation and immigrants.
##
## Really though, all we want is a function that will generate new
## phenotypes from the current system.  So this set of three functions
## can be swapped out with something less silly later.  That would
## allow for things like sampling from the fitness landscape, or
## picking points estimated to be right at the top of the landscape,
## etc.
##
## For now, we have gaussian mutation of existing phenotypes (given a
## mutantal VCV), plus "uniform" sampling from the current trait space
## (using latin hypercube sampling).
make_new_phenotypes <- function(n_mutants, vcv, n_immigrants, bounds) {
  mutation    <- make_mutation(n_mutants, vcv)
  immigration <- make_immigration(n_immigrants, bounds)
  function(sys) {
    traits <- sys[["traits"]]
    weights <- sys[["seed_rain"]]
    rbind(mutation(traits, weights), immigration())
  }
}

## Mutation: Draw (on average) n_mutants from the population with a
## mutational variance of vcv on the log scale.
make_mutation <- function(n_mutants, vcv) {
  n_traits <- ncol(vcv)
  function(traits, weights) {
    assert_that(is.matrix(traits)             &&
                ncol(traits) == n_traits      &&
                nrow(traits) == length(weights))

    n <- rpois(1, n_mutants)
    if (n == 0) {
      m <- matrix(nrow=0, ncol=n_traits)
    } else {
      i <- sample(length(weights), n, replace=TRUE, prob=weights)
      m <- exp(log(traits[i,,drop=FALSE]) + rmvnorm(n, sigma=vcv))
    }
    m
  }
}

## Immigration: Same from the bounds, in log space.
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

######################################################################

## There is some tension here between what we want to operate in (a
## couple of simple matrices) and what we *need* to operate in
## (Parameters and CohortSchedule).  This is going to become
## particularly apparent when it comes time to pull delete
## strategies.
setup_parameters <- function(traits, seed_rain, p0) {
  p <- p0$copy()
  strategy <- p$strategy_default
  for (i in seq_len(nrow(traits))){
    new_strategy <- strategy$copy()
    new_strategy$set_parameters(as.list(traits[i,]))
    p$add_strategy(new_strategy)
  }
  p$seed_rain <- seed_rain
  p
}


setup_schedule <- function(times, max_t) {
  schedule <- new(CohortSchedule, length(times))
  schedule$max_time <- max_t
  schedule$all_times <- times
  schedule
}

drop_strategies_parameters <- function(p, drop) {
  assert_that(is.logical(drop) && length(drop) == p$size)
  ret <- p$copy()
  if (any(drop)) {
    ret$clear()
    for (i in which(!drop)) # keep these
      ret$add_strategy(p[[i]])
  }
  ret
}

drop_strategies_schedule <- function(schedule, drop) {
  assert_that(is.logical(drop) && length(drop) == schedule$n_species)
  if (any(drop)) {
    keep <- !drop
    ret <- new(CohortSchedule, sum(keep))
    ret$max_time <- schedule$max_time
    ret$all_times <- schedule$all_times[keep]
  } else {
    ret <- schedule$copy()
  }
  ret
}

## Do we need to add functions for deaths, new_phenotypes, seed_rain0
## and times0 here?  Probably not.
make_run_simulation <- function(p0, max_t) {
  force(p0)
  force(max_t)
  function(sys) {
    p <- setup_parameters(sys[["traits"]], sys[["seed_rain"]], p0)
    schedule <- setup_schedule(sys[["times"]], max_t)
    res <- build_schedule(p)
    rain_out <- unname(attr(res, "seed_rain", exact=TRUE)[,"out"])
    list(traits=sys[["traits"]], seed_rain=rain_out, times=res$all_times)
  }
}

make_deaths <- function(eps) {
  force(eps)
  function(sys) {
    keep <- sys[["seed_rain"]] >= eps
    list(traits    = sys[["traits"]][keep,,drop=FALSE],
         seed_rain = sys[["seed_rain"]][keep],
         times     = sys[["times"]][keep])
  }
}

make_births <- function(new_phenotypes, seed_rain0, times0) {
  force(new_phenotypes)
  force(seed_rain0)
  force(times0)
  function(sys) {
    traits.new  <- new_phenotypes(sys[["traits"]], sys[["seed_rain"]])
    n_new <- nrow(traits.new)
    list(traits    = rbind(sys[["traits"]], traits.new),
         seed_rain = c(sys[["seed_rain"]], rep(seed_rain0,   n_new)),
         times     = c(sys[["times"]],     rep(list(times0), n_new)))
  }
}
