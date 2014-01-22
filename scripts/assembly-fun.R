## Generator function for the "mutation" step.  There are two things
## here; mutation and immigrants.
make.new.phenotypes <- function(n.mutants, vcv, n.immigrants, bounds) {
  mutation    <- make.mutation(n.mutants, vcv)
  immigration <- make.immigration(n.immigrants, bounds)
  function(traits, weights)
    rbind(mutation(traits, weights), immigration())
}

## Mutation: Draw (on average) n.mutants from the population with a
## mutational variance of vcv on the log scale.
make.mutation <- function(n.mutants, vcv) {
  n.traits <- ncol(vcv)
  function(traits, weights) {
    assert_that(is.matrix(traits)             &&
                ncol(traits) == n.traits      &&
                nrow(traits) == length(weights))

    n <- rpois(1, n.mutants)
    if (n == 0)
      return(matrix(nrow=0, ncol=n.traits))
    i <- sample(length(weights), n, replace=TRUE, prob=weights)
    exp(log(traits[i,,drop=FALSE]) + rmvnorm(n, sigma=vcv))
  }
}

## Immigration: Same from the bounds, in log sace.
make.immigration <- function(n.immigrants, bounds) {
  lower <- log(bounds[,1])
  range <- log(bounds[,2]) - lower
  n.traits <- nrow(bounds)
  function() {
    n <- rpois(1, n.immigrants)
    if (n == 0)
      return(matrix(nrow=0, ncol=n.traits))
    u <- t(randomLHS(n, n.traits))
    exp(t(lower + range * u))
  }
}

######################################################################

## There is some tension here between what we want to operate in (a
## couple of simple matrices) and what we *need* to operate in
## (Parameters and CohortSchedule).  This is going to become
## particularly apparent when it comes time to pull delete
## strategies.
setup.parameters <- function(traits, seed_rain, p0) {
  p <- p0$copy()
  for (i in seq_len(nrow(traits)))
    p$add_strategy(new(Strategy, as.list(traits[i,])))
  p$seed_rain <- seed_rain

  p
}

setup.schedule <- function(times, max.t) {
  schedule <- new(CohortSchedule, length(times))
  schedule$max_time <- max.t
  schedule$all_times <- times
  schedule
}

drop.strategies.parameters <- function(p, drop) {
  assert_that(is.logical(drop) && length(drop) == p$size)
  ret <- p$copy()
  if (any(drop)) {
    ret$clear()
    for (i in which(!drop)) # keep these
      ret$add_strategy(p[[i]])
  }
  ret
}

drop.strategies.schedule <- function(schedule, drop) {
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

make.run <- function(p0, max.t, build.args=list()) {
  build.args <- modifyList(list(nsteps=10, eps=1e-3, verbose=TRUE),
                           build.args)
  force(p0)
  force(max.t)
  function(sys) {
    p <- setup.parameters(sys[["traits"]], sys[["seed_rain"]], p0)
    schedule <- setup.schedule(sys[["times"]], max.t)

    res <- build.schedule(p, schedule, build.args$nsteps, build.args$eps,
                          progress=FALSE, verbose=build.args$verbose)
    rain.out <- unname(attr(res, "seed_rain", exact=TRUE)[,"out"])

    list(traits=sys[["traits"]], seed_rain=rain.out, times=res$all_times)
  }
}

make.deaths <- function(eps) {
  force(eps)
  function(sys) {
    keep <- sys[["seed_rain"]] >= eps
    list(traits    = sys[["traits"]][keep,,drop=FALSE],
         seed_rain = sys[["seed_rain"]][keep],
         times     = sys[["times"]][keep])
  }
}

make.births <- function(new.phenotypes, seed_rain0, times0) {
  force(new.phenotypes)
  force(seed_rain0)
  force(times0)
  function(sys) {
    traits.new  <- new.phenotypes(sys[["traits"]], sys[["seed_rain"]])
    n.new <- nrow(traits.new)
    list(traits    = rbind(sys[["traits"]], traits.new),
         seed_rain = c(sys[["seed_rain"]], rep(seed_rain0,   n.new)),
         times     = c(sys[["times"]],     rep(list(times0), n.new)))
  }
}
