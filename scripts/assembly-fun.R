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
