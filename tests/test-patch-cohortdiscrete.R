source("helper-tree.R")
options(error=traceback)

context("Patch [CohortDiscrete]")

p <- new(Parameters)
p$add_strategy(new(Strategy))

## A plant that will be the same in terms of strategy (and initial
## mass).
cmp <- new(Plant, p$get_strategy(0))

patch.p <- new(Patch,  p)
patch.c <- new(PatchC, p)

expect_that(patch.p$height_max, is_identical_to(0.0))
expect_that(patch.c$height_max, is_identical_to(0.0))

## Add a single seed.
patch.p$add_seeds(1)
patch.c$add_seeds(1)

plants.p <- patch.p$get_plants()
plants.c <- patch.c$get_plants()

expect_that(length(plants.c), equals(1))
expect_that(length(plants.c[[1]]), equals(1))
expect_that(length(plants.c[[c(1,1)]]), equals(1))
expect_that(class(plants.c[[c(1,1)]]),
            equals("Rcpp_CohortDiscrete", check.attr=FALSE))
expect_that(class(plants.p[[c(1,1)]]),
            equals("Rcpp_Plant", check.attr=FALSE))

## Then clear both 
patch.p$clear()
patch.c$clear()

## Add a *pair* of seeds.
patch.p$add_seeds(2)
patch.c$add_seeds(2)

## And expect that the size of the system is 6 with no cohorts, and 3
## with cohorts.
expect_that(patch.p$ode_size, equals(6))
expect_that(patch.c$ode_size, equals(3))

## Now, grow this pair of seeds deterministically for a bit.
f <- function(patch, t) {
  res <- list()
  while ( patch$age < t ) {
    patch$step_deterministic()
    res <- c(res, list(c(patch$age, patch$ode_values)))
  }
  do.call(rbind, res)
}

res.p <- f(patch.p, 10)
res.c <- f(patch.c, 10)

expect_that(res.p[,1:4], equals(res.c[,1:4]))
expect_that(res.p[,5:7], equals(res.p[,2:4]))

