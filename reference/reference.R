library(tree)
library(testthat)
library(digest)

p <- new(Parameters)
p$add_strategy(new(Strategy, list(lma=0.1)))
p$add_strategy(new(Strategy, list(lma=0.3)))
p$seed_rain <- c(1.1, 2.2)

path <- "ref"
## Evolve path hard coded to my machine for now.
## Running with version: d649534d6d92fd7343a57c9d26d51eec8e155f2d
path.evolve <- "~/Documents/Projects/veg/falster-traitdiversity/src"
## run.reference(path, p, path.evolve=path.evolve, verbose=TRUE)

output <- load.reference.output(path)
test_that("Output contains correct parameters", {
  compare.parameters <- function(p1, p2)
    isTRUE(all.equal(tree:::reference.from.parameters(p1),
                     tree:::reference.from.parameters(p2)))
  expect_that(compare.parameters(output$parameters, p), is_true())
})

test_that("Output is as expected (harsh test)", {
  expected.sha <- "be8be2963d91c85d9bec07e4eb43751f827e4266"
  expect_that(digest(output[names(output) != "parameters"], "sha1"),
              is_identical_to(expected.sha))
})

## 1: Disturbance:
age <- output$patch_age

test_that("Output disturbance calculations match", {
  ## tree version:
  d <- new(Disturbance, p$parameters[["mean_disturbance_interval"]])
  f.survival <- Vectorize(function(x) d$survival_probability(0, x))
  f.density  <- Vectorize(function(x) d$density(x))

  scal <- integrate(f.survival, 0, Inf)$value
  expect_that(f.density(age$age), equals(age$density))
  expect_that(f.survival(age$age) / scal, equals(age$density))
})

## 2: Find out when the cohorts were introduced.

## We're going to have to do some tweaking to get the model run at the
## same time points.  This is not yet implemented.
