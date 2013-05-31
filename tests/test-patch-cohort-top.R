source("helper-tree.R")

context("Patch [CohortTop]")

p <- new(Parameters)
p$add_strategy(new(Strategy))

## An individual CohortTop that will be the same in terms of strategy
## (and initial mass).
cmp <- new(CohortTop, p[[1]])

patch.c <- new(PatchCohortTop, p)

expect_that(patch.c$height_max, is_identical_to(cmp$height))

test_that("Initial patch strcture is sensible", {
  expect_that(patch.c$size, equals(1))
  expect_that(patch.c$ode_size, equals(0))
  expect_that(patch.c$time, equals(0))
})

r <- pi/2
patch.c$set_seed_rain(seed_rain(r))
expect_that(patch.c$environment$seed_rain$seed_rain,
            is_identical_to(r))

## This should not work, because we can't set environment from within
## patch.
expect_that(patch.c$environment$seed_rain <- seed_rain(r * 2),
            throws_error())

## At this point, doing this should fail -- with only one individual
## the light environment is not defined.  At the same time, we should
## probably recover more gracefully and agree that it is also an empty
## light environment.
patch.c$compute_light_environment()
expect_that(patch.c$leaf_area_above(0), is_identical_to(0))

patch.c$add_seedling(1)
expect_that(patch.c$ode_size, equals(4))

## Now that we've got started, we should not be able to set the seed
## rain:
expect_that(patch.c$set_seed_rain(seed_rain(1.0)),
            throws_error())

cmp$compute_initial_conditions(patch.c$environment)

expect_that(patch.c$ode_values,
            is_identical_to(cmp$ode_values))
expect_that(patch.c$ode_rates,
            is_identical_to(cmp$ode_rates))

y <- patch.c$ode_values
patch.c$set_ode_values(0, y)
expect_that(patch.c$ode_values, is_identical_to(y))

## NOTE: These should be identical, but are merely equal...
expect_that(patch.c$derivs(0, y),
            equals(cmp$ode_rates))

patch.c$step()
patch.c$add_seedling(1)
expect_that(patch.c$ode_size,
            equals(cmp$ode_size * patch.c$n_individuals))

patch.c$clear()
test_that("clear was successful", {
  expect_that(patch.c$n_individuals, equals(0))
  expect_that(patch.c$ode_size, equals(0))
  expect_that(patch.c$time, equals(0))
})

t <- patch.c$time

patch.c$add_seedling(1)
h <- patch.c$height[[1]]

while (patch.c$time < 25) {
  patch.c$step()
  t <- c(t, patch.c$time)
  h <- c(h, patch.c$height[[1]])
}

## TODO: This is not really a test, but we need to look at this and
## see if it makes any sense at all.
if (FALSE) {
  plot(t, h, type="l")
  plot(patch.c$environment$light_environment$xy, type="l")
}

if (FALSE) {
  ## Add a cohort every '1'.
  dt <- 1
  t.next <- 1
  patch.c$clear()
  patch.c$add_seedling(1)
  t <- patch.c$time
  h <- list(patch.c$height[[1]])

  ## This is quite a bit slower than Daniel's version!  Probably much of
  ## that is coming from things like the adaptive light refinement and
  ## he integration for computing photosynthesis.
  while (patch.c$time < 20) {
    patch.c$step()
    if (patch.c$time > t.next) {
      patch.c$add_seedling(1)
      t.next <- t.next + dt
    }
    t <- c(t, patch.c$time)
    h <- c(h, list(patch.c$height[[1]]))
  }

  n <- length(h[[length(h)]])
  h <- t(sapply(h, function(x) c(x, rep(NA, n-length(x)))))
  matplot(t, h, type="l", col="black", lty=1)
}

rm(patch.c)
gc()
