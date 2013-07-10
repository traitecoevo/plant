source("helper-tree.R")

context("Patch [CohortTop]")

p <- new(Parameters)
p$add_strategy(new(Strategy))
p$seed_rain <- pi/2

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

## At this point, doing this should fail -- with only one individual
## the light environment is not defined.  At the same time, we should
## probably recover more gracefully and agree that it is also an empty
## light environment.
patch.c$compute_light_environment()
expect_that(patch.c$leaf_area_above(0), is_identical_to(0))

patch.c$add_seedling(1)
expect_that(patch.c$ode_size, equals(4))

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

solver <- solver.from.odetarget(patch.c, p$control$ode_control)
solver$step()
patch.c$add_seedling(1)
expect_that(patch.c$ode_size,
            equals(cmp$ode_size * patch.c$n_individuals))

patch.c$reset()
test_that("clear was successful", {
  expect_that(patch.c$n_individuals, equals(0))
  expect_that(patch.c$ode_size, equals(0))
  expect_that(patch.c$time, equals(0))
})

t <- patch.c$time

patch.c$add_seedling(1)
h <- patch.c$height[[1]]

while (patch.c$time < 25) {
  solver$step()
  t <- c(t, patch.c$time)
  h <- c(h, patch.c$height[[1]])
}

## TODO: This is not really a test, but we need to look at this and
## see if it makes any sense at all.
if (interactive()) {
  plot(t, h, type="l")
  plot(patch.c$environment$light_environment$xy, type="l")
}

patch.c$reset()
patch.c$add_seedling(1)
solver <- solver.from.odetarget(patch.c, p$control$ode_control)

tt <- seq(0, 25, length=26)
hh <- patch.c$height[[1]]
for (ti in tt[-1]) {
  solver$advance(ti)
  hh <- c(hh, patch.c$height[[1]])
}

if (interactive()) {
  plot(t, h, type="l")
  points(tt, hh)
}

expect_that(hh, equals(spline(t, h, xout=tt)$y, tolerance=1e-7))

test_that("OK at end of sequence", {
  expect_that(patch.c$time, is_identical_to(tt[[length(tt)]]))
  solver$advance(tt[[length(tt)]])
  expect_that(patch.c$time, is_identical_to(tt[[length(tt)]]))
  expect_that(solver$advance(tt[[length(tt)]] - 1e-8),
              throws_error())
})

rm(patch.c)
gc()
