source("helper-tree.R")

context("CohortDiscrete")

s <- new(Strategy)
pars.s <- s$parameters

plant <- new(Plant, s)
cohort0 <- new(CohortDiscrete, s)
cohort1 <- new(CohortDiscrete, s, 1)
cohort2 <- new(CohortDiscrete, s, 2)

expect_that(new(CohortDiscrete, s, 0),
            throws_error())
expect_that(new(CohortDiscrete, s, -1),
            throws_error())

expect_that(cohort0$n_individuals, equals(1))
expect_that(cohort1$n_individuals, equals(1))
expect_that(cohort2$n_individuals, equals(2))

h0 <- 10
plant$height <- h0
cohort0$height <- h0
cohort1$height <- h0
cohort2$height <- h0

expect_that(cohort0$vars_size, is_identical_to(plant$vars_size))
expect_that(cohort1$vars_size, is_identical_to(plant$vars_size))
expect_that(cohort2$vars_size, is_identical_to(plant$vars_size))

## These have not actually been set/computed yet.
expect_that(cohort0$vars_phys, is_identical_to(plant$vars_phys))
expect_that(cohort1$vars_phys, is_identical_to(plant$vars_phys))
expect_that(cohort2$vars_phys, is_identical_to(plant$vars_phys))

expect_that(cohort0$ode_rates, is_identical_to(plant$ode_rates))
expect_that(cohort1$ode_rates, is_identical_to(plant$ode_rates))
expect_that(cohort2$ode_rates, is_identical_to(plant$ode_rates))

expect_that(cohort0$ode_values, is_identical_to(plant$ode_values))
expect_that(cohort1$ode_values, is_identical_to(plant$ode_values))
expect_that(cohort2$ode_values, is_identical_to(plant$ode_values))

hh <- seq(0, plant$height, length=101)
cmp <- sapply(hh, plant$leaf_area_above)
expect_that(sapply(hh, cohort0$leaf_area_above),
            is_identical_to(cmp))
expect_that(sapply(hh, cohort1$leaf_area_above),
            is_identical_to(cmp))
expect_that(sapply(hh, cohort2$leaf_area_above),
            is_identical_to(2L*cmp))

env <- test.environment(plant$height)
light.env <- attr(env, "light.env")

plant$compute_vars_phys(env)
cohort0$compute_vars_phys(env)
cohort1$compute_vars_phys(env)
cohort2$compute_vars_phys(env)

expect_that(cohort0$vars_phys, is_identical_to(plant$vars_phys))
expect_that(cohort1$vars_phys, is_identical_to(plant$vars_phys))
expect_that(cohort2$vars_phys, is_identical_to(plant$vars_phys))

expect_that(cohort0$leaf_area, is_identical_to(plant$leaf_area))
expect_that(cohort1$leaf_area, is_identical_to(plant$leaf_area))
expect_that(cohort2$leaf_area, is_identical_to(2L*plant$leaf_area))

## Grow the plants for a bit so that we can test the offspring
## calculations.
derivs <- function(t, y, pars) {
  plant <- pars$plant
  light.env <- pars$light.env

  plant$set_ode_values(t, y)
  plant$compute_vars_phys(light.env)
  plant$ode_rates
}

## Make a bigger light environment, so the plant has room to grow
env <- test.environment(plant$height)
light.env <- attr(env, "light.env")

h0 <- 10
plant$height <- h0
cohort0$height <- h0
cohort1$height <- h0
cohort2$height <- h0

obj.p <- new(OdeR, derivs, new.env(), list(plant=plant,   light.env=env))
obj.0 <- new(OdeR, derivs, new.env(), list(plant=cohort0, light.env=env))
obj.1 <- new(OdeR, derivs, new.env(), list(plant=cohort1, light.env=env))
obj.2 <- new(OdeR, derivs, new.env(), list(plant=cohort2, light.env=env))

t0 <- 0.0
y0 <- plant$ode_values
tmp <- obj.p$derivs(t0, y0)
expect_that(obj.0$derivs(t0, y0), is_identical_to(tmp))
expect_that(obj.1$derivs(t0, y0), is_identical_to(tmp))
expect_that(obj.2$derivs(t0, y0), is_identical_to(tmp))

tt <- seq(0, 8, length=51)
res.p <- obj.p$run(tt, y0)
res.0 <- obj.0$run(tt, y0)
res.1 <- obj.1$run(tt, y0)
res.2 <- obj.2$run(tt, y0)

expect_that(res.0, is_identical_to(res.p))
expect_that(res.1, is_identical_to(res.p))
expect_that(res.2, is_identical_to(res.p))

n.p <- plant$offspring()
expect_that(cohort0$offspring(), is_identical_to(n.p))
expect_that(cohort1$offspring(), is_identical_to(n.p))
expect_that(cohort2$offspring(), is_identical_to(2L * n.p))

expect_that(plant$offspring(), is_identical_to(0L))
expect_that(cohort0$offspring(), is_identical_to(0L))
expect_that(cohort1$offspring(), is_identical_to(0L))
expect_that(cohort2$offspring(), is_identical_to(0L))

test_that("System state get/set works", {
  vals.p <- plant$state
  vals.0 <- cohort0$state
  vals.1 <- cohort1$state
  vals.2 <- cohort2$state

  expect_that(cohort0$state_size, equals(plant$ode_size+1))
  expect_that(vals.p, is_identical_to(plant$ode_values))
  expect_that(vals.0, is_identical_to(c(cohort0$ode_values, 1)))
  expect_that(vals.1, is_identical_to(c(cohort1$ode_values, 1)))
  expect_that(vals.2, is_identical_to(c(cohort1$ode_values, 2)))

  expect_that(cohort0$state, is_identical_to(cohort1$state))

  x <- cohort0$ode_values + runif(plant$ode_size)
  expect_that(cohort0$state <- x,          throws_error())
  expect_that(cohort0$state <- c(x, 1, 2), throws_error())

  n <- cohort0$n_individuals + 1
  cohort0$state <- c(x, n)
  expect_that(cohort0$state, is_identical_to(c(x, n)))

  cohort0$state <- cohort1$state
})

y1 <- res.p[,ncol(res.p)]
y1[2] <- 0.3

cohort2$n_individuals <- 10
expect_that(cohort2$n_individuals, equals(10))
cohort2$n_individuals <- 2

f.p <- function(obj, y) {
  t <- 0.0 # arbitrary, ignored.
  obj$set_ode_values(t, y)
  obj$died()
}

f.c <- function(obj, y) {
  n <- obj$n_individuals
  t <- 0.0 # arbitrary, ignored.
  obj$set_ode_values(t, y)
  ret <- obj$died()
  ret <- c(ret, obj$n_individuals)
  obj$n_individuals <- n
  ret
}

nrep <- 1000
set.seed(1)
d.p <- replicate(nrep, f.p(plant, y1))
set.seed(1)
d.0 <- replicate(nrep, f.c(cohort0, y1))
set.seed(1)
d.1 <- replicate(nrep, f.c(cohort1, y1))
set.seed(1)
d.2 <- replicate(nrep, f.c(cohort2, y1))

expect_that(as.logical(d.0[1,]),
            is_identical_to(d.p))
expect_that(d.0, is_identical_to(d.1))

## Whenever we died, the population has decreased to 0
expect_that(all(as.logical(d.0[1,]) == (d.0[2,] == 0)),
            is_true())

## This also holds for the 2 individual case:
expect_that(all(as.logical(d.2[1,]) == (d.2[2,] == 0)),
            is_true())

## Quick check that probably belongs in test-plant
n.p <- sum(d.p)
r <- qbinom(c(1/10, 9/10), nrep, 1 - exp(-y1[2]))
expect_that(n.p,
            is_within_interval(r))

## Similar check for the 2-individual case
err <- max(abs(table(d.2[2,]) - table(rbinom(nrep, 2, exp(-y1[2])))))
expect_that(err/nrep, is_less_than(1/10))
