source("helper-tree.R")

context("Species [CohortTop]")

s <- new(Strategy)
seed <- new(CohortTop, s)
h0 <- seed$height

sp <- new(SpeciesCT, s)

## Note that this initialises with a special base cohort always.
expect_that(sp$size, equals(0))
expect_that(sp$n_individuals, equals(0))
expect_that(sp$height_max, is_identical_to(seed$height))

env <- test.environment(seed$height * 10, seed.rain=1.0)

## Causes initial conditions to be estimated:
sp$compute_vars_phys(env)
seed$compute_initial_conditions(env)

sp$add_seeds(1)

expect_that(sp[[1]]$vars_phys,
            is_identical_to(seed$vars_phys))
expect_that(sp$ode_values,
            is_identical_to(seed$ode_values))

## Internal and test seed report same values:
expect_that(sp$seed$vars_phys,
            is_identical_to(seed$vars_phys))
expect_that(sp$seed$ode_values,
            is_identical_to(seed$ode_values))
## Finished with that seed now, so delete.
rm(seed)

plant.seed <- new(Plant, s)
expect_that(sp$germination_probability(env),
            is_identical_to(plant.seed$germination_probability(env)))

expect_that(sp$leaf_area_above(0), equals(0))

sp$height <- 1
h <- 0
x <- c(sp$seed$height, sp$height)
y <- c(sp$seed$leaf_area_above(h), sp[[1]]$leaf_area_above(h))

expect_that(sp$leaf_area_above(h),
            is_identical_to(trapezium(x, y)))

## Better tests: I want cases where:
## 1. empty: throws error

sp$clear()

## Re-set up the initial conditions
sp$compute_vars_phys(env)

expect_that(sp[[1]], throws_error())

## TODO: What to do about negative height?

## 1: empty species (no cohorts) has no leaf area above any height:
test_that("Empty species has no leaf area", {
  expect_that(sp$leaf_area_above(0), equals(0))
  expect_that(sp$leaf_area_above(10), equals(0))
  expect_that(sp$leaf_area_above(Inf), equals(0))
})

## 2: Cohort up against boundary has no leaf area:
sp$add_seeds(1)
test_that("Species with only boundary cohort no leaf area", {
  expect_that(sp$leaf_area_above(0), equals(0))
  expect_that(sp$leaf_area_above(10), equals(0))
  expect_that(sp$leaf_area_above(Inf), equals(0))
})

## 3: Single cohort; one round of trapezium:
h.top <- h0 * 4
sp$height <- h.top
cmp <- function(h, sp) {
  x <- c(sp$height, sp$seed$height)
  y <- c(sapply(sp$plants, function(p) p$leaf_area_above(h)),
         sp$seed$leaf_area_above(h))
  trapezium(rev(x), rev(y))
}

test_that("Leaf area sensible with one cohort", {
  ## At base and top
  expect_that(sp$leaf_area_above(0), equals(cmp(0, sp)))
  expect_that(sp$leaf_area_above(h.top), equals(0))
  ## Part way up (and above bottom seed boundary condition)
  expect_that(sp$leaf_area_above(h.top * .5),
              equals(cmp(h.top * .5, sp)))
})

sp$add_seeds(1)
sp$height <- h.top * c(1, .6)

test_that("Leaf area sensible with two cohorts", {
  ## At base and top
  expect_that(sp$leaf_area_above(0), equals(cmp(0, sp)))
  expect_that(sp$leaf_area_above(h.top), equals(0))
  ## Part way up (below bottom cohort, above boundarty condition)
  expect_that(sp$leaf_area_above(h.top * .5),
              equals(cmp(h.top * .5, sp)))
  ## Within the top pair (excluding the seed)
  expect_that(sp$leaf_area_above(h.top * .8),
              equals(cmp(h.top * .8, sp)))
})

sp$add_seeds(1)
sp$height <- h.top * c(1, .75, .6)
test_that("Leaf area sensible with two cohorts", {
  ## At base and top
  expect_that(sp$leaf_area_above(0), equals(cmp(0, sp)))
  expect_that(sp$leaf_area_above(h.top), equals(0))
  ## Part way up (below bottom cohort, above boundarty condition)
  expect_that(sp$leaf_area_above(h.top * .5),
              equals(cmp(h.top * .5, sp)))
  ## Within the top pair (excluding the seed)
  expect_that(sp$leaf_area_above(h.top * .8),
              equals(cmp(h.top * .8, sp)))
})

rm(sp)
gc()
