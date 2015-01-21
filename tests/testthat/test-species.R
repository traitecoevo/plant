if (interactive()) {
  devtools::load_all("../../")
  library(testthat)
  source("helper-tree2.R")
}

## TODO: The tests here really warrant splitting into different chunks
## - this was ported over from tree1 where the tests were loose in the
## file.

## TODO: Should 'height' become 'heights' to better indicate it's a
## plural?

context("Species")

env <- test_environment(3, seed_rain=1.0)

test_that("Basics", {
  s <- Strategy()
  sp <- Species(s)
  seed <- Cohort(s)
  plant <- Plant(s)
  h0 <- seed$height

  expect_that(sp$size, equals(0))
  expect_that(sp$height_max, is_identical_to(h0))
  expect_that(sp$species, is_identical_to(NULL))
  expect_that(sp$height, is_identical_to(numeric(0)))
  expect_that(sp$leaf_area, is_identical_to(numeric(0)))
  expect_that(sp$leaf_area_error, is_identical_to(numeric(0)))
  expect_that(sp$ode_size, equals(0))
  expect_that(sp$ode_state, is_identical_to(numeric(0)))
  expect_that(sp$ode_rates, is_identical_to(numeric(0)))

  ## Causes initial conditions to be estimated:
  sp$compute_vars_phys(env)
  seed$compute_initial_conditions(env)

  ## Internal and test seed report same values:
  expect_that(sp$seed$vars_phys,
              is_identical_to(seed$vars_phys))
  expect_that(sp$seed$ode_state,
              is_identical_to(seed$ode_state))

  sp$add_seed()
  expect_that(sp$size, equals(1))

  plants <- sp$plants
  expect_that(plants, is_a("list"))
  expect_that(length(plants), equals(1))
  expect_that(plants[[1]]$vars_phys, is_identical_to(seed$vars_phys))
  expect_that(sp$height, equals(seed$height))
  expect_that(sp$leaf_area, equals(seed$leaf_area))
  ## NOTE: Didn't check ode values

  ## Internal and test seed report same values:
  expect_that(sp$seed$vars_phys,
              is_identical_to(seed$vars_phys))

  expect_that(sp$plant_at(1), is_a("Cohort"))
  expect_that(sp$plant_at(1)$vars_phys,
              is_identical_to(plants[[1]]$vars_phys))

  ## Not sure about this -- do we need more immediate access?
  expect_that(sp$seed$plant$germination_probability(env),
              is_identical_to(plant$germination_probability(env)))

  expect_that(sp$leaf_area_above(0), equals(0))

  sp$height <- 1

  h <- 0
  x <- c(sp$seed$height, sp$height)
  y <- c(sp$seed$leaf_area_above(h),
         sp$plant_at(1)$leaf_area_above(h))

  expect_that(sp$leaf_area_above(h),
              is_identical_to(trapezium(x, y)))

  ## Better tests: I want cases where:
  ## 1. empty: throws error
  sp$clear()

  ## Re-set up the initial conditions
  sp$compute_vars_phys(env)

  expect_that(sp$plant_at(1), throws_error("Index 1 out of bounds"))
  expect_that(sp$plant_at(0), throws_error("Invalid value for index"))
})

## 1: empty species (no cohorts) has no leaf area above any height:
test_that("Empty species has no leaf area", {
  sp <- Species(Strategy())
  expect_that(sp$leaf_area_above(0), equals(0))
  expect_that(sp$leaf_area_above(10), equals(0))
  expect_that(sp$leaf_area_above(Inf), equals(0))
})

## 2: Cohort up against boundary has no leaf area:
test_that("Species with only boundary cohort no leaf area", {
  sp <- Species(Strategy())
  sp$add_seed()
  sp$compute_vars_phys(env)
  expect_that(sp$leaf_area_above(0), equals(0))
  expect_that(sp$leaf_area_above(10), equals(0))
  expect_that(sp$leaf_area_above(Inf), equals(0))
})

cmp_leaf_area_above <- function(h, sp) {
  x <- c(sp$height, sp$seed$height)
  y <- c(sapply(sp$plants, function(p) p$leaf_area_above(h)),
         sp$seed$leaf_area_above(h))
  trapezium(rev(x), rev(y))
}

## 3: Single cohort; one round of trapezium:
test_that("Leaf area sensible with one cohort", {
  sp <- Species(Strategy())
  sp$compute_vars_phys(env)
  sp$add_seed()
  h_top <- sp$height_max * 4
  sp$height <- h_top

  ## At base and top
  expect_that(sp$leaf_area_above(0), is_more_than(0))
  expect_that(sp$leaf_area_above(0), equals(cmp_leaf_area_above(0, sp)))

  expect_that(sp$leaf_area_above(h_top), is_identical_to(0.0))

  ## Part way up (and above bottom seed boundary condition)
  expect_that(sp$leaf_area_above(h_top * .5),
              equals(cmp_leaf_area_above(h_top * .5, sp)))

  ode_state <- sp$ode_state
  p <- sp$plant_at(1)
  expect_that(length(ode_state), equals(4))
  expect_that(ode_state, is_identical_to(p$ode_state))
})

test_that("Leaf area sensible with two cohorts", {
  sp <- Species(Strategy())
  sp$compute_vars_phys(env)
  sp$add_seed()
  h_top <- sp$height_max * 4
  sp$add_seed()
  sp$height <- h_top * c(1, .6)

  ## At base and top
  expect_that(sp$leaf_area_above(0), is_more_than(0))
  expect_that(sp$leaf_area_above(0), equals(cmp_leaf_area_above(0, sp)))
  expect_that(sp$leaf_area_above(h_top), equals(0))
  ## Part way up (below bottom cohort, above boundarty condition)
  expect_that(sp$leaf_area_above(h_top * .5),
              equals(cmp_leaf_area_above(h_top * .5, sp)))
  ## Within the top pair (excluding the seed)
  expect_that(sp$leaf_area_above(h_top * .8),
              equals(cmp_leaf_area_above(h_top * .8, sp)))

  ode_state <- sp$ode_state
  plants <- sp$plants
  expect_that(length(ode_state), equals(4 * sp$size))
  expect_that(ode_state,
              is_identical_to(unlist(lapply(plants, function(p) p$ode_state))))
})

test_that("Leaf area sensible with three cohorts", {
  sp <- Species(Strategy())
  sp$compute_vars_phys(env)
  sp$add_seed()
  h_top <- sp$height_max * 4
  sp$add_seed()
  sp$add_seed()
  sp$height <- h_top * c(1, .75, .6)

  ## At base and top
  expect_that(sp$leaf_area_above(0), is_more_than(0))
  expect_that(sp$leaf_area_above(0), equals(cmp_leaf_area_above(0, sp)))
  expect_that(sp$leaf_area_above(h_top), equals(0))
  ## Part way up (below bottom cohort, above boundarty condition)
  expect_that(sp$leaf_area_above(h_top * .5),
              equals(cmp_leaf_area_above(h_top * .5, sp)))
  ## Within the top pair (excluding the seed)
  expect_that(sp$leaf_area_above(h_top * .8),
              equals(cmp_leaf_area_above(h_top * .8, sp)))

  cmp_leaf_area <- sapply(seq_len(sp$size),
                          function(i) sp$plant_at(i)$leaf_area)
  expect_that(sp$leaf_area,
              is_identical_to(cmp_leaf_area))

  cmp <- local_error_integration(sp$height, cmp_leaf_area, 1)
  expect_that(sp$leaf_area_error,
              is_identical_to(cmp))

  ode_state <- sp$ode_state
  plants <- sp$plants
  expect_that(length(ode_state), equals(4 * sp$size))
  expect_that(ode_state,
              is_identical_to(unlist(lapply(plants, function(p) p$ode_state))))
})
