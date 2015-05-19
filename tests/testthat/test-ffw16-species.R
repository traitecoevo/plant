## TODO: The tests here really warrant splitting into different chunks
## - this was ported over from tree1 where the tests were loose in the
## file.

context("FFW16_Species")

test_that("Basics", {
  env <- test_environment(3, seed_rain=1.0)
  s <- FFW16_Strategy()
  sp <- FFW16_Species(s)
  seed <- FFW16_Cohort(s)
  plant <- FFW16_PlantPlus(s)
  h0 <- seed$height

  expect_that(sp$size, equals(0))
  expect_that(sp$height_max, is_identical_to(h0))
  expect_that(sp$species, is_identical_to(NULL))
  expect_that(sp$heights, is_identical_to(numeric(0)))
  expect_that(sp$log_densities, is_identical_to(numeric(0)))
  expect_that(sp$area_leafs, is_identical_to(numeric(0)))
  expect_that(sp$area_leafs_error(1.0), is_identical_to(numeric(0)))
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
  expect_that(sp$heights, equals(seed$height))
  expect_that(sp$log_densities, equals(seed$log_density))
  expect_that(sp$area_leafs, equals(seed$area_leaf))
  ## NOTE: Didn't check ode values

  ## Internal and test seed report same values:
  expect_that(sp$seed$vars_phys,
              is_identical_to(seed$vars_phys))

  expect_that(sp$plant_at(1), is_a("FFW16_Cohort"))
  expect_that(sp$plant_at(1)$vars_phys,
              is_identical_to(plants[[1]]$vars_phys))

  ## Not sure about this -- do we need more immediate access?
  expect_that(sp$seed$plant$germination_probability(env),
              is_identical_to(plant$germination_probability(env)))

  expect_that(sp$area_leaf_above(0), equals(0))

  sp$heights <- 1

  h <- 0
  x <- c(sp$seed$height, sp$heights)
  y <- c(sp$seed$area_leaf_above(h),
         sp$plant_at(1)$area_leaf_above(h))

  expect_that(sp$area_leaf_above(h),
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
  sp <- FFW16_Species(FFW16_Strategy())
  expect_that(sp$area_leaf_above(0), equals(0))
  expect_that(sp$area_leaf_above(10), equals(0))
  expect_that(sp$area_leaf_above(Inf), equals(0))
})

## 2: Cohort up against boundary has no leaf area:
test_that("FFW16_Species with only boundary cohort no leaf area", {
  env <- test_environment(3, seed_rain=1.0)
  sp <- FFW16_Species(FFW16_Strategy())
  sp$add_seed()
  sp$compute_vars_phys(env)
  expect_that(sp$area_leaf_above(0), equals(0))
  expect_that(sp$area_leaf_above(10), equals(0))
  expect_that(sp$area_leaf_above(Inf), equals(0))
})

cmp_area_leaf_above <- function(h, sp) {
  x <- c(sp$heights, sp$seed$height)
  y <- c(sapply(sp$plants, function(p) p$area_leaf_above(h)),
         sp$seed$area_leaf_above(h))
  trapezium(rev(x), rev(y))
}

## 3: Single cohort; one round of trapezium:
test_that("Leaf area sensible with one cohort", {
  env <- test_environment(3, seed_rain=1.0)
  sp <- FFW16_Species(FFW16_Strategy())
  sp$compute_vars_phys(env)
  sp$add_seed()
  h_top <- sp$height_max * 4
  sp$heights <- h_top

  ## At base and top
  expect_that(sp$area_leaf_above(0), is_more_than(0))
  expect_that(sp$area_leaf_above(0), equals(cmp_area_leaf_above(0, sp)))

  expect_that(sp$area_leaf_above(h_top), is_identical_to(0.0))

  ## Part way up (and above bottom seed boundary condition)
  expect_that(sp$area_leaf_above(h_top * .5),
              equals(cmp_area_leaf_above(h_top * .5, sp)))

  ode_state <- sp$ode_state
  p <- sp$plant_at(1)
  expect_that(length(ode_state), equals(4))
  expect_that(ode_state, is_identical_to(p$ode_state))
})

test_that("Leaf area sensible with two cohorts", {
  env <- test_environment(3, seed_rain=1.0)
  sp <- FFW16_Species(FFW16_Strategy())
  sp$compute_vars_phys(env)
  sp$add_seed()
  h_top <- sp$height_max * 4
  sp$add_seed()
  sp$heights <- h_top * c(1, .6)

  ## At base and top
  expect_that(sp$area_leaf_above(0), is_more_than(0))
  expect_that(sp$area_leaf_above(0), equals(cmp_area_leaf_above(0, sp)))
  expect_that(sp$area_leaf_above(h_top), equals(0))
  ## Part way up (below bottom cohort, above boundarty condition)
  expect_that(sp$area_leaf_above(h_top * .5),
              equals(cmp_area_leaf_above(h_top * .5, sp)))
  ## Within the top pair (excluding the seed)
  expect_that(sp$area_leaf_above(h_top * .8),
              equals(cmp_area_leaf_above(h_top * .8, sp)))

  ode_state <- sp$ode_state
  plants <- sp$plants
  expect_that(length(ode_state), equals(4 * sp$size))
  expect_that(ode_state,
              is_identical_to(unlist(lapply(plants, function(p) p$ode_state))))
})

test_that("Leaf area sensible with three cohorts", {
  env <- test_environment(3, seed_rain=1.0)
  sp <- FFW16_Species(FFW16_Strategy())
  sp$compute_vars_phys(env)
  sp$add_seed()
  h_top <- sp$height_max * 4
  sp$add_seed()
  sp$add_seed()
  sp$heights <- h_top * c(1, .75, .6)

  ## At base and top
  expect_that(sp$area_leaf_above(0), is_more_than(0))
  expect_that(sp$area_leaf_above(0), equals(cmp_area_leaf_above(0, sp)))
  expect_that(sp$area_leaf_above(h_top), equals(0))
  ## Part way up (below bottom cohort, above boundarty condition)
  expect_that(sp$area_leaf_above(h_top * .5),
              equals(cmp_area_leaf_above(h_top * .5, sp)))
  ## Within the top pair (excluding the seed)
  expect_that(sp$area_leaf_above(h_top * .8),
              equals(cmp_area_leaf_above(h_top * .8, sp)))

  cmp_area_leaf <- sapply(seq_len(sp$size),
                          function(i) sp$plant_at(i)$area_leaf)
  expect_that(sp$area_leafs,
              is_identical_to(cmp_area_leaf))

  cmp    <- local_error_integration(sp$heights, cmp_area_leaf, 1.0)
  cmp_pi <- local_error_integration(sp$heights, cmp_area_leaf, pi)

  expect_that(sp$area_leafs_error(),    is_identical_to(cmp))
  expect_that(sp$area_leafs_error(1.0), is_identical_to(cmp))
  expect_that(sp$area_leafs_error(pi),  is_identical_to(cmp_pi))

  ode_state <- sp$ode_state
  plants <- sp$plants
  expect_that(length(ode_state), equals(4 * sp$size))
  expect_that(ode_state,
              is_identical_to(unlist(lapply(plants, function(p) p$ode_state))))
})
