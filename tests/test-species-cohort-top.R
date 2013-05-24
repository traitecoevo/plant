source("helper-tree.R")

context("Species [CohortTop]")

s <- new(Strategy)
cmp <- new(Plant, s)

sp <- new(SpeciesCT, s)

## Note that this initialises with a special base cohort always.
expect_that(sp$size, equals(0))
expect_that(sp$n_individuals, equals(0))
expect_that(sp$height_max, is_identical_to(cmp$height))

env <- test.environment(cmp$height * 10)
env$seed_rain <- new(SeedRain, 1)

sp$add_seeds(1)
sp$compute_vars_phys(env)
cmp$compute_vars_phys(env)

expect_that(sp[[1]]$vars_phys,
            is_identical_to(cmp$vars_phys))

seed <- new(Plant, s)
expect_that(sp$germination_probability(env),
            is_identical_to(seed$germination_probability(env)))

## Not sure if this is actually desirable at present:
expect_that(sp$leaf_area_above(0), equals(0))

## TODO: Need a test for leaf_area_above() with non-empty species.

rm(sp)
gc()
