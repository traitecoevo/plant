source("helper-tree.R")

context("Species [CohortTop]")

s <- new(Strategy)
cmp <- new(CohortTop, s)

sp <- new(SpeciesCT, s)

## Note that this initialises with a special base cohort always.
expect_that(sp$size, equals(0))
expect_that(sp$n_individuals, equals(0))
expect_that(sp$height_max, is_identical_to(cmp$height))

env <- test.environment(cmp$height * 10)
env$seed_rain <- seed_rain(1.0)

## Causes initial conditions to be estimated:
sp$compute_vars_phys(env)
cmp$compute_initial_conditions(env)

sp$add_seeds(1)

expect_that(sp[[1]]$vars_phys,
            is_identical_to(cmp$vars_phys))
expect_that(sp$ode_values,
            is_identical_to(cmp$ode_values))

seed <- new(Plant, s)
expect_that(sp$germination_probability(env),
            is_identical_to(seed$germination_probability(env)))

expect_that(sp$leaf_area_above(0), equals(0))

sp$height <- 1
h <- 0
x <- c(cmp$height, sp$height)
y <- c(cmp$leaf_area_above(h), sp[[1]]$leaf_area_above(h))

expect_that(sp$leaf_area_above(h),
            is_identical_to(tree_module$trapezium(x, y)))


rm(sp)
gc()
