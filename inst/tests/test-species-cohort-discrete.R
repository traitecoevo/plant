source("helper-tree.R")
options(error=traceback)

context("Species [CohortDiscrete]")

s <- new(Strategy)
cmp <- new(Plant, s)

sp <- new(SpeciesC, s)
expect_that(sp$size, equals(0))
expect_that(sp$n_individuals, equals(0))
expect_that(sp$height_max, is_identical_to(cmp$height))

sp$add_seeds(2)
expect_that(sp$size, equals(1))
expect_that(sp$n_individuals, equals(2))
expect_that(sp$ode_size, equals(3))

expect_that(sp$size, equals(1))
expect_that(sp[[1]]$vars_size,
            is_identical_to(cmp$vars_size))
expect_that(length(sp$plants), equals(1))

h0 <- 10
cmp$height <- h0

expect_that(sp$height <- numeric(0),
            throws_error())
expect_that(sp$height <- c(h0, h0),
            throws_error())
sp$height <- h0

expect_that(sp[[1]]$vars_size,
            is_identical_to(cmp$vars_size))

env <- test.environment(sp$height_max)

sp$compute_vars_phys(env)
cmp$compute_vars_phys(env)

expect_that(sp[[1]]$vars_phys,
            is_identical_to(cmp$vars_phys))

seed <- new(Plant, s)
expect_that(sp$germination_probability(env),
            is_identical_to(seed$germination_probability(env)))

rm(sp)
gc()
