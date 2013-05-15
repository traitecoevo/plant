source("helper-tree.R")

context("Species [CohortTop]")

s <- new(Strategy)
cmp <- new(Plant, s)

sp <- new(SpeciesCT, s)

## Note that this initialises with a special base cohort always.
expect_that(sp$size, equals(1))
expect_that(sp$n_individuals, equals(1))
expect_that(sp$height_max, is_identical_to(cmp$height))

## Same light environment as test-plant:
hh <- seq(0, cmp$height * 10, length=101)
light.env <- function(x, hmax)
  exp(x/(hmax*2)) - 1 + (1 - (exp(.5) - 1))/2
ee <- light.env(hh, max(hh))
env <- new(Spline)
env$init(hh, ee)

sp$compute_vars_phys(env)
cmp$compute_vars_phys(env)

expect_that(sp[[1]]$vars_phys,
            is_identical_to(cmp$vars_phys))

seed <- new(Plant, s)
expect_that(sp$germination_probability(env),
            is_identical_to(seed$germination_probability(env)))

rm(sp)
gc()
