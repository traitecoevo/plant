source("helper-tree.R")
options(error=traceback)

context("Species [Plant]")

s <- new(Strategy)
cmp <- new(Plant, s)

sp <- new(Species, s)
expect_that(sp$size, equals(0))
expect_that(sp$n_individuals, equals(0))
expect_that(sp$height_max, is_identical_to(0.0))

sp$add_seeds(1)
expect_that(sp$size, equals(1))
expect_that(sp$n_individuals, equals(1))
expect_that(sp$ode_size, equals(3))

plants <- sp$plants
expect_that(length(plants), equals(1))
expect_that(plants[[1]]$vars_size,
            is_identical_to(cmp$vars_size))

cmp$set_mass_leaf(pi)

expect_that(sp$set_mass_leaf(numeric(0)),
            throws_error())
expect_that(sp$set_mass_leaf(c(pi, pi),),
            throws_error())
sp$set_mass_leaf(pi)

expect_that(sp$plants[[1]]$vars_size,
            is_identical_to(cmp$vars_size))

## Same light environment as test-plant:
hh <- seq(0, sp$height_max, length=101)
light.env <- function(x, hmax)
  exp(x/(hmax*2)) - 1 + (1 - (exp(.5) - 1))/2
ee <- light.env(hh, max(hh))
env <- new(Spline)
env$init(hh, ee)

sp$compute_vars_phys(env)
cmp$compute_vars_phys(env)

expect_that(sp$plants[[1]]$vars_phys,
            is_identical_to(cmp$vars_phys))

seed <- new(Plant, s)
expect_that(sp$germination_probability(env),
            is_identical_to(seed$germination_probability(env)))


rm(sp)
gc()

