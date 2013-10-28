source("helper-tree.R")
options(error=traceback)

context("Species [Plant]")

s <- new(Strategy)
cmp <- new(Plant, s)

sp <- new(Species, s)
expect_that(sp$size, equals(0))
expect_that(sp$n_individuals, equals(0))
expect_that(sp$height_max, is_identical_to(cmp$height))

sp$add_seeds(1)
expect_that(sp$size, equals(1))
expect_that(sp$n_individuals, equals(1))
expect_that(sp$ode_size, equals(3))

plants <- sp$plants
expect_that(length(plants), equals(1))
expect_that(plants[[1]]$vars_size,
            is_identical_to(cmp$vars_size))

## Check that access via "[[" works.
expect_that(sp[[1]]$vars_size,
            is_identical_to(cmp$vars_size))
expect_that(sp[[2]]$vars_size,
            throws_error())

## Height to set things to.
h0 <- 10

cmp$height <- h0

expect_that(sp$height <- numeric(0),
            throws_error())
expect_that(sp$height <- c(h0, h0),
            throws_error())
sp$height <- h0

expect_that(sp$plants[[1]]$vars_size,
            is_identical_to(cmp$vars_size))
expect_that(sp[[1]]$vars_size,
            is_identical_to(cmp$vars_size))

env <- test.environment(sp$height_max)

sp$compute_vars_phys(env)
cmp$compute_vars_phys(env)

expect_that(sp$plants[[1]]$vars_phys,
            is_identical_to(cmp$vars_phys))

seed <- new(Plant, s)
expect_that(sp$germination_probability(env),
            is_identical_to(seed$germination_probability(env)))

## Bizarrely, this includes about 100 points; more than I'd have
## thought, and not sure why.  Could be the error in assimilation
## calculation?  Could just be that we're being a bit enthusiastic
## about error checking...
sp$compute_assimilation_spline(env)
tmp <- sp$assimilation_spline
xy <- tmp$xy

## Test approximate plant:

ctrl <- new(Control,
            list(plant_assimilation_approximate_use=TRUE))
s.approx <- new(Strategy)
s.approx$control <- ctrl
sp.approx <- new(Species, s.approx)

sp.approx$add_seeds(1)
sp.approx$height <- sp$height

## Spline is empty on initialisation
expect_that(sp.approx$assimilation_spline$size, equals(0))

sp.approx$compute_vars_phys(env)

## Check that we did compute the spline:
expect_that(sp.approx$assimilation_spline$xy, is_identical_to(tmp$xy))

## The rates should be equal, but not identical, to the fully computed
## rates.
expect_that(sp.approx$ode_rates, equals(sp$ode_rates))
expect_that(identical(sp.approx$ode_rates, sp$ode_rates),
            is_false())

rm(sp)
gc()
