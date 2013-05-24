source("helper-tree.R")

context("PlantApprox")

s <- new(Strategy)

plant <- new(Plant, s)
n.spline <- 50
height.max <- 10
plant.spline <- new(PlantSpline, s, height.max, n.spline)

approx <- new(PlantApprox, s, plant.spline)
plant <- new(Plant, s)

expect_that(approx$height,    is_identical_to(plant$height))
expect_that(approx$vars_size, is_identical_to(plant$vars_size))

expect_that(approx$ode_size, equals(3))

h0 <- 10
plant$height <- h0
approx$height <- h0

expect_that(approx$height,    is_identical_to(plant$height))
expect_that(approx$vars_size, is_identical_to(plant$vars_size))

expect_that(approx$vars_phys, is_identical_to(plant$vars_phys))

plant$ode_values_set(c(h0, .1, .2))
approx$ode_values_set(c(h0, .1, .2))
expect_that(approx$vars_phys, is_identical_to(plant$vars_phys))

## Generate a light environment:
env <- test.environment(height.max * 1.1)

plant$compute_vars_phys(env)
## Update the underlying spline (a controlling class would normally do
## this).
approx$compute_vars_phys_spline(env)
## Then perhaps update the physiology for this individual.
approx$compute_vars_phys(env)
plant.spline$compute_vars_phys(env)

expect_that(all(approx$vars_phys == 0.0), is_true())

## The values should be unchanged.
expect_that(approx$ode_values, is_identical_to(plant$ode_values))

## Check against the plant spline....
expect_that(plant.spline$ode_rates(plant$height),
            equals(plant$ode_rates))
## This one does not work?
expect_that(approx$ode_rates,
            is_identical_to(plant.spline$ode_rates(plant$height)))

## Now, outside of the range of plants:
h.large <- height.max * 1.05
plant$height <- h.large
approx$height <- h.large

expect_that(approx$vars_size, is_identical_to(plant$vars_size))

plant$compute_vars_phys(env)
approx$compute_vars_phys(env)
plant.spline$compute_vars_phys(env)

expect_that(approx$vars_phys, is_identical_to(plant$vars_phys))
expect_that(approx$ode_rates, is_identical_to(plant$ode_rates))
expect_that(plant.spline$ode_rates(m), throws_error())

f.p <- function(obj, y) {
  obj$ode_values_set(y)
  obj$died()
}

##      Size, death, birth
y1 <- c(h0,   0.3,   1.5)
nrep <- 100
set.seed(1)
d.p <- replicate(nrep, f.p(plant, y1))
set.seed(1)
d.a <- replicate(nrep, f.p(approx, y1))
expect_that(d.a, is_identical_to(d.p))

plant$ode_values_set(y1)
expect_that(plant$offspring(), equals(y1[3] %/% 1))
expect_that(plant$ode_values[3], equals(y1[3] %% 1))
approx$ode_values_set(y1)
expect_that(approx$offspring(), equals(y1[3] %/% 1))
expect_that(approx$ode_values[3], equals(y1[3] %% 1))

rm(approx)
gc()
