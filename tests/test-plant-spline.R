source("helper-tree.R")

context("PlantSpline")

s <- new(Strategy)

plant <- new(Plant, s)
n.spline <- 50
height.max <- 10
plant.spline <- new(PlantSpline, s, height.max, n.spline)

## Expected leaf mass:
expect_that(plant.spline$height_max,
            is_identical_to(height.max))

height <- seq(plant$height, plant.spline$height_max,
              length=n.spline)

## Check that the plants are correctly spaced:
expect_that(sapply(plant.spline$plants, function(x) x$height),
            equals(height))

env <- test.environment(height.max)

## Get all the plants to estimate their physiological variables:
plant.spline$compute_vars_phys(env)

## Then extract the plants:
plants <- plant.spline$plants

## And the full spline:
spline <- plant.spline$plants_approx

## Now, the work starts.  First, check that the light environment was
## set correctly in each plant.  We'll create a *new* vector of plants
## with the same

## Check that the rates look OK between the spline and the underlying
## plants.
expect_that(t(sapply(plants, function(x) x$ode_rates)),
            is_identical_to(spline$y))

## And check against a manual solution:
f <- function(h, env, plant) {
  plant$height <- h
  plant$compute_vars_phys(env)
  plant$ode_rates
}

expect_that(spline$y,
            equals(t(sapply(height, f, env, plant))))

## Now, let's see how this responds to getting rates for intermediate
## values:
hmid <- (height[-length(height)] + height[-1])/2

rates.exact  <- t(sapply(hmid, f, env, plant))
rates.spline <- spline$eval(hmid)

## This is not that exact, because the splines are poorly estimated.
## However, this is the sort of thing we could build a spline up on at
## the beginning with full sunlight and hope that that does not change
## over the simulation?
expect_that(rates.spline,
            equals(rates.exact, tolerance=2e-3))

## TODO: Temporary until ode basis changes...
mass <- sapply(height, plant$mass_leaf_given_height)
mmid <- sapply(hmid, plant$mass_leaf_given_height)

expect_that(t(sapply(mass, function(m) plant.spline$ode_rates(m))),
            equals(spline$y))
expect_that(t(sapply(mmid, function(m) plant.spline$ode_rates(m))),
            equals(rates.spline))

expect_that(plant.spline$ode_rates(plant$mass_leaf_given_height(plant.spline$height_max) + 1e-8),
            throws_error())
