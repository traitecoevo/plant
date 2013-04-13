source("helper-tree.R")

context("PlantSpline")

s <- new(Strategy)

plant <- new(Plant, s)
n.spline <- 50
plant.spline <- new(PlantSpline, s, n.spline)

## Expected leaf mass:
mass.leaf <- seq(plant$mass_leaf, 5, length=n.spline)

## Check that the plants are correctly spaced:
expect_that(sapply(plant.spline$plants, function(x) x$mass_leaf),
            equals(mass.leaf))

## Same light environment as test-plant.R:
last <- function(x) x[[length(x)]]
hmax <- last(plant.spline$plants)$height
hh <- seq(0, hmax, length=101)
light.env <- function(x, hmax)
  exp(x/(max(hmax)*2)) - 1 + (1 - (exp(.5) - 1))/2
ee <- light.env(hh, max(hh))
env <- new(Spline)
env$init(hh, ee)

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
f <- function(m, env, plant) {
  plant$set_mass_leaf(m)
  plant$compute_vars_phys(env)
  plant$ode_rates
}

expect_that(spline$y,
            equals(t(sapply(mass.leaf, f, env, plant))))

## Now, let's see how this responds to getting rates for intermediate
## values:
m <- (mass.leaf[-length(mass.leaf)] + mass.leaf[-1])/2

rates.exact  <- t(sapply(m, f, env, plant))
rates.spline <- spline$eval(m)

## This is not that exact, because the splines are poorly estimated.
## However, this is the sort of thing we could build a spline up on at
## the beginning with full sunlight and hope that that does not change
## over the simulation?
expect_that(rates.spline,
            equals(rates.exact, tolerance=2e-3))
