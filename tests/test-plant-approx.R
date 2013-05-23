source("helper-tree.R")

context("PlantApprox")

s <- new(Strategy)

plant <- new(Plant, s)
n.spline <- 50
mass.leaf.max <- 5
plant.spline <- new(PlantSpline, s, mass.leaf.max, n.spline)

approx <- new(PlantApprox, s, plant.spline)
plant <- new(Plant, s)

expect_that(approx$mass_leaf, is_identical_to(plant$mass_leaf))
expect_that(approx$vars_size, is_identical_to(plant$vars_size))

expect_that(approx$ode_size, equals(3))

h0 <- 10
plant$height <- h0
approx$height <- h0

expect_that(approx$mass_leaf, is_identical_to(plant$mass_leaf))
expect_that(approx$vars_size, is_identical_to(plant$vars_size))

expect_that(approx$vars_phys, is_identical_to(plant$vars_phys))

m0 <- plant$mass_leaf_given_height(h0)
plant$ode_values_set(c(m0, .1, .2))
approx$ode_values_set(c(m0, .1, .2))
expect_that(approx$vars_phys, is_identical_to(plant$vars_phys))

## Generate a light environment:
last <- function(x) x[[length(x)]]
hmax <- last(plant.spline$plants)$height * 1.1
env <- test.environment(hmax * 1.1)

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
expect_that(plant.spline$ode_rates(plant$mass_leaf),
            equals(plant$ode_rates, tolerance=1e-5))
expect_that(approx$ode_rates,
            is_identical_to(plant.spline$ode_rates(plant$mass_leaf)))

## Now, outside of the range of plants:
h.large <- hmax / 1.1 * 1.05
f <- function(x) {
  m <- plant$mass_leaf
  on.exit(plant$set_mass_leaf(m))
  plant$set_mass_leaf(x)
  plant$height - h.large
}

m <- uniroot(f, c(plant$mass_leaf, mass.leaf.max * 2))$root

plant$set_mass_leaf(m)
approx$set_mass_leaf(m)
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
y1 <- c(m0,   0.3,   1.5)
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
