source("helper-tree.R")
options(error=traceback)

context("Patch")

p <- new(Parameters)
p$add_strategy(list())

## A plant that will be the same in terms of strategy (and initial
## mass).
cmp <- new(Plant, p$get_strategy(0))

patch <- new(Patch, p)

expect_that(patch$size, equals(p$size))

## We've not added any seeds yet, so this should be the empty list
plants <- patch$get_plants(0L)
expect_that(plants, is_identical_to(list()))

## And the height must be zero
expect_that(patch$height_max, is_identical_to(0.0))

## Add a single seed to this
patch$add_seed(0)

plants <- patch$get_plants(0L)
expect_that(length(plants), equals(1))
expect_that(length(plants[[1]]), equals(1))

expect_that(plants[[1]]$vars_size,
            is_identical_to(cmp$vars_size))

expect_that(patch$height_max,
            is_identical_to(cmp$vars_size[["height"]]))

cmp$set_mass_leaf(pi)
patch$set_mass_leaf(pi, 0L)
expect_that(patch$get_mass_leaf(0L), is_identical_to(pi))

plants <- patch$get_plants(0L)
expect_that(plants[[1]]$vars_size,
            is_identical_to(cmp$vars_size))

## Compute the light environment
patch$compute_light_environment()
env <- patch$get_light_environment()

## And compare that with the manually computed light environment
target <- function(x) sapply(x, patch$canopy_openness)
xx.eval <- env$xy[,1]
expect_that(env$xy[,2],
            is_identical_to(target(xx.eval)))

## And then check that the error is under control
xx.mid <- (xx.eval[-1] + xx.eval[-length(xx.eval)]) / 2
yy.mid <- target(xx.mid)
zz.mid <- env$eval(xx.mid)
expect_that(zz.mid, equals(yy.mid, tolerance=2e-8))

err <- pmax(abs(zz.mid - yy.mid), abs(1 - zz.mid / yy.mid))
expect_that(all(err < 1e-6), is_true())

## Compute the physiological parameters:
cmp$compute_vars_phys(env)
patch$compute_vars_phys()

## And compare against the single plant.
plants <- patch$get_plants(0L)
expect_that(plants[[1]]$vars_phys,
            equals(cmp$vars_phys))

## One species, one individual
expect_that(patch$ode_size, equals(3))

y <- patch$ode_values()
expect_that(y, is_identical_to(c(pi, 0, 0)))
dydt <- patch$ode_rates()

cmp.dydt <- unname(cmp$vars_phys[c("growth_rate", "mortality_rate",
                                   "fecundity_rate")])
expect_that(dydt, is_identical_to(cmp.dydt))

expect_that(patch$derivs(y), is_identical_to(dydt))

patch$step_deterministic()
y.new <- patch$ode_values()
expect_that(all(y.new > y), is_true())

rm(patch)
gc()
