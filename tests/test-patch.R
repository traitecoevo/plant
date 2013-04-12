source("helper-tree.R")
options(error=traceback)

context("Patch [Plant]")

p <- new(Parameters)
p$add_strategy(new(Strategy))

## A plant that will be the same in terms of strategy (and initial
## mass).
cmp <- new(Plant, p$get_strategy(0))

patch <- new(Patch, p)

expect_that(patch$size, equals(p$size))

## We've not added any seeds yet, so this should be a list of length 1
## containing the empty list.
plants <- patch$get_plants()
expect_that(plants, is_identical_to(list(list())))

## And the height must be zero
expect_that(patch$height_max, is_identical_to(0.0))

## Add a single seed to this
patch$add_seeds(1)

plants <- patch$get_plants()
expect_that(length(plants), equals(1))
expect_that(length(plants[[1]]), equals(1))
expect_that(length(plants[[c(1,1)]]), equals(1))
expect_that(patch$n_individuals, equals(1))

expect_that(plants[[c(1,1)]]$vars_size,
            is_identical_to(cmp$vars_size))

expect_that(patch$height_max,
            is_identical_to(cmp$vars_size[["height"]]))

cmp$set_mass_leaf(pi)
patch$set_mass_leaf(pi, 0L)
expect_that(patch$get_mass_leaf(0L), is_identical_to(pi))

plants <- patch$get_plants()
expect_that(plants[[c(1,1)]]$vars_size,
            is_identical_to(cmp$vars_size))

## Compute the light environment
patch$compute_light_environment()
env <- patch$light_environment

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
expect_that(max(err), is_less_than(1e-6))

## Compute the physiological parameters:
cmp$compute_vars_phys(env)
patch$compute_vars_phys()

## And compare against the single plant.
plants <- patch$get_plants()
expect_that(plants[[c(1,1)]]$vars_phys,
            equals(cmp$vars_phys))

## One species, one individual
expect_that(patch$n_individuals, equals(1))
expect_that(patch$ode_size,      equals(3))

y <- patch$ode_values
expect_that(y, is_identical_to(c(pi, 0, 0)))
dydt <- patch$ode_rates

cmp.dydt <- unname(cmp$vars_phys[c("growth_rate", "mortality_rate",
                                   "fecundity_rate")])
expect_that(dydt, is_identical_to(cmp.dydt))

expect_that(patch$derivs(0.0, y), is_identical_to(dydt))

patch$step_deterministic()
y.new <- patch$ode_values
expect_that(all(y.new > y), is_true())

## Now, wrap up the plant as a new derivatives function:
derivs <- function(t, y, pars)
  pars$derivs(t, y)

tt <- seq(0, 15, length=51)
obj <- new(OdeR, derivs, new.env(), patch)
expect_that(obj$derivs(0.0, y),
            is_identical_to(dydt))

library(deSolve)
derivs.d <- function(...)
  list(derivs(...))
cmp.run <- unname(rk(y, tt, derivs.d, patch,
                     method=rkMethod("rk45ck"), hini=1e-8, rtol=1e-8,
                     atol=1e-8)[-1,-1])

## TODO: This is quite slow (0.1s).  This is about 20% slower than the
## rk45ck in deSolve (0.08s), so given the overhead involved in going
## to and from R, there is some serious room for improvement in the
## ODE solver.  This probably should not be too surprising as it's
## based on GSL, and that is meant to be fairly slow.
##
## It actually looks like some of the overhead could be coming from
## the R<->C communication?  Time will tell.
expect_that(t(obj$run(tt, y)),
            equals(cmp.run))

## Reset so that we are starting from the "correct" starting point for
## a single individual.
patch$clear()
patch$add_seeds(1)

set.seed(1)
while ( patch$n_individuals == 1 && patch$age < 15.0 ) {
  patch$step_deterministic()
  patch$step_stochastic()
}
expect_that(patch$n_individuals == 2 && patch$age < 15,
            is_true())

## Check that the germination numbers look correct, by comparing fit
## to a binomial distribution with the appropriate parameters.  Being
## stochastic, this may give false positives sometimes.
seed <- new(Plant, p$get_strategy(0))
pr.cmp <- p$get_parameters()[["Pi_0"]] *
  seed$germination_probability(patch$light_environment)

set.seed(1)
n <- 100000
tmp <- replicate(100, patch$germination(n))
test <- suppressWarnings(ks.test(tmp, pbinom, n, pr.cmp))
expect_that(test$p.value, is_greater_than(1/5))

rm(patch)
gc()
