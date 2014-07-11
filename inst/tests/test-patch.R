source("helper-tree.R")

context("Patch [Plant]")

p <- new(Parameters)
p$add_strategy(new(Strategy))

## A plant that will be the same in terms of strategy (and initial
## mass).
cmp <- new(Plant, p[[1]])

patch <- new(Patch, p)

test_that("Parameters can be pulled from Patch", {
  p.cmp <- patch$parameters
  expect_that(p.cmp$size, is_identical_to(p$size))
  expect_that(p.cmp$parameters, is_identical_to(p$parameters))
  expect_that(p.cmp[[1]]$parameters, is_identical_to(p[[1]]$parameters))
})

expect_that(patch$size, equals(p$size))

## We've not added any seeds yet, so the only species should have no
## plants.
expect_that(patch$size, equals(1))
expect_that(patch[[1]]$plants, is_identical_to(list()))

## Check that alternative access approach via $species works:
expect_that(length(patch$species), equals(1))
expect_that(patch$species[[1]]$plants, is_identical_to(list()))

## And the height must be zero
expect_that(patch$height_max, is_identical_to(cmp$height))

## Add a single seed to this
patch$add_seedlings(1)

expect_that(patch$size, equals(1)) # no change

expect_that(patch[[1]]$size, equals(1))
expect_that(patch[[1]]$n_individuals, equals(1))
expect_that(patch[[1]][[1]]$vars_size,
            is_identical_to(cmp$vars_size))
expect_that(patch$height_max,
            is_identical_to(cmp$height))

h0 <- 10
cmp$height <- h0
patch$height <- list(h0)
expect_that(patch$height, is_identical_to(list(h0)))

plants <- patch[[1]]$plants
expect_that(plants[[1]]$vars_size,
            is_identical_to(cmp$vars_size))

## Compute the light environment
patch$compute_light_environment()
env <- patch$environment

## And compare that with the manually computed light environment
target <- function(x) sapply(x, patch$canopy_openness)
xx.eval <- env$light_environment$xy[,1]
expect_that(env$light_environment$xy[,2],
            is_identical_to(target(xx.eval)))

## And then check that the error is under control
xx.mid <- (xx.eval[-1] + xx.eval[-length(xx.eval)]) / 2
yy.mid <- target(xx.mid)
zz.mid <- env$light_environment$eval(xx.mid)
expect_that(zz.mid, equals(yy.mid, tolerance=2e-8))

err <- pmax(abs(zz.mid - yy.mid), abs(1 - zz.mid / yy.mid))
expect_that(max(err), is_less_than(1e-6))

## Compute the physiological parameters:
cmp$compute_vars_phys(env)
patch$compute_vars_phys()

## And compare against the single plant.
plants <- patch[[1]]$plants
expect_that(plants[[1]]$vars_phys,
            equals(cmp$vars_phys))

## One species, one individual
ode_size_plant <- new(Plant,new(Strategy))$ode_size

expect_that(patch$n_individuals, equals(1))
expect_that(patch$ode_size,      equals(ode_size_plant))

y <- patch$ode_values
expect_that(y, equals(c(h0, 0, 0,0,0)))
dydt <- patch$ode_rates

cmp.dydt <- unname(c(cmp$vars_phys[c("height_growth_rate",
                                   "mortality_rate", "fecundity_rate")],
                    cmp$vars_growth_decomp[c("dheartwood_area_dt","dheartwood_mass_dt")]))
expect_that(dydt, is_identical_to(cmp.dydt))

expect_that(patch$derivs(0.0, y), is_identical_to(dydt))

solver <- solver_from_ode_target(patch, p$control$ode_control)

expect_that(solver$derivs(patch$time, patch$ode_values),
            is_identical_to(dydt))

solver$step()
y.new <- patch$ode_values
expect_that(all(y.new > y), is_true())

## Compare running with the main OdeR solver and one from deSolve
## (why?)
tt <- seq(0, 15, length=51)

library(deSolve)
derivs.d <- function(t, y, pars)
  list(pars$derivs(t, y))

## Note that this is twice as slow as the deSolve rk method!  Some of
## this (perhaps 20%) is due to overhead with the R->C->R->C
## convoluted calls, especially with the interactions with S4 methods.
## The rest is probably tied up with the extra parameter set that we
## always do with the GSL ode solver.
cmp.run <- unname(rk(y, tt, derivs.d, patch,
                     method=rkMethod("rk45ck"), hini=1e-8, rtol=1e-8,
                     atol=1e-8)[-1,-1])
expect_that(t(solver$run(tt, y)),
            equals(cmp.run))

## Check a single individual added via add_seedling:
patch$reset()
expect_that(patch$n_individuals, equals(0))
expect_that(patch$add_seedling(patch$size() + 1), throws_error())
expect_that(patch$add_seedling(0), throws_error())
patch$add_seedling(1)
expect_that(patch$n_individuals, equals(1))

## Reset so that we are starting from the "correct" starting point for
## a single individual.
patch$reset()
patch$add_seedlings(1)
solver$reset()
solver$set_state(patch$ode_values, patch$time)

## Helper function to move through the stochastic parts of the life
## cycle (deaths and births)
step.stochastic <- function(patch) {
  patch$deaths()
  patch$add_seeds(patch$births())
}

enough.time <- 30
set.seed(1)
while (patch$n_individuals == 1 && patch$time < enough.time) {
  solver$step()
  step.stochastic(patch)
  solver$set_state(patch$ode_values, patch$time)
}

test_that("Patch growth matched expectation", {
  expect_that(patch$n_individuals, equals(2))
  expect_that(patch$time, is_within_interval(0.0, enough.time))
})

## Test the patch rescaling when changing height:
ctrl.new <- list(environment_light_rescale_usually=TRUE)
p$set_control_parameters(ctrl.new)
patch.scale <- new(Patch, p)
patch.scale$add_seedlings(1)
patch.scale$add_seedlings(1)
patch.scale$set_ode_values(patch$time, patch$ode_values)
patch.scale$compute_light_environment()

## Then, add new plants and force a rescale of the light environment
env1 <- patch.scale$environment$light_environment
h <- patch.scale$height
r <- 1.1
h[] <- lapply(h, function(x) x * r)

## Setting new heights will cause the light environment to be scaled
## (but not recalculated by default).
patch.scale$height <- h
env2 <- patch.scale$environment$light_environment

test_that("Light environment is scaled", {
  expect_that(nrow(env2), is_identical_to(nrow(env1)))
  expect_that(range(env2$xy[,1]), is_identical_to(c(0, max(unlist(h)))))
  expect_that(env2$xy[,1], equals(env1$xy[,1] * r))
  expect_that(sapply(env2$xy[,1], patch.scale$canopy_openness),
              is_identical_to(env2$xy[,2]))
})

patch.scale$compute_light_environment()
env3 <- patch.scale$environment$light_environment
test_that("Recomputed light environment would be different", {
  expect_that(nrow(env3$xy), is_greater_than(nrow(env2$xy)))
})

## Check that the germination numbers look correct, by comparing fit
## to a binomial distribution with the appropriate parameters.  Being
## stochastic, this may give false positives sometimes.
seed <- new(Plant, p[[1]])
pr.cmp <- p$parameters[["Pi_0"]] *
  seed$germination_probability(patch$environment)

set.seed(1)
n <- 100000
tmp <- replicate(100, patch$germination(n))
test <- suppressWarnings(ks.test(tmp, pbinom, n, pr.cmp))
expect_that(test$p.value, is_greater_than(1/5))

test_that("State get/set works", {
  patch2 <- new(Patch, p)
  env <- patch2$environment
  state <- patch2$state

  state.cmp <- list(species=list(new(Species, p[[1]])$state),
                    environment=env$state)
  expect_that(state, is_identical_to(state.cmp))

  ## Add a bunch of seedlings:
  set.seed(1)
  n <- 10
  for (i in seq_len(n))
    patch2$add_seedlings(1)
  h <- sort(runif(n, seed$height, 15), decreasing=TRUE)
  t <- runif(1)

  patch2$height <- list(h)
  patch2$time   <- t

  state <- patch2$state
  state.cmp <- list(species=list(patch2[[1]]$state),
                    environment=patch2$environment$state)
  expect_that(state, is_identical_to(state.cmp))

  state$species[[1]][2:3,] <- runif(n * 2)
  patch2$state <- state

  ## Check that this forces the light environment to be built:
  patch3 <- new(Patch, patch2$parameters)
  patch3$force_state(state)
  expect_that(patch3$state, is_identical_to(state))
  expect_that(patch3$environment$light_environment$xy,
              is_identical_to(patch2$environment$light_environment$xy))

  ## Check some invalid states fail:
  expect_that(patch2$state <- list(), throws_error())
  expect_that(patch2$state <- unname(state), throws_error())
  expect_that(patch2$state <- state[1], throws_error())
  expect_that(patch2$state <- state[2], throws_error())
  state.broken <- state
  state.broken$species[[1]] <- state.broken$species[[1]][,-4]
  expect_that(patch2$state <- state.broken, throws_error())

  ## Second check with completely different state (redundant?)
  patch3$force_state(patch$state)
  expect_that(patch3$state, is_identical_to(patch$state))
  expect_that(patch3$environment$light_environment$xy,
              is_identical_to(patch$environment$light_environment$xy))
})


expect_that(patch$disturbance_regime, is_a("Rcpp_Disturbance"))

test_that("Non-resident individuals do not contribute to environment", {
  s1 <- new(Strategy, list(lma=0.06))
  s2 <- new(Strategy, list(lma=0.20))

  set.seed(1)
  h <- list(sort(runif(15, 1, 20), decreasing=TRUE),
            sort(runif(18, 1, 20), decreasing=TRUE))

  p <- new(Parameters)
  p$add_strategy(s1)
  p$add_strategy(s2)
  expect_that(p$is_resident, is_identical_to(rep(TRUE, 2)))

  p1 <- new(Parameters)
  p1$add_strategy(s1)

  p2 <- new(Parameters)
  p2$add_strategy(s2)

  setup <- function(p, h) {
    patch <- new(Patch, p)
    for (i in seq_along(h))
      for (j in seq_along(h[[i]]))
        patch$add_seedling(i)
    patch$height <- h
    patch
  }

  patch.both <- setup(p, h)
  patch.1 <- setup(p1, h[1])
  patch.2 <- setup(p2, h[2])

  p$is_resident <- c(TRUE, FALSE) # Both species, but 1 as resident
  patch.1.mut <- setup(p, h)

  p$is_resident <- c(FALSE, TRUE) # Both species but 2 as resident
  patch.2.mut <- setup(p, h)

  expect_that(patch.1$height_max,
              is_identical_to(patch.1.mut$height_max))
  expect_that(patch.2$height_max,
              is_identical_to(patch.2.mut$height_max))

  light.env <- function(patch)
    patch$environment$light_environment$xy
  expect_that(light.env(patch.1),
              is_identical_to(light.env(patch.1.mut)))
  expect_that(light.env(patch.2),
              is_identical_to(light.env(patch.2.mut)))
})

rm(patch)
gc()
