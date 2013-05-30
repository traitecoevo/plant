source("helper-tree.R")
library(numDeriv)

context("CohortTop")

s <- new(Strategy)
pars.s <- s$parameters
plant <- new(Plant, s)

coh <- new(CohortTop, s)

## Set up cohort with empty light environment
env <- test.environment(2 * plant$height,
                        light.env=function(x) rep(1, length(x)))
env$seed_rain <- seed_rain(1.0)

## Initial conditions:
coh$compute_initial_conditions(env)

## Different ode size to plant
expect_that(coh$ode_size, equals(4))

## Set up plant too:
plant$compute_vars_phys(env)
p.germ <- plant$germination_probability(env)

y <- plant$ode_values
g <- plant$vars_phys[["height_growth_rate"]]
expect_that(coh$ode_values,
            equals(c(y[1], -log(p.germ), y[3], env$seed_rain_rate/g)))

## First, need to compute the gradient of height growth rate with
## respect to height.  This check can be removed once everything seems
## to be working as it uses mostly internal code.
p2 <- new(Plant, s)
growth.rate.given.height <- function(height, p, env) {
  p$height <- height
  p$compute_vars_phys(env)
  p$vars_phys[["height_growth_rate"]]
}

grad.forward <- function(f, x, dx, ...) {
  (f(x + dx, ...) - f(x, ...)) / dx
}

## Quick sanity check:
expect_that(growth.rate.given.height(plant$height, p2, env),
            equals(plant$vars_phys[["height_growth_rate"]]))

## With height:
ctrl <- coh$control
method.args <- list(d=ctrl$parameters$cohort_gradient_eps,
                    eps=ctrl$parameters$cohort_gradient_eps)
dgdm.accurate <- grad(growth.rate.given.height, plant$height,
                      p=p2, env=env, method.args=method.args)
dgdm.simple <- grad(growth.rate.given.height, plant$height,
                    "simple", method.args=method.args,
                    p=p2, env=env)
dgdm.forward <- grad.forward(growth.rate.given.height, plant$height,
                             method.args$eps, p=p2, env=env)

dgdm <- coh$growth_rate_gradient(env)
expect_that(dgdm, is_identical_to(dgdm.forward))
expect_that(dgdm, equals(dgdm.simple))
expect_that(dgdm, equals(dgdm.accurate, tolerance=0.002))

## Check with Richardson extrapolation (requires making a whole new
## Cohort object, as the Control slot is readonly).
ctrl$set_parameters(list(cohort_gradient_richardson=1))
s2 <- new(Strategy)
s2$control <- ctrl
c2 <- new(CohortTop, s2)
dgdm2 <- c2$growth_rate_gradient(env)
expect_that(dgdm2, equals(dgdm.forward, tolerance=0.002))
expect_that(dgdm2, equals(dgdm.simple, tolerance=0.002))
expect_that(dgdm2, equals(dgdm.accurate))

## Then, start comparing ODE rates.
dydt <- plant$ode_rates
coh$compute_vars_phys(env)

env$time <- 10
patch.survival <- env$patch_survival(0)

expect_that(coh$ode_rates,
            equals(c(dydt[1:2],
                     dydt[3] * patch.survival * coh$survival_probability,
                     dydt[2] + dgdm)))

## Delete to make sure we don't crash on cleanup
rm(coh)
gc()
