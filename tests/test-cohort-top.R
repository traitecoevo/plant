source("helper-tree.R")
library(numDeriv)

context("CohortTop")

s <- new(Strategy)
pars.s <- s$parameters
plant <- new(Plant, s)

coh <- new(CohortTop, s)

## Set up cohort with empty light environment
h <- seq(0, 2 * plant$height, length=10)
env <- new(Spline)
env$init(h, rep(1.0, length(h)))

## Initial conditions:
seed.input <- 1
coh$compute_initial_conditions(env, seed.input)

## Different ode size to plant
expect_that(coh$ode_size, equals(4))

## Set up plant too:
plant$compute_vars_phys(env)
p.germ <- plant$germination_probability(env)

y <- plant$ode_values
g <- plant$vars_phys[["growth_rate"]]
expect_that(coh$ode_values,
            equals(c(y[1], -log(p.germ), y[3], seed.input/g)))

## First, need to compute the gradient of growth rate with respect to
## leaf mass.  This check can be removed once everything seems to be
## working as it uses mostly internal code.
p2 <- new(Plant, s)
growth.rate.given.mass <- function(mass_leaf, p, env) {
  p$set_mass_leaf(mass_leaf)
  p$compute_vars_phys(env)
  p$vars_phys[["growth_rate"]]
}

grad.forward <- function(f, x, dx, ...) {
  (f(x + dx, ...) - f(x, ...)) / dx
}

## Quick sanity check:
expect_that(growth.rate.given.mass(plant$mass_leaf, p2, env),
            equals(plant$vars_phys[["growth_rate"]]))

dgdm.accurate <- grad(growth.rate.given.mass, plant$mass_leaf,
                      p=p2, env=env) #, method.args=list(show.details=TRUE))
dgdm.simple <- grad(growth.rate.given.mass, plant$mass_leaf,
                    "simple", method.args=list(eps=1e-6),
                    p=p2, env=env)
dgdm.forward <- grad.forward(growth.rate.given.mass, plant$mass_leaf,
                             1e-6, p=p2, env=env)

dgdm <- coh$growth_rate_gradient(env)
expect_that(dgdm, is_identical_to(dgdm.forward))
expect_that(dgdm, equals(dgdm.simple))
expect_that(dgdm, equals(dgdm.accurate, tolerance=0.002))

## Then, start comparing ODE rates.
dydt <- plant$ode_rates
patch.survival <- 0.9
coh$compute_vars_phys_surv(env, patch.survival)

expect_that(coh$ode_rates,
            equals(c(dydt[1:2],
                     dydt[3] * patch.survival * coh$survival_probability,
                     dydt[2] + dgdm)))

## Delete to make sure we don't crash on cleanup
rm(coh)
gc()
