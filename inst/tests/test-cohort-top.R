source("helper-tree.R")

context("CohortTop")

s <- new(Strategy)
pars.s <- s$parameters
plant <- new(Plant, s)

coh <- new(CohortTop, s)

## Set up cohort with empty light environment
env <- test.environment(2 * plant$height,
                        light.env=function(x) rep(1, length(x)),
                        seed.rain=1.0)

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
            equals(c(y[1], -log(p.germ), y[3],
                     log(p.germ * env$seed_rain_rate/g))))

test_that("State get/set works", {
  c2 <- new(CohortTop, s)
  c2$compute_initial_conditions(env) # time = 0, pr_surv = 1
  expect_that(c2$state_size, equals(5))
  expect_that(c2$state, is_identical_to(c(c2$ode_values, 1)))

  env$time <- 10
  c2$compute_initial_conditions(env) # time > 0, pr_surv < 1
  expect_that(c2$state, is_identical_to(c(c2$ode_values, env$patch_survival)))

  x <- c2$state + runif(5)
  c2$state <- x
  expect_that(c2$state, is_identical_to(x))
  expect_that(c2$state, is_identical_to(c(c2$ode_values, x[5])))

  expect_that(c2$state <- x[-1],   throws_error())
  expect_that(c2$state <- c(x, 1), throws_error())

  expect_that(c2$fecundity, is_identical_to(x[[3]]))
})

expect_that(coh$leaf_area,
            equals(exp(coh$state[[4]]) *
                   coh$vars_size[["leaf_area"]]))

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
ctrl <- s$control
method.args <- list(d=ctrl$parameters$cohort_gradient_eps,
                    eps=ctrl$parameters$cohort_gradient_eps)
dgdh.accurate <- numDeriv::grad(growth.rate.given.height, plant$height,
                                p=p2, env=env, method.args=method.args)
dgdh.simple <- numDeriv::grad(growth.rate.given.height, plant$height,
                              "simple", method.args=method.args,
                              p=p2, env=env)
dgdh.forward <- grad.forward(growth.rate.given.height, plant$height,
                             method.args$eps, p=p2, env=env)

dgdh <- coh$growth_rate_gradient(env)
expect_that(dgdh, is_identical_to(dgdh.forward))
expect_that(dgdh, equals(dgdh.simple))
expect_that(dgdh, equals(dgdh.accurate, tolerance=0.002))

## Not necessary (and not run as part of the tests), but it's nice to
## see how this actually lines up with the data:
if (interactive()) {
  hh <- seq(plant$height*.5, plant$height*1.5, length=101)
  gr <- sapply(hh, growth.rate.given.height, p2, env)
  p2$height <- plant$height
  h.focus <- plant$height
  g.focus <- growth.rate.given.height(plant$height, p2, env)
  plot(gr ~ hh, xlab="Height", ylab="Growth rate")
  points(g.focus ~ h.focus, col="red", pch=19)
  ## Intercept by solving y = m*x + c for c => (c = y - m * x).
  abline(g.focus - dgdh * plant$height, dgdh)
}

## Check with Richardson extrapolation (requires making a whole new
## Cohort object, as the Control slot is readonly).
ctrl$set_parameters(list(cohort_gradient_richardson=1))
s2 <- new(Strategy)
s2$control <- ctrl
c2 <- new(CohortTop, s2)
dgdh2 <- c2$growth_rate_gradient(env)
expect_that(dgdh2, equals(dgdh.forward, tolerance=0.002))
expect_that(dgdh2, equals(dgdh.simple, tolerance=0.002))
expect_that(dgdh2, equals(dgdh.accurate))

## Then, start comparing ODE rates.
dydt <- plant$ode_rates
coh$compute_vars_phys(env)

env$time <- 10
patch.survival <- env$patch_survival

expect_that(coh$ode_rates,
            equals(c(dydt[1:2],
                     dydt[3] * patch.survival * coh$survival_probability,
                     -dydt[2] - dgdh)))

test_that("Approximate assimilation works", {
  ctrl.full <- new(Control)
  ctrl.full$set_parameters(list(cohort_gradient_richardson=1))
  s.full <- new(Strategy)
  s.full$control <- ctrl.full

  ctrl.approx <- new(Control)
  ctrl.approx$set_parameters(list(plant_assimilation_approximate_use=TRUE))
  s.approx <- new(Strategy)
  s.approx$control <- ctrl.approx

  assim <- compute_assimilation_fn(s.approx, plant$height,
                                   env$light_environment$max, env)
  expect_that(assim$max, is_identical_to(env$light_environment$max))
  s.approx$assimilation_fn <- assim

  ## Then set this within the strategy so that it would be available to
  ## a plant.
  c.full   <- new(CohortTop, s.full)
  c.approx <- new(CohortTop, s.approx)

  env$time <- 0
  c.full$compute_initial_conditions(env)
  c.approx$compute_initial_conditions(env)

  env$time <- 10
  h <- (plant$height + env$light_environment$max) / 2
  c.full$height <- h
  c.approx$height <- h

  c.full$compute_vars_phys(env)
  c.approx$compute_vars_phys(env)

  ## Size should be exactly the same:
  expect_that(c.approx$vars_size, is_identical_to(c.full$vars_size))

  ## Expect that the physiological variables will be similar, but not
  ## exactly the same.
  expect_that(c.approx$vars_phys, equals(c.full$vars_phys))
  expect_that(identical(c.approx$vars_phys, equals(c.full$vars_phys)),
              is_false())

  r.full <- c.full$ode_rates
  r.approx <- c.approx$ode_rates
  expect_that(r.full[1:3], equals(r.approx[1:3]))
  expect_that(r.full[4],   equals(r.approx[4], tolerance=1e-6))
})


## Delete to make sure we don't crash on cleanup
rm(coh)
gc()
