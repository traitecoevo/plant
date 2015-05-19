context("FFW16_Cohort")

## TODO: This is all just ported over from tree1 and needs splitting
## into units.

test_that("Ported from tree1", {
  s <- FFW16_Strategy()

  plant <- FFW16_PlantPlus(s)
  cohort <- FFW16_Cohort(s)

  expect_that(cohort, is_a("FFW16_Cohort"))
  expect_that(cohort$plant, is_a("Plant<FFW16>"))

  env <- test_environment(2 * plant$height,
                          light_env=function(x) rep(1, length(x)),
                          seed_rain=1.0)

  ## The big unknown is the growth rate gradient calculation; that is,
  ## the derivative d(dh/dt)/dh.
  growth_rate_given_height <- function(height, plant, env) {
    plant$height <- height
    plant$compute_vars_phys(env)
    plant$internals[["height_dt"]]
  }
  grad_forward <- function(f, x, dx, ...) {
    (f(x + dx, ...) - f(x, ...)) / dx
  }

  plant$compute_vars_phys(env)
  p2 <- FFW16_PlantPlus(s)

  ## First, a quick sanity check that our little function behaves as
  ## expected:
  expect_that(growth_rate_given_height(plant$height, p2, env),
              equals(plant$internals[["height_dt"]]))

  ## With height:
  ctrl <- s$control
  method_args <- list(d=ctrl$cohort_gradient_eps,
                      eps=ctrl$cohort_gradient_eps)

  ## With a plant, manually compute the growth rate gradient using
  ## Richarson extrapolation:
  dgdh_richardson <- numDeriv::grad(growth_rate_given_height, plant$height,
                                    plant=p2, env=env,
                                    method.args=method_args)
  ## And also using plain forward differencing:
  dgdh_forward <- grad_forward(growth_rate_given_height, plant$height,
                               method_args$eps, plant=p2, env=env)

  ## These agree, but not that much:
  expect_that(dgdh_forward, equals(dgdh_richardson, tolerance=1e-6))

  ## Now, do this with the cohort:
  dgdh <- cohort$growth_rate_gradient(env)

  expect_that(dgdh, is_identical_to(dgdh_forward))

  ## Again with Richarson extrapolation:
  cohort2 <- FFW16_Cohort(FFW16_Strategy(control=Control(cohort_gradient_richardson=TRUE)))
  expect_that(cohort2$plant$strategy$control$cohort_gradient_richardson,
              is_true())

  ## NOTE: Not sure why this is not identical: it's either a bug
  ## somewhere or an issue due to reusing intervals.
  dgdh2 <- cohort2$growth_rate_gradient(env)
  ## expect_that(dgdh2, is_identical_to(dgdh_richardson))
  expect_that(dgdh2, not(is_identical_to(dgdh)))

  ## p <- cohort2$plant
  ## p$compute_vars_phys(env)
  ## f <- function(x) {
  ##   growth_rate_given_height(x, p, env)
  ## }
  ## dgdh3 <- test_gradient_richardson(f, p$height, ctrl$cohort_gradient_eps,
  ##                                   ctrl$cohort_gradient_richardson_depth)
  ## dgdh3 - dgdh2

  ## This is entirely optional, but kind of nice to see.
  if (interactive()) {
    hh <- seq(plant$height * 0.5, plant$height * 1.5, length.out=101)
    gr <- sapply(hh, growth_rate_given_height, p2, env)
    p2$height <- plant$height
    h_focus <- plant$height
    g_focus <- growth_rate_given_height(plant$height, p2, env)
    plot(gr ~ hh, xlab="Height", ylab="Growth rate")
    points(g_focus ~ h_focus, col="red", pch=19)
    ## Intercept by solving y = m*x + c for c => (c = y - m * x).
    abline(g_focus - dgdh * plant$height, dgdh)
  }
})

## TODO: Not done yet:
##   * Check that the initial conditions are actually correct
##   * Check that the rates computed are actually correct
test_that("ODE interface", {
  s <- FFW16_Strategy()
  plant <- FFW16_PlantPlus(s)
  cohort <- FFW16_Cohort(s)
  env <- test_environment(2 * plant$height,
                          light_env=function(x) rep(1, length(x)),
                          seed_rain=1.0)

  cohort$compute_initial_conditions(env)
  plant$compute_vars_phys(env)

  ## Different ode size to Plant, but that's not tested here (TODO)
  expect_that(cohort$ode_size, equals(4))

  ## Mortality is different because that's what Cohorts track
  v <- setdiff(names(cohort$plant$internals), "mortality")
  expect_that(cohort$plant$internals[v],
              is_identical_to(plant$internals[v]))

  ## Set up plant too:
  pr_germ <- plant$germination_probability(env)

  y <- plant$ode_state
  g <- plant$internals[["height_dt"]]

  ## Ode *values*:
  cmp <- c(plant$height,
           -log(pr_germ),
           0.0, # fecundity
           log(pr_germ * env$seed_rain_dt / g))
  expect_that(cohort$ode_state, equals(cmp))

  expect_that(cohort$fecundity, is_identical_to(0.0));

  ## Ode *rates*:
  env$time <- 10
  patch_survival <- env$patch_survival

  rates <- plant$internals

  cmp <- c(rates[["height_dt"]],
           rates[["mortality_dt"]],
           ## This is different to the approach in tree1?
           rates[["fecundity_dt"]] *
             patch_survival * exp(-cohort$plant$mortality),
           -rates[["mortality_dt"]] - cohort$growth_rate_gradient(env))
  expect_that(cohort$ode_rates, equals(cmp))
})

test_that("leaf area calculations", {
  s <- FFW16_Strategy()
  env <- test_environment(10,
                          light_env=function(x) rep(1, length(x)),
                          seed_rain=1.0)

  plant <- FFW16_PlantPlus(s)
  cohort <- FFW16_Cohort(s)
  h <- cohort$height

  expect_that(cohort$log_density,      equals(-Inf)) # zero
  expect_that(exp(cohort$log_density), equals(0.0)) # zero
  expect_that(cohort$area_leaf, equals(0)) # zero density
  cohort$compute_initial_conditions(env)
  expect_that(cohort$ode_state[[4]], equals(cohort$log_density))
  density <- exp(cohort$log_density)
  expect_that(cohort$area_leaf,
              equals(plant$area_leaf * density))
  expect_that(cohort$area_leaf_above(h / 2),
              equals(plant$area_leaf_above(h / 2) * density))

  h <- 8.0
  plant$height <- h
  v <- cohort$ode_state
  v[[1]] <- h
  cohort$ode_state <- v
  expect_that(plant$height, is_identical_to(h))
  expect_that(cohort$height, is_identical_to(h))

  expect_that(cohort$area_leaf,
              equals(plant$area_leaf * density))
  expect_that(cohort$area_leaf_above(h / 2),
              equals(plant$area_leaf_above(h / 2) * density))
})
