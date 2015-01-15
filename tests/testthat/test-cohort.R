if (interactive()) {
  devtools::load_all("../../")
  library(testthat)
  source("helper-tree2.R")
}

context("Cohort")

test_that("Creation", {
  s <- Strategy()

  plant <- Plant(s)
  cohort <- Cohort(s)

  expect_that(cohort, is_a("Cohort"))
  expect_that(cohort$plant, is_a("Plant"))

  env <- test_environment(2 * plant$height,
                          light_env=function(x) rep(1, length(x)),
                          seed_rain=1.0)

  ## The big unknown is the growth rate gradient calculation; that is,
  ## the derivative d(dh/dt)/dh.
  growth_rate_given_height <- function(height, plant, env) {
    plant$height <- height
    plant$compute_vars_phys(env)
    plant$vars_phys[["height_growth_rate"]]
  }
  grad_forward <- function(f, x, dx, ...) {
    (f(x + dx, ...) - f(x, ...)) / dx
  }

  plant$compute_vars_phys(env)
  p2 <- Plant(s)

  ## First, a quick sanity check that our little function behaves as
  ## expected:
  expect_that(growth_rate_given_height(plant$height, p2, env),
              equals(plant$vars_phys[["height_growth_rate"]]))

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
  cohort2 <- Cohort(Strategy(control=Control(cohort_gradient_richardson=TRUE)))
  expect_that(cohort2$plant$strategy$control$cohort_gradient_richardson,
              is_true())

  dgdh2 <- cohort2$growth_rate_gradient(env)
  expect_that(dgdh2, is_identical_to(dgdh_richardson))
  expect_that(dgdh2, not(is_identical_to(dgdh)))

  ## This is entirely optional, but kind of nice to see.
  if (interactive()) {
    hh <- seq(plant$height * 0.5, plant$height * 1.5, length=101)
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
