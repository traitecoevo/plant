context("ODE fixed-step (forward Euler)")

## Forward Euler is the alternative to the adaptive Cash-Karp RKCK solver: a
## single derivative evaluation per step on a uniform grid (the way many
## industry-standard DGVMs integrate). These tests cover (1) the bare solver via
## the OdeRunner, (2) the SCM fixed_time_step path converging on the adaptive
## result, and (3) the guards against combining it with the mutant-replay path.

derivs_lorenz <- function(y, pars) {
  c(pars[[1]] * (y[[2]] - y[[1]]),
    pars[[2]] * y[[1]] - y[[2]] - y[[1]] * y[[3]],
    -pars[[3]] * y[[3]] + y[[1]] * y[[2]])
}

test_that("advance_euler matches a hand-rolled forward Euler", {
  pars <- c(sigma = 10.0, R = 28.0, b = 8.0 / 3.0)
  y0 <- c(21, 21, 21)
  times <- seq(0, 1, by = 0.001)

  lo <- Lorenz(pars[[1]], pars[[2]], pars[[3]])
  lo$ode_state <- y0
  sys <- OdeRunner("Lorenz")(lo)
  sys$advance_euler(times)

  ## Reference: explicit forward Euler in R over the same grid.
  y <- y0
  for (i in 2:length(times)) {
    h <- times[[i]] - times[[i - 1]]
    y <- y + h * derivs_lorenz(y, pars)
  }

  expect_equal(sys$state, y, tolerance = 1e-12)
  expect_equal(sys$times, times)
})

test_that("advance_euler does plain Euler, not the RKCK step", {
  ## advance_fixed drives the full 6-stage RKCK step; advance_euler does one
  ## derivative evaluation. Over a coarse grid the two must therefore differ,
  ## and advance_euler must equal a single explicit Euler update.
  pars <- c(sigma = 10.0, R = 28.0, b = 8.0 / 3.0)
  y0 <- c(21, 21, 21)
  times <- c(0, 0.5, 1.0)

  lo <- Lorenz(pars[[1]], pars[[2]], pars[[3]])
  lo$ode_state <- y0
  euler <- OdeRunner("Lorenz")(lo)
  euler$advance_euler(times)

  lo2 <- Lorenz(pars[[1]], pars[[2]], pars[[3]])
  lo2$ode_state <- y0
  rkck <- OdeRunner("Lorenz")(lo2)
  rkck$advance_fixed(times)

  ## Single explicit Euler update by hand.
  y <- y0
  for (i in 2:length(times)) {
    h <- times[[i]] - times[[i - 1]]
    y <- y + h * derivs_lorenz(y, pars)
  }

  expect_equal(euler$state, y, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(euler$state, rkck$state)))
})

test_that("advance_euler validates its time grid", {
  lo <- Lorenz(10.0, 28.0, 8.0 / 3.0)
  sys <- OdeRunner("Lorenz")(lo)
  expect_error(sys$advance_euler(numeric(0)),
               "must be vector of at least length 1")
  expect_error(sys$advance_euler(c(1.0, 2.0)),
               "First element in 'times' must be same as current time")
})

test_that("SCM fixed-step Euler converges on the adaptive result", {
  x <- "FF16"
  p <- expand_parameters(trait_matrix(0.08, "lma"),
                         scm_base_parameters(x), birth_rate_list = 1.0)
  env <- Environment(x)

  ref <- run_scm(p, env, control_accurate())$net_reproduction_ratios

  dts <- c(0.5, 0.25, 0.125)
  err <- vapply(dts, function(dt) {
    r0 <- run_scm(p, env, Control(fixed_time_step = dt))$net_reproduction_ratios
    abs(r0 - ref)
  }, numeric(1))

  ## Error shrinks monotonically as the step shrinks (first-order convergence).
  expect_true(all(diff(err) < 0))
  ## Roughly halving the step roughly halves the error (allow generous slack).
  ratios <- err[-length(err)] / err[-1]
  expect_true(all(ratios > 1.4 & ratios < 3.0))

  ## A fine step lands close to the adaptive reference.
  fine <- run_scm(p, env, Control(fixed_time_step = 0.02))$net_reproduction_ratios
  expect_equal(fine, ref, tolerance = 0.03)
})

test_that("SCM fixed-step walks the expected uniform grid", {
  x <- "FF16"
  p <- expand_parameters(trait_matrix(0.08, "lma"),
                         scm_base_parameters(x), birth_rate_list = 1.0)
  env <- Environment(x)

  dt <- 0.25
  scm <- run_scm(p, env, Control(fixed_time_step = dt))
  ode_times <- scm$ode_times

  ## Introductions (schedule event times) are the only points where the grid
  ## restarts; between them every interval is <= dt (the last sub-interval of
  ## each window may be shorter).
  steps <- diff(ode_times)
  steps <- steps[steps > 0] # drop the exact-coincidence joins between windows
  expect_true(all(steps <= dt + 1e-8))

  ## Finer steps mean strictly more evaluations.
  scm_fine <- run_scm(p, env, Control(fixed_time_step = dt / 2))
  expect_gt(length(scm_fine$ode_times), length(ode_times))
})

test_that("fixed_time_step is rejected on the mutant-replay paths", {
  x <- "FF16"
  p <- expand_parameters(trait_matrix(0.08, "lma"),
                         scm_base_parameters(x), birth_rate_list = 1.0)
  env <- Environment(x)

  ## save_RK45_cache (the resident pass that feeds mutant fitness) has no Euler
  ## analogue: refused at construction.
  expect_error(
    run_scm(p, env, Control(fixed_time_step = 0.5, save_RK45_cache = TRUE)),
    "incompatible with save_RK45_cache")

  ## Pinned ode-time replay (use_ode_times) is refused mid-run.
  p$ode_times <- run_scm(p, env, Control())$ode_times
  expect_error(
    run_scm(p, env, Control(fixed_time_step = 0.5), use_ode_times = TRUE),
    "not supported for ode-time replay")
})
