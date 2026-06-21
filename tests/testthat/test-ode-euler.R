context("ODE fixed-step (forward Euler)")

## Forward Euler is the alternative to the adaptive Cash-Karp RKCK solver: a
## single derivative evaluation per step on a uniform grid (the way many
## industry-standard DGVMs integrate). The bare-solver behaviour (advance_euler
## on the Lorenz fixture) now lives with the solver itself in the odelia
## package; here we cover plant's use of it: (1) the SCM fixed_time_step path
## converging on the adaptive result, and (2) the guards against combining it
## with the mutant-replay path.

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
