# Initialising a patch from an exported state (revival of #304).
#
# The core guarantee: a patch run can be exported mid/late-run and re-imported
# to seed a fresh run, which then reproduces the original trajectory. The seed
# is loaded bit-for-bit; only forward adaptive-solver stepping (whose step-size
# state is not checkpointed) diverges, hence the tolerance on the resumed
# fitness rather than exact equality.

context("Initial patch state (export / re-import)")

## Run an SCM, export the patch at a mid-history step, and run a fresh SCM
## seeded from that export.
run_and_resume <- function(p, x, e, frac = 0.5) {
  env <- Environment(x)
  ctrl <- Control()

  scm <- SCM(x, e)(p, env, ctrl)
  scm$collect <- TRUE
  scm$run()

  k <- max(2L, as.integer(floor(length(scm$history) * frac)))
  state <- export_patch_state(scm, step = k)

  p2 <- set_initial_state(p, state)
  scm2 <- SCM(x, e)(p2, env, ctrl)
  scm2$collect <- TRUE
  scm2$run()

  list(scm = scm, scm2 = scm2, state = state, p2 = p2, k = k)
}

test_that("FF16 single-species round-trip reproduces the run", {
  x <- "FF16"; e <- "FF16_Env"
  p0 <- scm_base_parameters(x)
  p1 <- expand_parameters(trait_matrix(0.08, "lma"), p0, birth_rate_list = 1.0)

  r <- run_and_resume(p1, x, e)
  state <- r$state

  # The seeded patch (first collected snapshot of the resumed run) must equal
  # the exported state exactly.
  seed <- r$scm2$history[[1]]
  expect_equal(seed$time, state$time)
  expect_equal(seed$ode_state, state$ode_state)
  expect_equal(lapply(seed$species, function(s) s$node_times), state$node_times)
  expect_equal(lapply(seed$species, function(s) s$pr_patch_survival_at_birth),
               state$pr_patch_survival)
  expect_equal(as.integer(vapply(seed$species, function(s) s$size, numeric(1))), state$n)

  # The resumed run reproduces the original fitness (small residual is solver
  # step-size divergence across the checkpoint, not a state error).
  expect_equal(r$scm2$net_reproduction_ratios,
               r$scm$net_reproduction_ratios, tolerance = 1e-3)

  # Re-importing the same state is deterministic.
  env <- Environment(x); ctrl <- Control()
  scm3 <- SCM(x, e)(r$p2, env, ctrl); scm3$run()
  expect_identical(scm3$net_reproduction_ratios, r$scm2$net_reproduction_ratios)
})

test_that("Resumed schedule excludes already-imported introductions", {
  x <- "FF16"; e <- "FF16_Env"
  p0 <- scm_base_parameters(x)
  p1 <- expand_parameters(trait_matrix(0.08, "lma"), p0, birth_rate_list = 1.0)

  r <- run_and_resume(p1, x, e)
  # nodes seeded + future introductions == all of the original nodes
  n_seeded <- sum(r$state$n)
  n_future <- sum(vapply(r$state$node_schedule_times, length, integer(1)))
  n_total_orig <- sum(vapply(r$scm$patch$species, function(s) s$size, numeric(1)))
  expect_equal(n_seeded + n_future, n_total_orig)

  # all residual introductions are at or after the export time (the node due at
  # exactly the export time has not been seeded yet, so it stays in the schedule)
  expect_true(all(unlist(r$state$node_schedule_times) >= r$state$time - 1e-8))
  # ...and none coincides with an already-seeded node's introduction time
  seeded_times <- unlist(lapply(r$scm2$history[[1]]$species,
                                function(s) s$node_times))
  expect_true(all(unlist(r$state$node_schedule_times) > max(seeded_times)))
})

test_that("K93 multi-species round-trip reproduces the run", {
  x <- "K93"; e <- "K93_Env"
  p0 <- scm_base_parameters(x)
  p0$max_patch_lifetime <- 35.10667

  # three species varying b_0 (growth scaling); other K93 traits at defaults
  sp <- trait_matrix(c(0.059, 0.063, 0.052), "b_0")
  p2 <- expand_parameters(sp, p0, birth_rate_list = list(20, 20, 20))

  r <- run_and_resume(p2, x, e, frac = 0.5)

  expect_equal(length(r$state$n), 3)
  expect_equal(r$scm2$net_reproduction_ratios,
               r$scm$net_reproduction_ratios, tolerance = 1e-3)
})

test_that("make_initial_state seeds a patch from a size distribution", {
  x <- "FF16"; e <- "FF16_Env"
  p0 <- scm_base_parameters(x)
  p1 <- expand_parameters(trait_matrix(0.08, "lma"), p0, birth_rate_list = 1.0)
  env <- Environment(x); ctrl <- Control()

  heights <- seq(1, 8, length.out = 10)
  state <- make_initial_state(p1, heights = heights, densities = rep(0.2, 10))

  # 10 nodes, flat ODE state = 10 * node_ode_size (FF16 = 7), no environment ODE
  expect_equal(state$n, 10L)
  expect_equal(length(state$ode_state), 10L * 7L)
  expect_equal(state$time, 0)
  # all introduced at age 0, with a non-zero pr_patch_survival at birth (so the
  # fecundity rate does not divide by zero)
  expect_true(all(unlist(state$node_times) == 0))
  expect_true(all(unlist(state$pr_patch_survival) > 0))

  p2 <- set_initial_state(p1, state)
  seeded <- SCM(x, e)(p2, env, ctrl)
  # the seeded patch holds the requested nodes, ordered by decreasing height
  expect_equal(seeded$patch$species[[1]]$size, 10)
  expect_equal(sort(seeded$patch$species[[1]]$heights, decreasing = TRUE),
               seeded$patch$species[[1]]$heights)

  out <- run_scm(p2, env, ctrl, collect = TRUE)
  expect_true(is.finite(out$net_reproduction_ratios))
  expect_gt(out$net_reproduction_ratios, 0)
})

test_that("implausibly dense initial conditions are rejected", {
  x <- "FF16"; e <- "FF16_Env"
  p0 <- scm_base_parameters(x)
  p1 <- expand_parameters(trait_matrix(0.08, "lma"), p0, birth_rate_list = 1.0)
  env <- Environment(x); ctrl <- Control()

  # many large plants at very high density -> exploding density rates
  state <- make_initial_state(p1, heights = seq(8, 12, length.out = 20),
                              densities = rep(50, 20))
  expect_error(SCM(x, e)(set_initial_state(p1, state), env, ctrl),
               "non-finite densities")
})

test_that("seeded runs work with tidy outputs and plot_size_distribution", {
  x <- "FF16"
  p0 <- scm_base_parameters(x)
  p1 <- expand_parameters(trait_matrix(0.08, "lma"), p0, birth_rate_list = 1.0)

  st  <- make_initial_state(p1, heights = seq(1, 8, length.out = 10),
                            densities = rep(0.2, 10))
  res <- run_scm(set_initial_state(p1, st), collect = TRUE)

  # the standard tidy `species` table, with the seeded nodes present at t = 0
  expect_true(all(c("time", "height", "density", "log_density", "node",
                    "species") %in% names(res$species)))
  expect_true(any(res$species$time == 0 & !is.na(res$species$density)))

  # the standard plotting and aggregation helpers consume it unchanged
  expect_s3_class(plot_size_distribution(res$species), "ggplot")
  totals <- integrate_over_size_distribution(FF16_expand_state(res)$species)
  expect_true(all(is.finite(totals$area_leaf)))
  expect_true(all(is.finite(totals$density)))
})
