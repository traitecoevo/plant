## Driving a run from outside (#519): run_to() between legs, the thinned
## snapshot record, and the burn-in that stops when the stand has settled.

ff16 <- function(lifetime) {
  p <- scm_base_parameters("FF16")
  p$max_patch_lifetime <- lifetime
  p$node_schedule_times <- list()
  add_strategies(p, trait_matrix(0.0825, "lma"))
}

new_scm <- function(p, ctrl = Control()) {
  SCM("FF16", "FF16_Env")(p, Environment("FF16"), empty_events(), ctrl)
}

test_that("a sequence of run_to() calls is a run, leg for leg", {
  p <- ff16(20)
  a <- new_scm(p)
  a$run()
  b <- new_scm(p)
  for (t in c(5, 10, 15, 20)) {
    b$run_to(t)
  }
  expect_true(b$complete)
  expect_identical(b$patch$ode_state, a$patch$ode_state)
  expect_identical(b$ode_times, a$ode_times)
  expect_identical(b$net_reproduction_ratios, a$net_reproduction_ratios)
  ## It stops at the first leg boundary at or past the time asked for, and the
  ## patch it exposes is the live one.
  c <- new_scm(p)
  c$run_to(7.3)
  expect_gte(c$time, 7.3)
  expect_false(c$complete)
  expect_equal(c$patch$time, c$time)
  expect_true(c$time %in% p$node_schedule_times[[1]])
  ## A finished run is left where it is.
  c$run_to(20)
  t_end <- c$time
  c$run_to(30)
  expect_equal(c$time, t_end)
  ## And run() starts over.
  c$run()
  expect_identical(c$patch$ode_state, a$patch$ode_state)
})

test_that("snapshots thin the record and keep every node's lifetime output", {
  p <- ff16(20)
  scm <- new_scm(p)
  scm$snapshot_interval <- 5
  scm$run()
  s <- scm$snapshots
  times <- vapply(s, `[[`, numeric(1), "time")
  ## One at the start, then the first leg boundary at or past 5, 10, 15, 20.
  expect_equal(floor(times / 5), 0:4)
  expect_true(all(times >= c(0, 5, 10, 15, 20)))
  last <- s[[length(s)]]$species[[1]]
  sp <- scm$patch$species[[1]]
  expect_equal(last$n_nodes, sp$size)
  expect_equal(last$fecundity, sp$net_reproduction_ratio_by_node)
  trapezium <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  expect_equal(last$density, trapezium(sp$node_times, exp(sp$log_densities_state)))
  expect_equal(last$height_max, max(sp$heights[exp(sp$log_densities_state) > 0]))
  expect_equal(last$n_retired, 0L)
  ## FF16 carries no environment state.
  expect_length(last$environment, 0)
  tidy <- tidy_snapshots(scm)
  expect_equal(nrow(tidy$species), 5)
  expect_equal(tidy$species$time, times)
  expect_equal(names(tidy$environment), "time")
  ## A run with no interval keeps no snapshots, and run() starts the record over.
  scm$snapshot_interval <- 0
  scm$run()
  expect_length(scm$snapshots, 0)
})

test_that("run_until_converged stops when the block means settle", {
  p <- ff16(400)
  p$patch_type <- "fixed"
  p$node_schedule_times <- list(node_schedule_times_periodic(0, 400, 1))
  ctrl <- Control(node_density_in_birth_date = TRUE, retire_density = 1e-8)
  scm <- new_scm(p, ctrl)
  r <- run_until_converged(scm, block = 20, tol = 0.02)
  expect_true(r$converged)
  expect_gt(r$time, 40)
  expect_lt(r$time, 400)
  expect_false(scm$complete)
  expect_equal(scm$patch$time, r$time)
  expect_true(all(c("density", "height_max", "rel_density") %in% names(r$trace)))
  ## The last two blocks are the ones that settled.
  n <- nrow(r$trace)
  expect_true(all(r$trace$rel_density[c(n - 1, n)] < 0.02))
  ## And the run can be exported where it stopped.
  state <- export_patch_state(scm)
  expect_equal(state$time, r$time)
  ## A tolerance nothing meets runs to the cap and says so.
  scm2 <- new_scm(p, ctrl)
  r2 <- run_until_converged(scm2, block = 20, tol = 1e-12, t_max = 100)
  expect_false(r2$converged)
  expect_gte(r2$time, 100)
})
