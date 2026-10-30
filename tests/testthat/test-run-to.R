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
  ctrl <- Control(node_density_coordinate = "birth_date", absorb_epsilon = 0.02)
  scm <- new_scm(p, ctrl)
  r <- run_until_converged(scm, block = 20, tol = 0.02)
  expect_true(r$converged)
  expect_gt(r$time, 40)
  expect_lt(r$time, 400)
  expect_false(scm$complete)
  expect_equal(scm$patch$time, r$time)
  expect_true(all(c("density", "rel_density") %in% names(r$trace)))
  ## The layer pins the tallest cohort, so by default it is not judged.
  expect_false("height_max" %in% names(r$trace))
  ## The last window of blocks is the one that settled, and the ones before it
  ## had not.
  n <- nrow(r$trace)
  expect_lt(r$trace$rel_density[[n]], 0.02)
  expect_true(all(is.na(r$trace$rel_density[1:2])))
  expect_true(all(r$trace$rel_density[3:(n - 1)] >= 0.02))
  ## And the run can be exported where it stopped.
  state <- export_patch_state(scm)
  expect_equal(state$time, r$time)
  ## A tolerance nothing meets runs to the cap and says so.
  scm2 <- new_scm(p, ctrl)
  r2 <- run_until_converged(scm2, block = 20, tol = 1e-12, t_max = 100)
  expect_false(r2$converged)
  expect_gte(r2$time, 100)
})

test_that("run_to() with collect on keeps the same history as run()", {
  p <- ff16(15)
  ev <- events(node_introductions(p), harvest(time = 7.5, fraction = 0.3))
  a <- run_scm(p, events = ev, run = FALSE)
  a$collect <- TRUE
  a$run()
  b <- run_scm(p, events = ev, run = FALSE)
  b$collect <- TRUE
  for (t in c(3, 7.5, 12, 15)) {
    b$run_to(t)
  }
  expect_identical(length(b$history), length(a$history))
  expect_identical(b$history[[length(b$history)]]$ode_state,
                   a$history[[length(a$history)]]$ode_state)
  expect_identical(b$event_log, a$event_log)
})

test_that("an export where run_to() stopped keeps the action due there", {
  ## run_to() stops on arrival at an instant, before the instant's actions run,
  ## so a harvest due at the stop is still to come in the export.
  p <- ff16(15)
  ev <- events(node_introductions(p), harvest(time = 7.5, fraction = 0.5))
  straight <- run_scm(p, events = ev)
  s <- run_scm(p, events = ev, run = FALSE)
  s$run_to(7.2)
  expect_equal(s$time, 7.5)
  expect_length(s$event_log$time, 0)
  state <- export_patch_state(s)
  expect_true("harvest" %in% state$events$type)
  resumed <- run_scm(set_initial_state(p, state), events = state$events)
  expect_equal(resumed$event_log$type, "harvest")
  expect_equal(resumed$net_reproduction_ratios, straight$net_reproduction_ratios,
               tolerance = 1e-3)
})

test_that("run_scm(run = FALSE) returns an SCM to drive, and only that", {
  p <- ff16(10)
  s <- run_scm(p, run = FALSE)
  expect_false(s$complete)
  expect_equal(s$time, 0)
  expect_error(run_scm(p, run = FALSE, collect = TRUE), "unrun")
  expect_error(run_scm(p, run = FALSE, refine_schedule = TRUE), "unrun")
})

test_that("run_until_converged() refuses what it cannot judge", {
  p <- ff16(100)
  p$patch_type <- "fixed"
  p$node_schedule_times <- list(node_schedule_times_periodic(0, 100, 1))
  ctrl <- Control(node_density_coordinate = "birth_date", absorb_epsilon = 0.02)
  expect_error(run_until_converged(run_scm(p, ctrl = ctrl)), "already finished")
  s <- function() run_scm(p, ctrl = ctrl, run = FALSE)
  expect_error(run_until_converged(s(), block = 0), "block")
  expect_error(run_until_converged(s(), tol = -1), "tol")
  expect_error(run_until_converged(s(), t_max = 0), "after")
  ## FF16's environment carries no state: by default it is left out, asked for
  ## it is refused, and its entries are refused without it.
  expect_error(run_until_converged(s(), statistics = "environment"), "no state")
  expect_error(run_until_converged(s(), statistics = "density",
                                   environment_entries = 1), "environment")
  ## A last block shorter than `block` still runs the patch on to t_max.
  r <- run_until_converged(s(), block = 30, tol = 1e-12, t_max = 70,
                           statistics = "density")
  expect_false(r$converged)
  expect_gte(r$time, 70)
})

test_that("a statistic that is undefined in a block cannot settle, and says so", {
  ## A strategy whose seedlings make no net carbon never establishes, so no
  ## cohort is alive and the tallest living cohort is undefined in every block.
  p <- scm_base_parameters("FF16")
  p$patch_type <- "fixed"
  p$max_patch_lifetime <- 60
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 1)
  p$strategies[[1]]$pars$r_l <- 1e6
  p$node_schedule_times <- list(node_schedule_times_periodic(0, 60, 1))
  s <- run_scm(p, ctrl = Control(node_density_coordinate = "birth_date"), run = FALSE)
  expect_warning(r <- run_until_converged(s, block = 20, tol = 0.1,
                                          statistics = c("density", "height_max")),
                 "height_max")
  expect_false(r$converged)
})

test_that("a statistic is judged on its own scale, however small", {
  ## A stand too sparse to shade itself is the same stand at any birth rate,
  ## scaled, so whether it has settled cannot depend on the birth rate. Judged
  ## against a floor of machine epsilon, a stand at 1e-300 read as settled two
  ## blocks in, while its density was still rising by 15 % a block.
  converge <- function(b, ...) {
    p <- ff16(200)
    p$patch_type <- "fixed"
    p$strategies[[1]]$birth_rate_y <- b
    p$node_schedule_times <- list(node_schedule_times_periodic(0, 200, 1))
    s <- run_scm(p, ctrl = Control(node_density_coordinate = "birth_date"), run = FALSE)
    run_until_converged(s, block = 20, tol = 0.02, t_max = 120,
                        statistics = "density", ...)
  }
  tiny <- converge(1e-300)
  small <- converge(1e-150)
  expect_lt(max(tiny$trace$density), 1e-280)
  expect_equal(tiny$trace$rel_density, small$trace$rel_density, tolerance = 1e-10)
  expect_equal(tiny$time, small$time)
  expect_gt(tiny$trace$rel_density[[4]], 0.1)
  expect_gt(tiny$time, 80)
  ## A floor in absolute terms is the caller's to set, and then it counts.
  floored <- converge(1e-300, tol_abs = 1e-200)
  expect_true(floored$converged)
  expect_equal(floored$time, 60)
})

test_that("statistics that settle by construction are not judged by default", {
  ## The top layer pins the tallest cohort: asked for, it warns.
  p <- ff16(100)
  p$patch_type <- "fixed"
  p$node_schedule_times <- list(node_schedule_times_periodic(0, 100, 1))
  ctrl <- Control(node_density_coordinate = "birth_date", absorb_epsilon = 0.02)
  expect_warning(run_until_converged(run_scm(p, ctrl = ctrl, run = FALSE),
                                     block = 20, t_max = 40,
                                     statistics = c("density", "height_max")),
                 "pinned")
  expect_error(run_until_converged(run_scm(p, ctrl = ctrl, run = FALSE), window = 1),
               "window")
  ## TF24's cumulative flux accumulators grow without bound: by default only the
  ## soil water is judged.
  p <- scm_base_parameters("TF24")
  p$patch_type <- "fixed"
  p$max_patch_lifetime <- 4
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = TF24_hyperpar,
                      birth_rate = 1)
  p$node_schedule_times <- list(0:3)
  s <- run_scm(p, ctrl = Control(node_density_coordinate = "birth_date",
                                 ode_step_size_max = 0.05), run = FALSE)
  n_env <- s$patch$environment$ode_size
  n_acc <- length(s$patch$state$env$soil_moist_cumulative_flux)
  expect_gt(n_acc, 0)
  r <- run_until_converged(s, block = 1, tol = 1e-12, t_max = 3)
  env <- grep("^env_", names(r$trace), value = TRUE)
  expect_identical(env, paste0("env_", seq_len(n_env - n_acc)))
})

test_that("run_to() refuses a time it cannot compare", {
  s <- run_scm(ff16(10), run = FALSE)
  expect_error(s$run_to(NaN), "NaN")
})
