## A run resumed from a saved state supports what a run from empty does: a
## pinned run of its own recorded schedule, an invasion against it, and the events
## still ahead of it. The burn-in of a long run is paid once and exported, so
## every later run is one of these (#519).

ff16 <- function(lifetime) {
  p <- scm_base_parameters("FF16")
  p$max_patch_lifetime <- lifetime
  p$node_schedule_times <- list()
  add_strategies(p, trait_matrix(0.0825, "lma"))
}

## A resident run to 10 and exported, then the parameters that resume it to 20.
resumed_parameters <- function() {
  p1 <- ff16(10)
  scm1 <- run_scm(p1)
  state <- export_patch_state(scm1)
  p2 <- set_initial_state(ff16(20), state)
  tt <- ff16(20)$node_schedule_times[[1]]
  p2$node_schedule_times <- list(tt[tt > state$time])
  list(p2 = p2, state = state)
}

test_that("an empty patch starts at initial_time", {
  p <- ff16(20)
  p$patch_type <- "fixed"
  tt <- p$node_schedule_times[[1]]
  p$initial_time <- 5
  p$node_schedule_times <- list(tt[tt >= 5])
  scm <- SCM("FF16", "FF16_Env")(p, Environment("FF16"), empty_events(), Control())
  ## The clock is the parameters' before anything runs.
  expect_equal(scm$time, 5)
  scm$run()
  expect_equal(scm$patch$time, 20)
  expect_equal(scm$patch$species[[1]]$size, sum(tt >= 5))
  ## Under constant drivers and no disturbance weighting the same run shifted
  ## to start at zero gives the same stand, to the rounding of the shifted
  ## clock.
  p0 <- p
  p0$initial_time <- 0
  p0$max_patch_lifetime <- 15
  p0$node_schedule_times <- list(tt[tt >= 5] - 5)
  scm0 <- run_scm(p0)
  expect_equal(scm$net_reproduction_ratios, scm0$net_reproduction_ratios,
               tolerance = 1e-6)
})

test_that("a resumed run reruns its own recorded schedule", {
  r <- resumed_parameters()
  scm2 <- run_scm(r$p2)
  p3 <- r$p2
  p3$ode_times <- scm2$ode_times
  p3$ode_step_sizes <- scm2$ode_step_sizes
  ## The first recorded time is the resumed time, not zero, and the pinned run
  ## crosses the gap to the first residual introduction by the schedule.
  expect_equal(p3$ode_times[[1]], r$state$time)
  scm3 <- run_scm(p3)
  expect_identical(scm3$patch$ode_state, scm2$patch$ode_state)
  expect_identical(scm3$net_reproduction_ratios, scm2$net_reproduction_ratios)
})

test_that("a mutant invades a resumed resident", {
  r <- resumed_parameters()
  scm2 <- run_scm(r$p2)
  resident_times <- scm2$patch$species[[1]]$node_times
  resident_nodes <- scm2$patch$species[[1]]$net_reproduction_ratio_by_node
  ## The invader starts where the resident resumed, with no cohorts of its own
  ## and the resident's residual schedule, and is the resident's own strategy,
  ## so its cohorts must reproduce the resident's cohorts born after the resume.
  pm <- ff16(20)
  pm$initial_time <- r$state$time
  pm$node_schedule_times <- r$p2$node_schedule_times
  scm2$run_mutant(pm)
  mutant_times <- scm2$patch$species[[1]]$node_times
  mutant_nodes <- scm2$patch$species[[1]]$net_reproduction_ratio_by_node
  expect_equal(mutant_times, r$p2$node_schedule_times[[1]])
  after <- resident_times > r$state$time
  expect_equal(mutant_nodes, resident_nodes[after], tolerance = 1e-8)
})

test_that("an exported state carries the events still to come", {
  p <- ff16(20)
  ev <- events(events_default(p),
               harvest(time = 5, fraction = 0.3),
               harvest(time = 15, fraction = 0.2))
  scm <- SCM("FF16", "FF16_Env")(p, Environment("FF16"), ev, Control())
  scm$collect <- TRUE
  scm$run()
  times <- vapply(scm$history, function(h) h$time, numeric(1))
  k <- which(times >= 10)[[1]]
  state <- export_patch_state(scm, step = k)
  expect_true(all(state$events$time >= state$time - 1e-8))
  is_harvest <- state$events$type == "harvest"
  expect_equal(state$events$time[is_harvest], 15)
  ## The residual introductions in the events are the residual schedule.
  expect_equal(state$events$time[state$events$type == "node_introduction"],
               state$node_schedule_times[[1]])
  ## Resumed with its events the run reproduces the original to the solver's
  ## step-size restart; without them the harvest at 15 is lost.
  p2 <- set_initial_state(p, state)
  with_events <- run_scm(p2, events = state$events)
  expect_equal(with_events$net_reproduction_ratios, scm$net_reproduction_ratios,
               tolerance = 1e-3)
  without <- run_scm(p2)
  expect_gt(abs(without$net_reproduction_ratios - scm$net_reproduction_ratios),
            1e-2 * scm$net_reproduction_ratios)
})

test_that("a resume with no introductions left is refused, not stopped at its start", {
  ## A run walks its schedule's entries; with none left it would report the
  ## exported state as the run's end, ten years short of it.
  r <- resumed_parameters()
  p2 <- r$p2
  p2$node_schedule_times <- list(numeric(0))
  expect_error(run_scm(p2), "no introductions left")
})

test_that("an action applied at the export time is not exported again", {
  p <- ff16(15)
  ev <- events(events_default(p), harvest(time = 15, fraction = 0.2))
  scm <- run_scm(p, events = ev)
  expect_equal(scm$event_log$time, 15)
  state <- export_patch_state(scm)
  expect_equal(state$time, 15)
  expect_false(any(state$events$type == "harvest"))
  expect_equal(length(state$node_schedule_times[[1]]), 0)
})

test_that("an invader introduced at the resume instant stands in the resident's canopy", {
  ## The resident resumes at 60 and introduces its first cohort at 60 itself:
  ## that insertion is the recording's first row after its start, so no step
  ## has handed the invader the resident's field yet, and made in the invader's
  ## own empty canopy its first cohort would be far off the resident's.
  ctrl <- Control(node_density_coordinate = "birth_date")
  mk <- function(lifetime) {
    p <- scm_base_parameters("FF16")
    p$patch_type <- "fixed"
    p$max_patch_lifetime <- lifetime
    add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 1)
  }
  p <- mk(160)
  p$node_schedule_times <- list(node_schedule_times_periodic(0, 160, 1))
  burn <- SCM("FF16", "FF16_Env")(p, Environment("FF16"), empty_events(), ctrl)
  burn$run_to(60)
  state <- export_patch_state(burn)
  sched <- node_schedule_times_periodic(60, 160, 1)
  p2 <- set_initial_state(mk(160), state)
  p2$node_schedule_times <- list(sched)
  res <- SCM("FF16", "FF16_Env")(p2, Environment("FF16"), empty_events(), ctrl)
  res$run()
  sp <- res$patch$species[[1]]
  resident <- sp$net_reproduction_ratio_by_node[sp$node_times >= 60]
  pm <- mk(160)
  pm$initial_time <- 60
  pm$node_schedule_times <- list(sched)
  res$run_mutant(pm)
  invader <- res$patch$species[[1]]$net_reproduction_ratio_by_node
  expect_gt(resident[[1]], 1)
  expect_equal(invader, resident, tolerance = 1e-8)
})

test_that("an invader of a resumed resident starts empty, at the resume time", {
  ## add_mutant() on the resumed parameters keeps their initial state, which
  ## would seed the invader with the resident's exported stand under its own
  ## traits; and an invader left at time zero is not where the pinned run starts.
  ctrl <- Control(node_density_coordinate = "birth_date")
  mk <- function(lifetime) {
    p <- scm_base_parameters("FF16")
    p$patch_type <- "fixed"
    p$max_patch_lifetime <- lifetime
    add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 1)
  }
  p <- mk(40)
  p$node_schedule_times <- list(node_schedule_times_periodic(0, 40, 1))
  burn <- run_scm(p, ctrl = ctrl, run = FALSE)
  burn$run_to(20)
  state <- export_patch_state(burn)
  p2 <- set_initial_state(mk(40), state)
  res <- run_scm(p2, ctrl = ctrl)
  seeded <- add_mutant(p2, trait_matrix(0.1, "lma"))
  expect_gt(length(seeded$initial_state), 0)
  expect_error(res$run_mutant(seeded), "initial state")
  empty <- mk(40)
  empty$node_schedule_times <- p2$node_schedule_times
  expect_error(res$run_mutant(empty), "initial_time")
  empty$initial_time <- state$time
  res$run_mutant(empty)
  expect_true(is.finite(res$net_reproduction_ratios))
})
