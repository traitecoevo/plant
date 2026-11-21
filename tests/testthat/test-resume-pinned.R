## A run resumed from a saved state supports what a run from empty does: a
## pinned run of its own recorded schedule, an invasion against it, and the events
## still ahead of it. The burn-in of a long run is paid once and exported, so
## every later run is one of these (#519).

## A resident run to 10 and exported, then the parameters that resume it to 20.
resumed_parameters <- function() {
  p1 <- ff16_parameters(10)
  solver1 <- run(p1)$solver
  state <- export_patch_state(solver1)
  p2 <- set_initial_state(ff16_parameters(20), state)
  tt <- ff16_parameters(20)$node_schedule_times[[1]]
  p2$node_schedule_times <- list(tt[tt > state$time])
  list(p2 = p2, state = state)
}

test_that("an empty patch starts at initial_time", {
  p <- ff16_parameters(20)
  p$fitness_weighting <- "none"
  tt <- p$node_schedule_times[[1]]
  p$initial_time <- 5
  p$node_schedule_times <- list(tt[tt >= 5])
  solver <- DeterministicPatchRunner("FF16", "FF16_Env")(p, Environment("FF16"), empty_events(), Control())
  ## The clock is the parameters' before anything runs.
  expect_equal(solver$time, 5)
  solver$run()
  expect_equal(solver$patch$time, 20)
  expect_equal(solver$patch$species[[1]]$size, sum(tt >= 5))
  ## Under constant drivers and no disturbance weighting the same run shifted
  ## to start at zero gives the same stand, to the rounding of the shifted
  ## clock.
  p0 <- p
  p0$initial_time <- 0
  p0$max_patch_lifetime <- 15
  p0$node_schedule_times <- list(tt[tt >= 5] - 5)
  solver0 <- run(p0)$solver
  expect_equal(solver$patch$net_reproduction_ratios, solver0$patch$net_reproduction_ratios,
               tolerance = 1e-6)
})

test_that("a resumed run reruns its own recorded schedule", {
  r <- resumed_parameters()
  solver2 <- run(r$p2)$solver
  p3 <- r$p2
  p3$ode_times <- solver2$ode_times
  p3$ode_step_sizes <- solver2$ode_step_sizes
  ## The first recorded time is the resumed time, not zero, and the pinned run
  ## crosses the gap to the first residual introduction by the schedule.
  expect_equal(p3$ode_times[[1]], r$state$time)
  solver3 <- run(p3)$solver
  expect_identical(solver3$patch$ode_state, solver2$patch$ode_state)
  expect_identical(solver3$patch$net_reproduction_ratios, solver2$patch$net_reproduction_ratios)
})

test_that("a mutant invades a resumed resident", {
  r <- resumed_parameters()
  solver2 <- run(r$p2)$solver
  resident_times <- solver2$patch$species[[1]]$node_times
  resident_nodes <- solver2$patch$species[[1]]$net_reproduction_ratio_by_node
  ## The invader starts where the resident resumed, with no cohorts of its own
  ## and the resident's residual schedule, and is the resident's own strategy,
  ## so its cohorts must reproduce the resident's cohorts born after the resume.
  pm <- ff16_parameters(20)
  pm$initial_time <- r$state$time
  pm$node_schedule_times <- r$p2$node_schedule_times
  solver2$run_mutant(pm)
  mutant_times <- solver2$patch$species[[1]]$node_times
  mutant_nodes <- solver2$patch$species[[1]]$net_reproduction_ratio_by_node
  expect_equal(mutant_times, r$p2$node_schedule_times[[1]])
  after <- resident_times > r$state$time
  expect_equal(mutant_nodes, resident_nodes[after], tolerance = 1e-8)
})

test_that("an exported state carries the events still to come", {
  p <- ff16_parameters(20)
  ev <- events(events_default(p),
               harvest(time = 5, fraction = 0.3),
               harvest(time = 15, fraction = 0.2))
  solver <- DeterministicPatchRunner("FF16", "FF16_Env")(p, Environment("FF16"), ev, Control())
  solver$collect <- TRUE
  solver$run()
  times <- vapply(solver$history, function(h) h$time, numeric(1))
  k <- which(times >= 10)[[1]]
  state <- export_patch_state(solver, step = k)
  expect_true(all(state$events$time >= state$time - 1e-8))
  is_harvest <- state$events$type == "harvest"
  expect_equal(state$events$time[is_harvest], 15)
  ## The residual introductions in the events are the residual schedule.
  expect_equal(state$events$time[state$events$type == "node_introduction"],
               state$node_schedule_times[[1]])
  ## Resumed with its events the run reproduces the original to the solver's
  ## step-size restart; without them the harvest at 15 is lost.
  p2 <- set_initial_state(p, state)
  with_events <- run(p2, events = state$events)$solver
  expect_equal(with_events$patch$net_reproduction_ratios, solver$patch$net_reproduction_ratios,
               tolerance = 1e-3)
  without <- run(p2)$solver
  expect_gt(abs(without$patch$net_reproduction_ratios - solver$patch$net_reproduction_ratios),
            1e-2 * solver$patch$net_reproduction_ratios)
})

test_that("a resume with no introductions left is refused, not stopped at its start", {
  ## A run walks its schedule's entries; with none left it would report the
  ## exported state as the run's end, ten years short of it.
  r <- resumed_parameters()
  p2 <- r$p2
  p2$node_schedule_times <- list(numeric(0))
  expect_error(run(p2), "no introductions left")
})

test_that("an action applied at the export time is not exported again", {
  p <- ff16_parameters(15)
  ev <- events(events_default(p), harvest(time = 15, fraction = 0.2))
  solver <- run(p, events = ev)$solver
  expect_equal(solver$event_log$time, 15)
  state <- export_patch_state(solver)
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
    p <- testing_parameters("FF16")
    p$fitness_weighting <- "none"
    p$max_patch_lifetime <- lifetime
    add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 1)
  }
  p <- mk(160)
  p$node_schedule_times <- list(node_schedule_times_periodic(0, 160, 1))
  burn <- DeterministicPatchRunner("FF16", "FF16_Env")(p, Environment("FF16"), empty_events(), ctrl)
  burn$run_to(60)
  state <- export_patch_state(burn)
  sched <- node_schedule_times_periodic(60, 160, 1)
  p2 <- set_initial_state(mk(160), state)
  p2$node_schedule_times <- list(sched)
  res <- DeterministicPatchRunner("FF16", "FF16_Env")(p2, Environment("FF16"), empty_events(), ctrl)
  res$run()
  sp <- res$patch$species[[1]]
  resident <- sp$net_reproduction_ratio_by_node[sp$node_times >= 60]
  pm <- mk(160)
  pm$initial_time <- 60
  pm$node_schedule_times <- list(sched)
  res$run_mutant(pm)
  invader <- res$patch$species[[1]]$net_reproduction_ratio_by_node
  expect_gt(resident[[1]], 1)
  ## A cohort's output carries the establishment averaged over its cell. The
  ## resident's cohort born at 60 has the half interval back to its cohort born
  ## at 59 in its cell; the invader's first has no cohort before it, so its cell
  ## is the half interval after it alone. That moves its output in the sixth
  ## digit and says nothing about the canopy, which every later cohort shares.
  expect_equal(invader[-1], resident[-1], tolerance = 1e-8)
  expect_equal(invader[[1]], resident[[1]], tolerance = 1e-4)
})

test_that("an invader of a resumed resident starts empty, at the resume time", {
  ## add_mutant() on the resumed parameters keeps their initial state, which
  ## would seed the invader with the resident's exported stand under its own
  ## traits; and an invader left at time zero is not where the pinned run starts.
  ctrl <- Control(node_density_coordinate = "birth_date")
  mk <- function(lifetime) {
    p <- testing_parameters("FF16")
    p$fitness_weighting <- "none"
    p$max_patch_lifetime <- lifetime
    add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 1)
  }
  p <- mk(40)
  p$node_schedule_times <- list(node_schedule_times_periodic(0, 40, 1))
  burn <- run(p, control = ctrl, until = 0)$solver
  burn$run_to(20)
  state <- export_patch_state(burn)
  p2 <- set_initial_state(mk(40), state)
  res <- run(p2, control = ctrl)$solver
  seeded <- add_mutant(p2, trait_matrix(0.1, "lma"))
  expect_gt(length(seeded$initial_state), 0)
  expect_error(res$run_mutant(seeded), "initial state")
  empty <- mk(40)
  empty$node_schedule_times <- p2$node_schedule_times
  expect_error(res$run_mutant(empty), "initial_time")
  empty$initial_time <- state$time
  res$run_mutant(empty)
  expect_true(is.finite(res$patch$net_reproduction_ratios))
})

test_that("an invader of a resumed resident is scored in a window that opens after the resume", {
  ## The residents' first cohort after the resume takes establishment from the
  ## interval before it, from their export; an invader's first cohort cannot. A
  ## window weighting that cohort is refused. Opened at the next introduction,
  ## it weighs only cohorts the two runs share, and a copy of the resident
  ## returns the resident's ratio.
  mk <- function(start) {
    p <- fitness_window(ff16_parameters(30), start, 25)
    p$node_schedule_times <- list(node_schedule_times_periodic(0, 30, 1))
    p
  }
  ctrl <- Control(node_density_coordinate = "birth_date")
  burn <- run(mk(10), control = ctrl, until = 10)
  for (start in c(10, 11)) {
    res <- run(mk(start), control = ctrl, from = burn)
    pm <- mk(start)
    pm$node_schedule_times <- res$scenario$p$node_schedule_times
    pm$initial_time <- 10
    if (start == 10) {
      expect_error(run_mutant(res, pm), "resume time t = 10[.0]* and the window weights it")
    } else {
      expect_equal(run_mutant(res, pm)$net_reproduction_ratios, res$net_reproduction_ratios,
                   tolerance = 1e-12)
    }
  }
})
