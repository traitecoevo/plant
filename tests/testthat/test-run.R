## run(): the one verb (#428). Each path is checked against the run_scm() call
## that does the same thing, so the new surface adds no numerics of its own.

outputs <- function(x) {
  scm <- if (inherits(x, "plant_run")) x$solver else x
  list(ode_state = scm$patch$ode_state,
       offspring_production = scm$offspring_production,
       net_reproduction_ratios = scm$net_reproduction_ratios)
}

test_that("run(p) is run_scm(p), and reports what it ran", {
  p <- ff16_parameters(20)
  res <- run(p)
  expect_s3_class(res, "plant_run")
  expect_identical(outputs(res), outputs(run_scm(p)))
  expect_true(res$complete)
  expect_identical(res$offspring_production, res$solver$offspring_production)
  expect_s3_class(res$scenario, "plant_scenario")
  expect_identical(res$scenario$p$max_patch_lifetime, 20)
  expect_null(res$mutants)
})

test_that("a scenario and its pieces run the same, and a result reruns", {
  p <- ff16_parameters(20)
  ctrl <- control(node_density_coordinate = "birth_date")
  ev <- events(node_introductions(p), harvest(time = 10, fraction = 0.3))
  a <- run(scenario(p, control = ctrl, events = ev))
  b <- run(p, control = ctrl, events = ev)
  expect_identical(outputs(a), outputs(b))
  expect_identical(outputs(a), outputs(run_scm(p, ctrl = ctrl, events = ev)))
  expect_identical(outputs(run(a$scenario)), outputs(a))
  expect_equal(a$event_log$type, b$event_log$type)
  expect_true("harvest" %in% a$event_log$type)
})

test_that("refinement through run() is run_scm()'s, and its schedule is the result's", {
  p <- ff16_parameters(20)
  res <- run(p, refine_schedule = TRUE)
  scm <- run_scm(p, refine_schedule = TRUE)
  expect_identical(outputs(res), outputs(scm))
  expect_identical(res$scenario$p, scm$parameters)
})

test_that("record keeps what tidy() and a gradient read", {
  p <- ff16_parameters(10)
  res <- run(p, record = "steps")
  collected <- run_scm(p, collect = TRUE)
  expect_identical(tidy(res), collected[c("steps", "n_spp", "species", "env")])
  expect_error(tidy(run(p)), "record = \"steps\"")
  expect_true(run(p, record = "states")$solver$keep_states)
  expect_false(run(p)$solver$keep_states)
})

test_that("a run stopped at `until` and continued is the run taken whole", {
  p <- ff16_parameters(20)
  part <- run(p, until = 8)
  expect_false(part$complete)
  expect_gte(part$time, 8)
  expect_lt(part$time, 20)
  whole <- run(part)
  expect_true(whole$complete)
  expect_identical(outputs(whole), outputs(run(p)))
  ## `part` described the run at t = 8; its solver has moved on since.
  expect_false(part$complete)
  expect_error(run(part), "moved on")
  expect_error(export_patch_state(part), "moved on")
  expect_error(run(whole, until = 5), "past `until`")
  ## Continuing to where a run stands takes no step, so an unrun one stays unrun.
  z <- run(p, until = 0)
  z2 <- run(z, until = 0)
  expect_identical(z2$time, 0)
  expect_identical(run_until_converged(z2, block = 5)$trace,
                   run_until_converged(run_scm(p, run = FALSE), block = 5)$trace)
})

test_that("until at the start builds the run unrun, as run_scm(run = FALSE) does", {
  p <- ff16_parameters(20)
  res <- run(p, until = 0)
  expect_identical(res$time, 0)
  expect_false(res$complete)
  ## So a convergence run from it starts with the t = 0 snapshot.
  a <- run_until_converged(res, block = 5)
  scm <- run_scm(p, run = FALSE)
  b <- run_until_converged(scm, block = 5)
  expect_identical(a$trace, b$trace)
  expect_identical(tidy_snapshots(a$run), tidy_snapshots(scm))
  expect_s3_class(a$run, "plant_run")
  expect_error(export_patch_state(res), "moved on")
  expect_identical(export_patch_state(a$run), export_patch_state(res$solver))
})

test_that("a resume is set_initial_state() with the state's pending events", {
  p <- ff16_parameters(20)
  p$patch_type <- "fixed"
  ctrl <- control(node_density_coordinate = "birth_date")
  ev <- events(node_introductions(p), harvest(time = 15, fraction = 0.3))
  part <- run(p, control = ctrl, events = ev, until = 8)
  st <- export_patch_state(part)
  resumed <- run(p, control = ctrl, from = part)
  manual <- run_scm(set_initial_state(p, st), ctrl = ctrl, events = st$events)
  expect_identical(outputs(resumed), outputs(manual))
  expect_identical(resumed$scenario$p$initial_time, st$time)
  expect_true("harvest" %in% resumed$event_log$type)
  expect_identical(outputs(run(p, control = ctrl, from = st)), outputs(resumed))
  ## The events are the state's: a scenario that brings its own is refused.
  expect_error(run(p, control = ctrl, events = ev, from = st), "none of its own")
})

test_that("a resume into a longer scenario keeps its introductions", {
  short <- ff16_parameters(20)
  short$patch_type <- "fixed"
  long <- ff16_parameters(40)
  long$patch_type <- "fixed"
  ## A harvest pending at the export, so the introductions travel as events.
  ev <- events(node_introductions(short), harvest(time = 15, fraction = 0.3))
  part <- run(short, events = ev, until = 8)
  resumed <- run(long, from = part)
  expect_true(resumed$complete)
  expect_true("harvest" %in% resumed$event_log$type)
  ## What the solver introduced: the long scenario's times from the export on.
  made <- resumed$solver$patch$species[[1]]$node_times
  wanted <- long$node_schedule_times[[1]]
  expect_identical(made[made > part$time + 1e-8], wanted[wanted > part$time + 1e-8])
  expect_gt(max(made), 20)
})

test_that("the functions that took an SCM take a result", {
  p <- ff16_parameters(20)
  res <- run(p)
  expect_identical(export_patch_state(res), export_patch_state(res$solver))
})

test_that("run_mutant() is the solver's own invasion, and retires the resident's result", {
  p <- ff16_parameters(20)
  pm <- add_mutant(p, trait_matrix(0.1, "lma"))

  scm <- run_scm(p)
  scm$run_mutant(pm)
  res <- run(p)
  inv <- run_mutant(res, pm)
  expect_s3_class(inv, "plant_run")
  expect_identical(inv$offspring_production, scm$offspring_production)
  expect_identical(inv$mutants, pm)
  expect_identical(inv$scenario, res$scenario)
  expect_output(print(inv), "1 species invading")
  inv3 <- run_mutant(run(p), add_mutant(p, trait_matrix(c(0.09, 0.1, 0.11), "lma"),
                                           birth_rate = as.list(rep(1, 3))))
  expect_output(print(inv3), "3 species invading")
  expect_error(export_patch_state(res), "moved on")
  expect_error(tidy(res), "moved on")
  expect_error(run_mutant(inv, pm), "resident")

  scm <- run_scm(p)
  scm$run_mutant_constant(pm)
  inv <- run_mutant(run(p), pm, constant = TRUE)
  expect_identical(inv$offspring_production, scm$offspring_production)

  expect_error(run_mutant(run(p, until = 5), pm), "not finished")
  expect_error(run_mutant(scm, pm), "plant_run")
})

test_that("run() refuses what it cannot do", {
  p <- ff16_parameters(5)
  expect_error(run(p, refine_schedule = TRUE, until = 2), "until")
  expect_error(run(p, until = c(1, 2)), "one time")
  expect_error(run(p, record = "everything"), "takes")
  expect_error(run(scenario(p), env = Environment("FF16")), "does not take env")
  expect_error(run(scenario(p), 3), "unnamed")
  expect_error(run(run(p, until = 2), record = "steps"), "until")
  expect_error(scenario(list()), "Parameters")
  expect_error(scenario(p, events = list()), "Events")
  expect_error(run(p, from = run_scm(p)), "export_patch_state")
})

## The stochastic solver and a single plant under run(). Each is the function
## run() calls, so each is asserted identical to it.

stochastic_stand <- function(lifetime = 20) {
  p <- testing_parameters("FF16")
  p$max_patch_lifetime <- lifetime
  add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 2)
}

test_that("a stochastic run is run_stochastic_collect()'s, under the same seed", {
  p <- stochastic_stand()
  set.seed(42)
  collected <- run_stochastic_collect(p)
  set.seed(42)
  res <- run(p, method = "stochastic", record = "steps")
  expect_s3_class(res, "plant_run")
  expect_identical(res$method, "stochastic")
  expect_identical(tidy(res), collected[c("time", "species", "env")])
  expect_true(res$complete)
  expect_null(res$offspring_production)
  ## The scenario records the arrivals drawn, not p's deterministic schedule.
  expect_identical(res$scenario$p, res$solver$parameters)
  expect_false(identical(res$scenario$p$node_schedule_times, p$node_schedule_times))
  expect_output(print(res), "stochastic")
  ## Unrecorded, it draws the same individuals and ends in the same state.
  set.seed(42)
  quiet <- run(p, method = "stochastic")
  expect_identical(quiet$solver$patch$state, res$solver$patch$state)
  expect_error(tidy(quiet), "record = \"steps\"")
})

test_that("a stochastic run refuses what it cannot do", {
  p <- stochastic_stand(5)
  expect_error(run(p, method = "stochastic", until = 2), "until")
  expect_error(run(p, method = "stochastic", refine_schedule = TRUE), "refine_schedule")
  expect_error(run(p, method = "stochastic", record = "states"), "states")
  expect_error(run(p, method = "stochastic",
                   events = events(node_introductions(p), harvest(time = 2, fraction = 0.3))),
               "events")
  expect_error(run(p, method = "stochastic", from = run(p)), "from")
  res <- run(p, method = "stochastic")
  expect_error(export_patch_state(res), "deterministic")
  expect_error(run(res), "deterministic")
  expect_error(run_mutant(res, p), "deterministic")
  expect_error(run(p, method = "neither"), "should be one of")
})

test_that("a single plant's run is grow_individual_to_*()'s", {
  ind <- FF16_Individual()
  env <- Environment("FF16")
  at_times <- run(ind, times = c(1, 5, 10))
  expect_s3_class(at_times, "plant_individual_run")
  expect_identical(unclass(at_times)[c("time", "state")],
                   grow_individual_to_time(ind, c(1, 5, 10), env)[c("time", "state")])
  expect_identical(tidy(at_times), tidy_individual(grow_individual_to_time(ind, c(1, 5, 10), env)))
  at_sizes <- run(ind, sizes = c(1, 5), env = env)
  expect_identical(unclass(at_sizes)[c("time", "state")],
                   grow_individual_to_size(ind, c(1, 5), "height", env)[c("time", "state")])
  expect_error(run(ind), "one of")
  expect_error(run(ind, times = 1, sizes = 1), "one of")
  expect_error(run(ind, times = 1, warn = FALSE), "does not take warn")
})
