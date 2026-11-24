## run(): the one verb (#428). Each path is checked against the solver driven
## by hand, as the entry points run() replaced drove it, so run() adds no
## numerics of its own.

## The deterministic runner built and run by hand, step for step as the
## run_scm() entry point that run() replaced did.
solver_by_hand <- function(p, env = NULL, ctrl = control(), events = NULL,
                           refine_schedule = FALSE, collect = FALSE, run = TRUE) {
  types <- extract_RcppR6_template_types(p, "Parameters")
  if (is.null(env)) env <- Environment(types[[1]])
  if (is.null(events)) events <- empty_events()
  solver <- do.call("DeterministicPatchRunner", types)(p, env, events, ctrl)
  if (!run) return(solver)
  if (collect) solver$collect <- TRUE
  if (refine_schedule) solver$refine_schedule() else solver$run()
  solver
}

## What that entry point returned when asked to collect, from a collected
## deterministic runner, and the control it ran under, which tidy() adds.
collected_by_hand <- function(solver, ctrl = control()) {
  results <- tidy_patch(lapply(solver$history, "[[", "state"))
  results[["offspring_production"]] <- solver$patch$offspring_production
  results[["net_reproduction_ratios"]] <- solver$patch$net_reproduction_ratios
  results[["p"]] <- solver$parameters
  results[["events"]] <- solver$events
  results[["event_log"]] <- solver$event_log
  results[["control"]] <- ctrl
  results
}

outputs <- function(x) {
  solver <- if (inherits(x, "plant_run")) x$solver else x
  list(ode_state = solver$patch$ode_state,
       offspring_production = solver$patch$offspring_production,
       net_reproduction_ratios = solver$patch$net_reproduction_ratios)
}

test_that("run(p) is the deterministic run by hand, and reports what it ran", {
  p <- ff16_parameters(20)
  res <- run(p)
  expect_s3_class(res, "plant_run")
  expect_identical(outputs(res), outputs(solver_by_hand(p)))
  expect_true(res$complete)
  expect_identical(res$offspring_production, res$solver$patch$offspring_production)
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
  expect_identical(outputs(a), outputs(solver_by_hand(p, ctrl = ctrl, events = ev)))
  expect_identical(outputs(run(a$scenario)), outputs(a))
  expect_equal(a$event_log$type, b$event_log$type)
  expect_true("harvest" %in% a$event_log$type)
})

test_that("refinement through run() is the deterministic runner's, and its schedule is the result's", {
  p <- ff16_parameters(20)
  res <- run(p, refine_schedule = TRUE)
  solver <- solver_by_hand(p, refine_schedule = TRUE)
  expect_identical(outputs(res), outputs(solver))
  expect_identical(res$scenario$p, solver$parameters)
})

test_that("a refinement that runs out of passes reports a run of the schedule it leaves", {
  p <- ff16_parameters(5)
  ## One pass at this tolerance bisects and does not converge.
  ctrl <- control(schedule_eps = 1e-6, schedule_nsteps = 1)
  res <- run(p, control = ctrl, refine_schedule = TRUE)
  expect_gt(length(res$scenario$p$node_schedule_times[[1]]),
            length(node_schedule_times_default(5)))
  ## Its p runs again to the same outputs, and refines again.
  expect_identical(outputs(run(res$scenario)), outputs(res))
  expect_no_error(run(res$scenario, refine_schedule = TRUE))
  ## Keeping states changes nothing it reports.
  states <- run(p, control = ctrl, refine_schedule = TRUE, record = "states")
  expect_identical(outputs(states), outputs(res))
  expect_identical(states$scenario$p, res$scenario$p)
  expect_no_error(run_mutant(res, add_mutant(res$scenario$p, trait_matrix(0.09, "lma"))))
})

## A birth rate stepping up from near zero exactly at a scheduled node: the
## error at that node is half the step times the interval above it, which
## bisecting the interval below never reduces, so each pass flags the node
## again. The schedule starts with that interval eight doubles wide.
stepped_birth_rate_parameters <- function() {
  p <- ff16_parameters(2)
  times <- node_schedule_times_default(2)
  at <- times[length(times) - 2L]
  stopifnot(at >= 1, at < 2) # adjacent doubles are .Machine$double.eps apart
  below <- at - .Machine$double.eps
  s <- p$strategies[[1]]
  s$birth_rate_x <- c(0, below, at, 2)
  s$birth_rate_y <- c(1e-6, 1e-6, 1, 1)
  s$is_variable_birth_rate <- TRUE
  p$strategies <- list(s)
  p$node_schedule_times <- list(sort(c(times, at - 8 * .Machine$double.eps)))
  list(p = p, at = at, below = below)
}

test_that("a refinement stops splitting an interval it has halved to adjacent doubles", {
  x <- stepped_birth_rate_parameters()
  ## The interval is split in four passes. The rest are left: a pass that can
  ## split nothing it flags ends the refinement, and without that stop this
  ## test would rerun the same schedule a thousand times.
  ctrl <- control(schedule_eps = 0.1, schedule_nsteps = 1000)
  expect_warning(res <- run(x$p, control = ctrl, refine_schedule = TRUE),
                 sprintf("below 1 node \\(species 1 at t = %.17g\\).*adjacent doubles", x$at))
  ## It names the node above the interval it gave up on, and leaves that
  ## interval whole and every node introduced once.
  expect_identical(res$solver$unsplittable_node_times, list(x$at))
  times <- res$scenario$p$node_schedule_times[[1]]
  expect_identical(anyDuplicated(times), 0L)
  expect_identical(times[match(x$at, times) - 1L], x$below)
  ## The result is a run of the schedule it leaves, with or without states.
  expect_identical(outputs(run(res$scenario)), outputs(res))
  expect_warning(states <- run(x$p, control = ctrl, refine_schedule = TRUE, record = "states"),
                 "adjacent doubles")
  expect_identical(outputs(states), outputs(res))
  expect_identical(states$scenario$p, res$scenario$p)
  ## A refinement that bisects and gives up on nothing leaves an empty entry
  ## per species.
  p <- ff16_parameters(5)
  expect_no_warning(plain <- run(p, control = control(schedule_eps = 1e-3), refine_schedule = TRUE))
  expect_gt(length(plain$scenario$p$node_schedule_times[[1]]),
            length(node_schedule_times_default(5)))
  expect_identical(plain$solver$unsplittable_node_times, list(numeric(0)))
})

test_that("a refinement allowed no pass runs the schedule it was given", {
  p <- ff16_parameters(5)
  res <- run(p, control = control(schedule_nsteps = 0), refine_schedule = TRUE)
  expect_identical(outputs(res), outputs(run(p)))
  expect_identical(outputs(run(res$scenario)), outputs(res))
})

## A refined scenario whose events name its introductions, run again: the
## same introductions, and the same outputs to rounding, since with events the
## run steps to p's ODE times rather than replaying their step sizes.
expect_reruns_refined <- function(res) {
  again <- run(res$scenario)
  expect_identical(again$solver$patch$species[[1]]$node_times,
                   res$solver$patch$species[[1]]$node_times)
  expect_equal(again$offspring_production, res$offspring_production,
               tolerance = 1e-9)
  expect_equal(again$net_reproduction_ratios, res$net_reproduction_ratios,
               tolerance = 1e-9)
}

test_that("a refined run whose events name its introductions reruns the refined schedule", {
  p <- ff16_parameters(10)
  ctrl <- control(node_density_coordinate = "birth_date", schedule_eps = 5e-4)
  ev <- events(node_introductions(p), harvest(time = 5, fraction = 0.3))
  res <- run(p, control = ctrl, events = ev, refine_schedule = TRUE)
  refined <- res$scenario$p$node_schedule_times[[1]]
  expect_gt(length(refined), length(p$node_schedule_times[[1]]))
  ## The events carry the refined introductions, and keep their actions.
  ev_res <- res$scenario$events
  intro <- ev_res$type == "node_introduction"
  expect_identical(ev_res$time[intro], refined)
  expect_identical(ev_res$type[!intro], "harvest")
  expect_reruns_refined(res)
})

test_that("a species introduced at no time refines and resumes", {
  p <- testing_parameters("FF16")
  p$max_patch_lifetime <- 5
  p$node_schedule_times <- list()
  p <- add_strategies(p, trait_matrix(c(0.0825, 0.12), "lma"), birth_rate = c(1, 1))
  t1 <- p$node_schedule_times[[1]]
  only_first <- p
  only_first$node_schedule_times <- list(t1, numeric(0))
  expect_identical(node_introductions(only_first)$target_index, rep(1L, length(t1)))
  ev <- events(node_introductions(only_first), harvest(time = 2, fraction = 0.2))
  ctrl <- control(node_density_coordinate = "birth_date", schedule_eps = 5e-4)
  res <- run(p, control = ctrl, events = ev, refine_schedule = TRUE)
  expect_length(res$scenario$p$node_schedule_times[[2]], 0L)
  expect_false(2L %in% res$scenario$events$target_index[
    res$scenario$events$type == "node_introduction"])
  expect_reruns_refined(res)
  ## Species 2 stops recruiting before the export, so has none left to come.
  q <- p
  q$node_schedule_times[[2]] <- c(0, 0.5, 1)
  part <- run(q, until = 2)
  resumed <- run(q, from = part)
  expect_true(resumed$complete)
  expect_length(resumed$scenario$p$node_schedule_times[[2]], 0L)
})

test_that("a run's parameters carry the introductions its events name", {
  p <- ff16_parameters(10)
  few <- p$node_schedule_times[[1]][c(1, 5, 20, 40, 60)]
  q <- p
  q$node_schedule_times <- list(few)
  res <- run(p, events = events(node_introductions(q)))
  expect_identical(res$scenario$p$node_schedule_times, list(few))
  expect_identical(res$solver$patch$species[[1]]$node_times, few)
})

test_that("a resumed run refines the introductions still to come", {
  p <- ff16_parameters(10)
  part <- run(p, until = 4)
  residual <- run(p, from = part)$scenario$p$node_schedule_times[[1]]
  ## One pass at this tolerance bisects and does not converge.
  ctrl <- control(schedule_eps = 1e-6, schedule_nsteps = 1)
  res <- run(p, control = ctrl, from = part, refine_schedule = TRUE)
  refined <- res$scenario$p$node_schedule_times[[1]]
  expect_gt(length(refined), length(residual))
  expect_true(all(residual %in% refined))
  ## Nothing goes before the first time still to come, whose interval reaches
  ## back past the resume.
  expect_identical(refined[[1]], residual[[1]])
  expect_reruns_refined(res)
  ## Converged, as in the run that found it.
  conv <- run(p, control = control(schedule_eps = 5e-4), from = part,
              refine_schedule = TRUE)
  expect_reruns_refined(conv)
})

test_that("record keeps what tidy() and a gradient read", {
  p <- ff16_parameters(10)
  res <- run(p, record = "steps")
  ## tidy() of a deterministic run is what the old collected run returned,
  ## field for field, in the same order.
  expect_identical(tidy(res), collected_by_hand(solver_by_hand(p, collect = TRUE)))
  expect_named(tidy(res), c("steps", "n_spp", "species", "env",
                            "offspring_production", "net_reproduction_ratios",
                            "p", "events", "event_log", "control"))
  expect_identical(tidy(res)$p, res$scenario$p)
  ev <- events(node_introductions(p), harvest(time = 5, fraction = 0.3))
  harvested <- run(p, control = control(node_density_coordinate = "birth_date"),
                   events = ev, record = "steps")
  by_hand <- solver_by_hand(p, ctrl = control(node_density_coordinate = "birth_date"),
                            events = ev, collect = TRUE)
  expect_identical(tidy(harvested),
                   collected_by_hand(by_hand, control(node_density_coordinate = "birth_date")))
  expect_true("harvest" %in% tidy(harvested)$event_log$type)
  expect_error(tidy(run(p)), "record = \"steps\"")
  expect_error(tidy(run(p, until = 0, record = "steps")), "no step yet")
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
                   run_until_converged(solver_by_hand(p, run = FALSE), block = 5)$trace)
})

test_that("until at the start builds the run unrun, as the deterministic runner is built", {
  p <- ff16_parameters(20)
  res <- run(p, until = 0)
  expect_identical(res$time, 0)
  expect_false(res$complete)
  ## So a convergence run from it starts with the t = 0 snapshot.
  a <- run_until_converged(res, block = 5)
  solver <- solver_by_hand(p, run = FALSE)
  b <- run_until_converged(solver, block = 5)
  expect_identical(a$trace, b$trace)
  expect_identical(tidy_snapshots(a$run), tidy_snapshots(solver))
  expect_s3_class(a$run, "plant_run")
  expect_error(export_patch_state(res), "moved on")
  expect_identical(export_patch_state(a$run), export_patch_state(res$solver))
})

test_that("a resume is set_initial_state() with the state's pending events", {
  p <- ff16_parameters(20)
  p$fitness_weighting <- "none"
  ctrl <- control(node_density_coordinate = "birth_date")
  ev <- events(node_introductions(p), harvest(time = 15, fraction = 0.3))
  part <- run(p, control = ctrl, events = ev, until = 8)
  st <- export_patch_state(part)
  resumed <- run(p, control = ctrl, from = part)
  manual <- solver_by_hand(set_initial_state(p, st), ctrl = ctrl, events = st$events)
  expect_identical(outputs(resumed), outputs(manual))
  expect_identical(resumed$scenario$p$initial_time, st$time)
  expect_true("harvest" %in% resumed$event_log$type)
  expect_identical(outputs(run(p, control = ctrl, from = st)), outputs(resumed))
  ## The events are the state's: a scenario that brings its own is refused.
  expect_error(run(p, control = ctrl, events = ev, from = st), "none of its own")
})

test_that("a resume into a longer scenario keeps its introductions", {
  short <- ff16_parameters(20)
  short$fitness_weighting <- "none"
  long <- ff16_parameters(40)
  long$fitness_weighting <- "none"
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

test_that("the functions that took a deterministic runner take a result", {
  p <- ff16_parameters(20)
  res <- run(p)
  expect_identical(export_patch_state(res), export_patch_state(res$solver))
})

test_that("run_mutant() is the solver's own invasion, and retires the resident's result", {
  p <- ff16_parameters(20)
  pm <- add_mutant(p, trait_matrix(0.1, "lma"))

  solver <- solver_by_hand(p)
  solver$run_mutant(pm)
  res <- run(p)
  inv <- run_mutant(res, pm)
  expect_s3_class(inv, "plant_run")
  expect_identical(inv$offspring_production, solver$patch$offspring_production)
  expect_identical(inv$mutants, pm)
  expect_identical(inv$scenario, res$scenario)
  expect_output(print(inv), "1 species invading")
  inv3 <- run_mutant(run(p), add_mutant(p, trait_matrix(c(0.09, 0.1, 0.11), "lma"),
                                           birth_rate = as.list(rep(1, 3))))
  expect_output(print(inv3), "3 species invading")
  expect_error(export_patch_state(res), "moved on")
  expect_error(tidy(res), "moved on")
  expect_error(run_mutant(inv, pm), "resident")

  solver <- solver_by_hand(p)
  solver$run_mutant_constant(pm)
  inv <- run_mutant(run(p), pm, constant = TRUE)
  expect_identical(inv$offspring_production, solver$patch$offspring_production)

  expect_error(run_mutant(run(p, until = 5), pm), "not finished")
  expect_error(run_mutant(solver, pm), "plant_run")
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
  expect_error(run(p, from = solver_by_hand(p)), "export_patch_state")
})

## The stochastic solver and a single plant under run(), each asserted
## identical to its solver driven by hand.

stochastic_stand <- function(lifetime = 20) {
  p <- testing_parameters("FF16")
  p$max_patch_lifetime <- lifetime
  add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 2)
}

test_that("a stochastic run is the runner's by hand, under the same seed", {
  p <- stochastic_stand()
  set.seed(42)
  ctrl <- control(stochastic_patch_area = 1)
  obj <- StochasticPatchRunner("FF16", "FF16_Env")(p, Environment("FF16"), ctrl)
  obj$draw_arrivals()
  collected <- tidy_stochastic_steps(collect_stochastic_steps(obj))
  set.seed(42)
  res <- run(p, control = ctrl, method = "stochastic", record = "steps")
  expect_s3_class(res, "plant_run")
  expect_identical(res$method, "stochastic")
  collected[["control"]] <- ctrl
  expect_identical(tidy(res), collected)
  expect_true(res$complete)
  expect_null(res$offspring_production)
  ## The scenario records the arrivals drawn, not p's deterministic schedule.
  expect_identical(res$scenario$p, res$solver$parameters)
  expect_false(identical(res$scenario$p$node_schedule_times, p$node_schedule_times))
  expect_output(print(res), "stochastic")
  ## Unrecorded, it draws the same individuals and ends in the same state.
  set.seed(42)
  quiet <- run(p, control = ctrl, method = "stochastic")
  expect_identical(quiet$solver$patch$state, res$solver$patch$state)
  expect_error(tidy(quiet), "record = \"steps\"")
})

test_that("a stochastic run refuses what it cannot do", {
  p <- stochastic_stand(5)
  expect_error(run(p, method = "stochastic", until = 2), "until")
  expect_error(run(p, method = "stochastic", refine_schedule = TRUE), "refine_schedule")
  expect_error(run(p, method = "stochastic", record = "states"), "states")
  expect_error(run(p, method = "stochastic", from = run(p)), "from")
  res <- run(p, control = control(stochastic_patch_area = 1), method = "stochastic")
  expect_error(export_patch_state(res), "deterministic")
  expect_error(run(res), "deterministic")
  expect_error(run_mutant(res, p), "deterministic")
  expect_error(run(p, method = "neither"), "should be one of")
})

test_that("a single plant's run is the individual grower's", {
  ind <- FF16_Individual()
  env <- Environment("FF16")
  at_times <- run(ind, times = c(1, 5, 10))
  expect_s3_class(at_times, "plant_individual_run")
  expect_identical(unclass(at_times)[c("time", "state")],
                   individual_grow_to_times(ind, c(1, 5, 10), env)[c("time", "state")])
  expect_identical(tidy(at_times),
                   tidy_individual(individual_grow_to_times(ind, c(1, 5, 10), env)))
  at_sizes <- run(ind, sizes = c(1, 5), env = env)
  expect_identical(unclass(at_sizes)[c("time", "state")],
                   individual_grow_to_sizes(ind, c(1, 5), "height", env)[c("time", "state")])
  expect_error(run(ind), "one of")
  expect_error(run(ind, times = 1, sizes = 1), "one of")
  expect_error(run(ind, times = 1, warn = FALSE), "does not take warn")
})
