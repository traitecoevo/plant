# The event queue and its R-facing wire format (issues #628, #522).

test_that("events() builds an empty schedule", {
  ev <- events()
  expect_equal(length(ev$time), 0)
  expect_equal(length(ev$type), 0)
  expect_equal(length(ev$params), 0)
})

test_that("constructors are vectorised over time and parameters", {
  ev <- rainfall_pulse(time = c(1, 2, 3), depth = c(0.01, 0.02, 0.03))
  expect_equal(ev$time, c(1, 2, 3))
  expect_equal(ev$type, rep("resource_pulse", 3))
  expect_equal(ev$params, list(0.01, 0.02, 0.03))

  ## A scalar parameter is recycled against the times.
  ev <- rainfall_pulse(time = c(1, 2), depth = 0.01)
  expect_equal(ev$params, list(0.01, 0.01))

  ## Anything else is an error naming the offending argument.
  expect_error(rainfall_pulse(time = c(1, 2, 3), depth = c(0.01, 0.02)),
               "'depth' must be length 1 or 3")
})

test_that("events() concatenates and sorts by time", {
  ev <- events(
    harvest(time = 40, fraction = 0.3),
    rainfall_pulse(time = c(1.5, 60), depth = 0.01),
    harvest(time = 20, fraction = 0.5, size_min = 10)
  )
  expect_equal(ev$time, c(1.5, 20, 40, 60))
  ## rainfall_pulse() is the TF24-flavoured name for a resource pulse: the
  ## generic layer knows about pools, not about water.
  expect_equal(ev$type, c("resource_pulse", "harvest", "harvest",
                          "resource_pulse"))
  ## Each event carries its own type's parameters, in constructor order:
  ## fraction, size_min, size_max.
  expect_equal(ev$params[[2]], c(0.5, 10, Inf))
  expect_equal(ev$params[[3]], c(0.3, 0, Inf))
  ## A pulse acts on the environment; harvest defaults to the whole patch.
  expect_equal(ev$target, c("environment", "patch", "patch", "environment"))
})

test_that("a target is validated against what the type can accept", {
  ## A resource pool belongs to the environment, not to any one species, so a
  ## pulse cannot be narrowed to one.
  expect_error(Events(time = 1, type = "resource_pulse", target = "species",
                      target_index = 1L, params = list(0.01)),
               "cannot be aimed at a single species")
  ## An introduction must say which species is being introduced.
  expect_error(Events(time = 1, type = "node_introduction", target = "patch",
                      target_index = 1L, params = list(numeric(0))),
               "must name the species")
  expect_error(Events(time = 1, type = "harvest", target = "nowhere",
                      target_index = 1L, params = list(c(0.5, 0, Inf))),
               "Unknown event target 'nowhere'")
  ## Harvest may be aimed either way.
  expect_equal(harvest(time = 1, fraction = 0.5)$target, "patch")
  expect_equal(harvest(time = 1, fraction = 0.5, species = 2)$target, "species")
  expect_equal(harvest(time = 1, fraction = 0.5, species = 2)$target_index, 2L)
})

test_that("Events validation rejects malformed input", {
  expect_error(Events(time = 1, type = "not_a_type", target = "environment", target_index = 1L,
                      params = list(numeric(0))),
               "Unknown event type 'not_a_type'")
  ## The message names the types that do exist, so the fix is visible.
  expect_error(Events(time = 1, type = "not_a_type", target = "environment", target_index = 1L,
                      params = list(numeric(0))),
               "resource_pulse")
  expect_error(Events(time = 1, type = "resource_pulse", target = "environment", target_index = 1L,
                      params = list(c(1, 2))),
               "expects 1 parameters but has 2")
  expect_error(Events(time = -1, type = "resource_pulse", target = "environment", target_index = 1L,
                      params = list(1)),
               "non-finite or negative time")
  expect_error(Events(time = c(1, 2), type = "resource_pulse",
                      target = "environment", target_index = 1L, params = list(1)),
               "same length")
})

test_that("node_introductions() reproduces the parameters' schedule", {
  p <- add_strategies(scm_base_parameters("FF16"),
                      trait_matrix(c(0.08, 0.1), "lma"),
                      birth_rate = list(1, 1))
  ev <- node_introductions(p)
  expect_equal(sum(ev$type == "node_introduction"), length(ev$time))
  ## One event per (species, time), and each species keeps its own times.
  for (i in seq_along(p$node_schedule_times)) {
    expect_equal(sort(ev$time[ev$target_index == i]),
                 sort(p$node_schedule_times[[i]]))
  }
})

test_that("events_default() is the schedule a run gets with no events", {
  p <- add_strategies(scm_base_parameters("FF16"), trait_matrix(1, "lma"))
  ev <- events_default(p)
  expect_true(all(ev$type == "node_introduction"))
  expect_identical(ev$time, events(node_introductions(p))$time)

  ## It composes: an Events object can be fed straight back into events(), so
  ## adding to an ordinary run does not mean rebuilding its schedule by hand.
  combined <- events(ev, harvest(time = 10, fraction = 0.5))
  expect_equal(length(combined$time), length(ev$time) + 1)
  expect_equal(sum(combined$type == "harvest"), 1)
  ## Still sorted, with the new event in its place rather than appended.
  expect_false(is.unsorted(combined$time))
})

test_that("an out-of-range species index is rejected", {
  p <- add_strategies(scm_base_parameters("FF16"), trait_matrix(1, "lma"))
  ev <- events(node_introductions(p))
  ## One species, so index 2 has nowhere to go.
  bad <- Events(time = c(0, 1), type = rep("node_introduction", 2),
                target = rep("species", 2), target_index = c(1L, 2L),
                params = list(numeric(0), numeric(0)))
  expect_error(SCM("FF16", "FF16_Env")(p, Environment("FF16"), bad, control()),
               "outside 1\\.\\.1")
})

test_that("the events path reproduces the default path exactly", {
  ## Expressing the same node introductions as events changes neither the stop
  ## times nor the actions, so the run must be identical -- not merely close.
  ## This is the guardrail on migrating introductions onto the shared queue.
  for (x in c("FF16", "K93")) {
    e <- environment_type(x)
    ## Two species, so that the tied introduction times -- every species shares
    ## the default schedule -- are exercised too. The trait has to be one the
    ## strategy actually carries: `lma` is FF16's and K93 has no such parameter
    ## (#637 now refuses an unknown trait name rather than ignoring it). Perturb
    ## the strategy's own default so this keeps working if the defaults move.
    tr <- if (x == "K93") "eta" else "lma"
    base <- do.call(sprintf("%s_Strategy", x), list())$pars[[tr]]
    p <- add_strategies(scm_base_parameters(x),
                        trait_matrix(c(base * 0.95, base * 1.05), tr),
                        birth_rate = list(1, 1))

    run <- function(ev) {
      scm <- SCM(x, e)(p, Environment(x), ev, control())
      scm$run()
      list(ode_times = scm$ode_times,
           net_reproduction_ratios = scm$net_reproduction_ratios,
           ode_state = scm$patch$ode_state)
    }

    expect_identical(run(events(node_introductions(p))),
                     run(empty_events()))
  }
})

test_that("a schedule round-trips through the events format", {
  p <- add_strategies(scm_base_parameters("FF16"),
                      trait_matrix(c(0.08, 0.1), "lma"),
                      birth_rate = list(1, 1))
  ev <- events(node_introductions(p))
  scm <- SCM("FF16", "FF16_Env")(p, Environment("FF16"), ev, control())
  back <- scm$events

  expect_identical(back$time, ev$time)
  expect_identical(back$type, ev$type)

  ## Species order within a tied time is not preserved, and is not meant to be:
  ## all species share the default introduction times, and the queue inserts an
  ## equal-time event ahead of the incumbent (the behaviour introductions have
  ## always had, kept so that migrating them onto the shared queue moves
  ## nothing). Order within a batch cannot matter anyway -- introduce_new_nodes
  ## does one environment recompute for the whole batch, and the ODE state is
  ## laid out by species index, not by introduction order. So compare the
  ## schedule as a set of (time, species) pairs.
  canon <- function(x) sort(paste(x$time, x$target_index))
  expect_identical(canon(back), canon(ev))
})

test_that("a resource pulse conserves mass and respects pool capacity", {
  env <- Environment("TF24")
  dz <- env$depth / env$get_soil_number_of_depths()
  sat <- env$soil_moist_sat

  ## A pulse the surface layer can absorb: all of it goes into layer 0, and
  ## the depth added is exactly delta_theta * dz.
  theta0 <- env$get_soil_water_state()[[1]]
  depth <- 0.005
  env$add_water_pulse(depth)
  flux <- env$get_soil_water_state_cumulative_flux()
  expect_equal(env$get_soil_water_state()[[1]], theta0 + depth / dz)
  expect_equal(flux[[1]], depth)      # sum_rainfall
  expect_equal(flux[[2]], depth)      # sum_infiltration: all accepted
  expect_equal(flux[[5]], 0)          # sum_pulse_runoff: nothing rejected
  ## Deeper layers are untouched: a pulse enters at the surface only.
  expect_equal(env$get_soil_water_state()[-1],
               rep(sat / 2, env$get_soil_number_of_depths() - 1))

  ## Water in equals water stored plus water shed, whatever the pulse size.
  env2 <- Environment("TF24")
  theta0 <- env2$get_soil_water_state()[[1]]
  big <- 0.5
  env2$add_water_pulse(big)
  f2 <- env2$get_soil_water_state_cumulative_flux()
  stored <- (env2$get_soil_water_state()[[1]] - theta0) * dz
  expect_equal(stored + f2[[5]], big)
  expect_equal(f2[[2]] + f2[[5]], f2[[1]])

  ## And the layer stops exactly at saturation rather than running past it.
  expect_equal(env2$get_soil_water_state()[[1]], sat)
  expect_gt(f2[[5]], 0)
})

test_that("a pulse into a saturated layer is shed entirely", {
  env <- Environment("TF24")
  n <- env$get_soil_number_of_depths()
  env$set_soil_water_state(rep(env$soil_moist_sat, n))
  env$add_water_pulse(0.02)
  flux <- env$get_soil_water_state_cumulative_flux()
  expect_equal(env$get_soil_water_state()[[1]], env$soil_moist_sat)
  expect_equal(flux[[2]], 0)     # nothing infiltrates
  expect_equal(flux[[5]], 0.02)  # all of it runs off
})

test_that("a resource pulse is refused by environments with no pools", {
  ## FF16 carries no soil state, so a pulse aimed at it is a modelling mistake
  ## and should say so rather than be quietly dropped. Tested through a run,
  ## which is the path a user actually takes.
  p <- add_strategies(scm_base_parameters("FF16"), trait_matrix(1, "lma"))
  ev <- events(node_introductions(p), rainfall_pulse(time = 1, depth = 0.01))
  scm <- SCM("FF16", "FF16_Env")(p, Environment("FF16"), ev, control())
  ## FF16 declares no resources at all, so the patch catches it before the
  ## environment is even asked, and says how many there are.
  expect_error(scm$run(), "this environment has 0 resources")

  expect_error(Environment("TF24")$add_water_pulse(-1),
               "finite and non-negative")
})

test_that("pulses wet the soil during a run", {
  ## Shorten the run before adding strategies: clearing node_schedule_times
  ## only takes effect when Parameters next crosses into C++ and re-validates,
  ## which add_strategies() does.
  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- 5
  p$node_schedule_times <- list()
  p <- add_strategies(p, trait_matrix(1, "lma"))

  run <- function(ev) {
    scm <- SCM("TF24", "TF24_Env")(p, Environment("TF24"), ev, control())
    scm$run()
    scm
  }

  base <- run(events(node_introductions(p)))
  pulsed <- run(events(node_introductions(p),
                       rainfall_pulse(time = c(1, 2, 3), depth = 0.02)))

  ## Every pulse is accounted for, and the run reaches the end.
  flux <- pulsed$patch$environment$get_soil_water_state_cumulative_flux()
  base_flux <- base$patch$environment$get_soil_water_state_cumulative_flux()
  expect_equal(flux[[1]] - base_flux[[1]], 0.06)   # sum_rainfall
  expect_equal(flux[[5]], 0)                       # the pulses themselves fit
  expect_equal(pulsed$time, base$time)

  ## The pulse's own 0.06 all infiltrates, but the run does not gain a full
  ## 0.06 of infiltration: a wetter surface sheds more of the *continuous*
  ## rain, through the saturation-excess term in compute_rates(). So the two
  ## channels interact, and the gain is strictly between zero and the pulse.
  ## (That shed water is currently not accumulated anywhere -- see #522.)
  infil_gain <- flux[[2]] - base_flux[[2]]
  expect_gt(infil_gain, 0)
  expect_lt(infil_gain, 0.06)

  ## The column balances: what it stored is what came in, less what drained
  ## and what the plants took. This holds with pulses in it precisely because
  ## a pulse adds to storage and to sum_infiltration together.
  balance <- function(scm) {
    e <- scm$patch$environment
    n <- e$get_soil_number_of_depths()
    dz <- e$depth / n
    f <- e$get_soil_water_state_cumulative_flux()
    stored <- sum(e$get_soil_water_state() - e$soil_moist_sat / 2) * dz
    stored - (f[[2]] - f[[3]] - f[[4]])
  }
  expect_equal(balance(base), 0, tolerance = 1e-6)
  expect_equal(balance(pulsed), 0, tolerance = 1e-6)

  ## The pulses actually did something: the extra water has to leave, and on
  ## this soil it leaves fast -- K(theta) rises as theta^16 -- so by the end of
  ## the run the pulsed column has drained more than the unpulsed one.
  ## (Which is why the *final* storage is not the thing to test.)
  expect_gt(flux[[3]], base_flux[[3]])

  ## And the extra water is fully accounted between the three sinks: whatever
  ## infiltrated over and above the base run either drained, was taken up, or
  ## is still in the column.
  e <- pulsed$patch$environment
  dz <- e$depth / e$get_soil_number_of_depths()
  stored_gain <- sum(e$get_soil_water_state() -
                     base$patch$environment$get_soil_water_state()) * dz
  expect_equal(stored_gain + (flux[[3]] - base_flux[[3]]) +
                 (flux[[4]] - base_flux[[4]]),
               infil_gain, tolerance = 1e-6)
})

## A short FF16 run, used by the demographic-event tests below.
ff16_run <- function(ev = NULL, max_lifetime = 20) {
  p <- scm_base_parameters("FF16")
  p$max_patch_lifetime <- max_lifetime
  p$node_schedule_times <- list()
  p <- add_strategies(p, trait_matrix(1, "lma"))
  scm <- SCM("FF16", "FF16_Env")(
    p, Environment("FF16"),
    if (is.null(ev)) events(node_introductions(p)) else ev(p), control())
  scm$run()
  scm
}

total_density <- function(scm) {
  sum(exp(scm$patch$species[[1]]$log_densities))
}

test_that("harvest removes a fraction of the standing density", {
  base <- ff16_run()
  thinned <- ff16_run(function(p) {
    events(node_introductions(p), harvest(time = 10, fraction = 0.5))
  })

  ## Density is the thing harvest acts on, and nothing puts individuals back
  ## into an existing cohort, so the removal is still visible at the end of the
  ## run. (Fitness is not the thing to test: at these patch ages it is ~1e-12,
  ## and harvest also relieves competition, so the two effects fight.)
  expect_lt(total_density(thinned), total_density(base))

  ## And the removal is recorded, with what was asked and what was done.
  log <- thinned$event_log
  expect_equal(length(log$time), 1)
  expect_equal(log$time, 10)
  expect_equal(log$type, "harvest")
  expect_equal(log$target, "patch")
  expect_equal(log$requested[[1]], c(0.5, 0, Inf))
  expect_equal(log$applied[[1]][[1]], 0.5)
  expect_gt(log$applied[[1]][[2]], 0)   # some nodes were actually touched
  expect_gt(log$applied[[1]][[3]], 0)   # and density actually came out
})

test_that("harvest respects its size band", {
  ## A band above every individual present takes nothing, so the run is
  ## untouched -- but the event still happened, and is still recorded.
  base <- ff16_run()
  spared <- ff16_run(function(p) {
    events(node_introductions(p),
           harvest(time = 10, fraction = 0.9, size_min = 1e6))
  })
  expect_equal(spared$net_reproduction_ratios, base$net_reproduction_ratios)
  expect_equal(spared$event_log$applied[[1]][[2]], 0)

  ## A band that covers everything takes more than one that covers the top.
  all_sizes <- ff16_run(function(p) {
    events(node_introductions(p), harvest(time = 10, fraction = 0.5))
  })
  tops_only <- ff16_run(function(p) {
    events(node_introductions(p),
           harvest(time = 10, fraction = 0.5, size_min = 5))
  })
  expect_gt(all_sizes$event_log$applied[[1]][[2]],
            tops_only$event_log$applied[[1]][[2]])
  expect_lt(total_density(all_sizes), total_density(tops_only))
})

test_that("removing an entire cohort is refused rather than silently infinite", {
  ## Caught when the schedule is built, not partway through a run.
  expect_error(harvest(time = 10, fraction = 1) |> events(),
               "which must be in \\[0, 1\\)")
})

test_that("a climate extreme accrues dose only above its threshold", {
  base <- ff16_run()

  ## Below the threshold there is nothing to accrue, so the event is a no-op
  ## with a recorded zero rather than an unrecorded nothing.
  mild <- ff16_run(function(p) {
    events(node_introductions(p),
           climate_extreme(time = 10, intensity = 30, threshold = 40))
  })
  expect_equal(mild$event_log$applied[[1]][[1]], 0)
  expect_equal(mild$net_reproduction_ratios, base$net_reproduction_ratios)

  ## Above it, damage rises with both peak intensity and duration -- which is
  ## the sub-integration doing its work, since duration enters nowhere else.
  damage <- function(intensity, duration) {
    scm <- ff16_run(function(p) {
      events(node_introductions(p),
             climate_extreme(time = 10, intensity = intensity,
                             duration = duration, threshold = 40))
    })
    scm$event_log$applied[[1]][[1]]
  }
  expect_gt(damage(45, 14 / 365), 0)
  expect_gt(damage(50, 14 / 365), damage(45, 14 / 365))
  expect_gt(damage(45, 28 / 365), damage(45, 14 / 365))

  ## The clock does not move across the event, however long it nominally lasts.
  long <- ff16_run(function(p) {
    events(node_introductions(p),
           climate_extreme(time = 10, intensity = 45, duration = 1))
  })
  expect_equal(long$time, base$time)
})

test_that("the event log records every applied event, in order", {
  scm <- ff16_run(function(p) {
    events(node_introductions(p),
           harvest(time = 5, fraction = 0.2),
           climate_extreme(time = 12, intensity = 45),
           harvest(time = 15, fraction = 0.1, size_min = 2))
  })
  log <- scm$event_log
  expect_equal(log$time, c(5, 12, 15))
  expect_equal(log$type, c("harvest", "climate_extreme", "harvest"))
  expect_equal(length(log$requested), 3)
  expect_equal(length(log$applied), 3)

  ## Node introductions are not logged: they are the schedule, not an
  ## intervention, and logging 140 of them would bury the three that matter.
  expect_false(any(log$type == "node_introduction"))

  ## reset() clears it, so the log always describes the run in front of you.
  scm$reset()
  expect_equal(length(scm$event_log$time), 0)
})

# Domain hooks (odelia #55/#56). These make two previously-fatal conditions into
# rejected steps instead. They are insurance rather than a demonstrated fix:
# #599's dense-stochastic failure, the case they were written for, no longer
# reproduces on develop, so there is nothing left that they visibly rescue. What
# can be pinned is that the mechanism is wired up and that it costs nothing when
# the state is fine.

test_that("the patch declares its state domain to the solver", {
  p <- add_strategies(scm_base_parameters("TF24"), trait_matrix(1, "lma"))
  scm <- SCM("TF24", "TF24_Env")(p, Environment("TF24"), empty_events(), control())
  patch <- scm$patch
  y <- patch$ode_state

  ## A sane state is accepted.
  expect_true(patch$ode_state_valid(y))

  ## Only the environment block is checked. It is the trailing part of the ODE
  ## vector, and it is where integrator overshoot shows up.
  n_env <- Environment("TF24")$ode_size
  expect_gt(n_env, 0)
  for (i in seq(length(y) - n_env + 1, length(y))) {
    bad <- y; bad[[i]] <- NaN
    expect_false(patch$ode_state_valid(bad))
    bad[[i]] <- Inf
    expect_false(patch$ode_state_valid(bad))
  }

  ## A node's log_density is legitimately -Inf (a cohort that never
  ## established), so the species block must NOT be rejected for that -- a
  ## blanket finiteness test would stall the solver at its minimum step.
  ff <- add_strategies(scm_base_parameters("FF16"), trait_matrix(1, "lma"))
  scm_ff <- SCM("FF16", "FF16_Env")(ff, Environment("FF16"), empty_events(),
                                    control())
  expect_true(scm_ff$patch$ode_state_valid(c(-Inf, 1, 2)))
})

# --- Review of #632 (@elijahmagistrado / GPT 5.6 Sol) -------------------------

test_that("collected results carry the events and the log", {
  ## Events are supplied separately from `p`, so a collected result that drops
  ## them no longer records what was asked for or what was done -- including how
  ## much of a pulse the soil took and how much it shed.
  p <- scm_base_parameters("FF16")
  p$max_patch_lifetime <- 2
  p$node_schedule_times <- list()
  p <- add_strategies(p, trait_matrix(1, "lma"))
  ev <- events(events_default(p), harvest(time = 1, fraction = 0.2))

  out <- run_scm(p, events = ev, collect = TRUE)
  expect_true(all(c("events", "event_log") %in% names(out)))

  scm <- run_scm(p, events = ev)
  expect_identical(out$events$time, scm$events$time)
  expect_identical(out$events$type, scm$events$type)
  expect_identical(out$event_log$time, scm$event_log$time)
  expect_identical(out$event_log$type, scm$event_log$type)
  expect_identical(out$event_log$applied, scm$event_log$applied)

  ## The whole requested schedule comes back, not just the interventions.
  expect_equal(length(out$events$time), length(ev$time))
  ## And the log holds exactly the non-introduction events.
  expect_equal(out$event_log$type, "harvest")

  ## An event-free run returns a valid empty log rather than dropping the field.
  bare <- run_scm(p, collect = TRUE)
  expect_true(all(c("events", "event_log") %in% names(bare)))
  expect_equal(length(bare$event_log$time), 0)

  ## Re-running does not accumulate records from the previous run.
  again <- run_scm(p, events = ev, collect = TRUE)
  expect_identical(again$event_log$time, out$event_log$time)
})

test_that("an event after the horizon is refused before the run starts", {
  p <- scm_base_parameters("FF16")
  p$max_patch_lifetime <- 2
  p$node_schedule_times <- list()
  p <- add_strategies(p, trait_matrix(1, "lma"))
  at <- function(t) events(events_default(p), harvest(time = t, fraction = 0.5))

  ## Inside and exactly at the horizon are both fine, and the boundary event is
  ## applied rather than quietly dropped.
  expect_no_error(run_scm(p, events = at(1.9)))
  edge <- run_scm(p, events = at(2))
  expect_equal(edge$event_log$time, 2)

  ## Past it, refused at construction -- naming the event, its time and the
  ## horizon. Silently discarding it would hide a units slip or a truncated run.
  expect_error(run_scm(p, events = at(2.1)), "after max_patch_lifetime = 2")
  expect_error(run_scm(p, events = at(2.1)), "harvest")

  ## Every type, and the offending row is identified among many.
  expect_error(run_scm(p, events = events(events_default(p),
                                          climate_extreme(time = 5, intensity = 45))),
               "climate_extreme")
  many <- events(events_default(p),
                 harvest(time = c(0.5, 1.0, 9.0), fraction = 0.1))
  expect_error(run_scm(p, events = many), "occurs at time 9")
})

test_that("a type refuses a target it cannot act on", {
  ## Previously accepted and then silently reinterpreted.
  expect_error(Events(time = 1, type = "harvest", target = "environment",
                      target_index = 1L, params = list(c(0.5, 0, Inf))),
               "cannot act on the environment")
  expect_error(Events(time = 1, type = "resource_pulse", target = "patch",
                      target_index = 1L, params = list(0.01)),
               "cannot act on the patch")
  ## The message says what the type will take.
  expect_error(Events(time = 1, type = "harvest", target = "environment",
                      target_index = 1L, params = list(c(0.5, 0, Inf))),
               "accepts patch, species")
})

test_that("event parameters are checked before the run, not during it", {
  ## Each of these is otherwise silent: a NaN intensity is a zero-damage event,
  ## a negative sensitivity gives negative applied mortality, and an inverted
  ## band harvests nothing while looking like it did something.
  expect_error(events(climate_extreme(time = 1, intensity = NaN)),
               "non-finite intensity")
  expect_error(events(climate_extreme(time = 1, intensity = 45, sensitivity = -1)),
               "negative sensitivity")
  expect_error(events(climate_extreme(time = 1, intensity = 45, duration = -1)),
               "negative duration")
  expect_error(events(harvest(time = 1, fraction = 0.5, size_min = 5, size_max = 2)),
               "below size_min")
  expect_error(events(rainfall_pulse(time = 1, depth = -0.01)),
               "negative amount")
})

test_that("simultaneous events keep their type order and their input order", {
  ## Across types: environment, then removals, then introductions.
  p <- scm_base_parameters("FF16")
  p$max_patch_lifetime <- 20
  p$node_schedule_times <- list()
  p <- add_strategies(p, trait_matrix(1, "lma"))
  ev <- events(harvest(time = 5, fraction = 0.1),
               climate_extreme(time = 5, intensity = 45),
               events_default(p))
  at5 <- ev$type[ev$time == 5]
  expect_equal(at5[at5 != "node_introduction"], c("climate_extreme", "harvest"))

  ## And the object agrees with the queue that will run it. This is the guard
  ## against R's idea of the order drifting from the C++ enum's.
  scm <- SCM("FF16", "FF16_Env")(p, Environment("FF16"), ev, control())
  expect_identical(scm$events$type, ev$type)
  expect_identical(scm$events$time, ev$time)

  ## Within one type at one time, input order is preserved. It matters: two
  ## pulses at an instant are capped in sequence against the same pool, so the
  ## order decides which record is credited with the accepted water.
  env <- Environment("TF24")
  env$extrinsic_drivers_set_constant("rainfall", 0)
  tf <- scm_base_parameters("TF24")
  tf$max_patch_lifetime <- 5
  tf$node_schedule_times <- list()
  tf <- add_strategies(tf, trait_matrix(1, "lma"))
  big <- 0.5   # far more than layer 0 can hold, so the first one takes it all
  pulses <- events(events_default(tf),
                   rainfall_pulse(time = c(2, 2), depth = c(big, big)))
  scm <- run_scm(tf, env = env, events = pulses)
  log <- scm$event_log
  expect_equal(length(log$time), 2)
  ## The first-applied pulse is credited with the accepted water; the second
  ## finds the layer full and is shed. Reversed input order would swap these.
  expect_gt(log$applied[[1]][[1]], 0)
  expect_equal(log$applied[[2]][[1]], 0)
  expect_equal(log$applied[[2]][[2]], big)
})

test_that("the log reports the density a removal actually took out", {
  base <- ff16_run()
  cut <- ff16_run(function(p) {
    events(node_introductions(p), harvest(time = 10, fraction = 0.5))
  })
  removed <- cut$event_log$applied[[1]][[3]]
  expect_gt(removed, 0)
  ## It is a density, not a count of numerical cohorts: the two differ, and the
  ## count alone was what the log used to claim.
  expect_false(isTRUE(all.equal(removed, cut$event_log$applied[[1]][[2]])))
})
