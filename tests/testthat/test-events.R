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
  expect_equal(ev$type, rep("rainfall_pulse", 3))
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
    partial_disturbance(time = 40, fraction = 0.3),
    rainfall_pulse(time = c(1.5, 60), depth = 0.01),
    harvest(time = 20, fraction = 0.5, height_min = 10)
  )
  expect_equal(ev$time, c(1.5, 20, 40, 60))
  ## harvest() and partial_disturbance() are thinning() under names that read
  ## better at their own call sites; they are one action, not three.
  expect_equal(ev$type, c("rainfall_pulse", "thinning", "thinning",
                          "rainfall_pulse"))
  ## Each event carries its own type's parameters, in constructor order:
  ## fraction, height_min, height_max.
  expect_equal(ev$params[[2]], c(0.5, 10, Inf))
  expect_equal(ev$params[[3]], c(0.3, 0, Inf))
  ## A pulse acts on the abiotic state; thinning defaults to the whole patch.
  expect_equal(ev$target, c("environment", "patch", "patch", "environment"))
})

test_that("a target is validated against what the type can accept", {
  ## Soil water is not owned by any one species, so a pulse cannot be narrowed.
  expect_error(Events(time = 1, type = "rainfall_pulse", target = "species",
                      target_index = 1L, params = list(0.01)),
               "cannot be aimed at a single species")
  ## An introduction must say which species is being introduced.
  expect_error(Events(time = 1, type = "node_introduction", target = "patch",
                      target_index = 1L, params = list(numeric(0))),
               "must name the species")
  expect_error(Events(time = 1, type = "thinning", target = "nowhere",
                      target_index = 1L, params = list(c(0.5, 0, Inf))),
               "Unknown event target 'nowhere'")
  ## Thinning may be aimed either way.
  expect_equal(thinning(time = 1, fraction = 0.5)$target, "patch")
  expect_equal(thinning(time = 1, fraction = 0.5, species = 2)$target, "species")
  expect_equal(thinning(time = 1, fraction = 0.5, species = 2)$target_index, 2L)
})

test_that("Events validation rejects malformed input", {
  expect_error(Events(time = 1, type = "not_a_type", target = "environment", target_index = 1L,
                      params = list(numeric(0))),
               "Unknown event type 'not_a_type'")
  ## The message names the types that do exist, so the fix is visible.
  expect_error(Events(time = 1, type = "not_a_type", target = "environment", target_index = 1L,
                      params = list(numeric(0))),
               "rainfall_pulse")
  expect_error(Events(time = 1, type = "rainfall_pulse", target = "environment", target_index = 1L,
                      params = list(c(1, 2))),
               "expects 1 parameters but has 2")
  expect_error(Events(time = -1, type = "rainfall_pulse", target = "environment", target_index = 1L,
                      params = list(1)),
               "non-finite or negative time")
  expect_error(Events(time = c(1, 2), type = "rainfall_pulse",
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
  combined <- events(ev, thinning(time = 10, fraction = 0.5))
  expect_equal(length(combined$time), length(ev$time) + 1)
  expect_equal(sum(combined$type == "thinning"), 1)
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
    ## the default schedule -- are exercised too.
    p <- add_strategies(scm_base_parameters(x),
                        trait_matrix(c(0.08, 0.1), "lma"),
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

test_that("a rainfall pulse conserves water and respects layer capacity", {
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

test_that("a rainfall pulse is refused by environments without soil water", {
  ## FF16 carries no soil state, so a pulse aimed at it is a modelling mistake
  ## and should say so rather than be quietly dropped. Tested through a run,
  ## which is the path a user actually takes.
  p <- add_strategies(scm_base_parameters("FF16"), trait_matrix(1, "lma"))
  ev <- events(node_introductions(p), rainfall_pulse(time = 1, depth = 0.01))
  scm <- SCM("FF16", "FF16_Env")(p, Environment("FF16"), ev, control())
  expect_error(scm$run(), "no soil water state")

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

test_that("thinning removes a fraction of the standing density", {
  base <- ff16_run()
  thinned <- ff16_run(function(p) {
    events(node_introductions(p), thinning(time = 10, fraction = 0.5))
  })

  ## Density is the thing thinning acts on, and nothing puts individuals back
  ## into an existing cohort, so the removal is still visible at the end of the
  ## run. (Fitness is not the thing to test: at these patch ages it is ~1e-12,
  ## and thinning also relieves competition, so the two effects fight.)
  expect_lt(total_density(thinned), total_density(base))

  ## And the removal is recorded, with what was asked and what was done.
  log <- thinned$event_log
  expect_equal(length(log$time), 1)
  expect_equal(log$time, 10)
  expect_equal(log$type, "thinning")
  expect_equal(log$target, "patch")
  expect_equal(log$requested[[1]], c(0.5, 0, Inf))
  expect_equal(log$applied[[1]][[1]], 0.5)
  expect_gt(log$applied[[1]][[2]], 0)   # some nodes were actually touched
})

test_that("thinning respects its height band", {
  ## A band above every individual present takes nothing, so the run is
  ## untouched -- but the event still happened, and is still recorded.
  base <- ff16_run()
  spared <- ff16_run(function(p) {
    events(node_introductions(p),
           thinning(time = 10, fraction = 0.9, height_min = 1e6))
  })
  expect_equal(spared$net_reproduction_ratios, base$net_reproduction_ratios)
  expect_equal(spared$event_log$applied[[1]][[2]], 0)

  ## A band that covers everything takes more than one that covers the top.
  all_sizes <- ff16_run(function(p) {
    events(node_introductions(p), thinning(time = 10, fraction = 0.5))
  })
  tops_only <- ff16_run(function(p) {
    events(node_introductions(p),
           thinning(time = 10, fraction = 0.5, height_min = 5))
  })
  expect_gt(all_sizes$event_log$applied[[1]][[2]],
            tops_only$event_log$applied[[1]][[2]])
  expect_lt(total_density(all_sizes), total_density(tops_only))
})

test_that("removing an entire cohort is refused rather than silently infinite", {
  expect_error(
    ff16_run(function(p) {
      events(node_introductions(p), thinning(time = 10, fraction = 1))
    }),
    "Thinning fraction must be in")
})

test_that("heat damage accrues only above its critical temperature", {
  base <- ff16_run()

  ## Below the threshold there is nothing to accrue, so the event is a no-op
  ## with a recorded zero rather than an unrecorded nothing.
  mild <- ff16_run(function(p) {
    events(node_introductions(p),
           heat_damage(time = 10, temperature = 30, temperature_crit = 40))
  })
  expect_equal(mild$event_log$applied[[1]][[1]], 0)
  expect_equal(mild$net_reproduction_ratios, base$net_reproduction_ratios)

  ## Above it, damage rises with both peak temperature and duration -- which is
  ## the sub-integration doing its work, since duration enters nowhere else.
  damage <- function(temperature, duration) {
    scm <- ff16_run(function(p) {
      events(node_introductions(p),
             heat_damage(time = 10, temperature = temperature,
                         duration = duration, temperature_crit = 40))
    })
    scm$event_log$applied[[1]][[1]]
  }
  expect_gt(damage(45, 14 / 365), 0)
  expect_gt(damage(50, 14 / 365), damage(45, 14 / 365))
  expect_gt(damage(45, 28 / 365), damage(45, 14 / 365))

  ## The clock does not move across the event, however long it nominally lasts.
  long <- ff16_run(function(p) {
    events(node_introductions(p),
           heat_damage(time = 10, temperature = 45, duration = 1))
  })
  expect_equal(long$time, base$time)
})

test_that("the event log records every applied event, in order", {
  scm <- ff16_run(function(p) {
    events(node_introductions(p),
           thinning(time = 5, fraction = 0.2),
           heat_damage(time = 12, temperature = 45),
           thinning(time = 15, fraction = 0.1, height_min = 2))
  })
  log <- scm$event_log
  expect_equal(log$time, c(5, 12, 15))
  expect_equal(log$type, c("thinning", "heat_damage", "thinning"))
  expect_equal(length(log$requested), 3)
  expect_equal(length(log$applied), 3)

  ## Node introductions are not logged: they are the schedule, not an
  ## intervention, and logging 140 of them would bury the three that matter.
  expect_false(any(log$type == "node_introduction"))

  ## reset() clears it, so the log always describes the run in front of you.
  scm$reset()
  expect_equal(length(scm$event_log$time), 0)
})
