# The event queue and its R-facing wire format (issue #522).

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
  expect_equal(ev$type, c("rainfall_pulse", "harvest", "partial_disturbance",
                          "rainfall_pulse"))
  ## Each event carries its own type's parameters, in constructor order.
  expect_equal(ev$params[[2]], c(0.5, 10))
})

test_that("Events validation rejects malformed input", {
  expect_error(Events(time = 1, type = "not_a_type", species_index = 1L,
                      params = list(numeric(0))),
               "Unknown event type 'not_a_type'")
  ## The message names the types that do exist, so the fix is visible.
  expect_error(Events(time = 1, type = "not_a_type", species_index = 1L,
                      params = list(numeric(0))),
               "rainfall_pulse")
  expect_error(Events(time = 1, type = "rainfall_pulse", species_index = 1L,
                      params = list(c(1, 2))),
               "expects 1 parameters but has 2")
  expect_error(Events(time = -1, type = "rainfall_pulse", species_index = 1L,
                      params = list(1)),
               "non-finite or negative time")
  expect_error(Events(time = c(1, 2), type = "rainfall_pulse",
                      species_index = 1L, params = list(1)),
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
    expect_equal(sort(ev$time[ev$species_index == i]),
                 sort(p$node_schedule_times[[i]]))
  }
})

test_that("an out-of-range species index is rejected", {
  p <- add_strategies(scm_base_parameters("FF16"), trait_matrix(1, "lma"))
  ev <- events(node_introductions(p))
  ## One species, so index 2 has nowhere to go.
  bad <- Events(time = c(0, 1), type = rep("node_introduction", 2),
                species_index = c(1L, 2L), params = list(numeric(0), numeric(0)))
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
  canon <- function(x) sort(paste(x$time, x$species_index))
  expect_identical(canon(back), canon(ev))
})

test_that("an unimplemented event type fails where the user can see it", {
  ## Types 2-5 are declared but not yet dispatched; until they are, reaching one
  ## must be a clear error rather than a silent no-op.
  p <- add_strategies(scm_base_parameters("FF16"), trait_matrix(1, "lma"))
  ev <- events(node_introductions(p),
               partial_disturbance(time = 1, fraction = 0.5))
  scm <- SCM("FF16", "FF16_Env")(p, Environment("FF16"), ev, control())
  expect_error(scm$run(), "Event type not implemented")
})
