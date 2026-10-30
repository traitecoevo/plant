
strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

test_that("Run SCM", {
  ## TF24 is excluded here: this test covers strategy-agnostic SCM machinery
  ## (run_next/reset/scheduling/ode_times reproducibility), already exercised
  ## by FF16 & K93. TF24's full-SCM behaviour is covered in
  ## test-strategy-tf24.R, and a full TF24 run here costs ~4 min (its
  ## hydraulics make each schedule run ~59s, and this test runs four) for no
  ## unique coverage.
  for (x in setdiff(names(strategy_types), "TF24")) {

    s <- strategy_types[[x]]()
    e <- environment_types[[x]]
    plant <- Individual(x, e)(s)
    node <- Node(x, e)(s)

    p <- Parameters(x, e)(strategies=list(s),
                          patch_area=1)
    
    env <- Environment(x)
    ## This test hand-builds a fine uniform schedule and checks the solver
    ## lands exactly on each introduction time, which needs the small
    ## ode_step_size_max of the accurate preset (the fast Control() default
    ## permits steps that overshoot these closely-spaced introductions).
    ctrl <- control_accurate()
    scm <- SCM(x, e)(p, env, empty_events(), ctrl)
    expect_inherits(scm, sprintf("SCM<%s,%s>", x, e))

    expect_equal(scm$parameters, p)

    ## Check that the underlying Patch really is a Patch<NodeTop>:
    expect_inherits(scm$patch, sprintf("Patch<%s,%s>", x, e))
    expect_equal(length(scm$patch$species), 1)
    expect_inherits(scm$patch$species[[1]], sprintf("Species<%s,%s>",x,e))
    expect_inherits(scm$patch$species[[1]]$new_node, sprintf("Node<%s,%s>",x,e))
    expect_identical(scm$patch$time, 0.0)

    sched <- scm$node_schedule
    cmp_sched <- make_node_schedule(p)
    expect_equal(sched$size, cmp_sched$size)
    expect_equal(sched$all_times, cmp_sched$all_times)

    ## If the schedule is for the wrong number of species, it should cause
    ## an error...
    sched2 <- plant:::NodeSchedule(sched$n_species + 1)
    expect_error(scm$node_schedule <- sched2, "Incorrect length input; expected 1, received 2")

    ## Build a schedule for 14 introductions from t=0 to t=5
    t <- seq(0, 5, length.out=14)
    sched$set_times(t, 1)
    sched$max_time <- max(t) + diff(t)[[1]]
    scm$node_schedule <- sched

    expect_identical(scm$node_schedule$all_times, sched$all_times)
    ## Parameters has been updated:
    expect_identical(scm$parameters$node_schedule_times, sched$all_times)

    ## Will be helpful for checking that things worked:
    times <- data.frame(start=t, end=c(t[-1], sched$max_time))

    ## Before starting, check that the SCM is actually empty
    expect_equal(scm$time, 0.0)
    expect_equal(scm$patch$node_ode_size, 0)
    expect_equal(scm$node_schedule$remaining, length(t))

    ode_size <- Node(x, e)(strategy_types[[x]]())$ode_size

    ## Run the whole schedule, collecting a patch snapshot after each
    ## introduction so we can check per-step progression.
    scm$collect <- TRUE
    scm$run()
    expect_equal(scm$node_schedule$remaining, 0)

    ## history[[1]] is the initial (empty) patch; the remaining snapshots are
    ## the state after each of the scheduled introductions.
    hist <- scm$history
    expect_equal(length(hist), length(t) + 1)
    expect_equal(hist[[1]]$time, 0.0)
    expect_equal(hist[[1]]$node_ode_size, 0)

    ## Each snapshot lands on the corresponding introduction end-time
    ## (= the next introduction's start time) ...
    expect_equal(sapply(hist[-1], function(h) h$time), times$end)
    ## ... and adds one node's worth of ODE state per introduction.
    expect_equal(sapply(hist[-1], function(h) h$node_ode_size),
                 ode_size * seq_along(t))

    ## Trying to set a schedule without resetting first fails
    expect_error(scm$node_schedule <- sched, "Cannot set schedule without resetting first")

    ## Reset everything
    ## "SCM reset successful"
    scm$reset()
    expect_equal(scm$time, 0.0)
    expect_equal(scm$patch$node_ode_size, 0)
    expect_equal(scm$node_schedule$remaining, length(t))

    ## At this point, and possibly before scm$net_reproduction_ratio is corrupt.

    ## This is stalling really badly, but it's not totally clear why.
    ## It's *not* the ODE system thrashing (thankfully) because the
    ## number of ODE times reported are not that bad.
    ##
    ## 50.1% in growth_rate_gradient(), and 45.4% in compute_rates()
    ## and 2.8% in initial_conditions() (so that's 98.3%) total.
    ## growth_rate_gradient and initial_conditions spend *all* their
    ## time doing growth_rate_gradient(), in turn all in
    ## compute_assimilation.
    ##
    list_to_matrix <- function(x) {
      n <- max(sapply(x, length))
      t(sapply(x, function(i) c(i, rep(NA, n-length(i)))))
    }
    run_scm_test <- function(scm, t_max=Inf) {
      species_index <- 1L
      scm$reset()
      scm$collect <- TRUE
      scm$run()
      ## Drop the initial empty-patch snapshot; keep one per introduction.
      snaps <- Filter(function(h) h$time < t_max, scm$history[-1])
      tt <- sapply(snaps, function(h) h$time)
      hh <- lapply(snaps, function(h) h$species[[species_index]]$height)
      hh <- list_to_matrix(hh)
      list(t=tt, h=hh)
    }

    ## Next, Run the whole schedule using the SCM.
    res_e_1 <- run_scm_test(scm)

    ## Then, check that resetting the node allows rerunning easily:
    ## SCM can be rerun successfully:
    scm$reset()
    res_e_2 <- run_scm_test(scm)
    expect_identical(res_e_2, res_e_1)

    ## Pull the times out of the SCM and set them in the schedule:
    sched <- scm$node_schedule
    sched$set_ode_steps(scm$ode_times, scm$ode_step_sizes,
                        scm$ode_soil_input_slopes, scm$ode_soil_substep_ends)
    scm$reset() # must reset
    scm$node_schedule <- sched

    ## So; this does not actually produce *exactly* the same output, which
    ## is very surprising.  It's definitely "close enough" but not exactly
    ## the same (and differs to right around the same order as the patch
    ## vs scm case).  I suspect that this might come from the difference
    ## between stepping to a point (requiring calculating the step size)
    ## and the stepping a particular step size (requiring calculating the
    ## final time).
    ## SCM with fixed times agrees:
    res_e_3 <- run_scm_test(scm)
    expect_identical(res_e_3, res_e_1)

    ## SCM can be rerun successfully with fixed times:
    scm$reset()
    res_e_4 <- run_scm_test(scm)
    expect_identical(res_e_4, res_e_3)
  }
})

test_that("A pinned run on a run's own steps reproduces it", {
  ## Pinning only the times leaves the pinned run to recover each step size by
  ## differencing them, and fl(fl(t + h) - t) != h, so it steps by something
  ## other than what was recorded. Pinning the sizes too takes the run's own steps.
  traits <- list(FF16 = trait_matrix(0.0825, "lma"),
                 FF16r = trait_matrix(0.0825, "lma"),
                 K93 = trait_matrix(0.059, "b_0"))
  for (x in intersect(names(strategy_types), names(traits))) {
    p <- add_strategies(scm_base_parameters(x), traits[[x]])
    env <- Environment(x)

    free <- SCM(x, environment_types[[x]])(p, env, empty_events(), Control())
    free$run()

    rerun <- function(with_sizes) {
      scm <- SCM(x, environment_types[[x]])(p, env, empty_events(), Control())
      sched <- scm$node_schedule
      sched$all_times <- free$node_schedule$all_times
      if (with_sizes) {
        sched$set_ode_steps(free$ode_times, free$ode_step_sizes,
                            free$ode_soil_input_slopes, free$ode_soil_substep_ends)
      } else {
        sched$set_ode_steps(free$ode_times, numeric(0), list(), list())
      }
      scm$node_schedule <- sched
      scm$run()
      scm
    }

    pinned <- rerun(TRUE)
    expect_identical(pinned$ode_times, free$ode_times)
    expect_identical(pinned$ode_step_sizes, free$ode_step_sizes)
    expect_identical(pinned$net_reproduction_ratios,
                     free$net_reproduction_ratios)
    expect_identical(pinned$patch$ode_state, free$patch$ode_state)

    ## The sizes are what makes it exact: they are not recoverable from the
    ## times they are recorded with.
    expect_false(identical(diff(free$ode_times), free$ode_step_sizes[-1]))
  }
})

test_that("an ODE schedule is installed whole or refused", {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"))
  scm <- SCM("FF16", "FF16_Env")(p, Environment("FF16"), empty_events(), Control())
  scm$run()
  sched <- scm$node_schedule

  ## The times and the sizes are two halves of one schedule, so they go in
  ## together: there is no setter that takes one and keeps the other, and so no
  ## way to end up holding a size recorded against a different time.
  expect_error(sched$ode_step_sizes <- scm$ode_step_sizes, "read-only")
  expect_error(sched$ode_times <- scm$ode_times, "read-only")

  slopes <- scm$ode_soil_input_slopes
  ends <- scm$ode_soil_substep_ends
  expect_error(sched$set_ode_steps(scm$ode_times, scm$ode_step_sizes[-1],
                                   slopes[-1], ends[-1]),
               "same length as ode_times")
  expect_error(sched$set_ode_steps(scm$ode_times, c(0, scm$ode_step_sizes[-1]),
                                   slopes, ends),
               "First step size must be NaN")
  ## Each step's soil record goes with its size, or none is given at all.
  expect_error(sched$set_ode_steps(scm$ode_times, scm$ode_step_sizes,
                                   slopes, ends[-1]),
               "as long as ode_step_sizes")
  expect_error(sched$set_ode_steps(scm$ode_times, numeric(0), slopes, ends),
               "as long as ode_step_sizes")
  expect_error(sched$set_ode_steps(scm$ode_times, scm$ode_step_sizes,
                                   replace(slopes, 2, list(0.5)), ends),
               "step 2 has one")

  sched$set_ode_steps(scm$ode_times, scm$ode_step_sizes, slopes, ends)
  expect_true(sched$using_ode_steps)
  expect_identical(sched$ode_times, scm$ode_times)
  expect_identical(sched$ode_step_sizes, scm$ode_step_sizes)
  expect_identical(sched$ode_soil_input_slopes, slopes)
  expect_identical(sched$ode_soil_substep_ends, ends)

  ## Times with no sizes is a grid to stop at rather than a run to rerun, and a
  ## schedule holding either uses it.
  sched$set_ode_steps(scm$ode_times, numeric(0), list(), list())
  expect_true(sched$using_ode_steps)
  expect_true(all(is.na(sched$ode_step_sizes)))

  sched$clear_ode_steps()
  expect_false(sched$using_ode_steps)
  expect_identical(sched$ode_times, numeric(0))
})

test_that("TF24's soil is substepped where that lets a step grow", {
  ## Substepped, the soil is integrated in substeps of its own under the uptake
  ## extrapolated over each step, so it no longer holds the stand to the soil's
  ## stability limit or its error. A wet year, whose rain the soil drains within
  ## days while the young stand takes up a small fraction of it, is where that
  ## pays.
  p0 <- scm_base_parameters("TF24", "TF24_Env")
  p0$max_patch_lifetime <- 1
  p0$patch_type <- "fixed"
  p <- add_strategies(p0, trait_matrix(0.1978791, "lma"))
  wet <- function() {
    e <- Environment("TF24")
    e$extrinsic_drivers_set_constant("rainfall", 3)
    e
  }
  run <- function(substep, tol = 1e-5) {
    run_scm(p, wet(), Control(node_density_coordinate = "birth_date",
                              ode_tol_rel = tol, ode_tol_abs = 1e-4 * tol,
                              ode_soil_substep = substep))
  }
  fitness <- function(scm) sum(scm$offspring_production)

  coupled <- run(FALSE)
  substepped <- run(TRUE)
  taken <- lengths(substepped$ode_soil_substep_ends) > 0
  expect_false(any(lengths(coupled$ode_soil_substep_ends) > 0))
  expect_gt(mean(taken), 0.3)
  expect_lt(length(substepped$ode_times), 0.8 * length(coupled$ode_times))
  ## Each record: an input slope per layer, and substeps ending at the step's
  ## end.
  n_soil <- Environment("TF24")$get_soil_number_of_depths()
  expect_true(all(lengths(substepped$ode_soil_input_slopes[taken]) == n_soil))
  expect_true(all(vapply(substepped$ode_soil_substep_ends[taken],
                         function(e) e[length(e)] == 1, TRUE)))
  ## Its error falls with the tolerance, against a tight run with the soil
  ## coupled.
  tight <- fitness(run(FALSE, 1e-8))
  error <- function(scm) abs(fitness(scm) / tight - 1)
  ## Within twice the tolerance rather than once: the error control reads the
  ## RMS over every state, and TF24's three flexible-allometry states (#516) are
  ## held at zero error by default, which loosens the effective tolerance: a
  ## cohort carries eleven states where it carried eight (1.1e-5 measured).
  expect_lt(error(substepped), 2e-5)
  expect_lt(error(run(TRUE, 3e-6)), error(substepped))
})

test_that("A run's steps answer to the soil's tolerance scale, and at one to nothing", {
  ## TF24's environment holds the soil, so it is the model the scale reaches.
  ## One year is enough: what is checked is each step's decision.
  p0 <- scm_base_parameters("TF24", "TF24_Env")
  p0$max_patch_lifetime <- 1
  p <- add_strategies(p0, trait_matrix(0.1978791, "lma"))
  run <- function(...) run_scm(p, Environment("TF24"), Control(...))
  same_run <- function(a, b) {
    expect_identical(a$ode_times, b$ode_times)
    expect_identical(a$ode_step_sizes, b$ode_step_sizes)
    expect_identical(a$patch$ode_state, b$patch$ode_state)
  }

  ## A scale of one is the default exactly.
  base <- run()
  same_run(run(ode_tol_scale_soil = 1), base)

  ## The soil binds steps, so loosening it lengthens them. By a hundredfold:
  ## the error is an RMS over every state, and TF24's three flexible-allometry
  ## states (#516), held at zero error, dilute the soil's share of it, so a
  ## tenfold loosening no longer removes a step on this stand.
  soil <- run(ode_tol_scale_soil = 100)
  expect_lt(length(soil$ode_times), length(base$ode_times))

  ## The soil's substeps take their tolerance from the control: at its default
  ## the run is the same, and looser it takes fewer substeps.
  substeps <- function(scm) sum(lengths(scm$ode_soil_substep_ends))
  substepped <- run(ode_soil_substep = TRUE)
  same_run(run(ode_soil_substep = TRUE, ode_soil_substep_tol = 1e-9), substepped)
  expect_lt(substeps(run(ode_soil_substep = TRUE, ode_soil_substep_tol = 1e-6)),
            substeps(substepped))

  ## A model with no soil gives the scale nothing to loosen.
  ff <- add_strategies(scm_base_parameters("FF16"), trait_matrix(0.0825, "lma"))
  same_run(run_scm(ff, ctrl = Control(ode_tol_scale_soil = 100)),
           run_scm(ff, ctrl = Control()))
})

test_that("The stochastic runner warns that it ignores the soil's tolerance scale", {
  p <- add_strategies(scm_base_parameters("FF16"), trait_matrix(0.0825, "lma"))
  runner <- function(ctrl) {
    StochasticPatchRunner("FF16", "FF16_Env")(p, Environment("FF16"), ctrl)
  }
  expect_warning(runner(Control(ode_tol_scale_soil = 10)), "ignores")
  expect_no_warning(runner(Control()))
})

test_that("TF24 takes a daily rain record without being stopped at its knots", {
  ## Daily rain through a monotone cubic changes form at every day beside a wet
  ## one. A step across such a knot reads the rain only at its stages and misses
  ## what falls between, so a step that spans one substeps the soil, whatever the
  ## stand draws, and lands a substep on each. A run given no stops then matches
  ## one stopped at every knot, with ode_soil_substep on and off.
  p0 <- scm_base_parameters("TF24", "TF24_Env")
  p0$max_patch_lifetime <- 2
  p0$patch_type <- "fixed"
  p <- add_strategies(p0, trait_matrix(0.1978791, "lma"))
  set.seed(4)
  days <- seq(0, 2 * 365 + 5)
  rain <- ifelse(runif(length(days)) < 0.15, rgamma(length(days), 0.8, 0.2), 0)
  knots <- local({
    wet <- rain != 0
    n <- length(wet)
    x <- (days / 365)[wet | c(FALSE, wet[-n]) | c(wet[-1], FALSE)]
    x[x > 0 & x < 2]
  })
  env <- function() {
    e <- Environment("TF24")
    e$extrinsic_drivers_set_variable("rainfall", days / 365, rain)
    e
  }
  ## At 1e-6 the match is 3e-6 on macOS but 5e-5 on Windows; at 1e-7 it is
  ## about 1e-7 here, so the bound holds with a margin on every platform.
  for (setting in list(c(substep = TRUE, tol = 1e-7), c(substep = FALSE, tol = 1e-7))) {
    ctrl <- Control(node_density_coordinate = "birth_date",
                    ode_tol_rel = setting[["tol"]],
                    ode_tol_abs = 1e-4 * setting[["tol"]],
                    ode_soil_substep = as.logical(setting[["substep"]]))
    stops <- events(events_default(p),
                    rainfall_pulse(time = knots, depth = rep(0, length(knots))))
    stopped <- run_scm(p, env(), ctrl, events = stops)
    free <- run_scm(p, env(), ctrl)
    expect_lt(abs(sum(free$offspring_production) /
                  sum(stopped$offspring_production) - 1), 2e-5)
    t <- free$ode_times
    substepped <- lengths(free$ode_soil_substep_ends[-1]) > 0
    spans <- vapply(seq_along(substepped), function(k)
      any(knots > t[k] & knots < t[k + 1]), TRUE)
    expect_true(any(spans))
    expect_true(all(substepped[spans]))
    ends <- unlist(Map(function(t0, h, e) t0 + e * h, t[-length(t)],
                       free$ode_step_sizes[-1], free$ode_soil_substep_ends[-1]))
    landed <- vapply(knots, function(k) any(abs(c(t, ends) - k) < 1e-12), TRUE)
    expect_true(all(landed))
    if (!setting[["substep"]]) {
      ## Stopped at every knot, no step spans one, and with ode_soil_substep
      ## off none substeps.
      expect_true(all(lengths(stopped$ode_soil_substep_ends) == 0))
    }
  }
})

test_that("A run that substepped the soil is rerun, swept and invaded as it ran", {
  ## Three years of seasonal rain, on most of whose steps the young stand takes
  ## up a small fraction of the water flowing through the soil. Each step keeps
  ## its record, so a pinned run, the sweep and an invader take it as the run did.
  p0 <- scm_base_parameters("TF24", "TF24_Env")
  p0$max_patch_lifetime <- 3
  p0$patch_type <- "fixed"
  lma <- 0.1978791
  p <- add_strategies(p0, trait_matrix(lma, "lma"))
  env <- function() {
    e <- Environment("TF24")
    t <- seq(0, 3, length.out = 601)
    e$extrinsic_drivers_set_variable("rainfall", t,
                                     0.25 * (1 + 0.95 * sin(2 * pi * t)))
    e
  }
  ctrl <- Control(node_density_coordinate = "birth_date", ode_tol_abs = 1e-8,
                  ode_soil_substep = TRUE)
  substepped <- run_scm(p, env(), ctrl, keep_states = TRUE)
  expect_gt(mean(lengths(substepped$ode_soil_substep_ends) > 0), 0.5)

  ## A pinned run of the schedule, records and all, repeats the run bit for bit.
  rerun <- function(q) {
    scm <- SCM("TF24", "TF24_Env")(q, env(), empty_events(), ctrl)
    sched <- scm$node_schedule
    sched$all_times <- substepped$node_schedule$all_times
    sched$set_ode_steps(substepped$ode_times, substepped$ode_step_sizes,
                        substepped$ode_soil_input_slopes,
                        substepped$ode_soil_substep_ends)
    scm$node_schedule <- sched
    scm$run()
    scm
  }
  again <- rerun(p)
  expect_identical(again$patch$ode_state, substepped$patch$ode_state)
  expect_identical(again$ode_soil_substep_ends, substepped$ode_soil_substep_ends)
  ## So does the schedule carried in the parameters. Its times alone, a grid,
  ## are refused where a step would substep the soil, and its sizes without
  ## their records where a step spans a knot of the rain.
  q <- p
  q$ode_times <- substepped$ode_times
  expect_error(run_scm(q, env(), ctrl), "step the run's schedule")
  q$ode_step_sizes <- substepped$ode_step_sizes
  expect_error(run_scm(q, env(), ctrl), "carries no substeps")
  q$ode_soil_input_slopes <- substepped$ode_soil_input_slopes
  q$ode_soil_substep_ends <- substepped$ode_soil_substep_ends
  carried <- run_scm(q, env(), ctrl)
  expect_identical(carried$patch$ode_state, substepped$patch$ode_state)

  ## The sweep holds each step's substeps and input slopes, as a tangent of the
  ## pinned schedule does, and both agree with its central differences.
  g <- stand_gradient(substepped, traits = "1.lma", hyperpar = NULL)
  traits <- census_trait_names_tf24(substepped)
  tangent <- ladder_trajectory_tangent_tf24(substepped,
                                            as.numeric(traits == "1.lma"))
  ## The sweep's rows: the census, then each species' fitness.
  metrics <- function(scm) c(stand_census(scm), stand_fitness(scm))
  expect_identical(names(metrics(substepped)), rownames(g$gradient))
  ## The tangent is the run at another scalar, not the run, so its value agrees
  ## to the budget test-gradient-ladder-introductions.R explains, not bit for bit
  ## (3.6e-13 on macOS, 1.25e-12 on Windows).
  expect_equal(tangent$value, unname(metrics(substepped)), tolerance = 1e-9)
  expect_equal(unname(g$gradient[, "1.lma"]), tangent$tangent, tolerance = 1e-6)
  at_lma <- function(x) {
    q <- p
    s <- q$strategies[[1]]
    pars <- s$pars
    pars$lma <- x
    s$pars <- pars
    q$strategies[[1]] <- s
    metrics(rerun(q))
  }
  ## The difference's truncation error at a relative step of 1e-4 is 7% on the
  ## first row; at 1e-5 it is 0.4%, and at 1e-6 rounding takes over.
  d <- 1e-5 * lma
  differenced <- (at_lma(lma + d) - at_lma(lma - d)) / (2 * d)
  expect_equal(unname(g$gradient[, "1.lma"]), unname(differenced),
               tolerance = 1e-2)

  ## An invader with the stand's traits, integrated on the run's steps with the
  ## soil read from the run, has the stand's fitness. Not to the bit: on this
  ## stand it misses by 4e-14 with the soil coupled too.
  ratios <- substepped$net_reproduction_ratios
  substepped$run_mutant(p)
  expect_equal(substepped$net_reproduction_ratios, ratios, tolerance = 1e-12)
})

test_that("schedule setting", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    p <- Parameters(x, e)(
      strategies=list(strategy_types[[x]]()),
      max_patch_lifetime=5.0)
    env <- Environment(x)
    ctrl <- Control()
    scm <- SCM(x, e)(p, env, empty_events(), ctrl)

    ## Then set a node schedule:
    ## Build a schedule for 14 introductions from t=0 to t=5
    sched <- scm$node_schedule
    t <- seq(0, sched$max_time, length.out=14)
    scm$set_node_schedule_times(list(t))

    ## Did set in the SCM:
    expect_identical(scm$node_schedule$all_times, list(t))

    ## And updated in the parameters:
    p2 <- scm$parameters
    expect_identical(p2$max_patch_lifetime, sched$max_time)
    expect_identical(p2$node_schedule_times, list(t))

    ## Remake the schedule:
    sched2 <- make_node_schedule(p2)
    expect_identical(sched2$max_time, sched$max_time)
    expect_identical(sched2$all_times, list(t))

    scm2 <- SCM(x, e)(p2, env, empty_events(), ctrl)
    expect_identical(scm2$node_schedule$max_time, sched2$max_time)
    expect_identical(scm2$node_schedule$all_times, sched2$all_times)
  }
})

test_that("Offspring production & error calculations correct", {
  for (x in c("FF16")) {
    e <- environment_types[[x]]
    p0 <- scm_base_parameters(x)
    p1 <- add_strategies(p0, trait_matrix(0.08, "lma"), birth_rate = 1.0)
    
    env <- Environment(x)
    ctrl <- Control()

    scm <- run_scm(p1, env, ctrl)
    expect_inherits(scm, sprintf("SCM<%s,%s>", x, e))

    net_reproduction_ratio_R <- function(scm, error=FALSE) {
      a <- scm$node_schedule$times(1)
      density <- purrr::map_dbl(a, ~ scm$patch$density(.x))
      net_reproduction_ratio_by_node_weighted <- density *
        scm$patch$species[[1]]$net_reproduction_ratio_by_node *
        scm$parameters$strategies[[1]]$pars$S_D
      total <- trapezium(a, net_reproduction_ratio_by_node_weighted)
      if (error)
        local_error_integration(a, net_reproduction_ratio_by_node_weighted, total)
      else total
    }

    expect_equal(scm$net_reproduction_ratio_for_species(1), net_reproduction_ratio_R(scm))
    expect_equal(scm$net_reproduction_ratios, net_reproduction_ratio_R(scm))
    expect_equal(scm$net_reproduction_ratio_errors[[1]], net_reproduction_ratio_R(scm, error=TRUE))

    lae_cmp <-
      scm$patch$species[[1]]$compute_competition_effect_by_nodes_error(scm$patch$compute_competition(0))
    expect_identical(scm$compute_competition_effect_error_by_node_for_species_i(1), lae_cmp)

  }
})

test_that("refinement_error_by_node collected in C++ matches per-step assembly", {
  for (x in c("FF16")) {
    e <- environment_types[[x]]
    p0 <- scm_base_parameters(x)
    p1 <- add_strategies(p0, trait_matrix(0.08, "lma"), birth_rate = 1.0)
    env <- Environment(x)
    ctrl <- Control()
    n_spp <- length(p1$strategies)

    ## New path: a single run with error collection enabled.
    scm <- SCM(x, e)(p1, env, empty_events(), ctrl)
    scm$collect_refinement_errors <- TRUE
    scm$run()
    new_total <- scm$refinement_error_by_node

    ## Reference: sample the competition error per introduction step and take
    ## the column-wise max with the reproduction error (the logic the R-side
    ## refinement loop used to perform). Rather than stepping the SCM, we run it
    ## once collecting a patch snapshot after each introduction, then recompute
    ## each step's competition error from the snapshot. A species counts as
    ## "added" at a step when its node count grew from the previous snapshot.
    scm_ref <- SCM(x, e)(p1, env, empty_events(), ctrl)
    scm_ref$collect <- TRUE
    scm_ref$run()
    lai_error <- rep(list(NULL), n_spp)
    sizes_prev <- rep(0, n_spp)
    for (h in scm_ref$history) {
      sizes_now <- sapply(seq_len(n_spp), function(i) h$species[[i]]$size)
      for (idx in which(sizes_now > sizes_prev)) {
        lai_error[[idx]] <- c(
          lai_error[[idx]],
          list(h$species[[idx]]$compute_competition_effect_by_nodes_error(
            h$compute_competition(0)))
        )
      }
      sizes_prev <- sizes_now
    }
    rbind_list <- function(z) do.call("rbind", as.list(z))
    lai_error <- lapply(lai_error, function(z) rbind_list(pad_matrix(z)))
    repro <- scm_ref$net_reproduction_ratio_errors
    f <- function(m) suppressWarnings(apply(m, 2, max, na.rm = TRUE))
    ref_total <- lapply(seq_len(n_spp), function(idx)
      f(rbind(lai_error[[idx]], repro[[idx]])))

    expect_equal(new_total, ref_total)
  }
})

test_that("Can create empty SCM", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    p <- Parameters(x, e)()
    env <- Environment(x)
    ctrl <- Control()
    scm <- SCM(x, e)(p, env, empty_events(), ctrl)

    ## Check light environment is empty:
    env <- scm$patch$environment
    patch <- scm$patch

    ## Clearing leaves the flat open field rather than no field: a query reads the
    ## interpolant's bounds, so there has to be a domain to read.
    expect_equal(nrow(env$light_availability$state), 3)
    expect_equal(env$light_availability$state[, 2], rep(1, 3))
    expect_equal(env$get_environment_at_height(0), 1.0)
  }
})

test_that("A second run on one SCM reproduces the first", {
  ## TF24 is the only model with environment ODE state (FF16 and K93 report
  ## ode_size 0), so it is the one that can carry state across a reset.
  ##
  ## What this covers is that reset() reaches Environment::clear(), whose
  ## clear_state() puts the soil and the flux accumulators back where the run
  ## began. Without it a second run continues out of the state the first left.
  ##
  ## A lifetime of 2, because the cost is the trajectory and the guarantee is
  ## not: the environment moves off its starting state within the first
  ## introductions, so a longer run buys price and no signal. The schedule is NOT thinned with it
  ## -- its 81 introductions are what the reproducibility check below compares,
  ## and they are the cheap half.
  p0 <- scm_base_parameters("TF24", "TF24_Env")
  p0$max_patch_lifetime <- 2
  p <- add_strategies(p0, trait_matrix(0.1978791, "lma"))
  new_scm <- function() SCM("TF24", "TF24_Env")(p, Environment("TF24"), empty_events(), Control())

  ## The environment's states are the tail of the patch ODE state, and before
  ## any introduction and after a reset they are all of it.
  scm <- new_scm()
  n_env <- scm$patch$environment$ode_size
  unrun <- utils::tail(new_scm()$patch$ode_state, n_env)

  scm$run()
  first_offspring <- scm$offspring_production
  first_state <- scm$patch$ode_state

  ## The run moved the environment, which is the state the reset has to put
  ## back. Asserted rather than assumed: a fixture too short to move it would
  ## pass the two checks below while guarding nothing.
  expect_false(isTRUE(all.equal(utils::tail(first_state, n_env), unrun)))

  scm$reset()
  expect_identical(utils::tail(scm$patch$ode_state, n_env), unrun)

  scm$run()
  expect_identical(scm$offspring_production, first_offspring)
  expect_identical(scm$patch$ode_state, first_state)
})



test_that("store_recording records one state per schedule_row", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    s <- strategy_types[[x]]()
    p <- Parameters(x, e)(strategies = list(s), patch_area = 1)
    scm <- SCM(x, e)(p, Environment(x), empty_events(), Control())

    ## A short schedule: eight introductions from t = 0 to t = 2.
    sched <- scm$node_schedule
    t <- seq(0, 2, length.out = 8)
    sched$set_times(t, 1)
    sched$max_time <- max(t) + diff(t)[[1]]
    scm$node_schedule <- sched

    ## Recording moves nothing: the state the run ends on is the last recorded.
    scm$run()
    forward_state <- scm$patch$ode_state
    grid <- scm$ode_times

    traj <- scm$store_recording()
    expect_identical(traj[[length(traj)]]$state, forward_state)

    ## One record per schedule_row, which is a row per accepted step and a row per
    ## introduction. The steps are the resolved grid itself rather than any second list
    ## of times, and the schedule is the steps, so an introduction cannot reach it.
    is_junction <- vapply(traj, function(r) r$introduction, logical(1))
    steps <- traj[!is_junction]
    expect_equal(length(steps), length(grid))
    expect_identical(vapply(steps, function(r) r$time, numeric(1)), grid)

    ## No step reached the first time, and every later time is its predecessor
    ## advanced by the recorded size, which is how the stepper reached it. A
    ## introduction between two steps does not disturb this, because it takes no time.
    h <- vapply(steps, function(r) r$step_size, numeric(1))
    times <- vapply(steps, function(r) r$time, numeric(1))
    expect_true(is.na(h[[1]]))
    expect_true(all(h[-1] > 0))
    expect_identical(times[-1], times[-length(times)] + h[-1])

    ## And what an introduction row is: no size, because no step reached it, and the
    ## time of the row below it, which holds the state its map ran on. Never the
    ## first row, for the same reason.
    at <- which(is_junction)
    expect_gt(length(at), 0L)
    expect_true(all(at > 1L))
    expect_true(all(is.na(vapply(traj[at], function(r) r$step_size, numeric(1)))))
    expect_identical(vapply(traj[at], function(r) r$time, numeric(1)),
                     vapply(traj[at - 1L], function(r) r$time, numeric(1)))

    ## The state widens at introductions, so the records are ragged, and the
    ## last is as wide as the patch the run ended on.
    widths <- vapply(traj, function(r) length(r$state), integer(1))
    expect_true(all(diff(widths) >= 0))
    expect_gt(max(widths), min(widths))
    expect_equal(widths[[length(widths)]], scm$patch$ode_size)

    ## The store is a function of the run: recording it twice gives the same
    ## records, down to the last bit of every state.
    expect_identical(scm$store_recording(), traj)
  }
})
