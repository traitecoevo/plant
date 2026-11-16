
strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

## ODE states contributed by the environment; only TF24 carries any (five
## soil-water layers then five cumulative fluxes).
n_environment_ode_states <- c(FF16 = 0L, TF24 = 10L, K93 = 0L)

## Control for a stochastic run, with TF24's soil integrated at a step size it
## can actually take.
##
## ⚠️ This caps the step size in the TEST, not in the library. Running TF24
## stochastically at the `Control()` default of 5 yr throws out of the leaf's
## collar root-find on a substantial fraction of seeds -- 5 of 40 measured on
## macOS, and 17 of 40 on develop's model code, so the dense regime is broken
## independently of this branch. The soil balance is stiff (the conductivity
## exponent is 2*n_psi+3 ~ 16, and K_sat/dz is ~543/yr at the defaults) and RKCK
## diverges rather than losing accuracy: a failing run reaches
## theta = [-51.4, 52.0, 0.146, nan, nan], and once the rates are non-finite the
## step is not even rejected, because adjust_step_size derives rmax from a NaN
## error estimate. At 0.05 yr, 0 of 60 seeds fail; the relationship is not
## monotone (0.2 fails MORE often than 5 does), so the value is measured, not
## derived from a stability limit.
##
## So this makes the suite deterministic without fixing anything for callers: a
## user who runs TF24 stochastically with default Control still hits it. That is
## tracked in #599, which is where the real fix belongs -- a per-environment cap,
## or odelia's stiff RODAS stepper for environments carrying stiff state. Delete
## this helper when #599 lands.
##
## FF16 and K93 keep plain Control(): they carry no environment ODE state, are not
## affected, and their seeded baselines below were derived without the cap.
stochastic_control <- function(x) {
  ctrl <- Control()
  if (startsWith(x, "TF24")) {
    ctrl$ode_step_size_max <- 0.05
  }
  ctrl
}

test_that("empty", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    set.seed(1)
    p <- Parameters(x, e)(strategies=list(strategy_types[[x]]()))
    
    env <- Environment(x)
    ctrl <- Control()

    ## Built, the runner holds nothing and stands at 0.
    obj <- StochasticPatchRunner(x, e)(p, env, ctrl)
    expect_identical(obj$time, 0.0)
    expect_equal(obj$schedule$size, 0)

    ## draw_arrivals() draws them: ascending, inside the patch's life, and
    ## recorded on its parameters. Its death checks are instants of their own.
    obj$draw_arrivals()
    sched <- obj$schedule
    expect_equal(sched$max_time, p$max_patch_lifetime)
    arrivals <- sched$all_times[[1]]
    expect_gt(length(arrivals), 0)
    expect_false(is.unsorted(arrivals))
    expect_true(all(arrivals > 0 & arrivals < p$max_patch_lifetime))
    expect_identical(obj$parameters$node_schedule_times, list(arrivals))
    expect_gt(sched$size, length(arrivals))

    ## Importantly, it stands at the first arrival.
    expect_identical(obj$time, arrivals[[1]])

    ## We're empty though.... The ODE system still carries the environment's own
    ## states (TF24's soil water and cumulative fluxes; none for FF16 or K93).
    expect_equal(obj$patch$species[[1]]$size, 0)
    expect_equal(obj$patch$node_ode_size, 0)
    expect_length(obj$patch$ode_state, n_environment_ode_states[[x]])

    res <- obj$run_next()
    expect_equal(res, 1L)
    expect_identical(obj$time, arrivals[[2]])

    ode_size <- Individual(x, e)(strategy_types[[x]]())$ode_size
    expect_equal(obj$patch$node_ode_size, ode_size)
    expect_equal(length(obj$patch$ode_state),
                 ode_size + n_environment_ode_states[[x]])
    expect_equal(obj$patch$size, 1)

    expect_false(obj$complete)
  }
})

test_that("the expected number of arrivals scales with patch area", {
  ## Arrivals are a Poisson process at the birth rate times the area, drawn
  ## when the runner is built. This compares realised counts against
  ## lifetime x birth_rate x area over a horizon long enough that the relative
  ## standard error is small: sd/mean is 1/sqrt(1000) = 3.2% at the smallest
  ## count, so the 20% tolerances are more than 5 Poisson standard deviations
  ## wide, and the seed is fixed. Under the old defect, where the area never
  ## reached the draw, the 2 and 8 m^2 cases came in at ~1000, failing by 2x
  ## and 8x. Only the draw is exercised, so one strategy type is enough.
  lifetime <- 100
  birth_rate <- 10
  areas <- c(1, 2, 8)

  set.seed(1)
  n <- vapply(areas, function(area) {
    p <- Parameters("FF16", "FF16_Env")(strategies = list(FF16_Model()),
                                        patch_area = area,
                                        max_patch_lifetime = lifetime)
    p$strategies[[1]]$birth_rate_y <- birth_rate
    obj <- StochasticPatchRunner("FF16", "FF16_Env")(p, Environment("FF16"),
                                                     Control())
    obj$draw_arrivals()
    length(obj$schedule$all_times[[1]])
  }, numeric(1))

  ## Doubling the area doubles the arrivals; eight times the area, eight times.
  expect_equal(n[[2]] / n[[1]], 2, tolerance = 0.2)
  expect_equal(n[[3]] / n[[1]], 8, tolerance = 0.2)

  ## And the absolute rate is right, not just its scaling.
  for (i in seq_along(areas)) {
    expect_equal(n[[i]], lifetime * birth_rate * areas[[i]], tolerance = 0.2)
  }
})

arrivals_for <- function(x, y, lifetime, area = 1, seed = 1) {
  p <- Parameters("FF16", "FF16_Env")(strategies = list(FF16_Model()),
                                      patch_area = area,
                                      max_patch_lifetime = lifetime)
  p$strategies[[1]]$is_variable_birth_rate <- TRUE
  p$strategies[[1]]$birth_rate_x <- x
  p$strategies[[1]]$birth_rate_y <- y
  set.seed(seed)
  obj <- StochasticPatchRunner("FF16", "FF16_Env")(p, Environment("FF16"),
                                                   Control())
  obj$draw_arrivals()
  obj$schedule$all_times[[1]]
}

test_that("a varying birth rate is the deterministic solver's curve", {
  ## The monotone spline the deterministic solver's inflow reads: zero until
  ## t = 1 (a monotone spline through zeros stays at zero), so nobody arrives
  ## before it, and 50 per year from 1.1.
  arrivals <- arrivals_for(c(0, 1, 1.1, 21), c(0, 0, 50, 50), 21, seed = 2)
  expect_gt(min(arrivals), 1)
  ## About 50 x 19.95 = 998, within four Poisson standard deviations.
  expect_lt(abs(length(arrivals) - 50 * 19.95), 4 * sqrt(50 * 19.95))
})

test_that("arrivals under a rising rate are thinned to its shape", {
  ## A ramp from 0 to 200 per year over 10 years: through two points the
  ## spline is the line, so 1000 arrivals are expected, at a mean time of two
  ## thirds of the span. Keeping every candidate would put the mean at half.
  arrivals <- arrivals_for(c(0, 10), c(0, 200), 10)
  n <- length(arrivals)
  expect_lt(abs(n - 1000), 4 * sqrt(1000))
  ## The arrival times have sd 10 / sqrt(18) = 2.36, so the mean's is 0.075.
  expect_lt(abs(mean(arrivals) - 20 / 3), 4 * 10 / sqrt(18 * n))
  expect_false(is.unsorted(arrivals, strictly = TRUE))
})

test_that("arrivals are drawn exactly where the spline overshoots between knots", {
  ## odelia's monotone slopes can leave a span that overshoots its data: here
  ## [1, 2] peaks just above 7. The bound is the cubic's own largest value, so
  ## the draw neither stops nor clips the overshoot, on any seed.
  for (seed in 1:20) {
    arrivals <- arrivals_for(0:4, c(0, 6, 7, 8, 108), 4, seed = seed)
    expect_false(is.unsorted(arrivals, strictly = TRUE))
  }
})

test_that("each species draws its own arrivals", {
  p <- Parameters("FF16", "FF16_Env")(
    strategies = list(FF16_Model(), FF16_Model()), max_patch_lifetime = 50)
  p$strategies[[1]]$birth_rate_y <- 2
  p$strategies[[2]]$birth_rate_y <- 20
  set.seed(3)
  obj <- StochasticPatchRunner("FF16", "FF16_Env")(p, Environment("FF16"),
                                                   Control())
  obj$draw_arrivals()
  n <- lengths(obj$schedule$all_times)
  expect_lt(abs(n[[1]] - 100), 4 * sqrt(100))
  expect_lt(abs(n[[2]] - 1000), 4 * sqrt(1000))
  expect_length(intersect(obj$schedule$all_times[[1]],
                          obj$schedule$all_times[[2]]), 0)
})

test_that("drawing arrivals needs a positive, finite patch area", {
  for (area in c(0, -1, Inf)) {
    p <- Parameters("FF16", "FF16_Env")(strategies = list(FF16_Model()),
                                        patch_area = area)
    obj <- StochasticPatchRunner("FF16", "FF16_Env")(p, Environment("FF16"),
                                                     Control())
    expect_error(obj$draw_arrivals(), "patch_area must be positive and finite")
  }
})

test_that("collect returns a well-formed, non-empty trajectory (#498)", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    set.seed(1)
    ## patch_area sets the population size: arrivals scale with it, so the
    ## default 1 m^2 is ~105 individuals over the default 105.32 yr lifetime.
    ## Keep it at 1 -- this test asserts well-formedness, not a large stand,
    ## and the cost is linear in the number of individuals.
    p <- Parameters(x, e)(strategies=list(strategy_types[[x]]()),
                          patch_area=1)
    res <- tidy(run(p, Environment(x), stochastic_control(x),
                    method = "stochastic", record = "steps"))

    ## Regression guard for #498: the collector used to read a removed `state`
    ## accessor and silently returned empty output, which `expect_silent` could
    ## not catch. Assert the trajectory is actually populated.
    expect_setequal(names(res), c("time", "species", "env", "control"))
    expect_s3_class(res$control, "Control")
    expect_gt(length(res$time), 1)
    expect_false(is.unsorted(res$time))          # patch age is non-decreasing
    expect_length(res$species, 1)

    ## The environment leg was reported under a name nothing produced, so this
    ## used to be a list of NULLs and the name check above passed anyway. Assert
    ## content, not shape.
    expect_length(res$env, length(res$time))
    expect_false(any(vapply(res$env, is.null, logical(1))))
    ## Every environment reports light availability, whatever else it carries.
    expect_true(all(vapply(res$env,
                           function(z) "light_availability" %in% names(z),
                           logical(1))))

    sp <- res$species[[1]]
    expect_equal(length(dim(sp)), 3)             # [variable, time, plant]
    expect_equal(dim(sp)[2], length(res$time))
    expect_true("height" %in% dimnames(sp)[[1]])

    ia <- attr(res$species, "is_alive")[[1]]
    expect_equal(dim(ia), c(length(res$time), dim(sp)[3]))  # [time, plant]

    ## At least one individual was introduced and some survive to the end.
    expect_gt(dim(sp)[3], 0)
    expect_gt(sum(ia[nrow(ia), ], na.rm = TRUE), 0)

    ## Heights are finite wherever an individual is alive. Padded (not-yet-born)
    ## cells are NA in `is_alive`, so select with which() to drop them.
    heights <- sp["height", , ]
    expect_true(all(is.finite(heights[which(ia == 1)])))
  }
})

test_that("collect output is reproducible and matches a seeded baseline (#482)", {
  ## With set.seed(1) the arrival schedule (R RNG) and the deaths (R::unif_rand,
  ## in C++) are fully reproducible, so the number of individuals introduced and
  ## the number alive at the final step are fixed. These golden values guard
  ## against trajectory-changing regressions in the stochastic tower; update
  ## them deliberately if the model/RNG use changes.
  ##
  ## patch_area is 1 m^2, the default, and DO NOT raise it: arrivals scale with
  ## area, so 50 m^2 means ~5300 individuals, measured at ~670x the run time for
  ## FF16 alone. 1 m^2 keeps the stand to ~105 individuals over the full default
  ## 105.32 yr lifetime, which is what every value below was generated at.
  ##
  ## The seeded run draws 101 arrivals (mean 105.32 = lifetime x birth_rate x
  ## area). n_total is how many of those established: establishment is a
  ## Bernoulli draw against the birth environment, so it is model-dependent --
  ## K93 takes all 101, while FF16 and TF24 lose a share in the shade.
  ## n_alive_final counts survivors at the last step; at this density most of the
  ## stand is thinned out, which is why these are far below the old sparse-stand
  ## values. TF24's survivor count is also sensitive to the environment's own ODE
  ## state being integrated, since that moves the mortality probabilities without
  ## changing the number of draws; FF16 and K93 carry no environment state.
  ##
  ## Note what re-deriving these costs. At patch_area = 50 the TF24 pair was
  ## 103/28 both before and after the leaf model moved out to standalone
  ## `phylloptim`, and that exact match of two seeded integers was a sharper
  ## statement that the swap preserved TF24's behaviour than any tolerance check
  ## could be. The area fix changes the run, so that particular equivalence is no
  ## longer what this test pins; it was established at the time (see the atm_kpa
  ## entry under Breaking changes in NEWS.md) and is not re-checked here.
  ## ⚠️ BOTH SIDES MOVED THIS BASELINE and the merged value is neither, so
  ## every entry below is measured from a run of this tree. Upstream's move is
  ## the paragraph after this one; ours is K93's survivor count, 2 -> 3, from the
  ## light field holding its knot data once. K93 reads neither a storage pool nor
  ## a leaf, so the light field is the only thing this branch changed that it
  ## sees -- and FF16, which shares that field and moves for nothing else here,
  ## is unmoved at 83 and 5.
  ## Bounding the storage pool by the shape of its own flow (#609) moved TF24's
  ## established count from 81 to 79 and left its survivor count at 3. FF16 and
  ## K93 are untouched, which is the discriminator: the change is in TF24's
  ## storage block and nothing else reads it.
  ## TF24 v13's default light (SWdown 900 W m^-2, 2056.5 PPFD above the canopy
  ## where the old PPFD driver gave 1800) moved TF24 to 54 established and 2
  ## surviving. FF16 and K93 read no light driver and are unmoved. That is still
  ## the count at a zero storage relaxation offset, pinned below; the default
  ## week's offset gives 45 and 1. The default stem curve moving to
  ## the TF24 trade-off function's (#611) moved both to 48 and 1: the offset no longer
  ## changes this stand's counts.
  ## Drawing deaths every 0.05 yr instead of at each arrival, and keeping each
  ## survivor's hazard reset (lost before whenever nobody died in a draw), moved
  ## every model's survivors: FF16 to 80 established and 3 surviving, TF24 to 47
  ## and 3 (47 and 1 at a zero offset, so the offset counts again), and K93 to
  ## 117 and 3. The number of draws changes, so these are a new stream rather
  ## than the old one corrected (#601).
  ## Drawing the arrivals in C++, from the model's own birth-rate curve, rather
  ## than in R (#601, re-blessed by Daniel), moved every model again: 101
  ## arrivals, FF16 84 established and 1 surviving, TF24 53 and 1 (52 and 1 at
  ## a zero offset), K93 101 and 3. A new stream again, not the old draws
  ## corrected.
  baseline <- list(
    FF16 = list(n_total = 84L, n_alive_final = 1L),
    TF24 = list(n_total = 53L, n_alive_final = 1L),
    K93  = list(n_total = 101L, n_alive_final = 3L)
  )
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    run_once <- function() {
      set.seed(1)
      p <- Parameters(x, e)(strategies=list(strategy_types[[x]]()),
                            patch_area=1)
      tidy(run(p, Environment(x), stochastic_control(x),
               method = "stochastic", record = "steps"))
    }
    res <- run_once()
    sp <- res$species[[1]]
    ia <- attr(res$species, "is_alive")[[1]]
    expect_equal(dim(sp)[3], baseline[[x]]$n_total)
    expect_equal(sum(ia[nrow(ia), ], na.rm = TRUE), baseline[[x]]$n_alive_final)

    ## Determinism: same seed -> identical trajectory.
    res2 <- run_once()
    expect_equal(res$time, res2$time)
    expect_equal(res$species, res2$species)
  }

  set.seed(1)
  s <- TF24_Model()
  s$pars$storage_relaxation_offset <- 0
  p <- Parameters("TF24", "TF24_Env")(strategies = list(s), patch_area = 1)
  res <- tidy(run(p, Environment("TF24"), stochastic_control("TF24"),
                  method = "stochastic", record = "steps"))
  ia <- attr(res$species, "is_alive")[[1]]
  expect_equal(dim(res$species[[1]])[3], 52L)
  expect_equal(sum(ia[nrow(ia), ], na.rm = TRUE), 1L)
})

test_that("TF24's environment states are integrated over a stochastic run", {
  ## The stochastic patch used to run the species alone, so TF24's soil water sat
  ## at its initial 0.214 for a whole run while the deterministic runner
  ## integrated it up to the
  ## drainage equilibrium and then drew it down as leaf area grew.
  p <- Parameters("TF24", "TF24_Env")(strategies=list(TF24_Model()),
                                      patch_area=50)
  p$max_patch_lifetime <- 8
  env <- Environment("TF24")
  expect_equal(env$get_soil_water_state(), rep(0.214, 5))

  set.seed(1)
  obj <- StochasticPatchRunner("TF24", "TF24_Env")(p, env, Control())
  obj$set_node_schedule_times(list(seq(0.2, 7.8, by = 0.2)))
  while (!obj$complete) obj$run_next()
  expect_gt(obj$patch$species[[1]]$size, 0)

  theta <- obj$patch$environment$get_soil_water_state()
  ## Moved, in whichever direction the light sets: at TF24 v13's default light
  ## the stand draws the upper layers below their initial value (0.149 at the
  ## top), where the dimmer pre-v13 default left every layer above it. What this
  ## test guards is that the state is integrated at all.
  expect_true(all(abs(theta - 0.214) > 1e-3))
  expect_true(all(theta < obj$patch$environment$soil_moist_sat))

  ## Cumulative rainfall, infiltration, deep drainage and root uptake. A non-zero
  ## uptake is the individuals' consumption reaching the water balance, which the
  ## per-node and per-species consumption_rate forwarders carry; it cannot exceed
  ## what infiltrated.
  flux <- obj$patch$environment$get_soil_water_state_cumulative_flux()
  expect_length(flux, 5)
  ## The four continuous accumulators. The fifth, sum_pulse_runoff, is fed only
  ## by rainfall-pulse events, so it is exactly zero on a run with none.
  expect_true(all(flux[1:4] > 0))
  expect_identical(flux[[5]], 0)
  expect_lt(flux[[4]], flux[[2]])
})

## The collected output now carries the environment, which only became
## meaningful once the stochastic solver started integrating it. Before that the
## soil state was frozen at its initial value, so a visible trajectory would have
## been a flat line; now it recharges and is drawn down, and this asserts the
## collector actually reports that rather than a list of NULLs.
test_that("the collected trajectory reports a moving soil state", {
  x <- "TF24"
  e <- environment_types[[x]]
  set.seed(1)
  p <- Parameters(x, e)(strategies=list(strategy_types[[x]]()), patch_area=1)
  res <- tidy(run(p, Environment(x), stochastic_control(x),
                  method = "stochastic", record = "steps"))

  expect_true("env" %in% names(res))
  expect_length(res$env, length(res$time))

  moist <- lapply(res$env, "[[", "soil_moist")
  expect_false(any(vapply(moist, is.null, logical(1))))
  expect_true(all(vapply(moist, length, integer(1)) == p$soil_number_of_depths))

  layer1 <- vapply(moist, function(z) z[[1]], 0)
  expect_true(all(is.finite(layer1)))
  ## Not frozen: the initial state is 0.214 and rainfall recharges it.
  expect_gt(diff(range(layer1)), 1e-3)
  expect_gt(max(layer1), 0.214)

  ## The cumulative fluxes are accumulators, so they cannot decrease.
  flux1 <- vapply(res$env, function(z) z[["soil_moist_cumulative_flux"]][[1]], 0)
  expect_false(is.unsorted(flux1))
  expect_gt(max(flux1), 0)
})

test_that("each leg starts where the last one landed", {
  ## The runner compares its time with the next arrival exactly. It read the
  ## patch's time, which the solver's final stage can leave one ulp off the
  ## landing time; this seed stopped at its second arrival with "Start time not
  ## what was expected".
  p <- base_parameters("FF16")
  p$max_patch_lifetime <- 30
  p$patch_area <- 1
  p <- suppressMessages(add_strategies(p, trait_matrix(0.0825, "lma"),
                                       birth_rate = 1))
  ## A death-check interval past the lifetime leaves the legs running arrival
  ## to arrival, which is where this seed failed; on the default grid most legs
  ## start at a grid point and land exactly whatever the runner reads.
  set.seed(6)
  obj <- stochastic_runner(p, Environment("FF16"),
                           Control(stochastic_death_interval = 100))
  arrivals <- obj$schedule$all_times[[1]]
  landed <- numeric(0)
  while (!obj$complete) {
    obj$run_next()
    landed <- c(landed, obj$time)
  }
  ## Every leg but the last ends on the next arrival, bit for bit.
  expect_identical(head(landed, -1), arrivals[-1])
  expect_identical(obj$time, p$max_patch_lifetime)
})

test_that("a lone plant's survival does not depend on how often deaths are drawn", {
  ## Each draw reads the hazard since the last one, so for a plant with nobody
  ## to shade it the interval cannot matter. When a draw in which nobody died
  ## left the solver unsynced, the survivor's reset was lost and the next draw
  ## counted the hazard again: at 0.02 yr no plant lived to ten years, against
  ## about one in ten drawn once at the end.
  skip_on_cran()
  p <- base_parameters("FF16")
  p$max_patch_lifetime <- 10
  p <- suppressMessages(add_strategies(p, trait_matrix(0.0825, "lma"),
                                       birth_rate = 1))
  survival <- function(interval, reps = 200) {
    ctrl <- Control(stochastic_death_interval = interval)
    alive <- vapply(seq_len(reps), function(s) {
      set.seed(s)
      obj <- stochastic_runner(p, Environment("FF16"), ctrl)
      sched <- obj$schedule
      sched$set_times(0.1, 1)
      obj$schedule <- sched
      obj$run()
      sp <- obj$patch$species[[1]]
      c(sp$size_individuals, sp$size)
    }, numeric(2))
    sum(alive[2, ]) / sum(alive[1, ])
  }
  fine <- survival(0.02)
  coarse <- survival(5)
  expect_gt(fine, 0.05)
  expect_lt(abs(fine - coarse), 0.06)
})

test_that("the death-check interval must be positive and finite", {
  p <- base_parameters("FF16")
  p <- suppressMessages(add_strategies(p, trait_matrix(0.0825, "lma"),
                                       birth_rate = 1))
  for (bad in c(0, -1, Inf, NaN)) {
    expect_error(stochastic_runner(p, Environment("FF16"),
                                   Control(stochastic_death_interval = bad)),
                 "stochastic_death_interval must be positive and finite")
  }
  ## Every check is held on the schedule, so one far below the lifetime is
  ## refused rather than allocated.
  expect_error(stochastic_runner(p, Environment("FF16"),
                                 Control(stochastic_death_interval = 1e-300)),
               "too small for this patch's lifetime")
})

test_that("the runner places its death checks on its schedule", {
  ## A check at every multiple of the interval after the first arrival and
  ## before the end, and one at the end. With arrivals at 0.1 and 0.33, an
  ## interval of 0.25 and a lifetime of 1, the checks fall at 0.25, 0.5, 0.75
  ## and 1: six instants in all. None falls before the first arrival, where
  ## nobody is alive to draw.
  p <- base_parameters("FF16")
  p$max_patch_lifetime <- 1
  p <- suppressMessages(add_strategies(p, trait_matrix(0.0825, "lma"),
                                       birth_rate = 1))
  set.seed(1)
  obj <- stochastic_runner(p, Environment("FF16"),
                           Control(stochastic_death_interval = 0.25))
  obj$set_node_schedule_times(list(c(0.1, 0.33)))
  expect_equal(obj$schedule$size, 6)
  ## Setting the schedule again replaces the checks rather than adding to them.
  obj$schedule <- obj$schedule
  expect_equal(obj$schedule$size, 6)
  expect_identical(obj$schedule$all_times, list(c(0.1, 0.33)))

  ## Lowering a copy's max_time drops the checks past it, and the runner
  ## places its own again: checks at 0.25 and 0.5 around the two arrivals.
  s <- obj$schedule
  s$max_time <- 0.5
  obj$schedule <- s
  expect_equal(obj$schedule$size, 4)

  ## An arrival on a grid point shares its instant with the check there.
  set.seed(1)
  obj <- stochastic_runner(p, Environment("FF16"),
                           Control(stochastic_death_interval = 0.25))
  obj$set_node_schedule_times(list(c(0.1, 0.5)))
  expect_equal(obj$schedule$size, 5)
  ## run_next() still steps from one arrival to the next, drawing the checks
  ## in between. The row at the shared instant is taken before its check is
  ## drawn, so the first plant still carries the hazard it built up since the
  ## check at 0.25.
  obj$run_next()
  expect_identical(obj$time, 0.5)
  sp <- obj$patch$state$species[[1]]
  expect_identical(attr(sp, "is_alive")[[1]], TRUE)
  expect_gt(sp["mortality", 1], 0)
  ## The last arrival carries on to the end of the patch's life.
  obj$run_next()
  expect_true(obj$complete)
  expect_identical(obj$time, 1)
})

test_that("a schedule the stochastic runner cannot take leaves it as it was", {
  p <- base_parameters("FF16")
  p$max_patch_lifetime <- 1
  p <- suppressMessages(add_strategies(p, trait_matrix(0.0825, "lma"),
                                       birth_rate = 1))
  set.seed(1)
  obj <- stochastic_runner(p, Environment("FF16"),
                           Control(stochastic_death_interval = 0.25))
  obj$set_node_schedule_times(list(c(0.1, 0.33)))
  size <- obj$schedule$size
  unchanged <- function() {
    expect_equal(obj$schedule$size, size)
    expect_identical(obj$schedule$all_times, list(c(0.1, 0.33)))
    expect_identical(obj$parameters$node_schedule_times, list(c(0.1, 0.33)))
  }
  ## Arrivals before the start or after the end of the patch's life.
  expect_error(obj$set_node_schedule_times(list(c(-1, 0.5))), "outside")
  unchanged()
  expect_error(obj$set_node_schedule_times(list(c(0.5, 2))), "outside")
  unchanged()
  ## A schedule with no end has nowhere for its last check.
  expect_error(obj$schedule <- {
    s <- Schedule(1)
    s$set_times(0.1, 1)
    s
  }, "finite max_time")
  unchanged()
})

test_that("the deterministic runner refuses a stochastic runner's schedule", {
  ## A death check has nothing to draw in a solver that carries mortality on
  ## density, and would still split the leg it falls in.
  p <- testing_parameters("FF16")
  p$max_patch_lifetime <- 5
  p$fitness_weighting <- "none"
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 5)
  p$node_schedule_times <- list(c(0, 1, 2))
  obj <- StochasticPatchRunner("FF16", "FF16_Env")(p, Environment("FF16"),
                                                   Control())
  obj$set_node_schedule_times(p$node_schedule_times)
  det <- DeterministicPatchRunner("FF16", "FF16_Env")(
    p, Environment("FF16"), empty_events(), Control())
  expect_error(det$schedule <- obj$schedule, "death checks")
  ## A whole-cohort harvest is refused on this path too, before any run.
  s <- det$schedule
  s$add_events(events(harvest(time = 1, fraction = 1)))
  expect_error(det$schedule <- s, "for the deterministic solver")
})

## A stand dense enough that a removal draws for dozens of individuals.
events_stand <- function() {
  p <- base_parameters("FF16")
  p$max_patch_lifetime <- 3
  suppressMessages(add_strategies(p, trait_matrix(0.0825, "lma"),
                                  birth_rate = 30))
}

test_that("a harvest of fraction 1 removes everyone alive", {
  ## The deterministic solver refuses this, since a cohort's density would go
  ## to -Inf; removing the last individual is routine here.
  p <- events_stand()
  set.seed(1)
  res <- run(p, method = "stochastic", record = "steps",
             events = events(harvest(time = 2.5, fraction = 1)))
  log <- res$event_log
  expect_identical(log$type, "harvest")
  expect_identical(log$time, 2.5)
  applied <- log$applied[[1]]
  expect_identical(applied[[1]], 1)
  expect_gt(applied[[2]], 0)
  ## Everyone drawn for is removed, reported per unit area.
  expect_identical(applied[[3]], applied[[2]] / p$patch_area)

  ## The row at the first arrival after the harvest holds nobody alive.
  st <- tidy(res)
  ia <- attr(st$species, "is_alive")[[1]]
  after <- which(st$time > 2.5)[[1]]
  expect_gt(sum(ia[after - 1, ], na.rm = TRUE), 0)
  expect_identical(sum(ia[after, ], na.rm = TRUE), 0L)
})

test_that("a harvest draws for its size band, with probability fraction", {
  p <- events_stand()
  applied <- function(ev) {
    set.seed(3)
    run(p, method = "stochastic", events = ev)$event_log$applied[[1]]
  }
  all <- applied(events(harvest(time = 2.5, fraction = 0.5)))
  none <- applied(events(harvest(time = 2.5, fraction = 0.5, size_min = 1e3)))
  ## A band nobody is in draws for nobody.
  expect_identical(none[2:3], c(0, 0))
  ## Each individual drawn for is removed with probability 0.5: within four
  ## binomial standard deviations (seeded, so this cannot flake).
  drawn <- all[[2]]
  expect_gt(drawn, 20)
  expect_lt(abs(all[[3]] * p$patch_area - 0.5 * drawn),
            4 * sqrt(0.25 * drawn))
})

test_that("a climate extreme draws against the deterministic survival", {
  p <- events_stand()
  extreme <- function(...) events(climate_extreme(time = 2.5, ...))
  stochastic <- function(ev) {
    set.seed(4)
    run(p, method = "stochastic", events = ev)$event_log$applied[[1]]
  }
  ## The fraction killed is the deterministic solver's, from the same dose.
  det <- run(p, events = events(node_introductions(p),
                                climate_extreme(time = 2.5, intensity = 45,
                                                sensitivity = 50)))
  st <- stochastic(extreme(intensity = 45, sensitivity = 50))
  expect_identical(st[[1]], det$event_log$applied[[1]][[1]])
  expect_gt(st[[1]], 0)
  expect_lt(st[[1]], 1)
  ## Below its threshold an extreme kills nobody, and draws for nobody.
  expect_identical(stochastic(extreme(intensity = 30))[2:3], c(0, 0))
  ## A dose that leaves no survival removes everyone drawn for.
  lethal <- stochastic(extreme(intensity = 1000, sensitivity = 1e4))
  expect_gt(lethal[[2]], 0)
  expect_identical(lethal[[3]] * p$patch_area, lethal[[2]])
})

test_that("a resource pulse reaches the stochastic patch's soil, even before the first arrival", {
  p <- Parameters("TF24", "TF24_Env")(strategies = list(TF24_Model()),
                                      patch_area = 1)
  p$max_patch_lifetime <- 2
  ctrl <- stochastic_control("TF24")
  depth <- 0.01
  rainfall <- function(res) {
    res$solver$patch$environment$get_soil_water_state_cumulative_flux()[[1]]
  }
  set.seed(5)
  res <- run(p, Environment("TF24"), ctrl, method = "stochastic",
             events = events(rainfall_pulse(time = c(1e-3, 1.5), depth = depth)))
  set.seed(5)
  dry <- run(p, Environment("TF24"), ctrl, method = "stochastic")
  expect_lt(1e-3, res$scenario$p$node_schedule_times[[1]][[1]])
  log <- res$event_log
  expect_identical(log$time, c(1e-3, 1.5))
  ## Each reports {accepted, shed}, and together they make up the depth.
  for (a in log$applied) {
    expect_equal(sum(a), depth)
  }
  expect_gt(log$applied[[1]][[1]], 0)
  ## The pulses reach the state the solver carries on from, not only the log:
  ## the run's cumulative rainfall holds both. The continuous accumulator moves
  ## at the ODE tolerance, since the pulses also split the legs.
  expect_equal(rainfall(res) - rainfall(dry), 2 * depth, tolerance = 1e-4)
})

test_that("an event before the first arrival on an empty FF16 patch is applied", {
  ## FF16's environment carries no state, so the leg to it is the fixed step
  ## an empty system takes.
  p <- events_stand()
  set.seed(6)
  res <- run(p, method = "stochastic",
             events = events(climate_extreme(time = 1e-3, intensity = 1000)))
  expect_lt(1e-3, res$scenario$p$node_schedule_times[[1]][[1]])
  expect_identical(res$event_log$time, 1e-3)
  expect_identical(res$event_log$applied[[1]][2:3], c(0, 0))
  expect_true(res$complete)
})

test_that("at a shared instant, a harvest acts before the arrival there", {
  p <- base_parameters("FF16")
  p$max_patch_lifetime <- 2
  p <- suppressMessages(add_strategies(p, trait_matrix(0.0825, "lma"),
                                       birth_rate = 1))
  set.seed(2)
  obj <- stochastic_runner(p, Environment("FF16"), Control())
  obj$set_node_schedule_times(list(c(0.1, 0.2, 1)))
  s <- obj$schedule
  s$add_events(events(harvest(time = 1, fraction = 1)))
  obj$schedule <- s
  while (!obj$complete) obj$run_next()
  ## Everyone who arrived before is removed, and the arrival at the harvest's
  ## own instant is not.
  alive <- attr(obj$patch$state$species[[1]], "is_alive")
  expect_identical(alive, c(FALSE, FALSE, TRUE))
  expect_identical(obj$event_log$applied[[1]][[2]], 2)
})

test_that("a harvest aimed at one species leaves the others", {
  p <- base_parameters("FF16")
  p$max_patch_lifetime <- 3
  p <- suppressMessages(add_strategies(p, trait_matrix(c(0.0825, 0.1), "lma"),
                                       birth_rate = list(30, 30)))
  set.seed(7)
  res <- run(p, method = "stochastic", record = "steps",
             events = events(harvest(time = 2.5, fraction = 1, species = 2)))
  st <- tidy(res)
  after <- which(st$time > 2.5)[[1]]
  alive <- function(i, row) {
    sum(attr(st$species, "is_alive")[[i]][row, ], na.rm = TRUE)
  }
  expect_gt(alive(2, after - 1), 0)
  expect_identical(alive(2, after), 0L)
  expect_gt(alive(1, after), 0)
})

test_that("a pulse into a pool the environment does not have is refused before the run", {
  ## FF16 has no resource pools. Refused when the schedule is set, so the
  ## runner is left as it was, rather than when the run reaches the pulse.
  p <- events_stand()
  set.seed(8)
  obj <- stochastic_runner(p, Environment("FF16"), Control())
  size <- obj$schedule$size
  s <- obj$schedule
  s$add_events(events(resource_pulse(time = 1.5, resource = 1, amount = 1)))
  expect_error(obj$schedule <- s, "does not exist")
  expect_equal(obj$schedule$size, size)
  expect_error(run(p, method = "stochastic",
                   events = events(resource_pulse(time = 1.5, resource = 1,
                                                  amount = 1))),
               "does not exist")
})

test_that("introductions in a stochastic run's events replace the drawn arrivals", {
  ## As in the deterministic solver: one arrival per time named.
  p <- events_stand()
  q <- p
  q$node_schedule_times <- list(c(1, 2))
  set.seed(9)
  res <- run(p, method = "stochastic",
             events = events(node_introductions(q),
                             harvest(time = 2.5, fraction = 0.5)))
  expect_identical(res$solver$schedule$all_times, list(c(1, 2)))
  expect_identical(res$scenario$p$node_schedule_times, list(c(1, 2)))
  ## The events' harvest is applied beside them.
  expect_identical(res$event_log$type, "harvest")
  expect_identical(res$event_log$time, 2.5)
})
