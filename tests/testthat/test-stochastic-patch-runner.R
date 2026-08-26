
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

    obj <- StochasticPatchRunner(x, e)(p, env, ctrl)
    expect_identical(obj$time, 0.0)

    sched <- obj$node_schedule
    expect_equal(sched$size, 0)
    expect_equal(sched$max_time, p$max_patch_lifetime)

    ## Now, create a new set of times:
    sched2 <- stochastic_schedule(p)
    expect_gt(sched2$size, 0)

    ## Does this need to happen twice?
    obj$node_schedule <- sched2
    expect_equal(obj$node_schedule$size, sched2$size)

    ## Importantly, this moves time forward to where the first
    ## introduction will be!
    expect_identical(obj$time, sched2$next_event$time_introduction)

    ## We're empty though.... The ODE system still carries the environment's own
    ## states (TF24's soil water and cumulative fluxes; none for FF16 or K93).
    expect_equal(obj$patch$species[[1]]$size, 0)
    expect_equal(obj$patch$node_ode_size, 0)
    expect_length(obj$patch$ode_state, n_environment_ode_states[[x]])

    res <- obj$run_next()
    expect_equal(res, 1L)
    expect_identical(obj$time, sched2$all_times[[1]][[2]])

    ode_size <- Individual(x, e)(strategy_types[[x]]())$ode_size
    expect_equal(obj$patch$node_ode_size, ode_size)
    expect_equal(length(obj$patch$ode_state),
                 ode_size + n_environment_ode_states[[x]])
    expect_equal(obj$patch$size, 1)

    expect_false(obj$complete)
  }
})

test_that("the expected number of arrivals scales with patch area", {
  ## stochastic_schedule() used to pass patch_area into
  ## stochastic_arrival_times()'s *third* argument, which is delta_t, leaving
  ## patch_area at its default 1. Arrivals therefore never scaled with area --
  ## the expected count came out as lifetime x birth_rate for every patch size --
  ## and the binning interval was silently set to the area.
  ##
  ## Arrivals are Poisson, so this compares realised counts against
  ## lifetime x birth_rate x area over a horizon long enough that the relative
  ## standard error is small. birth_rate 10 over 100 yr is 1000 expected
  ## arrivals per m^2 for only 100 / delta_t = 1000 draws, so the counts are
  ## large but the schedule is cheap. sd/mean is 1/sqrt(1000) = 3.2% at the
  ## smallest count, so the 20% tolerances below are all more than 5 Poisson
  ## standard deviations wide and cannot flake; the seed is fixed as well.
  ## Under the defect the 2 and 8 m^2 cases both come in at ~1000, failing by
  ## 2x and 8x. Only the R-side draw is exercised, which reads nothing from the
  ## strategy but its birth rate, so one strategy type is enough.
  lifetime <- 100
  birth_rate <- 10
  areas <- c(1, 2, 8)

  set.seed(1)
  n <- vapply(areas, function(area) {
    p <- Parameters("FF16", "FF16_Env")(strategies = list(FF16_Strategy()),
                                        patch_area = area,
                                        max_patch_lifetime = lifetime)
    p$strategies[[1]]$birth_rate_y <- birth_rate
    stochastic_schedule(p)$size
  }, numeric(1))

  ## Doubling the area doubles the arrivals; eight times the area, eight times.
  expect_equal(n[[2]] / n[[1]], 2, tolerance = 0.2)
  expect_equal(n[[3]] / n[[1]], 8, tolerance = 0.2)

  ## And the absolute rate is right, not just its scaling.
  for (i in seq_along(areas)) {
    expect_equal(n[[i]], lifetime * birth_rate * areas[[i]], tolerance = 0.2)
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
    res <- run_stochastic_collect(p, Environment(x), stochastic_control(x))

    ## Regression guard for #498: the collector used to read a removed `state`
    ## accessor and silently returned empty output, which `expect_silent` could
    ## not catch. Assert the trajectory is actually populated.
    expect_setequal(names(res), c("time", "species", "env", "p"))
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
  ## patch_area is 1 m^2, the default. It used to read 50, but that was inert:
  ## stochastic_schedule() passed it into stochastic_arrival_times()'s delta_t
  ## slot, so arrivals never scaled with area and the stand came out at ~105
  ## individuals whatever the area. Now that area-scaling works, 50 m^2 would
  ## mean ~5300 individuals, measured at ~670x the run time for FF16 alone;
  ## 1 m^2 keeps the stand the size this baseline has always actually run at,
  ## over the full default 105.32 yr lifetime. Every value below was regenerated
  ## from a run after the fix -- the stand is the same size but 50x denser, so
  ## establishment and mortality both differ from the old numbers.
  ##
  ## The seeded schedule holds 117 arrivals (mean 105.32 = lifetime x birth_rate
  ## x area). n_total is how many of those established: establishment is a
  ## Bernoulli draw against the birth environment, so it is model-dependent --
  ## K93 takes all 117, while FF16 and TF24 lose about a third in the shade.
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
  ## Bounding the storage pool by the shape of its own flow (#609) moved TF24's
  ## established count from 81 to 79 and left its survivor count at 3. FF16 and
  ## K93 are untouched, which is the discriminator: the change is in TF24's
  ## storage block and nothing else reads it.
  baseline <- list(
    FF16 = list(n_total = 83L, n_alive_final = 5L),
    TF24 = list(n_total = 79L, n_alive_final = 3L),
    K93  = list(n_total = 117L, n_alive_final = 2L)
  )
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    run_once <- function() {
      set.seed(1)
      p <- Parameters(x, e)(strategies=list(strategy_types[[x]]()),
                            patch_area=1)
      run_stochastic_collect(p, Environment(x), stochastic_control(x))
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
})

test_that("TF24's environment states are integrated over a stochastic run", {
  ## The stochastic patch used to run the species alone, so TF24's soil water sat
  ## at its initial 0.214 for a whole run while the SCM integrated it up to the
  ## drainage equilibrium and then drew it down as leaf area grew.
  p <- Parameters("TF24", "TF24_Env")(strategies=list(TF24_Strategy()),
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
  expect_true(all(theta > 0.214))
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
  res <- run_stochastic_collect(p, Environment(x), stochastic_control(x))

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
