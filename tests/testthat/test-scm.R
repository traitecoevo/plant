context("SCM-general")

strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

test_that("Ported from tree1", {
  for (x in names(strategy_types)) {

    s <- strategy_types[[x]]()
    e <- environment_types[[x]]
    plant <- Individual(x, e)(s)
    cohort <- Cohort(x, e)(s)

    p <- Parameters(x, e)(strategies=list(s),
                          seed_rain=pi/2,
                          patch_area=10,
                          is_resident=TRUE)
    
    if(grepl("K93", x))
      p$k_I <- 1e-3

    expect_error(scm <- SCM(x, e)(p), "Patch area must be exactly 1 for the SCM")

    p$patch_area <- 1.0
    scm <- SCM(x, e)(p)
    expect_is(scm, sprintf("SCM<%s,%s>", x, e))

    ## NOTE: I'm not sure where these are only equal and not identical.
    expect_equal(scm$parameters, p)

    ## Check that the underlying Patch really is a Patch<CohortTop>:
    expect_is(scm$patch, sprintf("Patch<%s,%s>", x, e))
    expect_equal(length(scm$patch$species), 1)
    expect_is(scm$patch$species[[1]], sprintf("Species<%s,%s>",x,e))
    expect_is(scm$patch$species[[1]]$seed, sprintf("Cohort<%s,%s>",x,e))
    expect_identical(scm$patch$time, 0.0)

    sched <- scm$cohort_schedule
    cmp_sched <- make_cohort_schedule(p)
    expect_equal(sched$size, cmp_sched$size)
    expect_equal(sched$all_times, cmp_sched$all_times)

    ## If the schedule is for the wrong number of species, it should cause
    ## an error...
    sched2 <- CohortSchedule(sched$n_species + 1)
    expect_error(scm$cohort_schedule <- sched2, "Incorrect length input; expected 1, received 2")

    ## Build a schedule for 14 introductions from t=0 to t=5
    t <- seq(0, 5, length.out=14)
    sched$set_times(t, 1)
    sched$max_time <- max(t) + diff(t)[[1]]
    scm$cohort_schedule <- sched

    expect_identical(scm$cohort_schedule$all_times, sched$all_times)
    ## Parameters has been updated:
    expect_identical(scm$parameters$cohort_schedule_times, sched$all_times)

    ## Will be helpful for checking that things worked:
    times <- data.frame(start=t, end=c(t[-1], sched$max_time))

    ## Before starting, check that the SCM is actually empty
    expect_equal(scm$time, 0.0)
    expect_equal(scm$patch$ode_size, 0)
    expect_equal(scm$patch$ode_size, 0)
    expect_false(scm$complete)

    ## TODO: The unlist here is annoying...
    ## Can be resolved by passing off to another utility function I
    ## think, but probably better done via traits (but I don't see how
    ## to make those work).
    i <- scm$run_next()

    ode_size <- Cohort(x, e)(strategy_types[[x]]())$ode_size

    expect_equal(scm$cohort_schedule$remaining, length(t) - 1)
    expect_false(scm$complete)
    expect_equal(i, 1)
    ## Note that this is the *second* time; the time of the next
    ## introduction, and the end time of the first introduction.
    expect_identical(scm$time, times$end[1])
    expect_identical(scm$time, times$start[2])
    expect_equal(scm$patch$ode_size, ode_size)

    ## Trying to set schedule for partly run scm fails
    expect_error(scm$cohort_schedule <- sched, "Cannot set schedule without resetting first")

    i <- unlist(scm$run_next())
    expect_equal(i, 1)
    ## SCM ran successfully:
    expect_equal(scm$cohort_schedule$remaining, length(t) - 2)
    expect_false(scm$complete)
    expect_identical(scm$time, times$end[2])
    expect_identical(scm$time, times$start[3])
    expect_equal(scm$patch$ode_size, ode_size * 2)

    ## Reset everything
    ## "SCM reset successful"
    scm$reset()
    expect_equal(scm$time, 0.0)
    expect_equal(scm$time, 0.0)
    expect_equal(scm$patch$ode_size, 0)
    expect_equal(scm$cohort_schedule$remaining, length(t))

    ## At this point, and possibly before scm$seed_rains is corrupt.

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
      tt <- hh <- NULL
      species_index <- 1L
      scm$reset()
      while (!scm$complete > 0 && scm$time < t_max) {
        scm$run_next()
        tt <- c(tt, scm$time)
        hh <- c(hh, list(scm$patch$species[[species_index]]$height))
      }
      hh <- list_to_matrix(hh)
      list(t=tt, h=hh)
    }

    ## Next, Run the whole schedule using the SCM.
    res_e_1 <- run_scm_test(scm)

    ## Then, check that resetting the cohort allows rerunning easily:
    ## SCM can be rerun successfully:
    scm$reset()
    res_e_2 <- run_scm_test(scm)
    expect_identical(res_e_2, res_e_1)

    ## Pull the times out of the SCM and set them in the schedule:
    sched <- scm$cohort_schedule
    sched$ode_times <- scm$ode_times
    sched$use_ode_times <- TRUE
    scm$reset() # must reset
    scm$cohort_schedule <- sched

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

test_that("schedule setting", {
  for (x in names(strategy_types)) { 
    e <- environment_types[[x]]
    p <- Parameters(x, e)(
      strategies=list(strategy_types[[x]]()),
      seed_rain=pi/2,
      is_resident=TRUE,
      cohort_schedule_max_time=5.0)
    scm <- SCM(x, e)(p)

    ## Then set a cohort schedule:
    ## Build a schedule for 14 introductions from t=0 to t=5
    sched <- scm$cohort_schedule
    t <- seq(0, sched$max_time, length.out=14)
    scm$set_cohort_schedule_times(list(t))

    ## Did set in the SCM:
    expect_identical(scm$cohort_schedule$all_times, list(t))

    ## And updated in the parameters:
    p2 <- scm$parameters
    expect_identical(p2$cohort_schedule_max_time, sched$max_time)
    expect_identical(p2$cohort_schedule_times, list(t))

    ## Remake the schedule:
    sched2 <- make_cohort_schedule(p2)
    expect_identical(sched2$max_time, sched$max_time)
    expect_identical(sched2$all_times, list(t))

    scm2 <- SCM(x, e)(p2)
    expect_identical(scm2$cohort_schedule$max_time, sched2$max_time)
    expect_identical(scm2$cohort_schedule$all_times, sched2$all_times)
  }
})

  ## ## TODO: This is a fairly inadequate set of tests; none of the failure
  ## ## conditions are tested, and it's undefined what will happen if we
  ## ## set a cohort schedule that leaves us between introduction points.
  ## test_that("State get/set works", {
  ##   ## Next, try and partly run the SCM, grab its state and push it into a
  ##   ## second copy.
  ##   scm$reset()
  ##   tmp <- run_scm_test(scm, sched$max_time / 2)
  ##   state <- scm$state

  ##   scm2 <- new(SCM, scm$parameters)
  ##   scm2$state <- state

  ##   expect_equal(scm2$state, scm$state)
  ##   ## Emergent things:
  ##   expect_equal(scm2$patch$environment$environment_interpolator$xy,
  ##                scm$patch$environment$environment_interpolator$xy)
  ##   expect_equal(scm2$ode_state, scm$ode_state)
  ##   expect_equal(scm2$ode_rates, scm$ode_rates)
  ##   # TODO: This needs implementing; requires get/set of the ODE solver
  ##   # state.
  ##   # expect_equal(scm2$time, scm$time)
  ## })

  ## test_that("Can set times directly", {
  ##   scm$reset()
  ##   times <- scm$times(1)
  ##   times2 <- sort(c(times, 0.5*(times[-1] + times[-length(times)])))
  ##   scm$set_times(times2, 1)
  ##   expect_identical(scm$times(1), times2)
  ##   expect_identical(scm$cohort_schedule$times(1), times2)
  ##   scm$run_next()
  ##   expect_error(scm$set_times(times, 1))
  ##   scm$reset()
  ##   scm$set_times(times, 1)
  ##   expect_identical(scm$times(1), times)
  ## })

test_that("Seed rain & error calculations correct", {
  for (x in c("FF16")) {
    context(sprintf("SCM-%s", x))
    e <- environment_types[[x]]
    p0 <- scm_base_parameters(x)
    p1 <- expand_parameters(trait_matrix(0.08, "lma"), p0, mutant=FALSE)
    
    if(grepl("K93", x))
      p1$k_I <- 1e-3

    scm <- run_scm(p1)
    expect_is(scm, sprintf("SCM<%s,%s>", x, e))

    seed_rain_R <- function(scm, error=FALSE) {
      a <- scm$cohort_schedule$times(1)
      d <- scm$patch$environment$disturbance_regime
      pa <- d$density(a)
      p <- scm$parameters
      scale <- scm$parameters$strategies[[1]]$S_D * p$seed_rain
      seeds <- pa * scm$patch$species[[1]]$seeds * scale
      total <- trapezium(a, seeds)
      if (error) local_error_integration(a, seeds, total) else total
    }

    expect_equal(scm$seed_rain(1), seed_rain_R(scm))
    expect_equal(scm$seed_rains, seed_rain_R(scm))
    expect_equal(scm$seed_rain_error[[1]], seed_rain_R(scm, error=TRUE))

    lae_cmp <-
      scm$patch$species[[1]]$competition_effects_error(scm$patch$compute_competition(0))
    expect_identical(scm$competition_effect_error(1), lae_cmp)

    int <- make_scm_integrate(scm)
    S_D <- scm$parameters$strategies[[1]]$S_D
    expect_equal(int("seeds_survival_weighted") * S_D, scm$seed_rain(1))

    res <- run_scm_collect(p1)
    int2 <- make_scm_integrate(res)

    expect_equal(int2("seeds_survival_weighted"), int("seeds_survival_weighted"))
    expect_equal(int2("height"), int("height"))
    expect_equal(int2("mortality"), int("mortality"))
    expect_equal(int2("fecundity"), int("fecundity"))
  }
})

test_that("Can create empty SCM", {
  context("SCM-empty")
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    p <- Parameters(x, e)()
    scm <- SCM(x, e)(p)

    ## Check light environment is empty:
    env <- scm$patch$environment
    patch <- scm$patch
    # This is no longer zero:
    # expect_equal(env$canopy$canopy_interpolator$size, 0)
    expect_equal(env$canopy_openness(0), 1.0)
  }
})
