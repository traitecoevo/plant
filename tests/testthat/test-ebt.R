context("EBT")

strategy_types <- get_list_of_strategy_types()

test_that("Ported from tree1", {
  for (x in names(strategy_types)) {
    s <- strategy_types[[x]]()
    plant <- Plant(x)(s)
    cohort <- Cohort(x)(s)

    p <- Parameters(x)(strategies=list(s),
                          seed_rain=pi/2,
                          patch_area=10,
                          is_resident=TRUE)

    expect_that(ebt <- EBT(x)(p),
                throws_error("Patch area must be exactly 1 for the EBT"))

    p$patch_area <- 1.0
    ebt <- EBT(x)(p)
    expect_that(ebt, is_a(sprintf("EBT<%s>", x)))

    ## NOTE: I'm not sure where these are only equal and not identical.
    expect_that(ebt$parameters, equals(p))

    ## Check that the underlying Patch really is a Patch<CohortTop>:
    expect_that(ebt$patch, is_a(sprintf("Patch<%s>", x)))
    expect_that(length(ebt$patch$species), equals(1))
    expect_that(ebt$patch$species[[1]], is_a(sprintf("Species<%s>",x)))
    expect_that(ebt$patch$species[[1]]$seed, is_a(sprintf("Cohort<%s>",x)))
    expect_that(ebt$patch$time, is_identical_to(0.0))

    sched <- ebt$cohort_schedule
    cmp_sched <- make_cohort_schedule(p)
    expect_that(sched$size, equals(cmp_sched$size))
    expect_that(sched$all_times, equals(cmp_sched$all_times))

    ## If the schedule is for the wrong number of species, it should cause
    ## an error...
    sched2 <- CohortSchedule(sched$n_species + 1)
    expect_that(ebt$cohort_schedule <- sched2,
                throws_error("Incorrect length input; expected 1, received 2"))

    ## Build a schedule for 14 introductions from t=0 to t=5
    t <- seq(0, 5, length.out=14)
    sched$set_times(t, 1)
    sched$max_time <- max(t) + diff(t)[[1]]
    ebt$cohort_schedule <- sched

    expect_that(ebt$cohort_schedule$all_times,
                is_identical_to(sched$all_times))
    ## Parameters has been updated:
    expect_that(ebt$parameters$cohort_schedule_times,
                is_identical_to(sched$all_times))

    ## Will be helpful for checking that things worked:
    times <- data.frame(start=t, end=c(t[-1], sched$max_time))

    ## Before starting, check that the EBT is actually empty
    expect_that(ebt$time,                equals(0.0))
    expect_that(ebt$patch$ode_size,      equals(0))
    expect_that(ebt$patch$ode_size,      equals(0))
    expect_that(ebt$complete,            is_false())

    ## TODO: The unlist here is annoying...
    ## Can be resolved by passing off to another utility function I
    ## think, but probably better done via traits (but I don't see how
    ## to make those work).
    i <- ebt$run_next()

    ode_size <- Cohort(x)(strategy_types[[x]]())$ode_size

    expect_that(ebt$cohort_schedule$remaining, equals(length(t) - 1))
    expect_that(ebt$complete, is_false())
    expect_that(i, equals(1))
    ## Note that this is the *second* time; the time of the next
    ## introduction, and the end time of the first introduction.
    expect_that(ebt$time, is_identical_to(times$end[1]))
    expect_that(ebt$time, is_identical_to(times$start[2]))
    expect_that(ebt$patch$ode_size, equals(ode_size))

    ## Trying to set schedule for partly run ebt fails
    expect_that(ebt$cohort_schedule <- sched,
                throws_error("Cannot set schedule without resetting first"))

    i <- unlist(ebt$run_next())
    expect_that(i, equals(1))
    ## EBT ran successfully:
    expect_that(ebt$cohort_schedule$remaining,
                equals(length(t) - 2))
    expect_that(ebt$complete, is_false())
    expect_that(ebt$time, is_identical_to(times$end[2]))
    expect_that(ebt$time, is_identical_to(times$start[3]))
    expect_that(ebt$patch$ode_size, equals(ode_size * 2))

    ## Reset everything
    ## "EBT reset successful"
    ebt$reset()
    expect_that(ebt$time,                equals(0.0))
    expect_that(ebt$time,                equals(0.0))
    expect_that(ebt$patch$ode_size,      equals(0))
    expect_that(ebt$cohort_schedule$remaining, equals(length(t)))

    ## At this point, and possibly before ebt$seed_rains is corrupt.

    ## This is stalling really badly, but it's not totally clear why.
    ## It's *not* the ODE system thrashing (thankfully) because the
    ## number of ODE times reported are not that bad.
    ##
    ## 50.1% in growth_rate_gradient(), and 45.4% in compute_vars_phys()
    ## and 2.8% in initial_conditions() (so that's 98.3%) total.
    ## growth_rate_gradient and initial_conditions spend *all* their
    ## time doing growth_rate_gradient(), in turn all in
    ## compute_assimilation.
    ##
    list_to_matrix <- function(x) {
      n <- max(sapply(x, length))
      t(sapply(x, function(i) c(i, rep(NA, n-length(i)))))
    }
    run_ebt_test <- function(ebt, t_max=Inf) {
      tt <- hh <- NULL
      species_index <- 1L
      ebt$reset()
      while (!ebt$complete > 0 && ebt$time < t_max) {
        ebt$run_next()
        tt <- c(tt, ebt$time)
        hh <- c(hh, list(ebt$patch$species[[species_index]]$height))
      }
      hh <- list_to_matrix(hh)
      list(t=tt, h=hh)
    }

    ## Next, Run the whole schedule using the EBT.
    res_e_1 <- run_ebt_test(ebt)

    ## Then, check that resetting the cohort allows rerunning easily:
    ## EBT can be rerun successfully:
    ebt$reset()
    res_e_2 <- run_ebt_test(ebt)
    expect_that(res_e_2, is_identical_to(res_e_1))

    ## Pull the times out of the EBT and set them in the schedule:
    sched <- ebt$cohort_schedule
    sched$ode_times <- ebt$ode_times
    sched$use_ode_times <- TRUE
    ebt$reset() # must reset
    ebt$cohort_schedule <- sched

    ## So; this does not actually produce *exactly* the same output, which
    ## is very surprising.  It's definitely "close enough" but not exactly
    ## the same (and differs to right around the same order as the patch
    ## vs ebt case).  I suspect that this might come from the difference
    ## between stepping to a point (requiring calculating the step size)
    ## and the stepping a particular step size (requiring calculating the
    ## final time).
    ## EBT with fixed times agrees:
    res_e_3 <- run_ebt_test(ebt)
    expect_that(res_e_3, is_identical_to(res_e_1))

    ## EBT can be rerun successfully with fixed times:
    ebt$reset()
    res_e_4 <- run_ebt_test(ebt)
    expect_that(res_e_4, is_identical_to(res_e_3))
  }
})

test_that("schedule setting", {
  for (x in names(strategy_types)) {
    p <- Parameters(x)(
      strategies=list(strategy_types[[x]]()),
      seed_rain=pi/2,
      is_resident=TRUE,
      cohort_schedule_max_time=5.0)
    ebt <- EBT(x)(p)

    ## Then set a cohort schedule:
    ## Build a schedule for 14 introductions from t=0 to t=5
    sched <- ebt$cohort_schedule
    t <- seq(0, sched$max_time, length.out=14)
    ebt$set_cohort_schedule_times(list(t))

    ## Did set in the EBT:
    expect_that(ebt$cohort_schedule$all_times,
                is_identical_to(list(t)))

    ## And updated in the parameters:
    p2 <- ebt$parameters
    expect_that(p2$cohort_schedule_max_time,
                is_identical_to(sched$max_time))
    expect_that(p2$cohort_schedule_times,
                is_identical_to(list(t)))

    ## Remake the schedule:
    sched2 <- make_cohort_schedule(p2)
    expect_that(sched2$max_time,
                is_identical_to(sched$max_time))
    expect_that(sched2$all_times,
                is_identical_to(list(t)))

    ebt2 <- EBT(x)(p2)
    expect_that(ebt2$cohort_schedule$max_time,
                is_identical_to(sched2$max_time))
    expect_that(ebt2$cohort_schedule$all_times,
                is_identical_to(sched2$all_times))
  }
})

  ## ## TODO: This is a fairly inadequate set of tests; none of the failure
  ## ## conditions are tested, and it's undefined what will happen if we
  ## ## set a cohort schedule that leaves us between introduction points.
  ## test_that("State get/set works", {
  ##   ## Next, try and partly run the EBT, grab its state and push it into a
  ##   ## second copy.
  ##   ebt$reset()
  ##   tmp <- run_ebt_test(ebt, sched$max_time / 2)
  ##   state <- ebt$state

  ##   ebt2 <- new(EBT, ebt$parameters)
  ##   ebt2$state <- state

  ##   expect_that(ebt2$state, equals(ebt$state))
  ##   ## Emergent things:
  ##   expect_that(ebt2$patch$environment$light_environment$xy,
  ##               equals(ebt$patch$environment$light_environment$xy))
  ##   expect_that(ebt2$ode_state, equals(ebt$ode_state))
  ##   expect_that(ebt2$ode_rates,  equals(ebt$ode_rates))
  ##   # TODO: This needs implementing; requires get/set of the ODE solver
  ##   # state.
  ##   # expect_that(ebt2$time,       equals(ebt$time))
  ## })

  ## test_that("Can set times directly", {
  ##   ebt$reset()
  ##   times <- ebt$times(1)
  ##   times2 <- sort(c(times, 0.5*(times[-1] + times[-length(times)])))
  ##   ebt$set_times(times2, 1)
  ##   expect_that(ebt$times(1), is_identical_to(times2))
  ##   expect_that(ebt$cohort_schedule$times(1), is_identical_to(times2))
  ##   ebt$run_next()
  ##   expect_that(ebt$set_times(times, 1), throws_error())
  ##   ebt$reset()
  ##   ebt$set_times(times, 1)
  ##   expect_that(ebt$times(1), is_identical_to(times))
  ## })

test_that("Seed rain & error calculations correct", {
  for (x in names(strategy_types)) {
    p0 <- ebt_base_parameters(x)
    p1 <- expand_parameters(trait_matrix(0.08, "lma"), p0, FALSE)

    ebt <- run_ebt(p1)
    expect_that(ebt, is_a(sprintf("EBT<%s>", x)))

    seed_rain_R <- function(ebt, error=FALSE) {
      a <- ebt$cohort_schedule$times(1)
      d <- ebt$patch$environment$disturbance_regime
      pa <- d$density(a)
      p <- ebt$parameters
      scale <- ebt$parameters$strategies[[1]]$Pi_0 * p$seed_rain
      seeds <- pa * ebt$patch$species[[1]]$seeds * scale
      total <- trapezium(a, seeds)
      if (error) local_error_integration(a, seeds, total) else total
    }

    expect_that(ebt$seed_rain(1), equals(seed_rain_R(ebt)))
    expect_that(ebt$seed_rains, equals(seed_rain_R(ebt)))
    expect_that(ebt$seed_rain_error[[1]],
                equals(seed_rain_R(ebt, error=TRUE)))

    lae_cmp <-
      ebt$patch$species[[1]]$area_leafs_error(ebt$patch$area_leaf_above(0))
    expect_that(ebt$area_leaf_error(1),
                is_identical_to(lae_cmp))

    int <- make_ebt_integrate(ebt)
    Pi_0 <- ebt$parameters$strategies[[1]]$Pi_0
    expect_that(int("seeds_survival_weighted") *
                  Pi_0, equals(ebt$seed_rain(1)))

    res <- run_ebt_collect(p1)
    int2 <- make_ebt_integrate(res)

    expect_that(int2("seeds_survival_weighted"),
                equals(int("seeds_survival_weighted")))
    expect_that(int2("area_leaf"),
                equals(int("area_leaf")))
  }
})

test_that("Can create empty EBT", {
  for (x in names(strategy_types)) {
    p <- Parameters(x)()
    ebt <- EBT(x)(p)

    ## Check light environment is empty:
    env <- ebt$patch$environment
    expect_that(env$light_environment$size, equals(0))
    expect_that(env$canopy_openness(0), equals(1.0))
  }
})
