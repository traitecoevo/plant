if (interactive()) {
  devtools::load_all("../../")
  library(testthat)
  source("helper-tree2.R")
}

context("EBT")

test_that("Ported from tree1", {
  p <- Parameters(strategies=list(Strategy()),
                  seed_rain=pi/2,
                  patch_area=10,
                  is_resident=TRUE)

  expect_that(ebt <- EBT(p),
              throws_error("Patch area must be exactly 1 for the EBT"))

  p$patch_area <- 1.0
  ebt <- EBT(p)
  expect_that(ebt, is_a("EBT"))

  ## NOTE: I'm not sure where these are only equal and not identical.
  expect_that(ebt$parameters, equals(p))

  ## Check that the underlying Patch really is a Patch<CohortTop>:
  expect_that(ebt$patch, is_a("Patch"))
  expect_that(length(ebt$patch$species), equals(1))
  expect_that(ebt$patch$species[[1]], is_a("Species"))
  expect_that(ebt$patch$species[[1]]$seed, is_a("Cohort"))

  sched <- ebt$cohort_schedule
  expect_that(sched$size, equals(0))
  expect_that(sched$next_event, throws_error("All events completed"))

  ## If the schedule is for the wrong number of species, it should cause
  ## an error...
  sched2 <- CohortSchedule(sched$n_species + 1)
  expect_that(ebt$cohort_schedule <- sched2,
              throws_error("Incorrect length input; expected 1, received 2"))

  ## Build a schedule for 14 introductions from t=0 to t=5
  t <- seq(0, 5, length=14)
  sched$set_times(t, 1)
  sched$max_time <- max(t) + diff(t)[[1]]
  ebt$cohort_schedule <- sched

  expect_that(sched$times(1),
              is_identical_to(ebt$cohort_schedule$times(1)))

  ## Will be helpful for checking that things worked:
  times <- data.frame(start=t, end=c(t[-1], sched$max_time))

  ## Before starting, check that the EBT is actually empty
  expect_that(ebt$time,                equals(0.0))
  expect_that(ebt$patch$ode_size,      equals(0))
  expect_that(ebt$patch$ode_size,      equals(0))
  expect_that(ebt$complete,            is_false())

  ## TODO: The unlist here is annoying...
  i <- unlist(ebt$run_next())

  expect_that(ebt$cohort_schedule$remaining, equals(length(t) - 1))
  expect_that(ebt$complete, is_false())
  expect_that(i, equals(1))
  ## Note that this is the *second* time; the time of the next
  ## introduction, and the end time of the first introduction.
  expect_that(ebt$time, is_identical_to(times$end[1]))
  expect_that(ebt$time, is_identical_to(times$start[2]))
  expect_that(ebt$patch$ode_size, equals(4))

  ## Trying to set schedule for partly run ebt fails
  expect_that(ebt$cohort_schedule <- sched,
              throws_error("Cannot set schedule without resetting first"))

  ebt$run_next()
  ## EBT ran successfully:
  expect_that(ebt$cohort_schedule$remaining,
              equals(length(t) - 2))
  expect_that(ebt$complete, is_false())
  expect_that(ebt$time, is_identical_to(times$end[2]))
  expect_that(ebt$time, is_identical_to(times$start[3]))
  expect_that(ebt$patch$ode_size, equals(4 * 2))

  ## Reset everything
  ## "EBT reset successful"
  ebt$reset()
  expect_that(ebt$time,                equals(0.0))
  expect_that(ebt$patch$time,          equals(0.0))
  expect_that(ebt$patch$ode_size,      equals(0))
  expect_that(ebt$cohort_schedule$remaining, equals(length(t)))

  ## Run the whole schedule using a Patch<CohortTop>, manually moving
  ## things along the schedule.
  patch <- Patch(p)
  sched$reset()
  species_index <- 1 # for getting state out.

  tt_p_start <- hh_p_start <- tt_p_end <- hh_p_end <- NULL
  solver <- solver_from_ode_target(patch, p$control$ode_control)
  while (sched$remaining > 0) {
    e <- sched$next_event
    if (!identical(patch$time, e$time_introduction))
      stop("Something terrible has happened")
    patch$add_seedling(e$species_index)

    ## Harvest statistics at start of step
    tt_p_start <- c(tt_p_start, patch$time)
    hh_p_start <- c(hh_p_start, list(patch$height[[species_index]]))

    ## Advance the solution
    solver$set_state(patch$ode_values, patch$time)
    solver$advance(e$time_end)
    sched$pop()

    ## Harvest statistics and end of step
    tt_p_end <- c(tt_p_end, patch$time)
    hh_p_end <- c(hh_p_end, list(patch$height[[species_index]]))
  }

  list_to_matrix <- function(x) {
    n <- max(sapply(x, length))
    t(sapply(x, function(i) c(i, rep(NA, n-length(i)))))
  }

  hh_p_start <- list_to_matrix(hh_p_start)
  hh_p_end <- list_to_matrix(hh_p_end)


  expect_that(nrow(hh_p_start), equals(length(t)))
  ## Start point is the leaf plant height:
  expect_that(diag(hh_p_start),
              equals(rep(new(Plant, p[[1]])$height, length(t))))

  run_ebt_test <- function(ebt, t_max=Inf) {
    tt <- hh <- NULL
    ebt$reset()
    while (!ebt$complete > 0 && ebt$time < t_max) {
      ebt$run_next()
      tt <- c(tt, ebt$time)
      hh <- c(hh, list(ebt$patch$height[[species_index]]))
    }
    hh <- list_to_matrix(hh)
    list(t=tt, h=hh)
  }

  ## Next, Run the whole schedule using the EBT.
  res_e_1 <- run_ebt_test(ebt)

  ## I'm actually quite surprised that the objects aren't identical.
  ## I've left the tolerance super strict here.
  ## EBT and Patch agree:
  expect_that(res_e_1$t, is_identical_to(tt_p_end))
  expect_that(res_e_1$h, equals(hh_p_end, tolerance=3e-11))
  expect_that(ebt$ode_values,
              equals(patch$ode_values, tolerance=2e-9))
  expect_that(ebt$ode_rates,
              equals(patch$ode_rates,  tolerance=1e-8))

  ## Then, check that resetting the cohort allows rerunning easily:
  ## EBT can be rerun successfully:
  ebt$reset()
  res_e_2 <- run_ebt_test(ebt)
  expect_that(res_e_2, is_identical_to(res_e_1))

  ## Pull the times out of the EBT and set them in the schedule:
  sched <- ebt$cohort_schedule
  sched$ode_times <- ebt$ode_times
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
  expect_that(res_e_3, equals(res_e_1, tolerance=1e-10))

  ## EBT can be rerun successfully with fixed times:
  ebt$reset()
  res_e_4 <- run_ebt_test(ebt)
  expect_that(res_e_4, is_identical_to(res_e_3))
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
##   expect_that(ebt2$ode_values, equals(ebt$ode_values))
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

## test_that("Seed rain & error calculations correct", {
##   p <- new(Parameters)
##   p$add_strategy(new(Strategy))
##   p$set_control_parameters(fast_control())
##   p$set_parameters(list(patch_area=1.0))   # See issue #13
##   p$seed_rain <- 1.1

##   ebt <- new(EBT, p)
##   ebt$reset()
##   ebt$cohort_schedule <- default_cohort_schedule(p)
##   ebt$run()

##   seed_rain.R <- function(ebt, error=FALSE) {
##     a <- ebt$cohort_schedule$times(1)
##     d <- ebt$patch$disturbance_regime
##     pa <- sapply(a, function(ai) d$density(ai))
##     p <- ebt$parameters
##     scale <- p$parameters[["Pi_0"]] * p$seed_rain
##     seeds <- pa * ebt$patch$species[[1]]$seeds * scale
##     total <- trapezium(a, seeds)
##     if (error) local_error_integration(a, seeds, total) else total
##   }

##   expect_that(ebt$seed_rain(1), equals(seed_rain.R(ebt)))
##   expect_that(ebt$seed_rains, equals(seed_rain.R(ebt)))
##   expect_that(ebt$seed_rain_error(1),
##               equals(seed_rain.R(ebt, error=TRUE)))

##   expect_that(ebt$seed_rain(0), throws_error())
##   expect_that(ebt$seed_rain(2), throws_error())

##   lae.cmp <-
##     ebt$patch$species[[1]]$leaf_area_error(ebt$patch$leaf_area_above(0))
##   expect_that(ebt$leaf_area_error(1),
##               is_identical_to(lae.cmp))
## })

## test_that("Can create empty EBT", {
##   p <- new(Parameters)
##   p$set_parameters(list(patch_area=1.0))
##   ebt <- new(EBT, p)

##   ## Check light environment is empty:
##   env <- ebt$patch$environment
##   expect_that(env$light_environment$size, equals(0))
##   expect_that(env$canopy_openness(0), equals(1.0))
## })

## test_that("Can create empty EBT with mutants", {
##   p <- new(Parameters)
##   p$set_parameters(list(patch_area=1.0))   # See issue #13
##   p$set_control_parameters(fast_control()) # A bit faster
##   p$add_strategy_mutant(new(Strategy, list(lma=0.1)))

##   t_max <- 10
##   schedule0 <- default_cohort_schedule(p, t_max)

##   ebt <- run_ebt(p, schedule0)
##   expect_that(ebt$time, equals(t_max))

##   ## Check light environment is empty:
##   env <- ebt$patch$environment
##   expect_that(env$light_environment$size, equals(0))
##   expect_that(env$canopy_openness(0), equals(1.0))
## })
