source("helper-tree.R")

context("EBT")

p <- new(Parameters)
p$add_strategy(new(Strategy))
p$seed_rain <- pi/2

expect_that(ebt <- new(EBT, p), throws_error())

p$set_parameters(list(patch_area=1.0))
ebt <- new(EBT, p)

test_that("Parameters can be pulled from EBT", {
  p.cmp <- ebt$parameters
  expect_that(p.cmp$size, is_identical_to(p$size))
  expect_that(p.cmp$parameters, is_identical_to(p$parameters))
  expect_that(p.cmp[[1]]$parameters, is_identical_to(p[[1]]$parameters))
})

## Check that the underlying Patch really is a Patch<CohortTop>, not
## one of the other Patch types (which would just fail miserably below
## here, if they even compile).
expect_that(inherits(ebt$patch, "Rcpp_PatchCohortTop"), is_true())

sched <- ebt$cohort_schedule
test_that("Empty CohortSchedule", {
  expect_that(sched$size, equals(0))
  expect_that(sched$next_event, throws_error())
})

## If the schedule is for the wrong number of species, it should cause
## an error...
sched2 <- new(CohortSchedule, sched$n_species + 1)
expect_that(ebt$cohort_schedule <- sched2, throws_error())

## Build a schedule for 14 introductions from t=0 to t=5
t <- seq(0, 5, length=14)
sched$set_times(t, 1)
sched$max_time <- max(t) + diff(t)[[1]]
ebt$cohort_schedule <- sched

expect_that(sched$times(1), is_identical_to(ebt$times(1)))

## Will be helpful for checking that things worked:
times <- data.frame(start=t, end=c(t[-1], sched$max_time))

## Before starting, check that the EBT is actually empty
test_that("EBT starts empty", {
  expect_that(ebt$time,                equals(0.0))
  expect_that(ebt$ode_size,            equals(0))
  expect_that(ebt$patch$time,          equals(0.0))
  expect_that(ebt$patch$ode_size,      equals(0))
  expect_that(ebt$patch$n_individuals, equals(0))
})

## Check that we can advance through and add the cohort at time zero.
ebt$run_next()
test_that("EBT adds cohort successfully", {
  expect_that(ebt$cohort_schedule$remaining, equals(length(t) - 1))
  expect_that(ebt$complete, is_false())
  ## Note that this is the *second* time; the time of the next
  ## introduction, and the end time of the first introduction.
  expect_that(ebt$time, is_identical_to(times$end[1]))
  expect_that(ebt$time, is_identical_to(times$start[2]))
  expect_that(ebt$ode_size,       equals(4))
  expect_that(ebt$patch$ode_size, equals(4))
})

test_that("Trying to set schedule for partly run ebt fails", {
  expect_that(ebt$cohort_schedule <- sched, throws_error())
})

## Introduce the second cohort:
ebt$run_next()
test_that("EBT ran successfully", {
  expect_that(ebt$cohort_schedule$remaining,
              equals(length(t) - 2))
  expect_that(ebt$complete, is_false())
  expect_that(ebt$time, is_identical_to(times$end[2]))
  expect_that(ebt$time, is_identical_to(times$start[3]))
  expect_that(ebt$ode_size,       equals(4 * 2))
  expect_that(ebt$patch$ode_size, equals(4 * 2))
})

## Reset everything
ebt$reset()
test_that("EBT reset successful", {
  expect_that(ebt$time,                equals(0.0))
  expect_that(ebt$ode_size,            equals(0))
  expect_that(ebt$patch$time,          equals(0.0))
  expect_that(ebt$patch$ode_size,      equals(0))
  expect_that(ebt$patch$n_individuals, equals(0))
  expect_that(ebt$cohort_schedule$remaining, equals(length(t)))
})

## Run the whole schedule using a Patch<CohortTop>, manually moving
## things along the schedule.
patch <- new(PatchCohortTop, p)
sched$reset()

species.index <- 1 # for getting state out.

tt.p.start <- hh.p.start <- tt.p.end <- hh.p.end <- NULL
solver <- solver.from.odetarget(patch, p$control$ode_control)
while (sched$remaining > 0) {
  e <- sched$next_event
  if (!identical(patch$time, e$time_introduction))
    stop("Something terrible has happened")
  patch$add_seedling(e$species_index)

  ## Harvest statistics at start of step
  tt.p.start <- c(tt.p.start, patch$time)
  hh.p.start <- c(hh.p.start, list(patch$height[[species.index]]))

  ## Advance the solution
  solver$set_state(patch$ode_values, patch$time)
  solver$advance(e$time_end)
  sched$pop()

  ## Harvest statistics and end of step
  tt.p.end <- c(tt.p.end, patch$time)
  hh.p.end <- c(hh.p.end, list(patch$height[[species.index]]))
}

list.to.matrix <- function(x) {
  n <- max(sapply(x, length))
  t(sapply(x, function(i) c(i, rep(NA, n-length(i)))))
}

hh.p.start <- list.to.matrix(hh.p.start)
hh.p.end <- list.to.matrix(hh.p.end)

test_that("Run looks successful", {
  expect_that(nrow(hh.p.start), equals(length(t)))
  ## Start point is the leaf plant height:
  expect_that(diag(hh.p.start),
              equals(rep(new(Plant, p[[1]])$height, length(t))))
})

run.ebt <- function(ebt, t.max=Inf) {
  tt <- hh <- NULL
  ebt$reset()
  while (!ebt$complete > 0 && ebt$time < t.max) {
    ebt$run_next()
    tt <- c(tt, ebt$time)
    hh <- c(hh, list(ebt$patch$height[[species.index]]))
  }
  hh <- list.to.matrix(hh)
  list(t=tt, h=hh)
}

## Next, Run the whole schedule using the EBT.
res.e.1 <- run.ebt(ebt)

## I'm actually quite surprised that the objects aren't identical.
## I've left the tolerance super strict here.
test_that("EBT and Patch agree", {
  expect_that(res.e.1$t, is_identical_to(tt.p.end))
  expect_that(res.e.1$h, equals(hh.p.end, tolerance=3e-11))
  expect_that(ebt$ode_values,
              equals(patch$ode_values, tolerance=2e-9))
  expect_that(ebt$ode_rates,
              equals(patch$ode_rates,  tolerance=1e-8))
})

## Then, check that resetting the cohort allows rerunning easily:
ebt$reset()
res.e.2 <- run.ebt(ebt)
test_that("EBT can be rerun successfully", {
  expect_that(res.e.2, is_identical_to(res.e.1))
})

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
res.e.3 <- run.ebt(ebt)
test_that("EBT with fixed times agrees", {
  expect_that(res.e.3, equals(res.e.1, tolerance=4e-12))
})

ebt$reset()
res.e.4 <- run.ebt(ebt)
test_that("EBT can be rerun successfully with fixed times", {
  expect_that(res.e.4, is_identical_to(res.e.3))
})

## TODO: This is a fairly inadequate set of tests; none of the failure
## conditions are tested, and it's undefined what will happen if we
## set a cohort schedule that leaves us between introduction points.
test_that("State get/set works", {
  ## Next, try and partly run the EBT, grab its state and push it into a
  ## second copy.
  ebt$reset()
  tmp <- run.ebt(ebt, sched$max_time / 2)
  state <- ebt$state

  ebt2 <- new(EBT, ebt$parameters)
  ebt2$state <- state

  expect_that(ebt2$state, equals(ebt$state))
  ## Emergent things:
  expect_that(ebt2$patch$environment$light_environment$xy,
              equals(ebt$patch$environment$light_environment$xy))
  expect_that(ebt2$ode_values, equals(ebt$ode_values))
  expect_that(ebt2$ode_rates,  equals(ebt$ode_rates))
  # TODO: This needs implementing; requires get/set of the ODE solver
  # state.
  # expect_that(ebt2$time,       equals(ebt$time))
})

test_that("Can set times directly", {
  ebt$reset()
  times <- ebt$times(1)
  times2 <- sort(c(times, 0.5*(times[-1] + times[-length(times)])))
  ebt$set_times(times2, 1)
  expect_that(ebt$times(1), is_identical_to(times2))
  expect_that(ebt$cohort_schedule$times(1), is_identical_to(times2))
  ebt$run_next()
  expect_that(ebt$set_times(times, 1), throws_error())
  ebt$reset()
  ebt$set_times(times, 1)
  expect_that(ebt$times(1), is_identical_to(times))
})

gc() # hide the "signalCondition" Rcpp issue

test_that("Fitness & error calculations correct", {
  p <- new(Parameters)
  p$add_strategy(new(Strategy))
  p$set_control_parameters(fast.control())
  p$set_parameters(list(patch_area=1.0))   # See issue #13
  p$seed_rain <- 1.1
  t.max <- p$disturbance$cdf(tree:::reference.pr.survival.eps)

  ebt <- new(EBT, p)
  ebt$reset()
  ebt$cohort_schedule <- default.schedule(30, t.max)
  ebt$run()

  fitness.R <- function(ebt, error=FALSE) {
    a <- ebt$cohort_schedule$times(1)
    d <- ebt$patch$disturbance_regime
    pa <- sapply(a, function(ai) d$density(ai))
    p <- ebt$parameters
    scale <- p$parameters[["Pi_0"]] * p$seed_rain
    seeds <- pa * ebt$patch$species[[1]]$seeds * scale
    total <- trapezium(a, seeds)
    if (error) local_error_integration(a, seeds, total) else total
  }

  expect_that(ebt$fitness(1), equals(fitness.R(ebt)))
  expect_that(ebt$fitnesses, equals(fitness.R(ebt)))
  expect_that(ebt$fitness_error(1),
              equals(fitness.R(ebt, error=TRUE)))

  expect_that(ebt$fitness(0), throws_error())
  expect_that(ebt$fitness(2), throws_error())

  lae.cmp <-
    ebt$patch$species[[1]]$leaf_area_error(ebt$patch$leaf_area_above(0))
  expect_that(ebt$leaf_area_error(1),
              is_identical_to(lae.cmp))
})

rm(ebt)
gc()
