source("helper-tree.R")

context("EBT")

p <- new(Parameters)
p$add_strategy(new(Strategy))
p$seed_rain <- pi/2

ebt <- new(EBT, p)

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

## Will be helpful for checking that things worked:
times <- data.frame(start=t, end=c(t[-1], sched$max_time))

## Before starting, check that the EBT is actually empty
test_that("EBT starts empty", {
  expect_that(ebt$time,                equals(0.0))
  expect_that(ebt$patch$time,          equals(0.0))
  expect_that(ebt$patch$ode_size,      equals(0))
  expect_that(ebt$patch$n_individuals, equals(0))
})

## Check that we can advance through and add the cohort at time zero.
ebt$run_next()
test_that("EBT adds cohort successfully", {
  expect_that(ebt$cohort_schedule$remaining, equals(length(t) - 1))
  ## Note that this is the *second* time; the time of the next
  ## introduction, and the end time of the first introduction.
  expect_that(ebt$time, is_identical_to(times$end[1]))
  expect_that(ebt$time, is_identical_to(times$start[2]))
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
  expect_that(ebt$time, is_identical_to(times$end[2]))
  expect_that(ebt$time, is_identical_to(times$start[3]))
  expect_that(ebt$patch$ode_size, equals(4 * 2))
})

## Reset everything
ebt$reset()
test_that("EBT reset successful", {
  expect_that(ebt$time,                equals(0.0))
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

run.ebt <- function(ebt) {
  tt <- hh <- NULL
  ebt$reset()
  while (ebt$cohort_schedule$remaining > 0) {
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
  expect_that(res.e.1$h, equals(hh.p.end, tolerance=1e-11))
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
  expect_that(res.e.3, equals(res.e.1, tolerance=1e-12))
})

ebt$reset()
res.e.4 <- run.ebt(ebt)
test_that("EBT can be rerun successfully with fixed times", {
  expect_that(res.e.4, is_identical_to(res.e.3))
})

if (interactive()) {
  matplot(tt.p.end, hh.p.end, type="o", col="black", pch=1, lty=1)
  matpoints(tt.1.e, hh.1.e, col="red", cex=.5, pch=1, type="o", lty=1)
}

rm(ebt)
gc()
